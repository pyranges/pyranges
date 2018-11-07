import pandas as pd
from natsort import natsorted
import numpy as np
from collections import defaultdict, OrderedDict

import logging

import ray

try:
    ray.init(local_mode=True, logging_level=logging.CRITICAL) # logging_level=logging.CRITICAL # local_mode=True
except:
    pass # Ray was already initialized

from tabulate import tabulate


from pyranges.genomicfeatures import GenomicFeaturesMethods
from pyranges.subset import get_string, get_slice, get_tuple
# from pyranges.methods import _cluster, _subtraction, _set_union, _set_intersection, _intersection, _nearest, _coverage, _overlap_write_both, _overlap, _tss, _tes, _jaccard, _lengths, _slack
from pyranges.multithreaded import (_cluster, pyrange_apply_single,
                                    _intersection, pyrange_apply)

def return_copy_if_view(df, df2):
    # https://stackoverflow.com/questions/26879073/checking-whether-data-frame-is-copy-or-view-in-pandas

    if df.values.base is df2.values.base:
        return df.copy(deep=True)
    else:
        return df


def create_df_dict(df):

    if "Strand" in df:
        grpby_key = "Chromosome Strand".split()
    else:
        grpby_key = "Chromosome"

    # {k: ray.put(v) for k, v in natsorted(df.groupby(grpby_key))}
    return {k: ray.put(v) for k, v in df.groupby(grpby_key)}


def create_pyranges_df(seqnames, starts, ends, strands=None):

    if isinstance(seqnames, str):
        seqnames = pd.Series([seqnames] * len(starts), dtype="category")

    if strands != None and strands != False:

        if isinstance(strands, str):
            strands = pd.Series([strands] * len(starts), dtype="category")

        columns = [seqnames, starts, ends, strands]
        lengths = list(str(len(s)) for s in columns)
        assert len(set(lengths)) == 1, "seqnames, starts, ends and strands must be of equal length. But are {}".format(", ".join(lengths))
        colnames = "Chromosome Start End Strand".split()
    else:
        columns = [seqnames, starts, ends]
        lengths = list(str(len(s)) for s in columns)
        assert len(set(lengths)) == 1, "seqnames, starts and ends must be of equal length. But are {}".format(", ".join(lengths))
        colnames = "Chromosome Start End".split()

    idx = range(len(starts))
    series_to_concat = []
    for s in columns:
        s = pd.Series(s, index=idx)
        series_to_concat.append(s)

    df = pd.concat(series_to_concat, axis=1)
    df.columns = colnames

    return df


def return_empty_if_one_empty(func):

    def extended_func(self, other, *args, **kwargs):

        if len(self) == 0 or len(other) == 0:
            df = pd.DataFrame(columns="Chromosome Start End Strand".split())
        else:
            df = func(self, other, *args, **kwargs)

        return df

    return extended_func


def return_empty_if_both_empty(func):

    def extended_func(self, other, *args, **kwargs):

        if len(self) == 0 and len(other) == 0:
            df = pd.DataFrame(columns="Chromosome Start End Strand".split())
        else:
            df = func(self, other, *args, **kwargs)

        return df

    return extended_func


def pyrange_or_df(func):

    def extension(self, other, *args, **kwargs):
        df = func(self, other, *args, **kwargs)

        if kwargs.get("df"):
            return df

        return PyRanges(df)

    return extension


def pyrange_or_df_single(func):

    def extension(self, *args, **kwargs):
        df = func(self, *args, **kwargs)

        if kwargs.get("df"):
            return df

        return PyRanges(df)

    return extension



class PyRanges():

    dfs = None
    gf = None

    def __init__(self, df=None, seqnames=None, starts=None, ends=None, strands=None, copy_df=True):

        if df is False or df is None:
            df = create_pyranges_df(seqnames, starts, ends, strands)
            self.__dict__["dfs"] = create_df_dict(df)
        elif isinstance(df, pd.DataFrame):
            self.__dict__["dfs"] = create_df_dict(df)
        else:
            # df is actually dict of dfs
            self.__dict__["dfs"] = df

        self.__dict__["ft"] = GenomicFeaturesMethods(self)


    def __len__(self):
        return sum(len(ray.get(d)) for d in self.values)

    def __setattr__(self, column_name, column):

        if column_name in "Chromosome Start End Strand".split():
            raise Exception("The columns Chromosome, Start, End or Strand can not be reset.")
        if column_name == "stranded":
            raise Exception("The stranded attribute is read-only. Create a new PyRanges object instead.")

        if not isinstance(column, str):
            if not len(self) == len(column):
                raise Exception("DataFrame and column must be same length.")

            column_to_insert = pd.Series(column, index=self.df.index)
        else:
            column_to_insert = pd.Series(column, index=self.df.index)

        pos = self.df.shape[1]
        if column_name in self.df:
            pos = list(self.df.columns).index(column_name)
            self.df.drop(column_name, inplace=True, axis=1)

        self.df.insert(pos, column_name, column_to_insert)



    def __eq__(self, other):

        return self.df.equals(other.df)

    def __getitem__(self, val):

        if isinstance(val, str):
            df = get_string(self, val)
        elif isinstance(val, tuple):
            df = get_tuple(self, val)
        elif isinstance(val, slice):
            df = get_slice(self, val)
        else:
            raise Exception("Not valid subsetter: {}".format(str(val)))

        return df


    def __str__(self):

        if len(self) == 0:
            return "Empty PyRanges"

        # keys = natsorted(list(self.dfs.keys()))
        if len(self.keys) == 1:

            first_key = self.keys[0]
            # last_key = list(self.dfs.keys())[-1]
            first_df = ray.get(self.dfs[first_key])

            # last_df = ray.get(self.dfs[last_key]).tail(3)
            h = first_df.head(3).astype(object)
            m = first_df.head(1).astype(object)
            t = first_df.tail(3).astype(object)
            m.loc[:,:] = "..."

            if len(self) > 6:
                s = pd.concat([h, m, t])
            elif len(self) == 6:
                s = pd.concat([h, t])
            else:
                s = h
        else:
            keys = self.keys
            first_key = keys[0]
            last_key = keys[-1]
            # first_key = self.keys[0]
            # last_key = self.keys[-1]
            first_df = ray.get(self.dfs[first_key]).head(3)
            last_df = ray.get(self.dfs[last_key]).tail(3)
            # last_df = self.dfs[list(self.dfs.keys())[-1]].tail(3)

            h = first_df.head(3).astype(object)
            m = first_df.head(1).astype(object)
            t = last_df.head(3).astype(object)
            m.loc[:,:] = "..."
            # m.index = ["..."]
            if len(self) > 6:

                # iterate from front until have three
                heads = []
                hl = 0
                for k in keys:
                    h = ray.get(self.dfs[k]).head(3)
                    first_df = h
                    hl += len(h)
                    heads.append(h)
                    if hl >= 3:
                        break

                tails = []
                tl = 0
                for k in keys[::-1]:
                    t = ray.get(self.dfs[k]).tail(3)
                    tl += len(t)
                    tails.append(t)
                    if tl >= 3:
                        break
                # iterate from back until have three

                h = pd.concat(heads).head(3).astype(object)
                t = pd.concat(tails).tail(3).astype(object)
                m = h.head(1).astype(object)

                m.loc[:,:] = "..."
                s = pd.concat([h, m, t])

            else:

                h = pd.concat(ray.get(self.dfs.values()))
                first_df = h
                s = h.astype(object)

        h = [c + "\n(" + str(t) + ")" for c, t in  zip(h.columns, first_df)]

        str_repr = tabulate(s, headers=h, tablefmt='psql', showindex=False) + \
                                        "\nPyRanges object has {} sequences from {} chromosomes.".format(len(self), len(self.dfs.keys()))
        return str_repr


    def __repr__(self):

        return str(self)


    @pyrange_or_df
    @return_empty_if_one_empty
    def overlap(self, other, strandedness=False, invert=False, how=None, **kwargs):

        "Want all intervals in self that overlap with other."

        df = _overlap(self, other, strandedness, invert, how)

        return df

    @pyrange_or_df
    @return_empty_if_one_empty
    def nearest(self, other, strandedness=False, suffix="_b", how=None, overlap=True, **kwargs):

        "Find the nearest feature in other."

        df = _nearest(self, other, strandedness, suffix, how, overlap)

        return df

    @return_empty_if_one_empty
    def intersection(self, other, strandedness=False, how=None):


        dfs = pyrange_apply(_intersection, self, other, strandedness=strandedness, how=how)

        return PyRanges(dfs)

    @return_empty_if_one_empty
    def set_intersection(self, other, strandedness=False, how=None):

        from pyranges.multithreaded import _set_intersection, pyrange_apply

        self_clusters = self.cluster()
        other_clusters = other.cluster()
        dfs = pyrange_apply(_set_intersection, self_clusters, other_clusters, strandedness=strandedness, how=how)

        # si = _set_intersection(self, other, strandedness, how)

        return PyRanges(dfs)

    @pyrange_or_df
    @return_empty_if_both_empty
    def set_union(self, other, strand=False):

        si = _set_union(self, other, strand)

        return si


    @pyrange_or_df
    def subtraction(self, other, strandedness=False):


        return _subtraction(self, other, strandedness)


    @pyrange_or_df
    @return_empty_if_one_empty
    def join(self, other, strandedness=False, new_pos=None, suffixes=["_a", "_b"], suffix="_b", how=None):

        df = _overlap_write_both(self, other, strandedness, new_pos, suffixes, suffix, how)

        return df


    def cluster(self, strand=None, max_dist=0, min_nb=1, **kwargs):

        df = pyrange_apply_single(_cluster, self, strand)

        return PyRanges(df)


    def coverage(self, value_col=None, stranded=False):

        return _coverage(self, value_col, stranded=stranded)


    @pyrange_or_df_single
    def slack(self, slack):

        return _slack(self, slack)



    @pyrange_or_df_single
    def tssify(self, slack=0):

        if not self.stranded:
            raise Exception("Cannot compute TSSes without strand info. Perhaps use slack() instead?")

        return _tss(self, slack)


    @pyrange_or_df_single
    def tesify(self, slack=0):

        return _tes(self, slack)

    def pos(self, *val):

        # turn int-tuple into slice
        newval = []
        for v in val:
            if isinstance(v, tuple) and len(v) == 2 and isinstance(v[0], int) and isinstance(v[1], int):
                newval.append(slice(v[0], v[1]))
            else:
                newval.append(v)

        val = tuple(newval)

        if len(val) == 1:
            val = val[0]

        if isinstance(val, str):
            df = get_string(self, val)
        elif isinstance(val, tuple):
            df = get_tuple(self, val)
        elif isinstance(val, slice):
            df = get_slice(self, val)
        else:
            raise Exception("Invalid type for indexer. Must be str, slice, or 2/3-tuple.")

        return return_copy_if_view(df, self.df)


    # always returns df
    def get(self, column, values, **kwargs):

        df = self.df
        if not column in self.df:
            raise Exception("Column {} not in PyRanges".format(column))

        if isinstance(values, str):
            df = df.loc[df[column] == values]
        else:
            df = df.loc[df[column].isin(values)]

        return return_copy_if_view(df, self.df)


    def lengths(self, **kwargs):

        df = _lengths(self)

        return df


    def jaccard(self, other, strandedness):

        return _jaccard(self, other, strandedness)


    @pyrange_or_df_single
    def sort(self, strand=True):

        if strand:
            return self.df.sort_values("Chromosome Strand".split())
        else:
            return self.df.sort_values("Chromosome")

    @property
    def keys(self):
        return natsorted(self.dfs.keys())

    @property
    def stranded(self):
        return len(list(self.keys)[0]) == 2

    @property
    def strands(self):

        if not self.stranded:
            raise Exception("PyRanges not stranded!")

        return natsorted(set([k[1] for k in self.keys]))

    @property
    def chromosomes(self):

        return natsorted(set([k[0] for k in self.keys]))

    @property
    def items(self):

        return natsorted([(k, df) for (k, df) in self.dfs.items()])

    @property
    def values(self):

        return [df for k, df in self.items]

    def as_df(self):

        return pd.concat(self.values)
