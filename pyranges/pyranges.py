import pandas as pd
from natsort import natsorted
import numpy as np
from collections import defaultdict, OrderedDict

import logging

import pyranges as pr


from tabulate import tabulate


from pyranges.genomicfeatures import GenomicFeaturesMethods
from pyranges.subset import get_string, get_slice, get_tuple
from pyranges.multithreaded import (_cluster, pyrange_apply_single,
                                    _write_both, _coverage,
                                    _intersection, pyrange_apply, _nearest,
                                    _overlap, _first_df, _subtraction, _tss,
                                    _tes, _slack, _sort, merge_dfs, _concat,
                                    _index_as_col, _relative_distance)

def fill_kwargs(kwargs):

    if not "strandedness" in kwargs:
        kwargs["strandedness"] = None
    if not "suffix" in kwargs:
        kwargs["suffix"] = "_b"
    if not "overlap" in kwargs:
        kwargs["overlap"] = True
    if not "how" in kwargs:
        kwargs["how"] = None
    if not "invert" in kwargs:
        kwargs["invert"] = None
    if not "new_pos" in kwargs:
        kwargs["new_pos"] = None
    if not "suffixes" in kwargs:
        kwargs["suffixes"] = ["_a", "_b"]
    if not "suffix" in kwargs:
        kwargs["suffixes"] = ["_b"]
    if not "sparse" in kwargs:
        kwargs["sparse"] = {"self": False, "other": False}

    return kwargs

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

    return {k: v for k, v in df.groupby(grpby_key)}


def create_pyranges_df(seqnames, starts, ends, strands=None):

    if isinstance(seqnames, str):
        seqnames = pd.Series([seqnames] * len(starts), dtype="category")

    if strands is not None:

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
        if isinstance(s, pd.Series):
            s = pd.Series(s.values, index=idx)
        else:
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






def set_dtypes(df, extended):

    if not extended:
        dtypes = {"Start": np.int32, "End": np.int32, "Chromosome": "category", "Strand": "category"}
    else:
        dtypes = {"Start": np.int64, "End": np.int64, "Chromosome": "category", "Strand": "category"}

    if not "Strand" in df:
        del dtypes["Strand"]


    if "ExonNumber" in df:
        df.ExonNumber = df.ExonNumber.fillna(-1)
        if df.ExonNumber.max() < 128:
            df.ExonNumber = df.ExonNumber.astype(np.int8)
        else:
            df.ExonNumber = df.ExonNumber.astype(np.int16)

    df.columns = (c.replace("#Chromosome", "Chromosome") for c in df.columns)

    for col, dtype in dtypes.items():

        if df[col].dtype.name != dtype:
            df[col] = df[col].astype(dtype)

    return df

def read_path(p):

    p = p.lower()
    if p.endswith((".gtf", ".gff", ".gtf.gz", ".gff.gz")):
        df = pr.read_gtf(p, output_df=True)
    elif p.endswith((".bed", ".bed.gz")):
        df = pr.read_bed(p, output_df=True)
    else:
        df = None
        # df = pr.

    return df


class PyRanges():

    dfs = None
    gf = None

    def __init__(self, df=None, seqnames=None, starts=None, ends=None, strands=None, copy_df=False, extended=False):

        if isinstance(df, PyRanges):
            raise Exception("Object is already a PyRange.")

        if copy_df:
            df = df.copy()

        if df is False or df is None:
            df = create_pyranges_df(seqnames, starts, ends, strands)

        if isinstance(df, str):
            df = read_path(df)

        if isinstance(df, pd.DataFrame):
            df = set_dtypes(df, extended)

        if isinstance(df, pd.DataFrame):
            self.__dict__["dfs"] = create_df_dict(df)
        else:
            # df is actually dict of dfs
            self.__dict__["dfs"] = df

        self.__dict__["ft"] = GenomicFeaturesMethods(self)


    def __len__(self):
        return sum([ len(d) for d in self.values() ])

    def __call__(self, eval_str):
        return self.eval(eval_str)

    def __getattr__(self, name):

        if name in self.values()[0]:
            return pd.concat([df[name] for df in self.values()])
        else:
            raise Exception("PyRanges object has no attribute", name)


    def __setattr__(self, column_name, column):

        if column_name in "Chromosome Strand".split():
            raise Exception("The columns Chromosome and Strand can not be reset.")

        isiterable = isinstance(column, list) or isinstance(column, pd.Series) or isinstance(column, np.ndarray)
        if isiterable:
            if not len(self) == len(column):
                raise Exception("DataFrame and column must be same length.")

        already_exists = column_name in self.values()[0]

        if already_exists:
            pos = list(self.values()[0].columns).index(column_name)
        else:
            pos = self.values()[0].shape[1]

        start_length, end_length = 0, 0

        dfs = {}
        for k, df in self.items():

            end_length += len(df)

            if already_exists:
                df = df.drop(column_name, axis=1)

            if isiterable:
                df.insert(pos, column_name, column[start_length:end_length])
            else:
                df.insert(pos, column_name, column)

            start_length = end_length

            dfs[k] = df

        self.__dict__["dfs"] = dfs


    def __eq__(self, other):

        return self.df.equals(other.df)

    def __getitem__(self, val):

        # print("self, val")
        # print(self, val)
        if isinstance(val, str):
            df = get_string(self, val)
        elif isinstance(val, tuple):
            df = get_tuple(self, val)
        elif isinstance(val, slice):
            df = get_slice(self, val)
        else:
            raise Exception("Not valid subsetter: {}".format(str(val)))

        # print("getitem " * 100)
        # print(df)
        # print(ray.get(df))
        return PyRanges(df)


    def __str__(self):

        # print("in str")
        if len(self) == 0:
            return "Empty PyRanges"

        # keys = natsorted(list(self.dfs.keys()))
        if len(self.keys()) == 1:

            first_key = self.keys()[0]
            # last_key = list(self.dfs.keys())[-1]
            first_df = self.dfs[first_key]

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
            keys = self.keys()
            first_key = keys[0]
            last_key = keys[-1]
            # first_key = self.keys[0]
            # last_key = self.keys[-1]
            first_df = self.dfs[first_key].head(3)
            last_df = self.dfs[last_key].tail(3)
            # last_df = self.dfs[list(self.dfs.keys())[-1]].tail(3)

            h = first_df.head(3).astype(object)
            m = first_df.head(1).astype(object)
            t = last_df.head(3).astype(object)
            m.loc[:,:] = "..."
            # m.index = ["..."]
            # print((len(h) + len(t)) < 6, len(self) >= 6)
            if (len(h) + len(t)) < 6:

                keys_covered = set()
                # iterate from front until have three
                heads = []
                hl = 0
                for k in keys:
                    keys_covered.add(k)
                    h = self.dfs[k].head(3)
                    first_df = h
                    hl += len(h)
                    heads.append(h)
                    if hl >= 3:
                        break

                tails = []
                tl = 0
                for k in keys[::-1]:
                    if k in keys_covered:
                        continue

                    t = self.dfs[k].tail(3)
                    tl += len(t)
                    tails.append(t)
                    if tl >= 3:
                        break
                # iterate from back until have three

                h = pd.concat(heads).head(3).astype(object)
                if tails:
                    t = pd.concat(tails).tail(3).astype(object)
                    if len(h) + len(t) > 6:
                        m = h.head(1).astype(object)
                        m.loc[:,:] = "..."
                        s = pd.concat([h, m, t])
                    else:
                        s = pd.concat([h, t])
                else:
                    s = h

            elif len(h) + len(t) == 6:
                m.loc[:,:] = "..."
                s = pd.concat([h, m, t])
            else:
                s = pd.concat([h, t])

        if False: # make setting
            if self.stranded:
                pos = s.Chromosome.astype(str) + " " + s.Start.astype(str) + "-" + s.End.astype(str) + " " + s.Strand.astype(str)
                s = s.drop("Chromosome Start End Strand".split(), axis=1)
                first_df = first_df.drop("Chromosome Start End Strand".split(), axis=1)
            else:
                pos = s.Chromosome.astype(str) + " " + s.Start.astype(str) + "-" + s.End.astype(str)
                s = s.drop("Chromosome Start End".split(), axis=1)
                first_df = first_df.drop("Chromosome Start End".split(), axis=1)

            s.insert(0, "Position", pos)
            h = [c + "\n(" + str(t) + ")" for c, t in  zip(s.columns, ["multiple types"] + list(first_df.dtypes))]
        else:
            dtypes = []
            for col, dtype in zip(s.columns, first_df.dtypes):
                if str(dtype) == "category":
                    dtype = first_df[col].cat.codes.dtype

                # dtype = str(dtype).replace("float", "f_").replace("int", "i_")
                dtypes.append(dtype)

            h = [c + "\n(" + str(t) + ")" for c, t in  zip(s.columns, list(dtypes))]

        str_repr = tabulate(s, headers=h, tablefmt='psql', showindex=False) + \
                                        "\nPyRanges object has {} sequences from {} chromosomes.".format(len(self), len(self.chromosomes))
        return str_repr


    def __repr__(self):

        return str(self)


    # @pyrange_or_df
    # @return_empty_if_one_empty
    def overlap(self, other, **kwargs):

        "Want all intervals in self that overlap with other."

        # print(kwargs)
        kwargs["sparse"] = {"self": False, "other": True}
        kwargs = fill_kwargs(kwargs)

        dfs = pyrange_apply(_overlap, self, other, **kwargs)

        return PyRanges(dfs)


    def no_overlap(self, other, **kwargs):

        "Want all intervals in self that do not overlap with other."

        # print(kwargs)
        kwargs = fill_kwargs(kwargs)
        kwargs["invert"] = True
        kwargs["sparse"] = {"self": False, "other": True}
        # print(kwargs)

        dfs = pyrange_apply(_overlap, self, other, **kwargs)

        return PyRanges(dfs)

    # @pyrange_or_df
    # @return_empty_if_one_empty
    def nearest(self, other, **kwargs):

        "Find the nearest feature in other."

        kwargs = fill_kwargs(kwargs)

        dfs = pyrange_apply(_nearest, self, other, **kwargs)

        return PyRanges(dfs)

    # @return_empty_if_one_empty
    def intersect(self, other, **kwargs):

        kwargs = fill_kwargs(kwargs)
        kwargs["sparse"] = {"self": False, "other": True}

        dfs = pyrange_apply(_intersection, self, other, **kwargs)

        return PyRanges(dfs)

    # @return_empty_if_one_empty
    def set_intersect(self, other, **kwargs):

        kwargs = fill_kwargs(kwargs)
        strandedness = kwargs["strandedness"]
        strand = True if strandedness else False
        self_clusters = self.cluster(strand=strand, **kwargs)
        other_clusters = other.cluster(strand=strand, **kwargs)
        dfs = pyrange_apply(_intersection, self_clusters, other_clusters, **kwargs)

        return PyRanges(dfs)

    def set_union(self, other, **kwargs):

        kwargs = fill_kwargs(kwargs)
        strandedness = kwargs["strandedness"]
        strand = True if strandedness else False

        pr = self.concat(other)
        # print(dfs)
        # print(pr)
        pr = pr.cluster(strand=strand, **kwargs)

        return pr


    # @pyrange_or_df
    def subtract(self, other, **kwargs):

        kwargs["sparse"] = {"self": False, "other": True}
        kwargs = fill_kwargs(kwargs)
        strandedness = kwargs["strandedness"]

        strand = True if strandedness else False
        other_clusters = other.cluster(strand=strand, **kwargs)
        result = pyrange_apply(_subtraction, self, other_clusters, **kwargs)

        return PyRanges(result)


    # @pyrange_or_df
    # @return_empty_if_one_empty
    def join(self, other, **kwargs):

        kwargs = fill_kwargs(kwargs)
        dfs = pyrange_apply(_write_both, self, other, **kwargs)

        return PyRanges(dfs)


    def cluster(self, strand=None, **kwargs):

        kwargs["sparse"] = True
        df = pyrange_apply_single(_cluster, self, strand, kwargs)

        return PyRanges(df)


    def coverage(self, value_col=None, strand=True):

        return _coverage(self, value_col, strand=strand)


    def apply(self, f, strand=True, as_pyranges=True, kwargs=None):

        if kwargs is None:
            kwargs = {}

        f = ray.remote(f)

        result = pyrange_apply_single(f, self, strand, kwargs)

        if not as_pyranges:
            return result
        else:
            return PyRanges(result)


    def apply_pair(self, other, f, kwargs, strand=True, as_pyranges=True):

        f = ray.remote(f)

        result = pyrange_apply(f, self, other, strand, kwargs)

        if not as_pyranges:
            return result
        else:
            return PyRanges(result)


    def eval(self, eval_cmd, strand=True, as_pyranges=True, kwargs=None):

        f = lambda df: eval(eval_cmd)

        if kwargs is None:
            kwargs = {}

        f = ray.remote(f)

        result = pyrange_apply_single(f, self, strand, kwargs)

        if not as_pyranges:
            return result
        else:
            return PyRanges(result)


    def concat(self, other):

        return PyRanges(_concat(self, other))


    def jaccard(self, other, **kwargs):

        kwargs = fill_kwargs(kwargs)
        strand = True if kwargs["strandedness"] else False

        intersection_sum = sum(v.sum() for v in self.set_intersect(other, **kwargs).lengths().values())

        union_sum = 0
        for gr in [self, other]:
            union_sum += sum(v.sum() for v in gr.cluster(strand=strand).lengths().values())

        denominator = (union_sum - intersection_sum)
        if denominator == 0:
            return 1
        else:
            jc = intersection_sum / denominator

        return jc



    def relative_distance(self, other, **kwargs):

        kwargs = fill_kwargs(kwargs)
        result = pyrange_apply(_relative_distance, self, other, **kwargs)

        result = pd.Series(np.concatenate(list(result.values())))

        not_nan = ~np.isnan(result)
        result.loc[not_nan] = np.floor(result[not_nan] * 100) / 100
        vc = result.value_counts(dropna=False).to_frame().reset_index()
        vc.columns = "reldist count".split()
        vc.insert(vc.shape[1], "total", len(result))
        vc.insert(vc.shape[1], "fraction", vc["count"] / len(result))
        vc = vc.sort_values("reldist", ascending=True)
        vc = vc.reset_index(drop=True)
        return vc


    def slack(self, slack):

        kwargs = {"slack": slack}
        prg = PyRanges(pyrange_apply_single(_slack, self, self.stranded, kwargs))

        return prg

    def tssify(self, slack=0):

        kwargs = {"slack": slack}
        return PyRanges(pyrange_apply_single(_tss, self, self.stranded, kwargs))

    def sort(self, columns=["Start", "End"], **kwargs):

        kwargs = fill_kwargs(kwargs)
        return PyRanges(pyrange_apply_single(_sort, self, self.stranded, kwargs))

    def tesify(self, slack=0):

        kwargs = {"slack": slack}
        return PyRanges(pyrange_apply_single(_tes, self, self.stranded, kwargs))

    def index_as_col(self, **kwargs):
        kwargs = fill_kwargs(kwargs)
        return PyRanges(pyrange_apply_single(_index_as_col, self, self.stranded, kwargs))

    def drop_empty(self):

        empty = []
        for k in self.dfs.keys():
            if self.dfs[k].empty:
                empty.append(k)

        for k in empty:
            del self.dfs[k]

    def empty(self):

        if len(self.values()) == 0:
            return True

        return False

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


    # @property
    def keys(self):
        return natsorted(self.dfs.keys())

    @property
    def stranded(self):
        # print(self.keys())
        return len(list(self.keys())[0]) == 2

    @property
    def strands(self):

        if not self.stranded:
            raise Exception("PyRanges not stranded!")

        return natsorted(set([k[1] for k in self.keys()]))

    @property
    def chromosomes(self):

        if self.stranded:
            return natsorted(set([k[0] for k in self.keys()]))
        else:
            return natsorted(set([k for k in self.keys()]))


    # @property
    def items(self):

        return natsorted([(k, df) for (k, df) in self.dfs.items()])

    # @property
    def values(self):

        return [df for k, df in self.items() if not df.empty]

    # @property
    # def objids(self):

    #     return [objid for k, objid in natsorted(self.dfs.items())]

    @property
    def df(self):

        return self.as_df()

    def as_df(self):

        if len(self) == 0:
            return pd.DataFrame()
        elif len(self) == 1:
            return self.values()[0]
        else:
            return pd.concat(self.values())


    def lengths(self):

        lengths = {}
        for k, df in self.items():
            lengths[k] = df.End - df.Start

        return lengths

    
    def midpoints(self):

        midpoints = {}
        for k, df in self.items():
            midpoints[k] = (df.End + df.Start) / 2

        return midpoints


    def summary(self):

        lengths = OrderedDict()
        lengths["pyrange"] = self.lengths()

        if self.stranded:
            c = self.cluster(strand=True)
            lengths["coverage_stranded"] = c.lengths()

        c = self.cluster(strand=False)
        lengths["coverage_unstranded"] = c.lengths()

        summaries = OrderedDict()

        for summary, d in lengths.items():
            summaries[summary] = pd.concat(d.values()).describe()


        summary = pd.concat(summaries.values(), axis=1)
        summary.columns = list(summaries)

        str_repr = tabulate(summary, headers=summary.columns, tablefmt='psql')
        print(str_repr)
