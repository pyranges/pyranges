import pandas as pd
import numpy as np
from collections import defaultdict

from tabulate import tabulate


from pyranges.genomicfeatures import GenomicFeaturesMethods
# from pyranges.genomicfeatures import _tss, _tes
from pyranges.subset import get_string, get_slice, get_tuple
from pyranges.methods import _cluster, _subtraction, _set_union, _set_intersection, _intersection, _nearest, _coverage, _overlap_write_both, _overlap, _tss, _tes, _jaccard, _lengths, _slack
# from pyranges.multithreaded import _cluster, pyrange_apply_single, _subtraction, _set_union, _set_intersection, _intersection, pyrange_apply, _nearest, _coverage, _write_both, _first_df

from ncls import NCLS

def return_copy_if_view(df, df2):
    # https://stackoverflow.com/questions/26879073/checking-whether-data-frame-is-copy-or-view-in-pandas

    if df.values.base is df2.values.base:
        return df.copy(deep=True)
    else:
        return df

# def is_view(df, df2):
#     # https://stackoverflow.com/questions/26879073/checking-whether-data-frame-is-copy-or-view-in-pandas

#     return df.values.base is df2.values.base

def create_ncls(cdf):

    return NCLS(cdf.Start.values,
                cdf.End.values,
                cdf.index.values)


def create_ncls_dict(df, strand):

    if not strand:
        grpby_key = "Chromosome"
    else:
        grpby_key = "Chromosome Strand".split()

    grpby = df.groupby(grpby_key)
    nclses = {key: create_ncls(cdf) for (key, cdf) in grpby}
    dd = defaultdict(NCLS)
    dd.update(nclses)

    return dd


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

    df = None
    gf = None

    def __init__(self, df=None, seqnames=None, starts=None, ends=None, strands=None, copy_df=True):

        df_given = True if not df is None else False

        if df is False or df is None:
            df = create_pyranges_df(seqnames, starts, ends, strands)
        else:
            assert "Chromosome" in df and "Start" in df and "End" in df
            df.index = range(len(df))

        # using __dict__ to avoid invoking setattr
        if "Strand" not in df:
            self.__dict__["stranded"] = False
        else:
            self.__dict__["stranded"] = True

        if copy_df and df_given:
            df = df.copy()
            df.Chromosome = df.Chromosome.astype("category")
            if "Strand" in df:
                df.Strand = df.Strand.astype("category")


        self.__dict__["df"] = df

        self.__dict__["__ncls__"] = create_ncls_dict(df, self.stranded)

        self.__dict__["ft"] = GenomicFeaturesMethods(self)


    def __len__(self):
        return self.df.shape[0]

    def __setattr__(self, column_name, column):

        if column_name in "Chromosome Start End Strand".split():
            raise Exception("The columns Chromosome, Start, End or Strand can not be reset.")
        if column_name == "stranded":
            raise Exception("The stranded attribute is read-only. Create a new PyRanges object instead.")

        if column_name == "df":
            self.__dict__["df"] = column
            return

        if not isinstance(column, str):
            if not len(self.df) == len(column):
                raise Exception("DataFrame and column must be same length.")

            column_to_insert = pd.Series(column, index=self.df.index)
        else:
            column_to_insert = pd.Series(column, index=self.df.index)

        pos = self.df.shape[1]
        if column_name in self.df:
            pos = list(self.df.columns).index(column_name)
            self.df.drop(column_name, inplace=True, axis=1)

        self.df.insert(pos, column_name, column_to_insert)


    def __getattr__(self, name):

        if name in self.df:
            return self.df[name]
        else:
            self.__dict__[name]

    def __eq__(self, other):

        return self.df.equals(other.df)

    def __getitem__(self, val):

        if isinstance(val, str):
            df = get_string(self, val)

        elif isinstance(val, tuple):
            df = get_tuple(self, val)
        elif isinstance(val, slice):
            df = get_slice(self, val)

        if not df._is_view:
            return PyRanges(df)
        else:
            return PyRanges(df.copy(deep=True))


    def __str__(self):

        if len(self.df) > 6:
            h = self.df.head(3).astype(object)
            t = self.df.tail(3).astype(object)
            m = self.df.head(1).astype(object)
            m.loc[:,:] = "..."
            m.index = ["..."]
            s = pd.concat([h, m, t])
        else:
            s = self.df

        h = [c + "\n(" + str(t) + ")" for c, t in  zip(self.df.columns, self.df.dtypes)]

        str_repr = tabulate(s, headers=h, tablefmt='psql', showindex=False) + \
                                        "\nPyRanges object has {} sequences from {} chromosomes.".format(self.df.shape[0], len(set(self.df.Chromosome)))
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

    @pyrange_or_df
    @return_empty_if_one_empty
    def intersection(self, other, strandedness=False, how=None):

        "Want the parts of the intervals in self that overlap with other."

        df = _intersection(self, other, strandedness, how)

        return df

    @pyrange_or_df
    @return_empty_if_one_empty
    def set_intersection(self, other, strandedness=False, how=None):

        si = _set_intersection(self, other, strandedness, how)

        return si

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


    @pyrange_or_df_single
    def cluster(self, strand=None, max_dist=0, min_nb=1, **kwargs):

        df = _cluster(self, strand, max_dist, min_nb)

        return df

    # @pyrange_or_df_single
    # def tile(self, tile_size=50):

    #     df = _tile(self, tile_size)
    #     return df


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



    # @pyrange_or_df
    # @return_empty_if_one_empty
    # def overlap(self, other, strandedness=False, invert=False, how=None, nb_cpu=1, **kwargs):

    #     "Want all intervals in self that overlap with other."

    #     df = pyrange_apply(_first_df, self, other, strandedness=strandedness,
    #                        invert=invert, how=how, n_jobs=nb_cpu, **kwargs)

    #     return df

    # @pyrange_or_df
    # @return_empty_if_one_empty
    # def nearest(self, other, strandedness=False, suffix="_b", how=None, overlap=True, nb_cpu=1, **kwargs):
    #     "Find the nearest feature in other."

    #     df = pyrange_apply(_nearest, self, other, strandedness=strandedness, suffix=suffix, how=how, overlap=overlap, n_jobs=nb_cpu)

    #     return df

    # @pyrange_or_df
    # @return_empty_if_one_empty
    # def intersection(self, other, strandedness=False, how=None, nb_cpu=1):

    #     "Want the parts of the intervals in self that overlap with other."

    #     df = pyrange_apply(_intersection, self, other, strandedness=strandedness, how=how, n_jobs=nb_cpu)
    #     # df = _intersection(self, other, strandedness=strandedness, how=how)

    #     return df

    # @pyrange_or_df
    # @return_empty_if_one_empty
    # def set_intersection(self, other, strandedness=False, how=None, nb_cpu=1):

    #     si = pyrange_apply(_set_intersection, self, other, strandedness=strandedness, how=how, n_jobs=nb_cpu)

    #     return si

    # @pyrange_or_df
    # @return_empty_if_both_empty
    # def set_union(self, other, strandedness=False, nb_cpu=1):

    #     si = pyrange_apply(_set_union, self, other, strandedness=strandedness, n_jobs=nb_cpu)

    #     return si


    # @pyrange_or_df
    # def subtraction(self, other, strandedness=False, nb_cpu=1):


    #     return pyrange_apply(_subtraction, self, other, strandedness=strandedness, n_jobs=nb_cpu)


    # @pyrange_or_df
    # @return_empty_if_one_empty
    # def join(self, other, strandedness=False, new_pos=None, suffixes=["_a", "_b"], how=None, nb_cpu=1, **kwargs):

    #     df = pyrange_apply(_write_both, self, other, strandedness=strandedness, new_pos=new_pos,
    #                        suffixes=suffixes, how=how, n_jobs=nb_cpu, **kwargs)

    #     return df


    # @pyrange_or_df_single
    # def cluster(self, strand=None, nb_cpu=1):

    #     df = pyrange_apply_single(_cluster, self, strand=strand, n_jobs=nb_cpu)

    #     return df


    # @pyrange_or_df_single
    # def tile(self, tile_size=50, nb_cpu=1):

    #     df = _tile(self, tile_size)
    #     return df


    # def coverage(self, value_col=None, stranded=False, nb_cpu=1):

    #     return _coverage(self, value_col, stranded=stranded, n_jobs=nb_cpu)
