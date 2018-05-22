import pandas as pd
import numpy as np
from collections import defaultdict

from tabulate import tabulate

from pyranges.methods import (_overlap, _cluster, _tile, _intersection,
                              _coverage, _overlap_write_both,
                              _set_intersection, _set_union, _subtraction,
                              _nearest)

from ncls import NCLS

# try:
#     dummy = profile
# except:
#     profile = lambda x: x


def create_ncls(cdf):

    return NCLS(cdf.Start.values,
                cdf.End.values,
                cdf.index.values)


def create_ncls_dict(df, strand):

    if not strand:
        grpby = df.groupby("Chromosome")
    else:
        grpby = df.groupby("Chromosome Strand".split())

    nclses = {key: create_ncls(cdf) for (key, cdf) in grpby}
    dd = defaultdict(NCLS)
    dd.update(nclses)

    return dd


def create_pyranges_df(seqnames, starts, ends, strands=None):

    if isinstance(seqnames, str):
        seqnames = [seqnames] * len(starts)

    if isinstance(strands, str):
        strands = [strands] * len(starts)

    if strands != None and strands != False:
        columns = [seqnames, starts, ends, strands]
        lengths = list(str(len(s)) for s in columns)
        assert len(set(lengths)) == 1, "seqnames, starts, ends and strands must be of equal length. But are {}".format(", ".join(lengths))
    else:
        columns = [seqnames, starts, ends]
        lengths = list(str(len(s)) for s in columns)
        assert len(set(lengths)) == 1, "seqnames, starts and ends must be of equal length. But are {}".format(", ".join(lengths))

    idx = range(len(starts))
    series_to_concat = []
    for s in columns:
        s = pd.Series(s, index=idx)
        series_to_concat.append(s)

    df = pd.concat(series_to_concat, axis=1)
    df.columns = "Chromosome Start End Strand".split()

    return df



class PyRanges():

    df = None

    def __init__(self, df=None, seqnames=None, starts=None, ends=None, strands=None):

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

        self.__dict__["df"] = df

        self.__dict__["__ncls__"] = create_ncls_dict(df, self.stranded)


    def __len__(self):
        return self.df.shape[0]

    def __setattr__(self, column_name, column):

        if column_name in "Chromosome Start End Strand".split():
            raise Exception("The columns Chromosome, Start, End or Strand can not be reset.")
        if column_name == "stranded":
            raise Exception("The stranded attribute is read-only. Create a new PyRanges object instead.")

        if not isinstance(column, str):
            if not len(self.df) == len(column):
                raise Exception("DataFrame and column must be same length.")

            column_to_insert = pd.Series(column, index=self.df.index)
        else:
            column_to_insert = pd.Series(column, index=self.df.index)

        self.df.insert(self.df.shape[1], column_name, column_to_insert)


    def __getattr__(self, name):

        if name in self.df:
            return self.df[name]
        else:
            self.__dict__[name]


    def __getitem__(self, val):

        pd.options.mode.chained_assignment = None
        if isinstance(val, str):
            if val in set(self.df.Chromosome):
                return PyRanges(self.df.loc[self.df.Chromosome == val])
            elif val in "+ -".split():
                return PyRanges(self.df.loc[self.df.Strand == val])
            else:
                raise Exception("Invalid choice for string subsetting PyRanges: {}. Must be either strand or chromosome".format(val))

        elif isinstance(val, tuple):

            # "chr1", 5:10
            if len(val) == 2 and val[0] in self.df.Chromosome.values and isinstance(val[1], slice):
                chromosome, loc = val
                start = loc.start or 0
                stop = loc.stop or max(self.df.loc[self.df.Chromosome == chromosome].End.max(), start)
                idxes = [r[2] for r in self.__ncls__[chromosome, "+"].find_overlap(start, stop)] + \
                        [r[2] for r in self.__ncls__[chromosome, "-"].find_overlap(start, stop)]

                return PyRanges(self.df.loc[idxes])

            # "+", 5:10
            if len(val) == 2 and val[0] in "+ -".split() and isinstance(val[1], slice):
                strand, loc = val
                start = loc.start or 0
                stop = loc.stop or max(self.df.loc[self.df.Chromosome == chromosome].End.max(), start)
                idxes = []
                for chromosome in self.df.Chromosome.drop_duplicates():
                    idxes.extend([r[2] for r in self.__ncls__[chromosome, strand].find_overlap(start, stop)])

                return PyRanges(self.df.loc[idxes])

            # "chr1", "+"
            if len(val) == 2 and val[1] in "+ -".split():

                chromosome, strand = val

                return PyRanges(self.df.loc[(self.df.Chromosome == chromosome) & (self.df.Strand == strand)])

            # "chr1", "+", 5:10
            elif len(val) == 3 and val[0] in self.df.Chromosome.values and val[1] in "+ -".split():

                chromosome, strand, loc = val
                start = loc.start or 0
                stop = loc.stop or max(self.df.loc[self.df.Chromosome == chromosome].End.max(), start)
                idxes = [r[2] for r in self.__ncls__[chromosome, strand].find_overlap(start, stop)]

                return PyRanges(self.df.loc[idxes])

        elif isinstance(val, slice):

            start = val.start or 0
            stop = val.stop or max(self.df.End.max(), start)

            idxes = []
            for it in self.__ncls__.values():
                idxes.extend([r[2] for r in it.find_overlap(start, stop)])

            return PyRanges(self.df.loc[idxes])

        pd.options.mode.chained_assignment = "warn"


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

        str_repr = tabulate(s, headers='keys', tablefmt='psql', showindex=False) + \
                                        "\nPyRanges object has {} sequences from {} chromosomes.".format(self.df.shape[0], len(set(self.df.Chromosome)))
        return str_repr


    def __repr__(self):

        return str(self)


    def overlap(self, other, strandedness=False, invert=False):

        "Want all intervals in self that overlap with other."

        df = _overlap(self, other, strandedness, invert)
        return PyRanges(df)

    def nearest(self, other, strandedness=False):

        "Want all intervals in self that overlap with other."

        df = _nearest(self, other, strandedness)

        return PyRanges(df)


    def intersection(self, other, strandedness=False):

        "Want the parts of the intervals in self that overlap with other."

        df = _intersection(self, other, strandedness)

        df.index.name = ""
        return PyRanges(df)


    def set_intersection(self, other, strandedness=False):

        si = _set_intersection(self, other, strandedness)

        return PyRanges(si)


    def set_union(self, other, strand=False):

        si = _set_union(self, other, strand)

        return PyRanges(si)


    def subtraction(self, other, strandedness=False):

        return PyRanges(_subtraction(self, other, strandedness))


    def join(self, other, strandedness=False, new_pos=None, suffixes=["_a", "_b"]):

        df = _overlap_write_both(self, other, strandedness, new_pos, suffixes)

        return PyRanges(df)


    def cluster(self, strand=None, df_only=False, max_dist=0, min_nb=1):

        df = _cluster(self, strand, max_dist, min_nb)

        if not df_only:
            return PyRanges(df)
        else:
            return df

    # def subtraction(self, other, df_only=False, max_dist=0, min_nb=1):

    #     _subtra


    def tile(self, tile_size=50):

        df = _tile(self, tile_size)
        return PyRanges(df)


    def coverage(self, value_col=None):

        return _coverage(self, value_col)
