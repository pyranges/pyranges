import pandas as pd
import numpy as np
from collections import defaultdict

from tabulate import tabulate

from natsort import natsorted

from pyranges.methods import (_overlap, _cluster, _tile, _inverse_intersection,
                              _intersection, _coverage, _overlap_write_both,
                              _set_intersection, _set_union)

from ncls import NCLS

try:
    dummy = profile
except:
    profile = lambda x: x


def create_ncls(cdf):

    return NCLS(cdf.Start.values,
                cdf.End.values,
                cdf.index.values)


def create_ncls_dict(df, n):

    if "Strand" not in df:
        grpby = df.groupby("Chromosome")
    else:
        grpby = df.groupby("Chromosome Strand".split())

    nclses = {key: create_ncls(cdf) for (key, cdf) in grpby}
    dd = defaultdict(NCLS)
    dd.update(nclses)

    return dd




class GRanges():

    df = None

    def __init__(self, df, n=1):

        df.index = range(len(df))

        # using __dict__ to avoid invoking setattr

        self.__dict__["df"] = df

        self.__dict__["__ncls__"] = create_ncls_dict(df, n)

    def __len__(self):
        return self.df.shape[0]

    def __setattr__(self, column_name, column):

        if not len(self.df) == len(column):
            raise Exception("DataFrame and column must be same length.")

        column.index = self.df.index
        self.df.insert(self.df.shape[1], column_name, column)

    def __getattr__(self, column_name):

        if column_name in self.df:
            return self.df[column_name]

        raise IndexError("Column name {} not in df.".format(column_name))


    def __getitem__(self, val):

        pd.options.mode.chained_assignment = None
        if isinstance(val, str):
            if val in set(self.df.Chromosome):
                return GRanges(self.df.loc[self.df.Chromosome == val])
            elif val in "+ -".split():
                return GRanges(self.df.loc[self.df.Strand == val])
            else:
                raise Exception("Invalid choice for string subsetting GRanges: {}. Must be either strand or chromosome".format(val))

        elif isinstance(val, tuple):

            # "chr1", 5:10
            if len(val) == 2 and val[0] in self.df.Chromosome.values and isinstance(val[1], slice):
                chromosome, loc = val
                start = loc.start or 0
                stop = loc.stop or max(self.df.loc[self.df.Chromosome == chromosome].End.max(), start)
                idxes = [r[2] for r in self.__ncls__[chromosome, "+"].find_overlap(start, stop)] + \
                        [r[2] for r in self.__ncls__[chromosome, "-"].find_overlap(start, stop)]

                return GRanges(self.df.loc[idxes])

            # "+", 5:10
            if len(val) == 2 and val[0] in "+ -".split() and isinstance(val[1], slice):
                strand, loc = val
                start = loc.start or 0
                stop = loc.stop or max(self.df.loc[self.df.Chromosome == chromosome].End.max(), start)
                idxes = []
                for chromosome in self.df.Chromosome.drop_duplicates():
                    idxes.extend([r[2] for r in self.__ncls__[chromosome, strand].find_overlap(start, stop)])

                return GRanges(self.df.loc[idxes])

            # "chr1", "+"
            if len(val) == 2 and val[1] in "+ -".split():

                chromosome, strand = val

                return GRanges(self.df.loc[(self.df.Chromosome == chromosome) & (self.df.Strand == strand)])

            # "chr1", "+", 5:10
            elif len(val) == 3 and val[0] in self.df.Chromosome.values and val[1] in "+ -".split():

                chromosome, strand, loc = val
                start = loc.start or 0
                stop = loc.stop or max(self.df.loc[self.df.Chromosome == chromosome].End.max(), start)
                idxes = [r[2] for r in self.__ncls__[chromosome, strand].find_overlap(start, stop)]

                return GRanges(self.df.loc[idxes])

        elif isinstance(val, slice):

            start = val.start or 0
            stop = val.stop or max(self.df.End.max(), start)

            idxes = []
            for it in self.__ncls__.values():
                idxes.extend([r[2] for r in it.find_overlap(start, stop)])

            return GRanges(self.df.loc[idxes])

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
                                        "\nGRanges object has {} sequences from {} chromosomes.".format(self.df.shape[0], len(set(self.df.Chromosome)))
        return str_repr


    def __repr__(self):

        return str(self)


    def overlap(self, other, strandedness=False, invert=False):

        "Want all intervals in self that overlap with other."

        df = _overlap(self, other, strandedness, invert)
        return GRanges(df)


    def intersection(self, other, strandedness=False, invert=False):

        "Want the parts of the intervals in self that overlap with other."

        if invert:
            df = _inverse_intersection(self, other, strandedness)
        else:
            df = _intersection(self, other, strandedness)

        df.index.name = ""
        return GRanges(df)


    def set_intersection(self, other, strandedness=False):

        si = _set_intersection(self, other, strandedness)

        return GRanges(si)


    def set_union(self, other, strand=False):

        si = _set_union(self, other, strand)

        return GRanges(si)


    def overlap_join(self, other, strandedness=False, new_pos=None, suffixes=["_a", "_b"]):

        df = _overlap_write_both(self, other, strandedness, new_pos, suffixes)

        return GRanges(df)


    def cluster(self, strand=None, df_only=False, max_dist=0, min_nb=1):

        df = _cluster(self, strand, max_dist, min_nb)

        if not df_only:
            return GRanges(df)
        else:
            return df



    def tile(self, tile_size=50):

        df = _tile(self, tile_size)
        return GRanges(df)


    def coverage(self, value_col=None):

        return _coverage(self, value_col)
