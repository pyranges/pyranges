
import pandas as pd
from collections import defaultdict

from quicksect import IntervalTree

from tabulate import tabulate

class GRanges():

    def __init__(self, df):

        assert "Start" in df and "End" in df and "Chromosome" in df, \
        """DataFrame missing at least one of Start, End or Chromosome columns.
These are the columns it contains: {}""".format(" ".join(df.columns))

        df = df.reset_index(drop=True)

        self.df = df

        self.__intervaltrees__ = defaultdict(IntervalTree)

        for chromosome, cdf in df.groupby("Chromosome"):

            it = IntervalTree()
            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(),
                                       cdf.End.tolist()):
                it.add(start, end, idx)

            self.__intervaltrees__[chromosome] = it


    def __sub__(self, other):

        "Want all intervals in self that do not overlap with other."

        indexes_of_nonoverlapping = []
        for chromosome, cdf in self.df.groupby("Chromosome"):

            it = other.__intervaltrees__[chromosome]

            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(), cdf.End.tolist()):

                if not it.search(start, end):
                    indexes_of_nonoverlapping.append(idx)

        return GRanges(self.df.loc[indexes_of_nonoverlapping])


    def __xor__(self, other):

        "Want all intervals in self that do not overlap with other."

        indexes_of_nonoverlapping = []
        for chromosome, cdf in self.df.groupby("Chromosome"):

            it = other.__intervaltrees__[chromosome]

            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(), cdf.End.tolist()):

                if not it.search(start, end):
                    indexes_of_nonoverlapping.append(idx)

        return GRanges(self.df.loc[indexes_of_nonoverlapping])

    def __getitem__(self, val):

        if isinstance(val, str):
            return GRanges(self.df.loc[self.df.Chromosome == val])

        elif isinstance(val, tuple):
            chromosome, loc = val
            start = loc.start or 0
            stop = loc.stop or max(self.df.loc[self.df.Chromosome == chromosome].End.max(), start)
            idxes = [r.data for r in self.__intervaltrees__[chromosome].search(start, stop)]

            return GRanges(self.df.loc[idxes])

        elif isinstance(val, slice):

            start = val.start or 0
            stop = val.stop or max(self.df.End.max(), start)

            idxes = []
            for it in self.__intervaltrees__.values():
                idxes.extend([r.data for r in it.search(start, stop)])

            return GRanges(self.df.loc[idxes])


    def __getattr__(self, col):

        try:
            return self.df[col]
        except:
            raise Exception("Column {} not found.".format(col))


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

        str_repr = tabulate(s, headers='keys', tablefmt='psql') + "\nGRanges object with {} sequences from {} chromosomes.".format(self.df.shape[0], len(set(self.df.Chromosome)))
        return str_repr
