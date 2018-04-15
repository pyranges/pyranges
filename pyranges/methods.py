import pandas as pd

import itertools


def _intersection(self, other, strandedness="same", invert=False):

    indexes_of_overlapping = []
    df = self.df
    if "Strand" in df and "Strand" in other.df and strandedness == "same":

        for (chromosome, strand), cdf in df.groupby("Chromosome Strand".split()):
            it = other.__intervaltrees__[chromosome, strand]

            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(), cdf.End.tolist()):

                if it.search(start, end):
                    indexes_of_overlapping.append(idx)

    other_strand = {"+": "-", "-": "+"}

    if "Strand" in df and "Strand" in other.df and strandedness == "opposite":

        for (chromosome, strand), cdf in df.groupby("Chromosome Strand".split()):
            opposite_strand = other_strand[strand]
            it = other.__intervaltrees__[chromosome, opposite_strand]

            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(), cdf.End.tolist()):

                if it.search(start, end):
                    indexes_of_overlapping.append(idx)

    else:

        for chromosome, cdf in df.groupby("Chromosome"):

            it = other.__intervaltrees__[chromosome, "+"]
            it2 = other.__intervaltrees__[chromosome, "-"]

            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(), cdf.End.tolist()):

                if it.search(start, end) or it2.search(start, end):
                    indexes_of_overlapping.append(idx)

    if not invert:
        return df.loc[indexes_of_overlapping]
    else:
        indexes_of_nonoverlapping = df[~df.index.isin(indexes_of_overlapping)].index
        return df.loc[indexes_of_nonoverlapping]


def _cluster(self, strand=False):

    df = self.df.copy()

    if strand:
        grp = ["Chromosome", "Strand"]
    else:
        grp = "Chromosome"

    clusters = []
    i = 0
    for _, cdf in df.groupby(grp):
        for _, cluster_df in cdf.groupby((cdf.End.shift() - cdf.Start).lt(0).cumsum()):
            i += 1
            clusters.append(["{}".format(i)] * cluster_df.shape[0])

    clusters = pd.Series(list(itertools.chain.from_iterable(clusters)))
    df.insert(df.shape[1], "ClusterID", clusters)
    return df




def _tile(self, tile_size=50):

    df = self.df.copy()

    df.Start = df.Start - (df.Start % tile_size)
    df.End = df.End - (df.End % tile_size) + tile_size - 1

    rows = []
    for _, row in df.iterrows():
        d = row.to_dict()
        for tile in range(row.Start, row.End, tile_size):
            d2 = d.copy()
            d2["Start"], d2["End"] = tile, tile + tile_size - 1
            rows.append(d2)


    df = pd.DataFrame.from_dict(rows)[df.columns]

    return df


def invert(self):

    pass

    # elif agg_func == "sum":

    #     for (chromosome, strand), cdf in df.groupby("Chromosome Strand".split()):
    #         cluster_df = cdf.groupby((cdf.End.shift() - cdf.Start).lt(0).cumsum()).sum()
    #         return cluster_dfs

    # else:

    #     for (chromosome, strand), cdf in df.groupby("Chromosome Strand".split()):
    #         cluster_df = cdf.groupby((cdf.End.shift() - cdf.Start).lt(0).cumsum()).apply(agg_func)
    #         return cluster_df
