import pandas as pd

import itertools


def _overlap(self, other, strandedness=False, invert=False):

    indexes_of_overlapping = []
    df = self.df
    other_strand = {"+": "-", "-": "+"}

    if "Strand" in df and "Strand" in other.df and strandedness == "same":

        for (chromosome, strand), cdf in df.groupby("Chromosome Strand".split()):
            it = other.__intervaltrees__[chromosome, strand]

            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(),
                                       (cdf.End - 1).tolist()):

                if it.search(start, end):
                    indexes_of_overlapping.append(idx)


    elif "Strand" in df and "Strand" in other.df and strandedness == "opposite":

        for (chromosome, strand), cdf in df.groupby("Chromosome Strand".split()):
            opposite_strand = other_strand[strand]
            it = other.__intervaltrees__[chromosome, opposite_strand]

            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(),
                                       (cdf.End - 1).tolist()):

                if it.search(start, end):
                    indexes_of_overlapping.append(idx)

    else:

        for chromosome, cdf in df.groupby("Chromosome"):

            it = other.__intervaltrees__[chromosome, "+"]
            it2 = other.__intervaltrees__[chromosome, "-"]

            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(),
                                       (cdf.End - 1).tolist()):

                if it.search(start, end) or it2.search(start, end):
                    indexes_of_overlapping.append(idx)

    # https://stackoverflow.com/a/7961425/992687
    indexes_of_overlapping = list(dict.fromkeys(indexes_of_overlapping))
    if not invert:
        return df.loc[indexes_of_overlapping]
    else:
        indexes_of_nonoverlapping = df[~df.index.isin(indexes_of_overlapping)].index
        return df.loc[indexes_of_nonoverlapping]

# def _intersection_cut():

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


def _inverse_intersection(self, other, strandedness=False):

    indexes_of_overlapping = []
    df = self.df

    other_strand = {"+": "-", "-": "+"}
    same_strand = {"+": "+", "-": "-"}

    if strandedness == "same":
        strand_dict = same_strand
    elif strandedness == "opposite":
        strand_dict = other_strand

    if "Strand" in df and "Strand" in other.df and strandedness:

        rowdicts = []
        for (chromosome, strand), cdf in df.groupby("Chromosome Strand".split()):

            strand = strand_dict[strand]

            it = other.__intervaltrees__[chromosome, strand]

            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(),
                                       (cdf.End - 1).tolist()):

                hits = it.search(start, end)

                # print("idx, start, end")
                # print(idx, start, end)
                for i in hits:

                    ostart, oend, oindex = i.data
                    # print("ostart, oend, oindex")
                    # print(ostart, oend, oindex)
                    if start - ostart < 0:
                        # print("if start - ostart < 0")
                        # print(start, ostart)
                        rowdicts.append({"Chromosome": chromosome, "Start":
                                         start, "End": ostart, "Index": idx})
                    elif end - oend > 0:
                        # print("end - oend > 0")
                        # print(end, oend)
                        rowdicts.append({"Chromosome": chromosome, "Start":
                                         end, "End": oend, "Index": idx})

                if not hits:
                    # print(idx, start, end)
                    rowdicts.append({"Chromosome": chromosome, "Start": start,
                                     "End": end + 1, "Index": idx})

    else:

        rowdicts = []
        for chromosome, cdf in df.groupby("Chromosome"):

            it = other.__intervaltrees__[chromosome, "+"]
            it2 = other.__intervaltrees__[chromosome, "-"]

            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(),
                                       (cdf.End - 1).tolist()):

                hits = it.search(start, end) + it2.search(start, end)
                for i in hits:

                    ostart, oend, oindex = i.data
                    if start - ostart < 0:
                        # print("if start - ostart < 0")
                        # print(start, ostart)
                        rowdicts.append({"Chromosome": chromosome, "Start":
                                         start, "End": ostart, "Index": idx})
                    elif end - oend > 0:
                        # print("end - oend > 0")
                        # print(end, oend)
                        rowdicts.append({"Chromosome": chromosome, "Start":
                                         end, "End": oend, "Index": idx})

                if not hits:
                    # print(idx, start, end)
                    rowdicts.append({"Chromosome": chromosome, "Start": start,
                                     "End": end + 1, "Index": idx})


    s_e_df = pd.DataFrame.from_dict(rowdicts).set_index("Index")["Chromosome Start End".split()]
    outdf = s_e_df.join(df.drop("Chromosome Start End".split(), axis=1))

    return outdf


def _intersection(self, other, strandedness=False):

    indexes_of_overlapping = []
    df = self.df

    other_strand = {"+": "-", "-": "+"}
    same_strand = {"+": "+", "-": "-"}

    if strandedness == "same":
        strand_dict = same_strand
    elif strandedness == "opposite":
        strand_dict = other_strand

    if "Strand" in df and "Strand" in other.df and strandedness:

        rowdicts = []
        for (chromosome, strand), cdf in df.groupby("Chromosome Strand".split()):

            strand = strand_dict[strand]

            it = other.__intervaltrees__[chromosome, strand]

            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(),
                                       (cdf.End - 1).tolist()):

                hits = it.search(start, end)

                # print("idx, start, end")
                # print(idx, start, end)
                for i in hits:

                    ostart, oend, oindex = i.data
                    # print("ostart, oend, oindex")
                    # print(ostart, oend, oindex)
                    rowdicts.append({"Chromosome": chromosome, "Start":
                                     max(start, ostart), "End": min(end, oend) + 1,
                                     "Index": idx})

    else:

        rowdicts = []
        for chromosome, cdf in df.groupby("Chromosome"):

            it = other.__intervaltrees__[chromosome, "+"]
            it2 = other.__intervaltrees__[chromosome, "-"]

            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(),
                                       (cdf.End - 1).tolist()):

                hits = it.search(start, end) + it2.search(start, end)
                for i in hits:

                    ostart, oend, oindex = i.data
                    rowdicts.append({"Chromosome": chromosome, "Start": max(start, ostart),
                                     "End": min(end, oend) + 1, "Index": idx})


    s_e_df = pd.DataFrame.from_dict(rowdicts)

    if s_e_df.empty:
        return pd.DataFrame(columns="Chromosome Start End".split())

    s_e_df = s_e_df.set_index("Index")["Chromosome Start End".split()]
    outdf = s_e_df.join(df.drop("Chromosome Start End".split(), axis=1))

    return outdf


def _coverage(ranges, value_col=None):

    try:
        import pyrle as rle
    except ImportError:
        raise Exception("Using the coverage method requires that pyrle is installed.")

    return rle.coverage(ranges)
