import pandas as pd
import numpy as np

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

                for ostart, oend, _ in hits:

                    if start - ostart < 0:
                        rowdicts.append({"Chromosome": chromosome, "Start":
                                         start, "End": ostart, "Index": idx})
                    elif end - oend > 0:
                        rowdicts.append({"Chromosome": chromosome, "Start":
                                         end, "End": oend - 1, "Index": idx})

                if not hits:
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
                for ostart, oend, _ in hits:

                    if start - ostart < 0:
                        rowdicts.append({"Chromosome": chromosome, "Start":
                                         start, "End": ostart, "Index": idx})
                    elif end - oend > 0:
                        rowdicts.append({"Chromosome": chromosome, "Start":
                                         end, "End": oend - 1, "Index": idx})

                if not hits:
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

                for ostart, oend, _ in hits:
                    rowdicts.append({"Chromosome": chromosome, "Start":
                                     max(start, ostart), "End": min(end, oend - 1) + 1,
                                     "Index": idx})

    else:

        rowdicts = []
        for chromosome, cdf in df.groupby("Chromosome"):

            it = other.__intervaltrees__[chromosome, "+"]
            it2 = other.__intervaltrees__[chromosome, "-"]

            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(),
                                       (cdf.End - 1).tolist()):

                hits = it.search(start, end) + it2.search(start, end)
                for ostart, oend, _ in hits:
                    rowdicts.append({"Chromosome": chromosome, "Start": max(start, ostart),
                                     "End": min(end, oend - 1) + 1, "Index": idx})


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


def _keep_both(idx_df, other, suffixes, new_pos):

    self_suffix, other_suffix = suffixes

    if new_pos:
        assert self_suffix and other_suffix, "If new_pos is union or intersection, you must give two suffixes of length > 0."
    else:
        assert other_suffix, "Must give a suffix for the right dataframe."

    if not new_pos:
        suffixes = ["", other_suffix]
        to_drop = ["IndexOther", "Chromosome" + other_suffix]
        self_suffix = ""
        to_rename = {}
    else:
        suffixes = [self_suffix, other_suffix]
        to_drop = ["IndexOther", "Chromosome" + other_suffix]
        to_rename = {"Chromosome" + self_suffix: "Chromosome"}

    idx_df = idx_df.merge(other.df, left_on="IndexOther", right_index=True, suffixes=suffixes)
    idx_df = idx_df.drop(to_drop, axis=1)
    idx_df = idx_df.rename(to_rename, axis="columns")

    if new_pos == "intersect":
        new_starts = np.where(idx_df["Start" + self_suffix] > idx_df["Start" + other_suffix], idx_df["Start" + self_suffix], idx_df["Start" + other_suffix])
        new_ends = np.where(idx_df["End" + self_suffix] < idx_df["End" + other_suffix], idx_df["End" + self_suffix], idx_df["End" + other_suffix])
        idx_df.insert(1, "Start", new_starts)
        idx_df.insert(2, "End", new_ends)

    if new_pos == "union":
        new_starts = np.where(idx_df["Start" + self_suffix] < idx_df["Start" + other_suffix], idx_df["Start" + self_suffix], idx_df["Start" + other_suffix])
        new_ends = np.where(idx_df["End" + self_suffix] > idx_df["End" + other_suffix], idx_df["End" + self_suffix], idx_df["End" + other_suffix])
        idx_df.insert(1, "Start", new_starts)
        idx_df.insert(2, "End", new_ends)

    return idx_df


def _overlap_write_both(self, other, strandedness=False, new_pos=None, suffixes=["_a", "_b"]):

    indexes_of_overlapping = []
    df = self.df
    other_strand = {"+": "-", "-": "+"}

    if "Strand" in df and "Strand" in other.df and strandedness == "same":

        for (chromosome, strand), cdf in df.groupby("Chromosome Strand".split()):
            it = other.__intervaltrees__[chromosome, strand]

            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(),
                                       (cdf.End - 1).tolist()):

                for hit in it.search(start, end):
                    oidx = hit[2]
                    indexes_of_overlapping.append({"IndexSelf": idx, "IndexOther": oidx})


    elif "Strand" in df and "Strand" in other.df and strandedness == "opposite":

        for (chromosome, strand), cdf in df.groupby("Chromosome Strand".split()):
            opposite_strand = other_strand[strand]
            it = other.__intervaltrees__[chromosome, opposite_strand]

            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(),
                                       (cdf.End - 1).tolist()):

                for hit in it.search(start, end):
                    oidx = hit[2]
                    indexes_of_overlapping.append({"IndexSelf": idx, "IndexOther": oidx})

    else:

        for chromosome, cdf in df.groupby("Chromosome"):

            it = other.__intervaltrees__[chromosome, "+"]
            it2 = other.__intervaltrees__[chromosome, "-"]

            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(),
                                       (cdf.End - 1).tolist()):

                for hit in it.search(start, end) + it2.search(start, end):
                    oidx = hit[2]
                    indexes_of_overlapping.append({"IndexSelf": idx, "IndexOther": oidx})

    if not indexes_of_overlapping:
        return pd.DataFrame(columns="Chromosome Start End Strand".split())

    idx_df = pd.DataFrame.from_dict(indexes_of_overlapping).set_index("IndexSelf")

    idx_df = idx_df.join(self.df, how="right")

    idx_df = _keep_both(idx_df, other, suffixes, new_pos)

    return idx_df
