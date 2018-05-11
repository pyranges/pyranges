import pandas as pd
import numpy as np

import itertools

# import multiprocessing
# from multiprocessing import Process
# from joblib import Parallel, delayed


from collections import OrderedDict

try:
    dummy = profile
except:
    profile = lambda x: x

# from pyranges.src.cython_methods import c_overlap, both_indexes

from time import time
import datetime

from collections import defaultdict

@profile
def pick_out_indexes_possibly_nonunique(df, indexes, invert=False):

    if isinstance(indexes, list):
        concat = np.concatenate(indexes)
        indexes = np.unique(concat)

    if not invert:
        return df.loc[df.index.isin(indexes)]
    else:
        return df.loc[~df.index.isin(indexes)]


@profile
def _overlap(self, other, strandedness, invert):

    assert strandedness in ["same", "opposite", False, None]

    df = self.df
    other_strand = {"+": "-", "-": "+"}
    same_strand = {"+": "+", "-": "-"}

    if strandedness == "same":
        strand_dict = same_strand
    elif strandedness == "opposite":
        strand_dict = other_strand

    dfs = OrderedDict()

    if "Strand" in df and "Strand" in other.df and strandedness:

        indexes_of_overlapping = []
        for (chromosome, strand), cdf in df.groupby("Chromosome Strand".split()):

            starts = cdf.Start.values
            ends = cdf.End.values
            indexes = cdf.index.values

            strand = strand_dict[strand]
            it = other.__ncls__[chromosome, strand]
            indexes = it.has_overlaps(starts, ends, indexes)

            dfs[chromosome, strand] = pick_out_indexes_possibly_nonunique(cdf, indexes, invert)

    else:

        for chromosome, cdf in df.groupby("Chromosome"):

            starts = cdf.Start.values
            ends = cdf.End.values
            indexes = cdf.index.values

            it = other.__ncls__[chromosome, "+"]
            it2 = other.__ncls__[chromosome, "-"]

            idx1 = it.has_overlaps(starts, ends, indexes)
            idx2 = it2.has_overlaps(starts, ends, indexes)

            dfs[chromosome] = pick_out_indexes_possibly_nonunique(cdf, [idx1, idx2], invert)

    return pd.concat(dfs.values())


def both_indexes(self, other, strandedness):

    assert strandedness in ["same", "opposite", False, None]

    df = self.df

    other_strand = {"+": "-", "-": "+"}
    same_strand = {"+": "+", "-": "-"}

    if strandedness == "same":
        strand_dict = same_strand
    elif strandedness == "opposite":
        strand_dict = other_strand

    self_d = defaultdict(list)
    other_d = defaultdict(list)

    if "Strand" in df and "Strand" in other.df and strandedness:

        for (chromosome, strand), cdf in df.groupby("Chromosome Strand".split()):

            starts = cdf.Start.values
            ends = cdf.End.values
            indexes = cdf.index.values

            strand = strand_dict[strand]
            it = other.__ncls__[chromosome, strand]

            _self_indexes, _other_indexes = it.all_overlaps_both(starts, ends, indexes)
            self_d[chromosome, strand].extend(_self_indexes)
            other_d[chromosome, strand].extend(_other_indexes)

    else:

        for chromosome, cdf in df.groupby("Chromosome"):

            it = other.__ncls__[chromosome, "+"]
            it2 = other.__ncls__[chromosome, "-"]

            starts = cdf.Start.values
            ends = cdf.End.values
            indexes = cdf.index.values

            for _it in [it, it2]:
                _self_indexes, _other_indexes = _it.all_overlaps_both(starts, ends, indexes)
                self_d[chromosome].extend(_self_indexes)
                other_d[chromosome].extend(_other_indexes)


    return self_d, other_d



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
    columns = list(df.columns)

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

            it = other.__ncls__[chromosome, strand]

            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(),
                                       (cdf.End).tolist()):

                hits = it.find_overlap_list(start, end)

                for ostart, oend, _ in hits:

                    if start - ostart < 0:
                        rowdicts.append({"Chromosome": chromosome, "Start":
                                         start, "End": ostart, "Index": idx})
                    elif end - oend > 0:
                        rowdicts.append({"Chromosome": chromosome, "Start":
                                         end, "End": oend - 1, "Index": idx})

                if not hits:
                    rowdicts.append({"Chromosome": chromosome, "Start": start,
                                     "End": end, "Index": idx})

    else:

        rowdicts = []
        for chromosome, cdf in df.groupby("Chromosome"):

            it = other.__ncls__[chromosome, "+"]
            it2 = other.__ncls__[chromosome, "-"]

            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(),
                                       (cdf.End).tolist()):

                hits = it.find_overlap_list(start, end) + it2.find_overlap_list(start, end)

                for ostart, oend, _ in hits:

                    print("ostart, oend", ostart, oend)

                    if start - ostart < 0:
                        print("here!")
                        rowdicts.append({"Chromosome": chromosome, "Start":
                                         start, "End": ostart, "Index": idx})
                    elif end - oend > 0:
                        rowdicts.append({"Chromosome": chromosome, "Start":
                                         end, "End": oend - 1, "Index": idx})

                if not hits:
                    rowdicts.append({"Chromosome": chromosome, "Start": start,
                                     "End": end, "Index": idx})


    if not rowdicts:
        return pd.DataFrame(columns=columns)

    s_e_df = pd.DataFrame.from_dict(rowdicts).set_index("Index")["Chromosome Start End".split()]
    outdf = s_e_df.join(df.drop("Chromosome Start End".split(), axis=1))

    return outdf

@profile
def _intersection(self, other, strandedness=None, njobs=1):

    sidx, oidx = both_indexes(self, other, strandedness)

    other_strand = {"+": "-", "-": "+"}

    if len(list(sidx.keys())[0]) == 2: # chromosome and strand
        grpby_key = "Chromosome Strand".split()
    else:
        grpby_key = "Chromosome"

    other_dfs = {k: d for k, d in other.df.groupby(grpby_key)}

    dfs = []

    for key, scdf in self.df.groupby(grpby_key):

        if len(key) == 2 and strandedness == "opposite":
            key = key[0], other_strand[key[1]]

        if not key in other_dfs:
            continue

        ocdf = other_dfs[key]

        scidx = sidx[key]
        ocidx = oidx[key]

        scdf = scdf.loc[scidx]
        ocdf = ocdf.loc[ocidx, ["Start", "End"]]

        new_starts = pd.Series(
            np.where(scdf.Start.values > ocdf.Start.values, scdf.Start, ocdf.Start),
            index=scdf.index, dtype=np.long)

        new_ends = pd.Series(
            np.where(scdf.End.values < ocdf.End.values, scdf.End, ocdf.End),
            index=scdf.index, dtype=np.long)

        pd.options.mode.chained_assignment = None  # default='warn'
        scdf.loc[:, "Start"] = new_starts
        scdf.loc[:, "End"] = new_ends
        pd.options.mode.chained_assignment = 'warn'

        dfs.append(scdf)

    df = pd.concat(dfs)

    return df




    # if len(sidx) == 0:
    #     return pd.DataFrame(columns="Chromosome Start End".split())

    # sdf = self.df.loc[sidx]
    # odf = other.df.loc[oidx, ["Start", "End"]]

    # new_starts = pd.Series(np.where(sdf.Start.values > odf.Start.values, sdf.Start, odf.Start), index=sdf.index, dtype=np.long)
    # new_ends = pd.Series(np.where(sdf.End.values < odf.End.values, sdf.End, odf.End), index=sdf.index, dtype=np.long)

    # pd.options.mode.chained_assignment = None  # default='warn'
    # sdf.loc[:, "Start"] = new_starts
    # sdf.loc[:, "End"] = new_ends
    # pd.options.mode.chained_assignment = "warn"  # default='warn'

    # return sdf


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
            it = other.__ncls__[chromosome, strand]

            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(),
                                       (cdf.End).tolist()):

                for hit in it.find_overlap(start, end):
                    oidx = hit[2]
                    indexes_of_overlapping.append({"IndexSelf": idx, "IndexOther": oidx})


    elif "Strand" in df and "Strand" in other.df and strandedness == "opposite":

        for (chromosome, strand), cdf in df.groupby("Chromosome Strand".split()):
            opposite_strand = other_strand[strand]
            it = other.__ncls__[chromosome, opposite_strand]

            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(),
                                       (cdf.End).tolist()):

                for hit in it.find_overlap(start, end):
                    oidx = hit[2]
                    indexes_of_overlapping.append({"IndexSelf": idx, "IndexOther": oidx})

    else:

        for chromosome, cdf in df.groupby("Chromosome"):

            it = other.__ncls__[chromosome, "+"]
            it2 = other.__ncls__[chromosome, "-"]

            for idx, start, end in zip(cdf.index.tolist(), cdf.Start.tolist(),
                                       (cdf.End).tolist()):

                for hit in it.find_overlap_list(start, end) + it2.find_overlap_list(start, end):
                    oidx = hit[2]
                    indexes_of_overlapping.append({"IndexSelf": idx, "IndexOther": oidx})

    if not indexes_of_overlapping:
        return pd.DataFrame(columns="Chromosome Start End Strand".split())

    idx_df = pd.DataFrame.from_dict(indexes_of_overlapping).set_index("IndexSelf")

    idx_df = idx_df.join(self.df, how="right")

    idx_df = _keep_both(idx_df, other, suffixes, new_pos)

    return idx_df




# def __intersection(scdf, ocdf, scidx, ocidx):

#     scdf = scdf.loc[scidx]
#     ocdf = ocdf.loc[ocidx, ["Start", "End"]]


#     new_starts = pd.Series(
#         np.where(scdf.Start.values > ocdf.Start.values, scdf.Start, ocdf.Start),
#         index=scdf.index, dtype=np.long)

#     new_ends = pd.Series(
#         np.where(scdf.End.values < ocdf.End.values, scdf.End, ocdf.End),
#         index=scdf.index, dtype=np.long)

#     pd.options.mode.chained_assignment = None  # default='warn'
#     scdf.loc[:, "Start"] = new_starts
#     scdf.loc[:, "End"] = new_ends
#     pd.options.mode.chained_assignment = 'warn'

#     return scdf

# def _intersection(self, other, strandedness=None, njobs=1):


#     sidx, oidx = both_indexes(self, other, strandedness)

#     dfs = []
#     if len(list(sidx.keys())[0]) == 2: # chromosome and strand
#         grpby_key = ("Chromosome", "Strand")
#     else:
#         grpby_key = "Chromosome"

#     other_dfs = defaultdict(lambda: default, {k: d for k, d in other.df.groupby(grpby_key)})
#     other_keys = set(other_dfs)

#     if not grpby_key in other_keys and len(grpby_key) == 2:
#         default = pd.DataFrame(columns="Chromosome Start End Strand".split())
#     elif not grpby_key in other_keys:
#         default = pd.DataFrame(columns="Chromosome Start End".split())

#     dfs = Parallel(n_jobs=1)(
#         delayed(__intersection)(cdf, other_dfs[key], sidx[key], oidx[key])
#                                   for (key, cdf) in self.df.groupby(grpby_key))

#     df = pd.concat(dfs)

#     return df
