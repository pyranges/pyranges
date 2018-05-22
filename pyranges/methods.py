import pandas as pd
import numpy as np

import itertools

# import multiprocessing
# from multiprocessing import Process
# from joblib import Parallel, delayed


from collections import OrderedDict

# try:
#     dummy = profile
# except:
#     profile = lambda x: x

# from pyranges.src.cython_methods import c_overlap, both_indexes

from time import time
import datetime

from collections import defaultdict

def pick_out_indexes_possibly_nonunique(df, indexes, invert=False):

    if isinstance(indexes, list):
        concat = np.concatenate(indexes)
        indexes = np.unique(concat)

    if not invert:
        return df.loc[df.index.isin(indexes)]
    else:
        return df.loc[~df.index.isin(indexes)]


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

    if self.stranded and other.stranded and strandedness:

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

    if strandedness:
        assert self.stranded and other.stranded, \
            "Can only do stranded searches when both PyRanges contain strand info"

    df = self.df

    other_strand = {"+": "-", "-": "+"}
    same_strand = {"+": "+", "-": "-"}

    if strandedness == "same":
        strand_dict = same_strand
    elif strandedness == "opposite":
        strand_dict = other_strand

    self_d = defaultdict(list)
    other_d = defaultdict(list)

    if strandedness: # know from assertion above that both have strand info

        for (chromosome, strand), cdf in df.groupby("Chromosome Strand".split()):

            starts = cdf.Start.values
            ends = cdf.End.values
            indexes = cdf.index.values

            strand = strand_dict[strand]
            it = other.__ncls__[chromosome, strand]

            _self_indexes, _other_indexes = it.all_overlaps_both(starts, ends, indexes)
            self_d[chromosome, strand].extend(_self_indexes)
            other_d[chromosome, strand].extend(_other_indexes)

    # if not strandedness, and other df has strand, need to search both ncls
    elif other.stranded:

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

    # if not strandedness, and other df has no strand info, need to search one ncls
    elif not other.stranded:

        for chromosome, cdf in df.groupby("Chromosome"):

            it = other.__ncls__[chromosome]

            starts = cdf.Start.values
            ends = cdf.End.values
            indexes = cdf.index.values

            _self_indexes, _other_indexes = it.all_overlaps_both(starts, ends, indexes)
            self_d[chromosome].extend(_self_indexes)
            other_d[chromosome].extend(_other_indexes)

    return self_d, other_d



def _cluster(self, strand=False, maxdist=0, minnb=1):

    from clustertree import find_clusters

    dfs = []

    idx_start, idx_end = 0, 0
    if strand:
        for (c, s), cdf in self.df.groupby(["Chromosome", "Strand"]):
            starts, ends = find_clusters(maxdist, minnb, cdf.Start.values, cdf.End.values)
            idx_end += len(starts)
            _df = pd.DataFrame({"Chromosome": c, "Start": starts, "End": ends, "Strand": s})
            _df.index = range(idx_start, idx_end)
            dfs.append(_df)
            idx_start = idx_end

        df = pd.concat(dfs)["Chromosome Start End Strand".split()]

    else:
        for c, cdf in self.df.groupby(["Chromosome"]):
            starts, ends = find_clusters(maxdist, minnb, cdf.Start.values, cdf.End.values)
            idx_end += len(starts)
            _df = pd.DataFrame({"Chromosome": c, "Start": starts, "End": ends})
            _df.index = range(idx_start, idx_end)
            dfs.append(_df)
            idx_start = idx_end

        df = pd.concat(dfs)["Chromosome Start End".split()]

    return df



def _tile(self, tile_size=50):

    "No no no this is slow! Write in C."

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

def _intersection(self, other, strandedness=None):

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

    assert new_pos in ["intersection", "union", False, None]

    sidx, oidx = both_indexes(self, other, strandedness)

    other_strand = {"+": "-", "-": "+"}

    if len(list(sidx.keys())[0]) == 2: # chromosome and strand
        grpby_key = "Chromosome Strand".split()
    else:
        grpby_key = "Chromosome"

    other_dfs = {k: d for k, d in other.df.groupby(grpby_key)}

    dfs = []

    idx_counter, current_end = 0, 0
    for key, scdf in self.df.groupby(grpby_key):

        if len(key) == 2 and strandedness == "opposite":
            key = key[0], other_strand[key[1]]

        if not key in other_dfs:
            continue

        ocdf = other_dfs[key].drop("Chromosome", 1)

        scidx = sidx[key]
        ocidx = oidx[key]

        scdf = scdf.loc[scidx]
        ocdf = ocdf.loc[ocidx]

        pd.options.mode.chained_assignment = None  # default='warn'
        pd.options.mode.chained_assignment = 'warn'

        current_end += len(scdf)
        scdf.index = range(idx_counter, current_end)
        ocdf.index = range(idx_counter, current_end)
        idx_counter += len(scdf)

        if not new_pos:
            _df = scdf.join(ocdf, rsuffix=suffixes[1])

        elif new_pos == "intersection":

            new_starts = pd.Series(
                np.where(scdf.Start.values > ocdf.Start.values, scdf.Start, ocdf.Start),
                index=scdf.index, dtype=np.long)

            new_ends = pd.Series(
                np.where(scdf.End.values < ocdf.End.values, scdf.End, ocdf.End),
                index=scdf.index, dtype=np.long)
            _df = scdf.join(ocdf, lsuffix=suffixes[0], rsuffix=suffixes[1])
            _df.insert(1, "Start", new_starts)
            _df.insert(2, "End", new_ends)
            _df.rename(index=str, columns={"Chromosome" + suffixes[0]: "Chromosome", "Strand" + suffixes[0]: "Strand"}, inplace=True)

        elif new_pos == "union":

            new_starts = pd.Series(
                np.where(scdf.Start.values < ocdf.Start.values, scdf.Start, ocdf.Start),
                index=scdf.index, dtype=np.long)

            new_ends = pd.Series(
                np.where(scdf.End.values > ocdf.End.values, scdf.End, ocdf.End),
                index=scdf.index, dtype=np.long)
            _df = scdf.join(ocdf, lsuffix=suffixes[0], rsuffix=suffixes[1])
            _df.insert(1, "Start", new_starts)
            _df.insert(2, "End", new_ends)
            _df.rename(index=str, columns={"Chromosome" + suffixes[0]: "Chromosome", "Strand" + suffixes[0]: "Strand"}, inplace=True)

        dfs.append(_df)


    df = pd.concat(dfs)

    return df


def _set_intersection(self, other, strandedness=None):

    strand = True if strandedness else False
    s = self.cluster(strand=strand)
    o = other.cluster(strand=strand)

    return _intersection(s, o, strandedness)


def _set_union(self, other, strand):

    if strand:
        assert "Strand" in self.df and "Strand" in other.df, \
            "Can only do stranded set union when both PyRanges contain strand info."

    from clustertree import find_clusters

    dfs = []

    if strand and len(list(self.__ncls__.keys())[0]) == 2: # chromosome and strand
        grpby_key = "Chromosome Strand".split()
        columns = "Chromosome Start End Strand".split()
    else:
        grpby_key = "Chromosome"
        columns = "Chromosome Start End".split()


    self_dfs =  {k: d for k, d in self.df.groupby(grpby_key)}
    other_dfs = {k: d for k, d in other.df.groupby(grpby_key)}

    idx_start, idx_end = 0, 0
    for key in set(self_dfs).union(other_dfs):

        if key in other_dfs and key in self_dfs:
            _starts = np.concatenate([
                self_dfs[key].Start.values,
                other_dfs[key].Start.values])
            _ends = np.concatenate([
                self_dfs[key].End.values,
                other_dfs[key].End.values])
        elif key in self_dfs and not key in other_dfs:
            _starts = self_dfs[key].Start.values
            _ends = self_dfs[key].End.values
        elif key in other_dfs and not key in self_dfs:
            _starts = other_dfs[key].Start.values
            _ends = other_dfs[key].End.values

        starts, ends = find_clusters(0, 1, _starts, _ends)
        idx_end += len(starts)

        if strand:
            _df = pd.DataFrame({"Chromosome": key[0], "Start": starts, "End": ends, "Strand": key[1]})
        else:
            _df = pd.DataFrame({"Chromosome": key, "Start": starts, "End": ends})

        _df.index = range(idx_start, idx_end)
        dfs.append(_df)
        idx_start = idx_end

    df = pd.concat(dfs)[columns]

    return df


def _nearest(self, other, strandedness, suffix="_b", is_sorted=False):

    try:
        from sorted_nearest import nearest
    except:
        raise ImportError("Cannot find package sorted_nearest")

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

        if not is_sorted:
            ocdf = ocdf.sort_values("Start End".split())
            scdf = scdf.sort_values("Start End".split())

        ocdf.index = range(len(ocdf))

        r_idx, dist = nearest(scdf.Start.values, scdf.End.values, ocdf.Start.values, ocdf.End.values)

        ocdf = ocdf.loc[pd.Series(r_idx)]
        ocdf.index = scdf.index
        ocdf.insert(ocdf.shape[1], "Distance", pd.Series(dist, index=ocdf.index))

        result = scdf.join(ocdf, rsuffix=suffix)

        dfs.append(result)

    return pd.concat(dfs)


def _subtraction(self, other, strandedness):

    strand = False
    if strandedness:
        assert "Strand" in self.df and "Strand" in other.df, \
            "Can only do stranded set difference when both PyRanges contain strand info."
        strand = True

    other = other.cluster(strand=strand)

    self_stranded = len(list(self.__ncls__.keys())[0]) == 2
    other_stranded = len(list(other.__ncls__.keys())[0]) == 2

    # print("self", self.__ncls__.keys())
    # print("other", other.__ncls__.keys())

    # print("self")
    # print(self)
    # print("other")
    # print(other)

    other_strand = {"+": "-", "-": "+"}
    same_strand = {"+": "+", "-": "-"}

    if strandedness == "same":
        strand_dict = same_strand
    elif strandedness == "opposite":
        strand_dict = other_strand

    if strandedness and self_stranded: # chromosome and strand
        grpby_key = "Chromosome Strand".split()
        columns = "Chromosome Start End Strand".split()
    else:
        grpby_key = "Chromosome"
        columns = "Chromosome Start End".split()

    self_dfs = {k: d for k, d in self.df.groupby(grpby_key)}
    other_dfs = {k: d for k, d in other.df.groupby(grpby_key)}

    dfs = []

    # need to store indexes in dict since if it is opposite stranded
    # cannot do removal in same iteration
    _idxes = dict()

    for key, scdf in self_dfs.items():

        # print("key before switch", key)
        if len(key) == 2:
            old_key = key
            key = key[0], strand_dict[key[1]]


        # print("key after switch", key)
        if not key in other_dfs:
            dfs.append(scdf)
            continue

        # scdf = self_dfs[key]
        odf = other_dfs[key]

        # check which self reads have partial overlap
        if not self_stranded:
            idx_self, idx_other, overlap_type = self.__ncls__[key].set_difference_helper(
                odf.Start.values,
                odf.End.values,
                odf.index.values)

        elif self_stranded and other_stranded and strandedness:
            # print("We are here!")
            idx_self, idx_other, overlap_type = self.__ncls__[old_key].set_difference_helper(
                odf.Start.values,
                odf.End.values,
                odf.index.values)
        elif self_stranded and not strandedness:
            idx_self1, idx_other1, overlap_type1 = self.__ncls__[key, "+"].set_difference_helper(
                odf.Start.values,
                odf.End.values,
                odf.index.values)
            idx_self2, idx_other2, overlap_type2 = self.__ncls__[key, "-"].set_difference_helper(
                odf.Start.values,
                odf.End.values,
                odf.index.values)

            idx_self = np.concatenate([idx_self1, idx_self2])

            if not idx_self.dtype == np.int64:
                idx_self = np.array(idx_self, dtype=np.int64)

            idx_other = np.concatenate([idx_other1, idx_other2])
            overlap_type = np.concatenate([overlap_type1, overlap_type2])


        # check what self reads have overlap at all
        if not other_stranded and self_stranded:
            ix_has_overlaps = other.__ncls__[key].has_overlaps(scdf.Start.values, scdf.End.values, scdf.index.values)
        elif other_stranded and not self_stranded:
            ix1 = other.__ncls__[key, "+"].has_overlaps(scdf.Start.values, scdf.End.values, scdf.index.values)
            ix2 = other.__ncls__[key, "-"].has_overlaps(scdf.Start.values, scdf.End.values, scdf.index.values)
            ix_has_overlaps = np.concat([ix1, ix2])
        elif other_stranded and self_stranded and not strandedness:
            ix_has_overlaps1 = other.__ncls__[key, "+"].has_overlaps(scdf.Start.values, scdf.End.values, scdf.index.values)
            ix_has_overlaps2 = other.__ncls__[key, "-"].has_overlaps(scdf.Start.values, scdf.End.values, scdf.index.values)
            ix_has_overlaps = np.unique(np.concatenate([ix_has_overlaps1, ix_has_overlaps2]))
        else:
            # print("other_stranded", other_stranded, "self_stranded", self_stranded)
            ix_has_overlaps = other.__ncls__[key].has_overlaps(scdf.Start.values, scdf.End.values, scdf.index.values)

        # print("self\n", self)
        # print("other\n", other)
        # print("idx_self", idx_self)
        # print("scdf", scdf)
        # print("idx_other", idx_other)
        # print("odf", odf)

        if strandedness == "opposite":
            scdf = self_dfs[old_key]
        else:
            pass

        starts = np.where(overlap_type == 0, scdf.loc[idx_self, "Start"], odf.loc[idx_other, "End"])
        ends = np.where(overlap_type == 0, odf.loc[idx_other, "Start"], scdf.loc[idx_self, "End"])
        starts = pd.Series(starts, index=idx_self)
        ends = pd.Series(ends, index=idx_self)

        # print("starts\n", starts)
        # print("ends\n", ends)

        # if len(key) == 2:
        #     # need to swtich back in case opposite strand
        #     _idxes[(key[0], strand_dict[key[1]]), "idx_self"] = idx_self
        #     _idxes[(key[0], strand_dict[key[1]]), "idx_other"] = idx_other
        #     _idxes[(key[0], strand_dict[key[1]]), "overlap_type"] = overlap_type
        #     _idxes[(key[0], strand_dict[key[1]]), "ix_has_overlaps"] = ix_has_overlaps
        # else:
        #     _idxes[key, "idx_self"] = idx_self
        #     _idxes[key, "idx_other"] = idx_other
        #     _idxes[key, "overlap_type"] = overlap_type
        #     _idxes[key, "ix_has_overlaps"] = ix_has_overlaps


        ix_no_overlaps = np.setdiff1d(scdf.index.values, ix_has_overlaps)
        # print(scdf.index.values)
        # print(ix_no_overlaps)

        # no_overlaps = cdf.loc[ix_no_overlaps]
        idx = np.unique(np.concatenate([ix_no_overlaps, idx_self]))
        idx.sort(kind="mergesort")
        # print(idx)

        # print(scdf)
        scdf = scdf.loc[idx]
        scdf.loc[idx_self, "Start"] = starts
        scdf.loc[idx_self, "End"] = ends
        # print("scdf", scdf)
        # print(scdf)

        scdf = scdf.sort_values(["Start", "End"])
        dfs.append(scdf)

    return pd.concat(dfs)
