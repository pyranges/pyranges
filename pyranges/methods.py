import pandas as pd
import numpy as np

import itertools

# import multiprocessing
# from multiprocessing import Process
# from joblib import Parallel, delayed

#nearest_next,
# nearest_previous)(nearest, nearest_nonoverlapping,

from collections import OrderedDict

from sorted_nearest import (nearest_previous_nonoverlapping,
                            nearest_next_nonoverlapping,
                            nearest_nonoverlapping, find_clusters)
try:
    dummy = profile
except:
    profile = lambda x: x


# from clustertree import find_clusters

import pyranges as pr

from time import time
import datetime

from collections import defaultdict

def pick_out_indexes_possibly_nonunique(df, indexes, invert=False):

    if isinstance(indexes, list) and indexes:
        concat = np.concatenate(indexes)
        indexes = np.unique(concat)

    if not invert:
        return df.loc[df.index.isin(indexes)]
    else:
        return df.loc[~df.index.isin(indexes)]


def _overlap(self, other, strandedness, invert, how=None):

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
            if how == "containment":
                indexes = it.has_containment(starts, ends, indexes)
            else:
                indexes = it.has_overlaps(starts, ends, indexes)


            dfs[chromosome, strand] = pick_out_indexes_possibly_nonunique(cdf, indexes, invert)

    else:

        for chromosome, cdf in df.groupby("Chromosome"):

            starts = cdf.Start.values
            ends = cdf.End.values
            indexes = cdf.index.values

            it = other.__ncls__[chromosome, "+"]
            it2 = other.__ncls__[chromosome, "-"]

            if how == "containment":
                idx1 = it.has_containments(starts, ends, indexes)
                idx2 = it2.has_containments(starts, ends, indexes)
            else:
                idx1 = it.has_overlaps(starts, ends, indexes)
                idx2 = it2.has_overlaps(starts, ends, indexes)

            dfs[chromosome] = pick_out_indexes_possibly_nonunique(cdf, [idx1, idx2], invert)

    return pd.concat(dfs.values())


def both_indexes(self, other, strandedness, how=None):

    assert strandedness in ["same", "opposite", False, None]
    assert how in "containment first".split() + [False, None]

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

    self_d = OrderedDict()
    other_d = OrderedDict()

    if how == "first" and not strandedness and other.stranded:
        other_dfs = OrderedDict([(k, v) for k, v in other.df.groupby("Chromosome")])

    if strandedness: # know from assertion above that both have strand info

        for (chromosome, strand), cdf in df.groupby("Chromosome Strand".split()):

            starts = cdf.Start.values
            ends = cdf.End.values
            indexes = cdf.index.values

            strand = strand_dict[strand]
            it = other.__ncls__[chromosome, strand]

            if not how:
                _self_indexes, _other_indexes = it.all_overlaps_both(starts, ends, indexes)
            elif how == "containment":
                _self_indexes, _other_indexes = it.all_containments_both(starts, ends, indexes)
            else:
                _self_indexes, _other_indexes = it.first_overlap_both(starts, ends, indexes)


            self_d[chromosome, strand] = _self_indexes
            other_d[chromosome, strand] = _other_indexes

    # if not strandedness, and other df has strand, need to search both ncls
    elif other.stranded:

        # print("elif other stranded " * 10, )
        # print("self\n", self)
        # print("other\n", other)
        for chromosome, cdf in df.groupby("Chromosome"):

            it = other.__ncls__[chromosome, "+"]
            it2 = other.__ncls__[chromosome, "-"]

            starts = cdf.Start.values
            ends = cdf.End.values
            indexes = cdf.index.values

            if not how:
                # print("all overlaps " * 10)
                _self_indexes, _other_indexes = it.all_overlaps_both(starts, ends, indexes)
                _self_indexes2, _other_indexes2 = it2.all_overlaps_both(starts, ends, indexes)
            elif how == "containment":
                # print("all containments " * 10)
                _self_indexes, _other_indexes = it.all_containments_both(starts, ends, indexes)
                _self_indexes2, _other_indexes2 = it2.all_containments_both(starts, ends, indexes)
            else:
                _self_indexes, _other_indexes = it.first_overlap_both(starts, ends, indexes)
                _self_indexes2, _other_indexes2 = it2.first_overlap_both(starts, ends, indexes)

                if chromosome in other_dfs:
                    odf = other_dfs[chromosome]
                else:
                    continue

                if not len(_other_indexes2):
                    pass
                elif not len(_other_indexes):
                    other_indexes = _other_indexes2
                else:
                    s = pd.Series(np.concatenate([_self_indexes, _self_indexes2])).duplicated(keep=False)
                    _dup1, _dup2 = s[:len(_self_indexes)].values, s[len(_self_indexes):].values
                    # print("dup1", _dup1)
                    # print("dup2", _dup2)
                    # print(_other_indexes[_dup1])
                    # print(_other_indexes2[_dup2])
                    cond = odf.loc[_other_indexes[_dup1]].Start.values < odf.loc[_other_indexes2[_dup2]].Start.values
                    # print("os", odf.loc[_other_indexes[_dup1]].Start.values)
                    # print("os2", odf.loc[_other_indexes2[_dup2]].Start.values)
                    # print("cond", cond)

                    _other_indexes_tmp = np.where(cond, _other_indexes[_dup1], _other_indexes2[_dup2])
                    _self_indexes_tmp = np.where(cond, _self_indexes[_dup1], _self_indexes2[_dup2])
                    # print(_self_indexes, ~_dup1)
                    # print("_self_indexes[~_dup1]", _self_indexes[~_dup1])
                    # print("_self_indexes[~_dup2]", _self_indexes2[~_dup2])
                    # print("_other_indexes[~_dup1]", _other_indexes[~_dup1])
                    # print("_other_indexes[~_dup2]", _other_indexes2[~_dup2])
                    # _self_indexes2 = np.concatenate(, _self_indexes2[~_dup2])
                    _self_indexes2 = np.concatenate([_self_indexes[~_dup1], _self_indexes2[~_dup2]])
                    _other_indexes2 = np.concatenate([_other_indexes[~_dup1], _other_indexes2[~_dup2]])

                    _self_indexes = _self_indexes_tmp
                    _other_indexes = _other_indexes_tmp
                    # print("~cond", ~cond)
                    # print("~_dup1", ~_dup1, _dup1)
                    # print("~_dup2", ~_dup2, _dup2)
                    # print("~_dup1", _other_indexes[~_dup1], _self_indexes[~_dup1])
                    # print("~_dup2", _other_indexes[~_dup2], _self_indexes[~_dup2])
                    # _other_indexes = np.where(~cond, _other_indexes[~_dup1], _other_indexes2[~_dup2])
                    # _self_indexes = np.where(~cond, _self_indexes[~_dup1], _self_indexes2[~_dup2])

            self_d[chromosome] = np.concatenate([_self_indexes, _self_indexes2])
            other_d[chromosome] = np.concatenate([_other_indexes, _other_indexes2])

            if not len(_self_indexes) or not len(_self_indexes2):
                self_d[chromosome] = np.array(self_d[chromosome], dtype=int)
                other_d[chromosome] = np.array(other_d[chromosome], dtype=int)


    # if not strandedness, and other df has no strand info, need to search one ncls
    elif not other.stranded:

        for chromosome, cdf in df.groupby("Chromosome"):

            it = other.__ncls__[chromosome]

            starts = cdf.Start.values
            ends = cdf.End.values
            indexes = cdf.index.values

            if not how:
                _self_indexes, _other_indexes = it.all_overlaps_both(starts, ends, indexes)
            elif how == "containment":
                _self_indexes, _other_indexes = it.all_containments_both(starts, ends, indexes)
            else:
                _self_indexes, _other_indexes = it.first_overlap_both(starts, ends, indexes)

            self_d[chromosome] = _self_indexes
            other_d[chromosome] = _other_indexes

    return self_d, other_d


def sort_one_by_one(d, col1, col2):
    """
    Equivalent to pd.sort_values(by=[col1, col2]), but faster.
    """
    d = d.sort_values(by=[col2])
    return d.sort_values(by=[col1], kind='mergesort')


def _cluster(self, strand=False, maxdist=0, minnb=1):

    dfs = []

    idx_start, idx_end = 0, 0
    if strand:

        for (c, s), cdf in self.df.groupby(["Chromosome", "Strand"]):
            cdf = cdf.sort_values("Start")
            starts, ends = find_clusters(cdf.Start.values, cdf.End.values)
            df = pd.DataFrame({"Chromosome": c, "Start": starts, "End": ends, "Strand": s})
            dfs.append(df)

        df = pd.concat(dfs, ignore_index=True)["Chromosome Start End Strand".split()]

    else:

        for c, cdf in self.df.groupby(["Chromosome"]):
            cdf = cdf.sort_values("Start")
            starts, ends = find_clusters(cdf.Start.values, cdf.End.values)
            df = pd.DataFrame({"Chromosome": c, "Start": starts, "End": ends})
            dfs.append(df)

        df = pd.concat(dfs, ignore_index=True)["Chromosome Start End".split()]


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


# def _multithreaded_intersection(key, scdf, ocdf, strandedness=None, how=None):


#     # is there anything to gain?

#     # create ncls

#     # both indexes

#     # intersect op

#     # c, s = list(ocdf.head(0)["Chromosome Strand".split()].values)

#     starts = np.copy(scdf.Start.values)
#     ends = np.copy(scdf.End.values)
#     indexes = np.copy(scdf.index.values)

#     it = pr.PyRanges(ocdf).__ncls__[key]

#     _self_indexes, _other_indexes = it.all_overlaps_both(starts, ends, indexes)

#     scdf = scdf.loc[_self_indexes]
#     ocdf = ocdf.loc[_other_indexes, ["Start", "End"]]

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




def _intersection(self, other, strandedness=None, how=None):

    sidx, oidx = both_indexes(self, other, strandedness, how)

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

    if dfs:
        df = pd.concat(dfs)
    else:
        df = pd.DataFrame(columns="Chromosome Start End Strand".split())

    return df



def _coverage(ranges, value_col=None, stranded=False):

    try:
        from pyranges import PyRles
    except ImportError:
        raise Exception("Using the coverage method requires that pyrle is installed.")

    return PyRles(ranges, value_col=value_col, stranded=stranded)


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


def _overlap_write_both(self, other, strandedness=False, new_pos=None, suffixes=["_a", "_b"], suffix="_b", how=None):

    assert new_pos in ["intersection", "union", False, None]

    sidx, oidx = both_indexes(self, other, strandedness, how)

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
            _df = scdf.join(ocdf, rsuffix=suffix)

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


    if dfs:
        df = pd.concat(dfs)
    else:
        return pd.DataFrame(columns="Chromosome Start End Strand".split())

    return df


def _set_intersection(self, other, strandedness=None, how=None):

    strand = True if strandedness else False
    s = self.cluster(strand=strand)
    o = other.cluster(strand=strand)

    return _intersection(s, o, strandedness, how)


def _set_union(self, other, strand):

    if strand:
        assert self.stranded and other.stranded, \
            "Can only do stranded set union when both PyRanges contain strand info."

    if len(self) == 0:
        return _cluster(other, strand=strand)
    elif len(other) == 0:
        return _cluster(self, strand=strand)

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


        cdf = pd.DataFrame({"Start": _starts, "End": _ends})["Start End".split()]
        # print("cdf" * 10 + "\n", cdf)
        # cdf = sort_one_by_one(cdf, "Start", "End")
        cdf = cdf.sort_values("Start")
        # print("cdf" * 10 + "\n", cdf)
        # print(clusters)
        starts, ends = find_clusters(cdf.Start.values, cdf.End.values)
        # print(starts, ends)
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


def _subtraction(self, other, strandedness):

    #print(self.df.to_csv(sep=" "))
    strand = False
    if strandedness:
        assert "Strand" in self.df and "Strand" in other.df, \
            "Can only do stranded set difference when both PyRanges contain strand info."
        strand = True

    if len(self) == 0:
        return pd.DataFrame(columns="Chromosome Start End Strand".split())
    if len(other) == 0:
        return self.df

    other = other.cluster(strand=strand)

    other_strand = {"+": "-", "-": "+"}
    same_strand = {"+": "+", "-": "-"}

    if strandedness == "same":
        strand_dict = same_strand
    elif strandedness == "opposite":
        strand_dict = other_strand

    if strandedness and self.stranded: # chromosome and strand
        grpby_key = "Chromosome Strand".split()
    else:
        grpby_key = "Chromosome"

    # print("self.df", self.df.to_csv(sep=" "))


    self_dfs = {k: d for k, d in self.df.groupby(grpby_key)}
    other_dfs = {k: d for k, d in other.df.groupby(grpby_key)}

    dfs = []

    # need to store indexes in dict since if it is opposite stranded
    # cannot do removal in same iteration
    _idxes = dict()

    for key, scdf in self_dfs.items():

        # the downstream union of self/missing indexes sorts them
        # therefore, we need them to be in correct order here
        scdf.index = pd.Index(range(len(scdf)))
        print("-------")
        print("scdf", scdf.to_csv(sep=" "))


        if len(key) == 2:
            old_key = key
            key = key[0], strand_dict[key[1]]


        if not key in other_dfs:
            dfs.append(scdf)
            continue

        # scdf = self_dfs[key]
        ogr = pr.PyRanges(other_dfs[key])
        # ogr.index = pd.Index(range(len(ogr)))
        #print("scdf\n", scdf)
        #print("odf\n", ogr)

        # check which self reads have partial overlap
        if not ogr.stranded and not strandedness:
            #print("1" * 100)
            idx_self, new_starts, new_ends = ogr.__ncls__[key].set_difference_helper(
               scdf.Start.values,
                scdf.End.values,
                scdf.index.values)

        elif self.stranded and ogr.stranded and strandedness:

            idx_self, new_starts, new_ends = ogr.__ncls__[key].set_difference_helper(
                scdf.Start.values,
                scdf.End.values,
                scdf.index.values)

        elif ogr.stranded and not strandedness:

            #print("3" * 100)
            idx_self1, new_starts1, new_ends1 = ogr.__ncls__[key, "+"].set_difference_helper(
                scdf.Start.values,
                scdf.End.values,
                scdf.index.values)
            idx_self2, new_starts2, new_ends2 = ogr.__ncls__[key, "-"].set_difference_helper(
                scdf.Start.values,
                scdf.End.values,
                scdf.index.values)

            idx_self = np.concatenate([idx_self1, idx_self2])

            if not idx_self.dtype == np.int64:
                idx_self = np.array(idx_self, dtype=np.int64)

            new_starts = np.concatenate([new_starts1, new_starts2])
            new_ends = np.concatenate([new_ends1, new_ends2])

        #print("scdf.index", scdf.index)
        missing_idx = pd.Index(scdf.index).difference(idx_self)
        #print("missing_idx " * 10, missing_idx)
        idx_to_drop = new_starts != -1
        #print(idx_to_drop)
        new_starts = new_starts[idx_to_drop]
        new_ends = new_ends[idx_to_drop]

        idx_self = idx_self[idx_to_drop]
        new_starts = pd.Series(new_starts, index=idx_self).sort_index()
        new_ends = pd.Series(new_ends, index=idx_self).sort_index()
        # idx_self = np.sort(idx_self)

        print("scdf", scdf.to_csv(sep=" "))
        print("scdf.index", scdf.index)
        print("idx_self", idx_self)
        print("missing_idx", missing_idx)
        scdf = scdf.reindex(missing_idx.union(idx_self))
        print("scdf", scdf.to_csv(sep=" "))
        print("scdf.index", scdf.index)
        print("idx_self", idx_self)

        if len(idx_self):
            # why does boolean indexing work, but not regular?

            # pd.options.mode.chained_assignment = None  # default='warn'
            #print("idx_self", idx_self)
            scdf.loc[scdf.index.isin(idx_self), "Start"] = new_starts
            scdf.loc[scdf.index.isin(idx_self), "End"] = new_ends

        # find those in scdf that were not in idx_self

        #print("scdf", scdf)

        dfs.append(scdf)

    if dfs:
        df = pd.concat(dfs)
    else:
        return pd.DataFrame(columns="Chromosome Start End Strand".split())

    #print("dfs", dfs)
    #print("df", df)
    return df


# def _next_nonoverlapping(left_ends, right_starts):

#     left_ends = left_ends.sort_values()
#     right_starts = right_starts.sort_values()
#     l_idx, r_idx, dist = nearest_next_nonoverlapping(left_ends.values - 1, left_ends.index.values, right_starts.values, right_starts.index.values)

#     return l_idx, r_idx, dist


# def _previous_nonoverlapping(left_starts, right_ends):

#     left_starts = left_starts.sort_values()
#     right_ends = right_ends.sort_values()
#     l_idx, r_idx, dist = nearest_previous_nonoverlapping(left_starts.values, left_starts.index.values, right_ends.values - 1, right_ends.index.values)

#     return l_idx, r_idx, dist

def _next_nonoverlapping(left_ends, right_starts, right_indexes):

    left_ends = left_ends.sort_values()
    right_starts = right_starts.sort_values()
    r_idx, dist = nearest_next_nonoverlapping(left_ends.values - 1, right_starts.values, right_indexes)
    r_idx = pd.Series(r_idx, index=left_ends.index).sort_index().values
    dist = pd.Series(dist, index=left_ends.index).sort_index().values

    return r_idx, dist


def _previous_nonoverlapping(left_starts, right_ends, right_indexes):

    left_starts = left_starts.sort_values()
    right_ends = right_ends.sort_values()
    r_idx, dist = nearest_previous_nonoverlapping(left_starts.values, right_ends.values - 1, right_ends.index.values)
    # print("ridx before", r_idx)
    r_idx = pd.Series(r_idx, index=left_starts.index).sort_index().values
    dist = pd.Series(dist, index=left_starts.index).sort_index().values
    # print("ridx after", r_idx)

    return r_idx, dist

# from memory_profiler import profile


def _overlapping_for_nearest(self, other, strandedness, suffix):
    nearest_df = pd.DataFrame(columns="Chromosome Start End Strand".split())

    sidx, oidx = both_indexes(self, other, strandedness, how="first")

    if list(oidx.values()):
        idxs = np.concatenate(list(sidx.values()))
        missing_overlap = self.df.index[~self.df.index.isin(idxs)]
        # print("missing overlap", missing_overlap)
        df_to_find_nearest_in = self.df.reindex(missing_overlap)

        odf = other.df.reindex(np.concatenate(list(oidx.values())))
        odf.index = np.concatenate(list(sidx.values()))
        sdf = self.df.reindex(np.concatenate(list(sidx.values())))

        nearest_df = sdf.join(odf, rsuffix=suffix).drop("Chromosome" + suffix, axis=1)
        nearest_df.insert(nearest_df.shape[1], "Distance", 0)
    else:
        df_to_find_nearest_in = self.df

    return nearest_df, df_to_find_nearest_in


def _nearest(self, other, strandedness, suffix="_b", how=None, overlap=True):

    if overlap:
        nearest_df, df_to_find_nearest_in = _overlapping_for_nearest(self, other, strandedness, suffix)
    else:
        df_to_find_nearest_in = self.df

    other_strand = {"+": "-", "-": "+"}

    if self.stranded and strandedness: # chromosome and strand
        grpby_key = "Chromosome Strand".split()
    else:
        grpby_key = "Chromosome"

    other_dfs = {k: d for k, d in other.df.groupby(grpby_key)}

    dfs = []

    for key, scdf in df_to_find_nearest_in.groupby(grpby_key):

        if len(key) == 2 and strandedness == "opposite":
            other_key = key[0], other_strand[key[1]]
        else:
            other_key = key

        if not other_key in other_dfs:
            continue

        ocdf = other_dfs[other_key]

        scdf.index = pd.Index(range(len(scdf)))

        if how == "next":
            r_idx, dist = _next_nonoverlapping(scdf.End, ocdf.Start, ocdf.index.values)
        elif how == "previous":
            r_idx, dist = _previous_nonoverlapping(scdf.Start, ocdf.End, ocdf.index.values)
        else:
            previous_r_idx, previous_dist = _previous_nonoverlapping(scdf.Start, ocdf.End, ocdf.index.values)

            next_r_idx, next_dist = _next_nonoverlapping(scdf.End, ocdf.Start, ocdf.index.values)

            r_idx, dist = nearest_nonoverlapping(previous_r_idx,
                                                 previous_dist,
                                                 next_r_idx, next_dist)

        ocdf = ocdf.reindex(r_idx, fill_value=-1) # instead of np.nan, so ints are not promoted to float

        ocdf.index = scdf.index
        ocdf.insert(ocdf.shape[1], "Distance", pd.Series(dist, index=ocdf.index).fillna(-1).astype(int))
        ocdf.drop("Chromosome", axis=1, inplace=True)

        r_idx = pd.Series(r_idx, index=ocdf.index)
        scdf = scdf.drop(r_idx.loc[r_idx == -1].index)

        result = scdf.join(ocdf, rsuffix=suffix)

        dfs.append(result)

    if dfs:
        df = pd.concat(dfs)
    else:
        df = pd.DataFrame(columns="Chromosome Start End Strand".split())


    if overlap and not df.empty and not nearest_df.empty:
        df = pd.concat([nearest_df, df])
    elif overlap and not nearest_df.empty:
        df = nearest_df

    return df


def _slack(self, slack):

    df = self.df.copy()
    df.Start = df.Start - slack
    df.loc[df.Start < 0, "Start"] = 0
    df.End = df.End + slack

    return df


def _tss(self, slack=0):

    df = self.df

    tss_pos = df.loc[df.Strand == "+"]

    tss_neg = df.loc[df.Strand == "-"].copy()

    # pd.options.mode.chained_assignment = None
    tss_neg.loc[:, "Start"] = tss_neg.End

    # pd.options.mode.chained_assignment = "warn"
    tss = pd.concat([tss_pos, tss_neg], sort=False)
    tss["End"] = tss.Start
    tss.End = tss.End + 1 + slack
    tss.Start = tss.Start - slack
    tss.loc[tss.Start < 0, "Start"] = 0

    return tss.reindex(df.index)


def _tes(self, slack=0):

    df = self.df

    tes_pos = df.loc[df.Strand == "+"]

    tes_neg = df.loc[df.Strand == "-"].copy()

    # pd.options.mode.chained_assignment = None
    tes_neg.loc[:, "Start"] = tes_neg.End

    # pd.options.mode.chained_assignment = "warn"
    tes = pd.concat([tes_pos, tes_neg], sort=False)
    tes["Start"] = tes.End
    tes.End = tes.End + 1 + slack
    tes.Start = tes.Start - slack
    tes.loc[tes.Start < 0, "Start"] = 0

    return tes.reindex(df.index)


def _lengths(self):

    df = self.df
    lengths = df.End - df.Start

    return lengths


def _jaccard(self, other, strandedness):

    if strandedness:
        strand = True
    else:
        strand = False

    il = self.set_intersection(other, strandedness).lengths().sum()
    # ul = self.set_union(other, strandedness).lengths().sum()
    s = self.cluster(True).lengths().sum()
    o = other.cluster(True).lengths().sum()

    # print("o", o)
    # print("s", s)
    # print("il", il)
    # print("ul", ul)



    # if il == ul:
    #     return 1
    if s + o == il:
        return 1

    return il / (s + o - il)
