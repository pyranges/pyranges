import numpy as np
import pandas as pd
from ncls import NCLS

import pyranges as pr

from natsort import natsorted

# from joblib import Parallel, delayed

from sorted_nearest import (find_clusters, nearest_previous_nonoverlapping,
                            nearest_next_nonoverlapping, nearest_nonoverlapping, find_clusters)

from collections import defaultdict

from functools import wraps

def parse_grpby_key(grpby_key):

    if isinstance(grpby_key, str):
        return grpby_key, False
    else:
        return grpby_key[0], grpby_key[1]


def return_empty_if_one_empty(func):

    @wraps(func)
    def extended_func(self, other, **kwargs):

        if len(self) == 0 or len(other) == 0:
            df = pd.DataFrame(columns="Chromosome Start End".split())
        else:
            df = func(self, other, **kwargs)

        return df

    return extended_func


def _pyrange_apply(function, scdf, other_dfs, grpby_key, **kwargs):

    outdfs = []
    if function.__name__ == "_set_union":
        self_dfs =  {k: d for k, d in scdf.groupby(grpby_key)}
        self_dfs = defaultdict(lambda: pd.DataFrame(columns="Chromosome Start End".split()), self_dfs)
        keys_union = natsorted(list(set(self_dfs).union(other_dfs)))
        for key in keys_union:
            outdfs.append(function(self_dfs[key], other_dfs[key], key=key, **kwargs))

    else:
        for key, df in natsorted(scdf.groupby(grpby_key)):
            outdfs.append(function(df, other_dfs[key], key=key, **kwargs))

    outdfs = [df for df in outdfs if not df.empty]

    if outdfs:
        df = pd.concat(outdfs)
        return df
    else:
        return pd.DataFrame(columns="Chromosome Start End Strand".split())


def pyrange_apply_single(function, self, **kwargs):

    strand = kwargs["strand"]


    if strand:
        assert self.stranded, \
            "Can only do stranded operation when PyRange contains strand info"

    if self.stranded and strand:
        grpby_key = ["Chromosome", "Strand"]
    else:
        grpby_key = "Chromosome"

    outdfs = []
    for df in natsorted(self.df.groupby(grpby_key)):

        outdfs.append(function(df, **kwargs))

    if outdfs:
        df = pd.concat(outdfs)
        return df
    else:
        return pd.DataFrame(columns="Chromosome Start End Strand".split())


def pyrange_apply(function, self, other, **kwargs):

    strandedness = kwargs["strandedness"]

    assert strandedness in ["same", "opposite", False, None]

    if strandedness:
        assert self.stranded and other.stranded, \
            "Can only do stranded operations when both PyRanges contain strand info"

    if self.stranded and other.stranded and strandedness:
        grpby_key = ["Chromosome", "Strand"]
    else:
        grpby_key = "Chromosome"

    other_strand = {"+": "-", "-": "+"}
    if strandedness == "opposite":
        other_dfs = {(c, other_strand[s]): v for (c, s), v in other.df.groupby(grpby_key)}
        other_dfs = defaultdict(lambda: pd.DataFrame(columns="Chromosome Start End Strand".split()), other_dfs)
    else:
        other_dfs = {key: v for key, v in other.df.groupby(grpby_key)}
        other_dfs = defaultdict(lambda: pd.DataFrame(columns="Chromosome Start End".split()), other_dfs)

    return _pyrange_apply(function, self.df, other_dfs, grpby_key, **kwargs)





def pick_out_indexes_possibly_nonunique(df, indexes, invert=False):

    if isinstance(indexes, list) and indexes:
        concat = np.concatenate(indexes)
        indexes = np.unique(concat)

    if not invert:
        return df.loc[df.index.isin(indexes)]
    else:
        return df.loc[~df.index.isin(indexes)]



@return_empty_if_one_empty
def _first_df(scdf, ocdf, how=False, invert=False, **kwargs):

    assert how in "containment first".split() + [False, None]
    starts = scdf.Start.values
    ends = scdf.End.values
    indexes = scdf.index.values

    it = NCLS(ocdf.Start.values, ocdf.End.values, ocdf.index.values)

    if not how:
        _self_indexes = it.has_overlaps(starts, ends, indexes)
    else:
        _self_indexes = it.has_containment(starts, ends, indexes)

    idxs = scdf.index.isin(_self_indexes)
    if invert:
        return scdf.loc[~idxs]
    else:
        return scdf.loc[idxs]


def _both_dfs(scdf, ocdf, how=False, **kwargs):

    assert how in "containment first".split() + [False, None]
    starts = scdf.Start.values
    ends = scdf.End.values
    indexes = scdf.index.values

    it = NCLS(ocdf.Start.values, ocdf.End.values, ocdf.index.values)

    if not how:
        _self_indexes, _other_indexes = it.all_overlaps_both(starts, ends, indexes)
    elif how == "containment":
        _self_indexes, _other_indexes = it.all_containments_both(starts, ends, indexes)
    else:
        _self_indexes, _other_indexes = it.first_overlap_both(starts, ends, indexes)

    _self_indexes = _self_indexes
    _other_indexes = _other_indexes

    return scdf.loc[_self_indexes], ocdf.loc[_other_indexes]




@return_empty_if_one_empty
def _intersection(scdf, ocdf, **kwargs):

    scdf, ocdf = _both_dfs(scdf, ocdf, **kwargs)

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

    return scdf


def _create_df_from_starts_ends(starts, ends, chromosome, strand=None):

    nidx = pd.Index(range(len(starts)))
    if strand:
        cluster_df = pd.DataFrame({"Chromosome": pd.Series(chromosome, dtype="category", index=nidx),
                                   "Start": starts, "End": ends,
                                   "Strand": pd.Series(strand, dtype="category", index=nidx)})
    else:
        cluster_df = pd.DataFrame({"Chromosome": pd.Series(chromosome, dtype="category", index=nidx), "Start": starts, "End": ends})

    return cluster_df


def _cluster(df, chromosome, strand=False, **kwargs):

    cdf = df.sort_values("Start")
    starts, ends = find_clusters(cdf.Start.values, cdf.End.values)

    cluster_df = _create_df_from_starts_ends(starts, ends, chromosome, strand)

    return cluster_df


@return_empty_if_one_empty
def _set_intersection(scdf, ocdf, strandedness=None, how=None, **kwargs):

    chromosome, strand = parse_grpby_key(kwargs["key"])

    strand = True if strandedness else False
    s = _cluster(scdf, chromosome, strand=strand)
    o = _cluster(ocdf, chromosome, strand=strand)

    return _intersection(s, o, strandedness=strandedness, how=how, **kwargs)


def _overlapping_for_nearest(scdf, ocdf, suffix, **kwargs):

    nearest_df = pd.DataFrame(columns="Chromosome Start End Strand".split())

    scdf2, ocdf2 = _both_dfs(scdf, ocdf, how="first")

    if not ocdf2.empty:
        # only copying data because of the eternal source buffer array is read only problem
        original_idx = scdf.index

        idxs = scdf2.index
        original_idx = scdf.index.copy(deep=True)
        missing_idxs = ~original_idx.isin(idxs)
        missing_overlap = scdf.index[missing_idxs]

        df_to_find_nearest_in = scdf.reindex(missing_overlap)

        odf = ocdf.reindex(ocdf2.index)
        odf.index = idxs
        sdf = scdf.reindex(idxs)

        nearest_df = sdf.join(odf, rsuffix=suffix).drop("Chromosome" + suffix, axis=1)
        nearest_df.insert(nearest_df.shape[1], "Distance", 0)
    else:
        df_to_find_nearest_in = scdf

    return nearest_df, df_to_find_nearest_in


def _next_nonoverlapping(left_ends, right_starts, right_indexes):

    left_ends = left_ends.sort_values()
    right_starts = right_starts.sort_values()
    r_idx, dist = nearest_next_nonoverlapping(left_ends.values - 1, right_starts.values, right_indexes)
    r_idx = pd.Series(r_idx, index=left_ends.index).sort_index().values
    dist = pd.Series(dist, index=left_ends.index).sort_index().values

    return r_idx, dist


def _previous_nonoverlapping(left_starts, right_ends):

    left_starts = left_starts.sort_values()
    right_ends = right_ends.sort_values()
    r_idx, dist = nearest_previous_nonoverlapping(left_starts.values, right_ends.values - 1, right_ends.index.values)

    r_idx = pd.Series(r_idx, index=left_starts.index).sort_index().values
    dist = pd.Series(dist, index=left_starts.index).sort_index().values

    return r_idx, dist

def sort_one_by_one(d, col1, col2):
    """
    Equivalent to pd.sort_values(by=[col1, col2]), but faster.
    """
    d = d.sort_values(by=[col2])
    return d.sort_values(by=[col1], kind='mergesort')


@return_empty_if_one_empty
def _nearest(scdf, ocdf, suffix="_b", how=None, overlap=True, **kwargs):

    if overlap:
        nearest_df, df_to_find_nearest_in = _overlapping_for_nearest(scdf, ocdf, suffix, **kwargs)
    else:
        df_to_find_nearest_in = scdf

    df_to_find_nearest_in = sort_one_by_one(df_to_find_nearest_in, "Start", "End")
    ocdf = sort_one_by_one(ocdf, "Start", "End")
    df_to_find_nearest_in.index = pd.Index(range(len(df_to_find_nearest_in)))

    if how == "next":
        r_idx, dist = _next_nonoverlapping(df_to_find_nearest_in.End, ocdf.Start, ocdf.index.values)
    elif how == "previous":
        r_idx, dist = _previous_nonoverlapping(df_to_find_nearest_in.Start, ocdf.End)
    else:
        previous_r_idx, previous_dist = _previous_nonoverlapping(df_to_find_nearest_in.Start, ocdf.End)

        next_r_idx, next_dist = _next_nonoverlapping(df_to_find_nearest_in.End, ocdf.Start, ocdf.index.values)

        r_idx, dist = nearest_nonoverlapping(previous_r_idx,
                                             previous_dist,
                                             next_r_idx, next_dist)

    ocdf = ocdf.reindex(r_idx, fill_value=-1) # instead of np.nan, so ints are not promoted to float

    ocdf.index = df_to_find_nearest_in.index
    ocdf.insert(ocdf.shape[1], "Distance", pd.Series(dist, index=ocdf.index).fillna(-1).astype(int))
    ocdf.drop("Chromosome", axis=1, inplace=True)

    r_idx = pd.Series(r_idx, index=ocdf.index)
    df_to_find_nearest_in = df_to_find_nearest_in.drop(r_idx.loc[r_idx == -1].index)

    df = df_to_find_nearest_in.join(ocdf, rsuffix=suffix)

    if overlap and not df.empty and not nearest_df.empty:
        df = pd.concat([nearest_df, df])
    elif overlap and not nearest_df.empty:
        df = nearest_df

    return df


def _set_union(scdf, ocdf, **kwargs):

    chromosome, strand = parse_grpby_key(kwargs["key"])

    strandedness = kwargs["strandedness"]
    strand = True if strandedness == "same" else False

    if len(scdf) == 0:
        return _cluster(ocdf, chromosome, strand=strand)
    elif len(ocdf) == 0:
        return _cluster(scdf, chromosome, strand=strand)

    _starts = np.concatenate([
        scdf.Start.values,
        ocdf.Start.values])
    _ends = np.concatenate([
        scdf.End.values,
        ocdf.End.values])

    cdf = pd.DataFrame({"Start": _starts, "End": _ends})["Start End".split()]
    cdf = cdf.sort_values("Start")
    starts, ends = find_clusters(cdf.Start.values, cdf.End.values)

    chromosome = scdf.head(1)["Chromosome"].iloc[0]
    if strandedness == "same":
        _strand = scdf.head(1)["Strand"].iloc[0]
    else:
        _strand = False

    cluster_df = _create_df_from_starts_ends(starts, ends, chromosome, _strand)

    return cluster_df


def _subtraction(scdf, ocdf, **kwargs):

    chromosome, strand = parse_grpby_key(kwargs["key"])

    if ocdf.empty or scdf.empty:
        return scdf

    strandedness = kwargs["strandedness"]
    strand = True if strandedness else False

    oc = _cluster(ocdf, chromosome, strand)
    o = NCLS(oc.Start.values, oc.End.values, oc.index.values)

    idx_self, new_starts, new_ends = o.set_difference_helper(
        scdf.Start.values,
        scdf.End.values,
        scdf.index.values)

    missing_idx = pd.Index(scdf.index).difference(idx_self)

    idx_to_drop = new_starts != -1

    new_starts = new_starts[idx_to_drop]
    new_ends = new_ends[idx_to_drop]

    idx_self = idx_self[idx_to_drop]
    new_starts = pd.Series(new_starts, index=idx_self).sort_index()
    new_ends = pd.Series(new_ends, index=idx_self).sort_index()
    idx_self = np.sort(idx_self)

    scdf = scdf.reindex(missing_idx.union(idx_self))

    if len(idx_self):

        scdf.loc[scdf.index.isin(idx_self), "Start"] = new_starts
        scdf.loc[scdf.index.isin(idx_self), "End"] = new_ends

    return scdf


def _coverage(ranges, value_col=None, stranded=False, **coverage):

    try:
        from pyranges import PyRles
    except ImportError:
        raise Exception("Using the coverage method requires that pyrle is installed.")

    return PyRles(ranges, value_col=value_col, stranded=stranded)




@return_empty_if_one_empty
def _write_both(scdf, ocdf, new_pos=False, **kwargs):

    suffixes = kwargs["suffixes"]

    ocdf = ocdf.drop("Chromosome", 1)
    scdf, ocdf = _both_dfs(scdf, ocdf, **kwargs)
    nix = pd.Index(range(len(scdf)))
    scdf.index = nix
    ocdf.index = nix

    if not new_pos:
        df = scdf.join(ocdf, rsuffix=suffixes[1])

    elif new_pos == "intersection":

        new_starts = pd.Series(
            np.where(scdf.Start.values > ocdf.Start.values, scdf.Start, ocdf.Start),
            index=scdf.index, dtype=np.long)

        new_ends = pd.Series(
            np.where(scdf.End.values < ocdf.End.values, scdf.End, ocdf.End),
            index=scdf.index, dtype=np.long)
        df = scdf.join(ocdf, lsuffix=suffixes[0], rsuffix=suffixes[1])
        df.insert(1, "Start", new_starts)
        df.insert(2, "End", new_ends)
        df.rename(index=str, columns={"Chromosome" + suffixes[0]: "Chromosome", "Strand" + suffixes[0]: "Strand"}, inplace=True)

    elif new_pos == "union":

        new_starts = pd.Series(
            np.where(scdf.Start.values < ocdf.Start.values, scdf.Start, ocdf.Start),
            index=scdf.index, dtype=np.long)

        new_ends = pd.Series(
            np.where(scdf.End.values > ocdf.End.values, scdf.End, ocdf.End),
            index=scdf.index, dtype=np.long)
        df = scdf.join(ocdf, lsuffix=suffixes[0], rsuffix=suffixes[1])
        df.insert(1, "Start", new_starts)
        df.insert(2, "End", new_ends)
        df.rename(index=str, columns={"Chromosome" + suffixes[0]: "Chromosome", "Strand" + suffixes[0]: "Strand"}, inplace=True)

    return df


if __name__ == "__main__":

    chip_f = "/mnt/scratch/endrebak/pyranges_benchmark/data/download/chip_1000000.bed.gz"
    background_f = "/mnt/scratch/endrebak/pyranges_benchmark/data/download/input_1000000.bed.gz"


    chip = pd.read_table(chip_f, sep="\t", usecols=[0, 1, 2, 5], header=None,
                            names="Chromosome Start End Strand".split(),
                            dtype={"Chromosome": "category", "Strand": "category"} )

    cgr = pr.PyRanges(chip, copy_df=False)

    background = pd.read_table(background_f, sep="\t", usecols=[0, 1, 2, 5],
                                header=None, names="Chromosome Start End Strand".split(),
                                dtype={"Chromosome": "category", "Strand": "category"})

    bgr = pr.PyRanges(background, copy_df=False)


# c = """chr1	3	6	h	0	+
# chr2	4	7	h	0	-"""

# c2 = """chr1	1	2	f	0	+
# chr2	6	7	f	0	-"""

# import pandas as pd
# from io import StringIO

# names = "Chromosome Start End Name Score Strand".split()
# df1 = pd.read_table(StringIO(c), header=None, names=names)
# df2 = pd.read_table(StringIO(c2), header=None, names=names)

# import dask.dataframe as dd

# ddf1 = dd.from_pandas(df1, npartitions=1)
# ddf2 = dd.from_pandas(df2, npartitions=1)

# ddf1.groupby("Chromosome").apply(func_name, ddf2)
