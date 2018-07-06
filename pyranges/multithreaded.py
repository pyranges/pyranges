import numpy as np
import pandas as pd
from ncls import NCLS

import pyranges as pr

from natsort import natsorted

from joblib import Parallel, delayed

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


def _pyrange_apply(function, scdf, other_dfs, grpby_key, n_jobs=1, **kwargs):

    print("in pyrange apply " * 10, n_jobs)
    if function.__name__ == "_set_union":
        self_dfs =  {k: d for k, d in scdf.groupby(grpby_key)}
        self_dfs = defaultdict(lambda: pd.DataFrame(columns="Chromosome Start End".split()), self_dfs)
        keys_union = natsorted(list(set(self_dfs).union(other_dfs)))
        outdfs = Parallel(n_jobs=n_jobs)(delayed(function)(self_dfs[key], other_dfs[key], key=key, n_jobs=n_jobs, **kwargs) for key in keys_union)

    elif not function.__name__ in ["_write_both", "_nearest"]:
        # for most methods, we do not need to know anything about other except start, end
        # i.e. less data that needs to be sent to other processes
        outdfs = Parallel(n_jobs=n_jobs)(delayed(function)(scdf, other_dfs[key][["Start", "End"]], key=key, n_jobs=n_jobs, **kwargs) for key, scdf in natsorted(scdf.groupby(grpby_key)))

    else:
        outdfs = Parallel(n_jobs=n_jobs)(delayed(function)(scdf, other_dfs[key], key=key, n_jobs=n_jobs, **kwargs) for key, scdf in natsorted(scdf.groupby(grpby_key)))

    outdfs = [df for df in outdfs if not df.empty]

    if outdfs:
        df = pd.concat(outdfs)
        return df
    else:
        return pd.DataFrame(columns="Chromosome Start End Strand".split())


def pyrange_apply_single(function, self, **kwargs):

    strand = kwargs["strand"]
    n_jobs = kwargs.get("n_jobs", 1)

    print("Using {} cores and strand {}".format(n_jobs, strand))

    if strand:
        assert self.stranded, \
            "Can only do stranded operation when PyRange contains strand info"

    if self.stranded and strand:
        grpby_key = ["Chromosome", "Strand"]
    else:
        grpby_key = "Chromosome"

    outdfs = Parallel(n_jobs=n_jobs)(delayed(function)(scdf, **kwargs) for key, scdf in natsorted(self.df.groupby(grpby_key)))


    outdfs = [df for df in outdfs if not df.empty]

    if outdfs:
        df = pd.concat(outdfs)
        return df
    else:
        return pd.DataFrame(columns="Chromosome Start End Strand".split())


def pyrange_apply(function, self, other, n_jobs, **kwargs):

    strandedness = kwargs["strandedness"]

    print("Using {} cores and strandedness {}".format(n_jobs, strandedness))

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

    return _pyrange_apply(function, self.df, other_dfs, grpby_key, n_jobs=n_jobs, **kwargs)




@return_empty_if_one_empty
def _first_df(scdf, ocdf, how=False, invert=False, n_jobs=1, **kwargs):

    assert how in "containment first".split() + [False, None]
    starts = scdf.Start.values
    ends = scdf.End.values
    indexes = scdf.index.values

    print("n_jobs " * 10)
    print(n_jobs)

    if n_jobs > 1:
        print("deepcopy")
        scdf = scdf.copy(deep=True)

    it = NCLS(ocdf.Start.values, ocdf.End.values, ocdf.index.values)

    if not how:
        _indexes = it.has_overlaps(starts, ends, indexes)
    elif how == "containment":
        _indexes = it.has_containments(starts, ends, indexes)

    if not invert:
        return scdf.reindex(_indexes)
    else:
        return scdf.loc[~scdf.index.isin(_indexes)]


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


def _overlapping_for_nearest(scdf, ocdf, suffix, n_jobs=1, **kwargs):

    nearest_df = pd.DataFrame(columns="Chromosome Start End Strand".split())

    scdf2, ocdf2 = _both_dfs(scdf, ocdf, how="first")

    if not ocdf2.empty:
        # only copying data because of the eternal source buffer array is read only problem
        if n_jobs > 1:
            scdf = scdf.copy(deep=True)
            original_idx = scdf.index.copy(deep=True)
        else:
            original_idx = scdf.index

        idxs = scdf2.index
        original_idx = scdf.index.copy(deep=True)
        missing_idxs = ~original_idx.isin(idxs)
        missing_overlap = scdf.index[missing_idxs]

        df_to_find_nearest_in = scdf.reindex(missing_overlap)

        odf = ocdf.reindex(ocdf2.index)
        odf.index = idxs
        sdf = scdf.reindex(idxs)

        nearest_df = sdf.join(odf, rsuffix=suffix)
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

    r_idx = pd.Series(r_idx, index=ocdf.index)
    df_to_find_nearest_in = df_to_find_nearest_in.drop(r_idx.loc[r_idx == -1].index)

    df = df_to_find_nearest_in.join(ocdf, rsuffix=suffix)

    if overlap and not df.empty and not nearest_df.empty:
        df = pd.concat([nearest_df, df])
    elif overlap and not nearest_df.empty:
        df = nearest_df

    df = df.drop("Chromosome" + suffix, axis=1)
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


def _coverage(ranges, value_col=None, stranded=False, n_jobs=1, **coverage):

    try:
        from pyranges import PyRles
    except ImportError:
        raise Exception("Using the coverage method requires that pyrle is installed.")

    return PyRles(ranges, value_col=value_col, stranded=stranded, nb_cpu=n_jobs)



@return_empty_if_one_empty
def _write_both(scdf, ocdf, new_pos=False, **kwargs):

    suffixes = kwargs["suffixes"]

    scdf, ocdf = _both_dfs(scdf, ocdf, **kwargs)
    nix = pd.Index(range(len(scdf)))
    scdf.index = nix
    ocdf.index = nix

    ocdf = ocdf.drop("Chromosome", axis=1)

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
