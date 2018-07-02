import numpy as np
import pandas as pd
from ncls import NCLS

import pyranges as pr

from joblib import Parallel, delayed

from sorted_nearest import (find_clusters, nearest_previous_nonoverlapping,
                            nearest_next_nonoverlapping, nearest_nonoverlapping, find_clusters)

from collections import defaultdict

def pyrange_apply(function, self, other, **kwargs):

    strandedness = kwargs["strandedness"]
    n_jobs = kwargs.get("n_jobs", 1)

    print("Using {} cores and strandedness {}".format(n_jobs, strandedness))

    assert strandedness in ["same", "opposite", False, None]

    if strandedness:
        assert self.stranded and other.stranded, \
            "Can only do stranded searches when both PyRanges contain strand info"

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


    outdfs = Parallel(n_jobs=n_jobs)(delayed(function)(scdf, other_dfs[key], **kwargs) for key, scdf in self.df.groupby(grpby_key))
    # else:
    # outdfs = []
    # for key, df in self.df.groupby(grpby_key):
    #     print(key)
    #     print(df)
    #     print(other_dfs[key])
    #     print("--"*100)
    #     print(function)

    #     outdfs.append(function(df, other_dfs[key], **kwargs))
    outdfs = [df for df in outdfs if not df.empty]

    if outdfs:
        df = pd.concat(outdfs)
        return df
    else:
        return pd.DataFrame(columns="Chromosome Start End Strand".split())


def _both_indexes(scdf, ocdf, how=False, **kwargs):

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


def return_empty_if_one_empty(func):

    def extended_func(self, other, **kwargs):

        if len(self) == 0 or len(other) == 0:
            df = pd.DataFrame(columns="Chromosome Start End".split())
        else:
            df = func(self, other, **kwargs)

        return df

    return extended_func


@return_empty_if_one_empty
def _intersection(scdf, ocdf, **kwargs):

    scdf, ocdf = _both_indexes(scdf, ocdf, **kwargs)

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


def _cluster(df, strand=False):

    cdf = df.sort_values("Start")
    starts, ends = find_clusters(cdf.Start.values, cdf.End.values)

    nidx = pd.Index(range(len(starts)))
    if strand:
        chromosome, strand = df.head(1)[["Chromosome", "Strand"]].iloc[0]
        cluster_df = pd.DataFrame({"Chromosome": pd.Series(chromosome, dtype="category", index=nidx),
                                   "Start": starts, "End": ends,
                                   "Strand": pd.Series(strand, dtype="category", index=nidx)})
    else:
        chromosome = df.head(1)["Chromosome"].iloc[0]
        cluster_df = pd.DataFrame({"Chromosome": pd.Series(chromosome, dtype="category", index=nidx), "Start": starts, "End": ends})

    return cluster_df

@return_empty_if_one_empty
def _set_intersection(scdf, ocdf, strandedness=None, how=None, **kwargs):

    strand = True if strandedness else False
    s = _cluster(scdf, strand=strand)
    o = _cluster(ocdf, strand=strand)

    return _intersection(s, o, strandedness=strandedness, how=how, **kwargs)


def _overlapping_for_nearest(scdf, ocdf, suffix, n_jobs=1, **kwargs):

    # scdf, ocdf = scdf.copy(), ocdf.copy()

    nearest_df = pd.DataFrame(columns="Chromosome Start End Strand".split())

    scdf2, ocdf2 = _both_indexes(scdf, ocdf, how="first")

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
    # print("ridx before", r_idx)
    r_idx = pd.Series(r_idx, index=left_starts.index).sort_index().values
    dist = pd.Series(dist, index=left_starts.index).sort_index().values
    # print("ridx after", r_idx)

    return r_idx, dist

@return_empty_if_one_empty
def _nearest(scdf, ocdf, suffix="_b", how=None, overlap=True, **kwargs):

    if overlap:
        nearest_df, df_to_find_nearest_in = _overlapping_for_nearest(scdf, ocdf, suffix, **kwargs)
    else:
        df_to_find_nearest_in = scdf

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

    # result = cgr. bgr, _intersection, strandedness="same")
