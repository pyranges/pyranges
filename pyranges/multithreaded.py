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

import ray

# def parse_grpby_key(grpby_key):

#     if isinstance(grpby_key, str):
#         return grpby_key, False
#     else:
#         return grpby_key[0], grpby_key[1]


def return_empty_if_one_empty(func):

    @wraps(func)
    def extended_func(self, other, **kwargs):

        if len(self) == 0 or len(other) == 0:
            df = pd.DataFrame(columns="Chromosome Start End".split())
        else:
            df = func(self, other, **kwargs)

        return df

    return extended_func


@ray.remote
def merge_strands(df1, df2):

    return ray.put(pd.concat([df1, df2]))



def pyrange_apply(function, self, other, **kwargs):

    strandedness = kwargs["strandedness"]

    other_strand = {"+": "-", "-": "+"}
    same_strand = {"+": "+", "-": "-"}

    if strandedness == "opposite":
        strand_dict = other_strand
    else:
        strand_dict = same_strand



    assert strandedness in ["same", "opposite", False, None]

    if strandedness:
        assert self.stranded and other.stranded, \
            "Can only do stranded operations when both PyRanges contain strand info"

    results = []

    items = natsorted(self.dfs.items())
    keys = natsorted(self.dfs.keys())


    if strandedness:

        for (c, s), df in items:
            os = strand_dict[s]
            odf = other.dfs[c, os]
            result = function.remote(df, odf, kwargs)
            results.append(result)

    else:

        if self.stranded and not other.stranded:

            for (c, s), df in items:
                odf = other.dfs[c]

                result = function.remote(df, odf, kwargs)
                results.append(result)


        elif not self.stranded and other.stranded:

            for c, df in items:

                odf1 = other.dfs[c, "+"]
                odf2 = other.dfs[c, "-"]
                odf = merge_strands.remote(odf1, odf2)

                result = function.remote(df, odf, kwargs)
                results.append(result)

        elif self.stranded and other.stranded:

            for (c, s), df in items:
                odf = other.dfs[c, s]
                result = function.remote(df, odf, kwargs)
                results.append(result)

        else:
            for c, df in items:
                odf = other.dfs[c]
                result = function.remote(df, odf, kwargs)
                results.append(result)

    results = ray.get(results)

    return {k: r[0] for k, r in zip(keys, results) if r is not None}



def pyrange_apply_single(function, self, **kwargs):

    strand = kwargs["strand"]

    if strand:
        assert self.stranded, \
            "Can only do stranded operation when PyRange contains strand info"

    if self.stranded and strand:
        grpby_key = ["Chromosome", "Strand"]
    else:
        grpby_key = "Chromosome"

    items = self.items

    results = []

    if strand:

        for (c, s), df in items:

            _strand = s
            result = function.remote(df, c, _strand)
            results.append(result)

        keys = self.keys
        print(keys)
        # raise

    elif not self.stranded:

        keys = []
        for c, df in items:

            result = function.remote(df, c, strand)
            results.append(result)
            keys.append(c)

    else:

        keys = []
        for c in self.chromosomes:

            df = self[c]
            df = merge_strands.remote(*df.dfs.values())
            df = df["Chromosome Start End".split()]
            result = function.remote(df, c, strand)
            results.append(result)
            keys.append(c)

    results = ray.get(results)

    return {k: r[0] for k, r in zip(keys, results) if r is not None}




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

@ray.remote
def _intersection(scdf, ocdf, kwargs):

    how = kwargs["how"]

    if ocdf.empty: # just return empty df
        return None
    # scdf2, ocdfe#2 = ray.get(_both_dfs.remote(scdf, ocdf, how=False))
    assert how in "containment first".split() + [False, None]
    starts = scdf.Start.values
    ends = scdf.End.values
    indexes = scdf.index.values

    oncls = NCLS(ocdf.Start.values, ocdf.End.values, ocdf.index.values)

    if not how or how is None:
        _self_indexes, _other_indexes = oncls.all_overlaps_both(starts, ends, indexes)
    elif how == "containment":
        _self_indexes, _other_indexes = oncls.all_containments_both(starts, ends, indexes)
    elif how == "both":
        _self_indexes, _other_indexes = oncls.first_overlap_both(starts, ends, indexes)

    _self_indexes = _self_indexes
    _other_indexes = _other_indexes

    scdf, ocdf = scdf.reindex(_self_indexes), ocdf.reindex(_other_indexes)

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

    if not scdf.empty:
        return [ray.put(scdf)]
    else:
        return None



def _create_df_from_starts_ends(starts, ends, chromosome, strand=None):

    nidx = pd.Index(range(len(starts)))
    if strand:
        cluster_df = pd.DataFrame({"Chromosome": pd.Series(chromosome, dtype="category", index=nidx),
                                   "Start": starts, "End": ends,
                                   "Strand": pd.Series(strand, dtype="category", index=nidx)})
    else:
        cluster_df = pd.DataFrame({"Chromosome": pd.Series(chromosome, dtype="category", index=nidx), "Start": starts, "End": ends})

    return cluster_df


@ray.remote
def _cluster(df, chromosome, strand=False):

    cdf = df.sort_values("Start")

    starts, ends = find_clusters(cdf.Start.values, cdf.End.values)

    cluster_df = _create_df_from_starts_ends(starts, ends, chromosome, strand)

    if not cluster_df.empty:
        return [ray.put(cluster_df)]
    else:
        return None


@ray.remote
def _set_intersection(scdf, ocdf, kwargs):

    chromosome = scdf.Chromosome.iloc[0]

    strand = True if kwargs["strandedness"] else False
    if strand:
        strand = scdf.Strand.iloc[0]

    s = _cluster.remote(scdf, chromosome, strand=strand)
    o = _cluster.remote(ocdf, chromosome, strand=strand)

    return _intersection.remote(s, o, kwargs)


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
