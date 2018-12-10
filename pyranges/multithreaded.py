import numpy as np
import pandas as pd
from ncls import NCLS

import pyranges as pr

from natsort import natsorted

from sorted_nearest import (find_clusters, nearest_previous_nonoverlapping,
                            nearest_next_nonoverlapping, nearest_nonoverlapping, find_clusters)

from collections import defaultdict

from functools import wraps

try:
    import ray
except:
    import pyranges.raymock as ray

import sys

@ray.remote
def _create_df_from_starts_ends(starts, ends, chromosome, strand=None):

    if not cluster_df.empty:
        return cluster_df
    else:
        return None


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
def merge_dfs(df1, df2, kwargs):

    if not df1.empty and not df2.empty:
        return pd.concat([df1, df2], sort=False)

    elif df1.empty and df2.empty:
        # can this happen?
        return None
    elif df1.empty:
        return df2
    else:
        return df1


def _concat(self, other):

    strand = False

    if self.stranded and not other.stranded:
        self = pr.PyRanges(pyrange_apply_single(merge_dfs, self, strand))
    elif not self.stranded and other.stranded:
        other = pr.PyRanges(pyrange_apply_single(merge_dfs, other, strand))

    keys = natsorted(set(self.keys() + other.keys()))

    concatted = []
    for k in keys:
        s = self.dfs.get(k, pd.DataFrame(columns="Chromosome Start End".split()))
        o = other.dfs.get(k, pd.DataFrame(columns="Chromosome Start End".split()))

        result = merge_dfs.remote(s, o, {})
        concatted.append(result)

    results = ray.get(concatted)

    return {k: v for (k, v) in zip(keys, results) if not v is None}


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
            # print("s", s, "os", os)
            # print(self.keys)
            # print(other.keys)

            # cannot do this with set subtraction
            # but need it for ...
            if not (c, os) in other.keys():
                # print(c, os, "not in ", other.keys)
                odf = pd.DataFrame(columns="Chromosome Start End".split())
            else:
                # print(other)
                odf = other[c, os].values()[0]


            # print(c, s)
            # print(df.head())
            # print(odf.head())
            # print(df.dtypes)
            # print(odf.dtypes)
            # print("other", other)
            # print("other[c, os]", other[c, os])
            # odf = other[c, os].values[0]
            result = function.remote(df, odf, kwargs)
            # print("successfully completed")
            # print(result)
            # print(result)
            # print(" --- " * 50)
            # print(result)
            results.append(result)

    else:

        if self.stranded and not other.stranded:

            for (c, s), df in items:

                if not c in other.chromosomes:
                    odf = pr.PyRanges(pd.DataFrame(columns="Chromosome Start End".split()))
                else:
                    odf = other.dfs[c]

                result = function.remote(df, odf, kwargs)
                results.append(result)


        elif not self.stranded and other.stranded:

            for c, df in items:

                if not c in other.chromosomes:
                    odf = pr.PyRanges(pd.DataFrame(columns="Chromosome Start End".split()))
                else:
                    odf1 = other[c, "+"]
                    odf2 = other[c, "-"]
                    # merge strands
                    odf = merge_dfs.remote(odf1, odf2)

                result = function.remote(df, odf, kwargs)
                results.append(result)

        elif self.stranded and other.stranded:

            for (c, s), df in self.items():

                if not c in other.chromosomes:
                    odfs = pr.PyRanges(pd.DataFrame(columns="Chromosome Start End".split()))
                else:
                    odfs = other[c].values()


                if len(odfs) == 2:
                    odf = merge_dfs.remote(*odfs, kwargs)
                elif len(odfs) == 1:
                    odf = odfs[0]
                else:
                    odf = pd.DataFrame(columns="Chromosome Start End".split())


                result = function.remote(df, odf, kwargs)
                results.append(result)

        else:

            for c, df in items:
                if not c in other.chromosomes:
                    odf = pd.DataFrame(columns="Chromosome Start End".split())
                else:
                    odf = other.dfs[c]

                result = function.remote(df, odf, kwargs)
                results.append(result)

    results = ray.get(results)


    return {k: r for k, r in zip(keys, results) if r is not None}



def pyrange_apply_single(function, self, strand, kwargs):

    if strand:
        assert self.stranded, \
            "Can only do stranded operation when PyRange contains strand info"

    if self.stranded and strand:
        grpby_key = ["Chromosome", "Strand"]
    else:
        grpby_key = "Chromosome"

    items = self.items()

    results = []

    if strand:

        for (c, s), df in items:

            kwargs["chromosome"] = c
            _strand = s
            kwargs["strand"] = _strand
            result = function.remote(df, kwargs)
            results.append(result)

        keys = self.keys()
        # raise

    elif not self.stranded:

        keys = []
        for c, df in items:

            kwargs["chromosome"] = c
            result = function.remote(df, kwargs)
            results.append(result)
            keys.append(c)

    else:

        keys = []
        for c in self.chromosomes:

            kwargs["chromosome"] = c

            dfs = self[c]
            # print(dfs.values)
            # print(type(dfs.values))
            # print(len(dfs.values))
            if len(dfs.keys()) == 2:
                df1, df2 = dfs.values()
                # merge strands
                df1 = merge_dfs.remote(df1, df2, kwargs)
            else:
                df1 = dfs.values()[0]
            # print(type( df1 ))
            # print(type( df2 ))
            # print(df1)
            # print(df2)
            result = function.remote(df1,kwargs)
            results.append(result)
            keys.append(c)

    results = ray.get(results)

    return {k: r for k, r in zip(keys, results) if r is not None}


@ray.remote
def _cluster(df, kwargs):

    chromosome, strand = kwargs["chromosome"], kwargs.get("strand", None)

    cdf = df.sort_values("Start")

    # print("ooooooool" * 10)
    # cdf = sort_one_by_one(df, "Start", "End")

    starts, ends = find_clusters(cdf.Start.values, cdf.End.values)

    nidx = pd.Index(range(len(starts)))
    if strand:
        cluster_df = pd.DataFrame({"Chromosome": pd.Series(chromosome, dtype="category", index=nidx),
                                   "Start": starts, "End": ends,
                                   "Strand": pd.Series(strand, dtype="category", index=nidx)})
    else:
        cluster_df = pd.DataFrame({"Chromosome": pd.Series(chromosome, dtype="category", index=nidx), "Start": starts, "End": ends})

    return cluster_df


@ray.remote
def _first_df(scdf, ocdf, kwargs):

    if scdf.empty or ocdf.empty:
        return None

    how = kwargs["how"]

    assert how in "containment first".split() + [False, None]
    starts = scdf.Start.values
    ends = scdf.End.values
    indexes = scdf.index.values

    it = NCLS(ocdf.Start.values, ocdf.End.values, ocdf.index.values)

    if not how:
        _indexes = it.has_overlaps(starts, ends, indexes)
    elif how == "containment":
        _indexes = it.has_containment(starts, ends, indexes)

    return scdf.reindex(_indexes)

@ray.remote
def _overlap(scdf, ocdf, kwargs):

    if scdf.empty or ocdf.empty:
        return None

    how = kwargs["how"]

    assert how in "containment first".split() + [False, None]
    starts = scdf.Start.values
    ends = scdf.End.values
    indexes = scdf.index.values

    it = NCLS(ocdf.Start.values, ocdf.End.values, ocdf.index.values)

    _indexes = it.all_overlaps_self(starts, ends, indexes)

    return scdf.reindex(_indexes)


def _both_dfs(scdf, ocdf, how=False):

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

    return scdf.reindex(_self_indexes), ocdf.reindex(_other_indexes)

@ray.remote
def _intersection(scdf, ocdf, kwargs):

    how = kwargs["how"]

    if ocdf.empty or scdf.empty: # just return empty df
        return None
    # scdf2, ocdfe#2 = ray.get(_both_dfs.remote(scdf, ocdf, how=False))
    assert how in "containment first".split() + [False, None]
    starts = scdf.Start.values
    ends = scdf.End.values
    indexes = scdf.index.values

    in_dtype = ocdf.Start.dtype

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
        index=scdf.index, dtype=in_dtype)

    new_ends = pd.Series(
        np.where(scdf.End.values < ocdf.End.values, scdf.End, ocdf.End),
        index=scdf.index, dtype=in_dtype)

    pd.options.mode.chained_assignment = None  # default='warn'
    scdf.loc[:, "Start"] = new_starts
    scdf.loc[:, "End"] = new_ends
    pd.options.mode.chained_assignment = 'warn'

    if not scdf.empty:
        return scdf
    else:
        return None





# @ray.remote
# def _set_intersection(scdf, ocdf, kwargs):

#     if len(scdf) == 0 or len(ocdf) == 0:
#         return None

#     chromosome = scdf.Chromosome.iloc[0]

#     strand = True if kwargs["strandedness"] else False
#     if strand:
#         strand = scdf.Strand.iloc[0]

#     s = _cluster.remote(scdf, chromosome, strand=strand)
#     o = _cluster.remote(ocdf, chromosome, strand=strand)

#     return _intersection.remote(s, o, kwargs)



def sort_one_by_one(d, col1, col2):
    """
    Equivalent to pd.sort_values(by=[col1, col2]), but faster.
    """
    d = d.sort_values(by=[col2])
    return d.sort_values(by=[col1], kind='mergesort')

def _overlapping_for_nearest(scdf, ocdf, suffix):

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



@ray.remote
def _nearest(scdf, ocdf, kwargs):


    if scdf.empty or ocdf.empty:
        return None

    overlap = kwargs["overlap"]
    how = kwargs["how"]
    suffix = kwargs["suffix"]

    if overlap:
        nearest_df, df_to_find_nearest_in = _overlapping_for_nearest(scdf, ocdf, suffix)
    else:
        df_to_find_nearest_in = scdf

    if not df_to_find_nearest_in.empty:
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


        # sys.stderr.flush()

        ocdf = ocdf.reindex(r_idx, fill_value=-1) # instead of np.nan, so ints are not promoted to float

        ocdf.index = df_to_find_nearest_in.index
        ocdf.insert(ocdf.shape[1], "Distance", pd.Series(dist, index=ocdf.index).fillna(-1).astype(int))

        r_idx = pd.Series(r_idx, index=ocdf.index)
        df_to_find_nearest_in = df_to_find_nearest_in.drop(r_idx.loc[r_idx == -1].index)

        df = df_to_find_nearest_in.join(ocdf, rsuffix=suffix)

    if overlap and "df" in locals() and not df.empty and not nearest_df.empty:
        df = pd.concat([nearest_df, df])
    elif overlap and not nearest_df.empty:
        df = nearest_df

    df = df.drop("Chromosome" + suffix, axis=1)
    return df



# @ray.remote
# def _nearest(scdf, ocdf, kwargs):

#     if scdf.empty or ocdf.empty:
#         return None

#     # suffix="_b", how=None, overlap=True
#     overlap = kwargs["overlap"]
#     how = kwargs["how"]
#     suffix = kwargs["suffix"]
#     strandedness = kwargs["strandedness"]

#     if overlap:
#         nearest_df, df_to_find_nearest_in = _overlapping_for_nearest(scdf, ocdf, strandedness, suffix)
#     else:
#         df_to_find_nearest_in = scdf


#     if not df_to_find_nearest_in.empty:
#         df_to_find_nearest_in = sort_one_by_one(df_to_find_nearest_in, "Start", "End")
#         ocdf = sort_one_by_one(ocdf, "Start", "End")
#         df_to_find_nearest_in.index = pd.Index(range(len(df_to_find_nearest_in)))

#         if how == "next":
#             l_idx, r_idx, dist = _next_nonoverlapping(df_to_find_nearest_in.End, ocdf.Start)
#         elif how == "previous":
#             l_idx, r_idx, dist = _previous_nonoverlapping(df_to_find_nearest_in.Start, ocdf.End)
#         else:
#             previous_l_idx, previous_r_idx, previous_dist = _previous_nonoverlapping(df_to_find_nearest_in.Start, ocdf.End)

#             next_l_idx, next_r_idx, next_dist = _next_nonoverlapping(df_to_find_nearest_in.End, ocdf.Start)

#             l_idx, r_idx, dist = nearest_nonoverlapping(previous_l_idx, previous_r_idx,
#                                                 previous_dist,
#                                                 next_l_idx, next_r_idx, next_dist)

#         ocdf = ocdf.reindex(r_idx, fill_value=-1) # instead of np.nan, so ints are not promoted to float

#         ocdf.index = df_to_find_nearest_in.index
#         ocdf.insert(ocdf.shape[1], "Distance", pd.Series(dist, index=ocdf.index).fillna(-1).astype(int))

#         r_idx = pd.Series(r_idx, index=ocdf.index)
#         df_to_find_nearest_in = df_to_find_nearest_in.drop(r_idx.loc[r_idx == -1].index)

#         df = df_to_find_nearest_in.join(ocdf, rsuffix=suffix)

#     if overlap and "df" in locals() and not df.empty and not nearest_df.empty:
#         df = pd.concat([nearest_df, df])
#     elif overlap and not nearest_df.empty:
#         df = nearest_df

#     # df = df.drop("Chromosome" + suffix, axis=1)
#     return df




# def _next_nonoverlapping(left_ends, right_starts):

#     # print("le", left_ends)
#     # print("rs", right_starts)

#     left_ends = left_ends.sort_values()
#     right_starts = right_starts.sort_values()
#     lidx = left_ends.index.values
#     ridx = right_starts.index.values
#     l_idx, r_idx, dist = nearest_next_nonoverlapping(left_ends.values - 1, lidx, right_starts.values, ridx)
#     r_idx = pd.Series(r_idx, index=left_ends.index).sort_index().values
#     dist = pd.Series(dist, index=left_ends.index).sort_index().values

#     return l_idx, r_idx, dist


# def _previous_nonoverlapping(left_starts, right_ends):

#     left_starts = left_starts.sort_values()
#     right_ends = right_ends.sort_values()
#     l_idx, r_idx, dist = nearest_previous_nonoverlapping(left_starts.values, left_starts.index.values, right_ends.values - 1, right_ends.index.values)
#     # print("ridx before", r_idx)
#     r_idx = pd.Series(r_idx, index=left_starts.index).sort_index().values
#     dist = pd.Series(dist, index=left_starts.index).sort_index().values
#     # print("ridx after", r_idx)

#     return l_idx, r_idx, dist

# # from memory_profiler import profile


# def _both_indexes(self, other, how):

#     # print("self", self)
#     # print("other", other)

#     starts = self.Start.values
#     ends = self.End.values
#     indexes = self.index.values

#     it = NCLS(other.Start.values, other.End.values, other.index.values)

#     if not how:
#         _self_indexes, _other_indexes = it.all_overlaps_both(starts, ends, indexes)
#     elif how == "containment":
#         _self_indexes, _other_indexes = it.all_containments_both(starts, ends, indexes)
#     else:
#         _self_indexes, _other_indexes = it.first_overlap_both(starts, ends, indexes)

#     return _self_indexes, _other_indexes


# def _overlapping_for_nearest(self, other, strandedness, suffix):
#     nearest_df = pd.DataFrame(columns="Chromosome Start End Strand".split())

#     sidx, oidx = _both_indexes(self, other, how="first")

#     if not len(oidx) == 0 or not len(sidx) == 0:

#         idxs = sidx
#         missing_overlap = self.index[~self.index.isin(idxs)]
#         # print("missing overlap", missing_overlap)
#         df_to_find_nearest_in = self.reindex(missing_overlap)

#         odf = other.reindex(oidx)
#         odf.index = sidx
#         sdf = self.reindex(sidx)

#         nearest_df = sdf.join(odf, rsuffix=suffix)
#         # print("nearest_df", nearest_df)
#         nearest_df = nearest_df.drop("Chromosome" + suffix, axis=1)
#         # print("nearest_df", nearest_df)
#         nearest_df.insert(nearest_df.shape[1], "Distance", 0)
#     else:
#         df_to_find_nearest_in = self

#     return nearest_df, df_to_find_nearest_in

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


@ray.remote
def _subtraction(scdf, ocdf, kwargs):

    if ocdf.empty or scdf.empty:
        return scdf

    strandedness = kwargs["strandedness"]
    strand = True if strandedness else False

    chromosome = scdf.Chromosome.head(1).iloc[0]
    kwargs["chromosome"] = chromosome

    if "Strand" in ocdf and strand:
        strand = scdf.Strand.head(1).iloc[0]
        kwargs["strand"] = strand

    oc = ray.get(_cluster.remote(ocdf, kwargs))
    # print("oc")
    # print(oc.head())
    # print(oc.dtypes)
    # print(oc.index.values.dtype)
    o = NCLS(oc.Start.values, oc.End.values, oc.index.values)
    # print(o)

    # print("s")
    # print(scdf.head())
    # print(scdf.dtypes)
    # print(scdf.index.values.dtype)
    idx_self, new_starts, new_ends = o.set_difference_helper(
        scdf.Start.values,
        scdf.End.values,
        scdf.index.values)

    missing_idx = pd.Index(scdf.index).difference(idx_self)

    idx_to_drop = new_starts != -1

    new_starts = new_starts[idx_to_drop]
    new_ends = new_ends[idx_to_drop]

    idx_self = idx_self[idx_to_drop]
    new_starts = pd.Series(new_starts, index=idx_self)#.sort_index()
    new_ends = pd.Series(new_ends, index=idx_self)#.sort_index()
    # idx_self = np.sort(idx_self)

    scdf = scdf.reindex(missing_idx.union(idx_self))

    if len(idx_self):
        scdf.loc[scdf.index.isin(idx_self), "Start"] = new_starts
        scdf.loc[scdf.index.isin(idx_self), "End"] = new_ends

    if not scdf.empty:
        return scdf
    else:
        return None


def _coverage(ranges, value_col=None, strand=True):

    try:
        from pyrle.methods import coverage
        from pyrle import PyRles
    except ImportError:
        raise Exception("Using the coverage method requires that pyrle is installed.")

    kwargs = {"value_col": value_col}
    return PyRles(pyrange_apply_single(coverage, ranges, strand, kwargs))



@ray.remote
def _write_both(scdf, ocdf, kwargs):

    if scdf.empty or ocdf.empty:
        return None

    suffixes = kwargs["suffixes"]
    suffix = kwargs["suffix"]
    how = kwargs["how"]
    new_pos = kwargs["new_pos"]

    scdf, ocdf = _both_dfs(scdf, ocdf, how=how)
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

def _lengths(df):

    lengths = df.End - df.Start

    return lengths

@ray.remote
def _jaccard(self, other, kwargs):

    strandedness = kwargs["strandedness"]

    if strandedness:
        strand = True
    else:
        strand = False

    s = _lengths(self).sum()
    o = _lengths(other).sum()

    res = ray.get(_intersection.remote(self, other, kwargs))
    if isinstance(res, pd.DataFrame):
        if not res.empty:
            il = _lengths(res).sum()
        else:
            il = 0
    else:
        il = 0

    return [ s, o, il ]


@ray.remote
def _tss(df, kwargs):

    slack = kwargs["slack"]

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


@ray.remote
def _tes(df, kwargs):

    slack = kwargs["slack"]

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


@ray.remote
def _slack(df, kwargs):

    slack = kwargs["slack"]
    # df = self.df.copy()
    df.Start = df.Start - slack
    df.loc[df.Start < 0, "Start"] = 0
    df.End = df.End + slack

    return df

@ray.remote
def _sort(df, kwargs):

    return sort_one_by_one(df, "Start", "End")


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
