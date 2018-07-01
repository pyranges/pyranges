import numpy as np
import pandas as pd
from ncls import NCLS

import pyranges as pr

from joblib import Parallel, delayed

from sorted_nearest import find_clusters

# (nearest_previous_nonoverlapping,
#                             nearest_next_nonoverlapping,
#                             nearest_nonoverlapping, find_clusters)
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
        print("returning outdfs " * 10)
        df = pd.concat(outdfs)
        return df
    else:
        print("returning empty " * 10)
        return pd.DataFrame(columns="Chromosome Start End Strand".split())


def _both_indexes(scdf, ocdf, strandedness=False, how=False, **kwargs):

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

    print("scdf")
    print(scdf)

    return scdf


def _cluster(df, strand=False, maxdist=0, minnb=1):

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

    print("in set intersection with df " * 10)
    print(scdf)
    strand = True if strandedness else False
    s = _cluster(scdf, strand=strand)
    o = _cluster(ocdf, strand=strand)

    print("s " * 10)
    print(s)
    print("o " * 10)
    print(o)

    return _intersection(s, o, strandedness=strandedness, how=how, **kwargs)

# def pyrange_or_df_single(func):

#     def extension(self, **kwargs):
#         df = func(self, **kwargs)

#         if kwargs.get("df_only"):
#             return df

#         return PyRanges(df)

#     return extension

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
