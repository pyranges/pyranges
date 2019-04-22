import numpy as np
import pandas as pd
from ncls import NCLS


def _both_dfs(scdf, ocdf, how=False):

    assert how in "containment first".split() + [False, None]
    starts = scdf.Start.values
    ends = scdf.End.values
    indexes = scdf.index.values

    ocdf = ocdf.reset_index(drop=True)
    it = NCLS(ocdf.Start.values, ocdf.End.values, ocdf.index.values)

    if not how:
        _self_indexes, _other_indexes = it.all_overlaps_both(
            starts, ends, indexes)
    elif how == "containment":
        _self_indexes, _other_indexes = it.all_containments_both(
            starts, ends, indexes)
    else:
        _self_indexes, _other_indexes = it.first_overlap_both(
            starts, ends, indexes)

    _self_indexes = _self_indexes
    _other_indexes = _other_indexes
    scdf = scdf.reindex(_self_indexes)
    ocdf = ocdf.reindex(_other_indexes)

    return scdf, ocdf


def _write_both(scdf, ocdf, kwargs):

    if scdf.empty or ocdf.empty:
        return None

    suffixes = kwargs["suffixes"]
    suffix = kwargs.get("suffix", "_b")
    how = kwargs["how"]
    new_pos = kwargs["new_pos"]

    scdf, ocdf = _both_dfs(scdf, ocdf, how=how)
    nix = pd.Index(range(len(scdf)))
    scdf.index = nix
    ocdf.index = nix

    ocdf = ocdf.drop("Chromosome", axis=1)

    if not new_pos:
        df = scdf.join(ocdf, rsuffix=suffix)

    elif new_pos == "intersection":

        new_starts = pd.Series(
            np.where(scdf.Start.values > ocdf.Start.values, scdf.Start,
                     ocdf.Start),
            index=scdf.index,
            dtype=np.long)

        new_ends = pd.Series(
            np.where(scdf.End.values < ocdf.End.values, scdf.End, ocdf.End),
            index=scdf.index,
            dtype=np.long)
        df = scdf.join(ocdf, lsuffix=suffixes[0], rsuffix=suffixes[1])
        df.insert(1, "Start", new_starts)
        df.insert(2, "End", new_ends)
        df.rename(
            index=str,
            columns={
                "Chromosome" + suffixes[0]: "Chromosome",
                "Strand" + suffixes[0]: "Strand"
            },
            inplace=True)

    elif new_pos == "union":

        new_starts = pd.Series(
            np.where(scdf.Start.values < ocdf.Start.values, scdf.Start,
                     ocdf.Start),
            index=scdf.index,
            dtype=np.long)

        new_ends = pd.Series(
            np.where(scdf.End.values > ocdf.End.values, scdf.End, ocdf.End),
            index=scdf.index,
            dtype=np.long)
        df = scdf.join(ocdf, lsuffix=suffixes[0], rsuffix=suffixes[1])
        df.insert(1, "Start", new_starts)
        df.insert(2, "End", new_ends)
        df.rename(
            index=str,
            columns={
                "Chromosome" + suffixes[0]: "Chromosome",
                "Strand" + suffixes[0]: "Strand"
            },
            inplace=True)

    return df
