import numpy as np
import pandas as pd
from ncls import NCLS


def _both_indexes(scdf, ocdf, how=False):

    assert (how in "containment first".split() + [False, None]) or isinstance(
        how, int)
    starts = scdf.Start.values
    ends = scdf.End.values
    indexes = scdf.index.values

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

    return _self_indexes, _other_indexes


def _both_dfs(scdf, ocdf, how=False):

    assert how in "containment first".split() + [False, None]

    _self_indexes, _other_indexes = _both_indexes(scdf, ocdf, how)

    scdf = scdf.reindex(_self_indexes)
    ocdf = ocdf.reindex(_other_indexes)

    return scdf, ocdf


def _write_both(scdf, ocdf, kwargs):

    if scdf.empty or ocdf.empty:
        return None

    if not kwargs.get("new_pos"):
        suffix = kwargs.get("suffix", "_b")
    else:
        suffix = kwargs.get("suffixes", "_a _b".split())[1]

    how = kwargs["how"]

    scdf, ocdf = _both_dfs(scdf, ocdf, how=how)
    nix = pd.Index(range(len(scdf)))
    scdf.index = nix
    ocdf.index = nix

    ocdf = ocdf.drop("Chromosome", axis=1)

    df = scdf.join(ocdf, rsuffix=suffix)

    return df
