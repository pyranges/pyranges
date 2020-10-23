import numpy as np
import pandas as pd
from ncls import NCLS


def _both_indexes(scdf, ocdf, how=False):

    assert (how in "containment first last outer right left".split() + [False, None]) or isinstance(
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
    elif how == "first":
        _self_indexes, _other_indexes = it.first_overlap_both(
            starts, ends, indexes)
    elif how == "last":
        _self_indexes, _other_indexes = it.last_overlap_both(
            starts, ends, indexes)
        six = scdf.index
        oix = ocdf.index
    elif how in ["outer", "left", "right"]:

        _self_indexes, _other_indexes = it.all_overlaps_both(
            starts, ends, indexes)

        missing_in_s = scdf.index.difference(_self_indexes)
        missing_in_o = ocdf.index.difference(_other_indexes)

        filler_s = np.ones(len(missing_in_o), dtype=int) * -1
        filler_o = np.ones(len(missing_in_s), dtype=int) * -1

        if how == "outer":
            _self_indexes = np.concatenate([_self_indexes, missing_in_s, filler_s])
            _other_indexes = np.concatenate([_other_indexes, filler_o, missing_in_o])
        elif how == "left":
            _self_indexes = np.concatenate([_self_indexes, missing_in_s])
            _other_indexes = np.concatenate([_other_indexes, filler_o])
        elif how == "right":
            _self_indexes = np.concatenate([_self_indexes, filler_s])
            _other_indexes = np.concatenate([_other_indexes, missing_in_o])

    return _self_indexes, _other_indexes

def null_types(h):
    h2 = h.copy()
    for n, d in zip(h, h.dtypes):
        if n in ["Chromosome", "Strand"]:
            continue

        d = str(d)

        if "int" in d or "float" in d:
            null = -1
        elif d == "str" or d == "object":
            null = "-1"
        elif d == "category":
            tmp_cat = h2[n].copy()
            tmp_cat = tmp_cat.cat.add_categories("-1")
            h2[n] = tmp_cat
            null = "-1"

        h2.loc[:, n] = null

    return h2


def _both_dfs(scdf, ocdf, how=False):

    _self_indexes, _other_indexes = _both_indexes(scdf, ocdf, how)

    if how in ["outer", "left", "right"]:

        sh = null_types(scdf.head(1))
        oh = null_types(ocdf.head(1))
        sh.index = [-1]
        oh.index = [-1]

        scdf = scdf.append(sh)
        ocdf = ocdf.append(oh)

        scdf = scdf.reindex(_self_indexes)
        ocdf = ocdf.reindex(_other_indexes)

        if "Strand" in scdf and "Strand" in ocdf:

            if how == "left":
                x = ocdf.index.values == -1
                ocdf.loc[x, "Strand"] = scdf[x].Strand.values
            elif how == "right":
                x = scdf.index.values == -1
                scdf.loc[x, "Strand"] = ocdf[x].Strand.values

    else:
        scdf = scdf.reindex(_self_indexes)
        ocdf = ocdf.reindex(_other_indexes)

    return scdf, ocdf


def _write_both(scdf, ocdf, **kwargs):

    how = kwargs["how"]

    if scdf.empty or ocdf.empty:
        if how in ["left", "outer"] and ocdf.empty:
            ocdf = null_types(kwargs["example_header_other"])
        elif how in ["right", "outer"] and scdf.empty:
            scdf = null_types(kwargs["example_header_self"])
        # TODO: add outer
        # elif how == "outer":
        #     if scdf.empty:
        else:
            return None

    if not kwargs.get("new_pos"):
        suffix = kwargs.get("suffix", "_b")
    else:
        suffix = kwargs.get("suffixes", "_a _b".split())[1]

    scdf, ocdf = _both_dfs(scdf, ocdf, how=how)
    nix = pd.Index(range(len(scdf)))
    scdf.index = nix
    ocdf.index = nix

    ocdf = ocdf.drop("Chromosome", axis=1)

    df = scdf.join(ocdf, rsuffix=suffix)

    if kwargs.get("report_overlap"):
        df["Overlap"] = df[["End", "End"+suffix]].min(axis=1) - df[["Start", "Start"+suffix]].max(axis=1)

    return df
