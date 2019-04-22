import numpy as np
import pandas as pd
from ncls import NCLS


def _intersection(scdf, ocdf, kwargs):

    how = kwargs["how"]

    if ocdf.empty or scdf.empty:
        return None

    assert how in "containment first".split() + [False, None]
    starts = scdf.Start.values
    ends = scdf.End.values
    indexes = scdf.index.values

    in_dtype = ocdf.Start.dtype

    oncls = NCLS(ocdf.Start.values, ocdf.End.values, ocdf.index.values)

    if not how or how is None:
        _self_indexes, _other_indexes = oncls.all_overlaps_both(
            starts, ends, indexes)
    elif how == "containment":
        _self_indexes, _other_indexes = oncls.all_containments_both(
            starts, ends, indexes)
    elif how == "first":
        _self_indexes, _other_indexes = oncls.first_overlap_both(
            starts, ends, indexes)

    _self_indexes = _self_indexes
    _other_indexes = _other_indexes

    scdf, ocdf = scdf.reindex(_self_indexes), ocdf.reindex(_other_indexes)

    new_starts = pd.Series(
        np.where(scdf.Start.values > ocdf.Start.values, scdf.Start,
                 ocdf.Start),
        index=scdf.index,
        dtype=in_dtype)

    new_ends = pd.Series(
        np.where(scdf.End.values < ocdf.End.values, scdf.End, ocdf.End),
        index=scdf.index,
        dtype=in_dtype)

    pd.options.mode.chained_assignment = None  # default='warn'
    scdf.loc[:, "Start"] = new_starts
    scdf.loc[:, "End"] = new_ends
    pd.options.mode.chained_assignment = 'warn'

    if not scdf.empty:
        return scdf
    else:
        return None


def _overlap(scdf, ocdf, kwargs):

    invert = kwargs["invert"]

    if scdf.empty or ocdf.empty:
        return None

    how = kwargs["how"]

    assert how in "containment first".split() + [False, None]
    starts = scdf.Start.values
    ends = scdf.End.values
    indexes = scdf.index.values

    it = NCLS(ocdf.Start.values, ocdf.End.values, ocdf.index.values)

    if not how:
        _indexes = it.all_overlaps_self(starts, ends, indexes)
    elif how == "containment":
        _indexes = it.has_containment(starts, ends, indexes)
    else:
        _indexes = it.has_overlaps(starts, ends, indexes)

    # dbg(_indexes)

    if invert:
        _indexes = scdf.index.difference(_indexes)
    # dbg(_indexes)
    # dbg(scdf)

    return scdf.reindex(_indexes)


# def _first_df(scdf, ocdf, kwargs):

#     if scdf.empty or ocdf.empty:
#         return None

#     how = kwargs["how"]

#     assert how in "containment first".split() + [False, None]
#     starts = scdf.Start.values
#     ends = scdf.End.values
#     indexes = scdf.index.values

#     it = NCLS(ocdf.Start.values, ocdf.End.values, ocdf.index.values)

#     if not how:
#         _indexes = it.has_overlaps(starts, ends, indexes)
#     elif how == "containment":
#         _indexes = it.has_containment(starts, ends, indexes)

#     return scdf.reindex(_indexes)
