import numpy as np
import pandas as pd
from ncls import NCLS


def _number_overlapping(scdf, ocdf, **kwargs):

    keep_nonoverlapping = kwargs.get("keep_nonoverlapping", True)
    column_name = kwargs.get("overlap_col", True)

    if scdf.empty:
        return None
    if ocdf.empty:
        if keep_nonoverlapping:
            df = scdf.copy()
            df.insert(df.shape[1], column_name, 0)
            return df
        else:
            return None

    oncls = NCLS(ocdf.Start.values, ocdf.End.values, ocdf.index.values)

    starts = scdf.Start.values
    ends = scdf.End.values
    indexes = scdf.index.values

    _self_indexes, _other_indexes = oncls.all_overlaps_both(
        starts, ends, indexes)

    s = pd.Series(_self_indexes)
    counts_per_read = s.value_counts()[s.unique()].reset_index()
    counts_per_read.columns = ["Index", "Count"]

    df = scdf.copy()

    if keep_nonoverlapping:
        _missing_indexes = np.setdiff1d(scdf.index, _self_indexes)
        missing = pd.DataFrame(data={"Index": _missing_indexes, "Count": 0}, index=_missing_indexes)
        counts_per_read = pd.concat([counts_per_read, missing])
    else:
        df = df.loc[_self_indexes]

    counts_per_read = counts_per_read.set_index("Index")

    df.insert(df.shape[1], column_name, counts_per_read)

    return df



def _coverage(scdf, ocdf, **kwargs):

    fraction_col = kwargs["fraction_col"]

    if scdf.empty:
        return None
    if ocdf.empty:
        df = scdf.copy()
        df.insert(df.shape[1], fraction_col, 0.0)
        return df

    oncls = NCLS(ocdf.Start.values, ocdf.End.values, ocdf.index.values)

    starts = scdf.Start.values
    ends = scdf.End.values
    indexes = scdf.index.values

    _lengths = oncls.coverage(starts, ends, indexes)
    _lengths = _lengths /  (ends - starts)
    _fractions = _lengths
    _fractions = _fractions.astype("float64")
    _fractions = np.nan_to_num(_fractions)

    scdf = scdf.copy()

    scdf.insert(scdf.shape[1], fraction_col, _fractions)

    return scdf
