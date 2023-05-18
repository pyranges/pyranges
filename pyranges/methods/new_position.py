import numpy as np
import pandas as pd


def _new_position(df, **kwargs):
    start_name1, end_name1, start_name2, end_name2 = kwargs["columns"]

    in_dtype = df.Start.dtype
    new_pos = kwargs.get("new_pos")

    start1 = df[start_name1]
    start2 = df[start_name2]
    end1 = df[end_name1]
    end2 = df[end_name2]

    if new_pos == "swap":
        df = df.copy()

        _from = [start_name1, end_name1, start_name2, end_name2]
        _to = [start_name2, end_name2, start_name1, end_name1]

        df[_from] = df[_to]

        return df

    if new_pos == "intersection":
        new_starts = pd.Series(
            np.where(start1.values > start2.values, start1, start2),
            index=df.index,
            dtype=in_dtype,
        )

        new_ends = pd.Series(
            np.where(end1.values < end2.values, end1, end2),
            index=df.index,
            dtype=in_dtype,
        )

    elif new_pos == "union":
        new_starts = pd.Series(
            np.where(start1.values < start2.values, start1, start2),
            index=df.index,
            dtype=in_dtype,
        )

        new_ends = pd.Series(
            np.where(end1.values > end2.values, end1, end2),
            index=df.index,
            dtype=in_dtype,
        )

    else:
        raise Exception("Invalid new pos: {}. Use False/None/union/intersection.".format(new_pos))

    df.loc[:, "Start"] = new_starts
    df.loc[:, "End"] = new_ends

    return df
