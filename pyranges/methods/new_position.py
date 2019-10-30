import pandas as pd
import numpy as np


def _new_position(df, kwargs):

    suffix1, suffix2 = kwargs["suffixes"]

    in_dtype = df.Start.dtype
    new_pos = kwargs.get("new_pos")

    if "Start" + suffix1 in df:
        start_name1 = "Start" + suffix1
        end_name1 = "End" + suffix1
    else:
        start_name1 = "Start"
        end_name1 = "End"

    strand_name1 = "Strand"
    strand_name2 = "Strand" + suffix2

    start_name2 = "Start" + suffix2
    end_name2 = "End" + suffix2

    start1 = df[start_name1]
    start2 = df[start_name2]
    end1 = df[end_name1]
    end2 = df[end_name2]

    if new_pos == "intersection":

        new_starts = pd.Series(
            np.where(start1.values > start2.values, start1, start2),
            index=df.index,
            dtype=in_dtype)

        new_ends = pd.Series(
            np.where(end1.values < end2.values, end1, end2),
            index=df.index,
            dtype=in_dtype)
        df = df.rename(columns={
            "Start": "Start" + suffix1,
            "End": "End" + suffix1
        })

        df.insert(1, "Start", new_starts)
        df.insert(2, "End", new_ends)
        df.rename(
            index=str,
            columns={
                "Chromosome" + suffix1: "Chromosome",
                "Strand" + suffix1: "Strand"
            },
            inplace=True)

    elif new_pos == "union":

        new_starts = pd.Series(
            np.where(start1.values < start2.values, start1, start2),
            index=df.index,
            dtype=in_dtype)

        new_ends = pd.Series(
            np.where(end1.values > end2.values, end1, end2),
            index=df.index,
            dtype=in_dtype)

        df = df.rename(columns={
            "Start": "Start" + suffix1,
            "End": "End" + suffix1
        })
        df.insert(1, "Start", new_starts)
        df.insert(2, "End", new_ends)
        df.rename(
            index=str,
            columns={
                "Chromosome" + suffix1: "Chromosome",
                "Strand" + suffix1: "Strand"
            },
            inplace=True)

    elif new_pos == "swap":

        df = df.copy()

        _from = [start_name1, end_name1, start_name2, end_name2]
        _to = [start_name2, end_name2, start_name1, end_name1]

        if "Strand" in df:
            _from.extend([strand_name1, strand_name2])
            _to.extend([strand_name2, strand_name1])

        df[_from] = df[_to]

    else:
        raise Exception(
            "Invalid new pos: {}. Use False/None/union/intersection.".format(
                new_pos))

    return df
