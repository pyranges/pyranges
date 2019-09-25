import pandas as pd
import numpy as np


def _new_position(df, kwargs):

    suffixes = kwargs.get("suffixes", ("_a", "_b"))
    in_dtype = df.Start.dtype
    new_pos = kwargs.get("new_pos")

    start1 = df["Start"]
    start2 = df["Start" + suffixes[1]]
    end1 = df["End"]
    end2 = df["End" + suffixes[1]]

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
            "Start": "Start" + suffixes[0],
            "End": "End" + suffixes[0]
        })

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
            np.where(start1.values < start2.values, start1, start2),
            index=df.index,
            dtype=in_dtype)

        new_ends = pd.Series(
            np.where(end1.values > end2.values, end1, end2),
            index=df.index,
            dtype=in_dtype)

        df = df.rename(columns={
            "Start": "Start" + suffixes[0],
            "End": "End" + suffixes[0]
        })
        df.insert(1, "Start", new_starts)
        df.insert(2, "End", new_ends)
        df.rename(
            index=str,
            columns={
                "Chromosome" + suffixes[0]: "Chromosome",
                "Strand" + suffixes[0]: "Strand"
            },
            inplace=True)
    else:
        raise Exception(
            "Invalid new pos: {}. Use False/None/union/intersection.".format(
                new_pos))

    return df
