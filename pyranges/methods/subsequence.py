#!/usr/bin/env python

import pandas as pd
import numpy as np

def _subseq(scdf, **kwargs):

    if scdf.empty:
        return None

    scdf = scdf.copy()
    scdf.insert(0, "__index__", np.arange(len(scdf)))

    by = kwargs.get("by") if kwargs.get("by") else "__index__"
    start = kwargs.get("start") if kwargs.get("start") else 0
    end = kwargs.get("end") if kwargs.get("end") else scdf.End.max()


    # no strand is treated the same way as positive strand
    strand = scdf.Strand.iloc[0] if "Strand" in scdf else "+"

    # print("start", start)
    # print("end", end)
    # print("by:", by)
    gs = scdf.sort_values("Start").groupby(by)
    gs = gs[["Start", "__index__"]].agg({"__index__": ["count", "first"], "Start": "first"})
    gs.columns = ['_'.join(col) if type(col) is tuple else col for col in gs.columns.values]
    gs = gs.rename(columns={"__index___count": "counts", "__index___first": "__index__", "Start_first": "Start"}).set_index("__index__")
    ge = scdf.sort_values("End").groupby(by)[["End", "__index__"]].agg({"__index__": "first", "End": "last"}).set_index("__index__")
    j = gs.join(ge).sort_index()
    if (strand == "+" and start >= 0) or (strand == "-" and start < 0):
        starts = j.Start + abs(start)
        # print("start 1")
    else:
        starts = j.End - abs(start)
        # print("start 2")

    if (strand == "+" and end >= 0) or (strand == "-" and end < 0):
        ends = j.Start + abs(end)
        # print("end 1")
    else:
        ends = j.End - abs(end)
        # print("end 2")

    starts = np.repeat(starts, j.counts).reset_index(drop=True)
    ends = np.repeat(ends, j.counts).reset_index(drop=True)

    if strand == "-":
        ends, starts = starts, ends

    scdf.insert(scdf.shape[1], "__min__", starts)
    scdf.insert(scdf.shape[1], "__max__", ends)

    # print(scdf)
    r = scdf[~((scdf.Start >= scdf.__max__) | (scdf.End <= scdf.__min__))].copy()
    # print(r)
    r.loc[:, "Start"] = np.maximum(r.Start, r.__min__)
    r.loc[:, "End"] = np.minimum(r.End, r.__max__)
    # print(r)
    r = r.drop(["__min__", "__max__", "__index__"], axis=1)

    return r
