#!/usr/bin/env python

import numpy as np

def _subseq(scdf, **kwargs):

    if scdf.empty:
        return None

    scdf = scdf.copy().reset_index(drop=True)

    by = kwargs.get("by", scdf.index.get_level_values())
    start = kwargs["start"] or 0
    end = kwargs["end"]

    # no strand is treated the same way as positive strand
    strand = scdf.Strand.iloc[0] if "Strand" in scdf else "+"

    if (strand == "+" and start < 0) or (strand == "-" and start >= 0):
        g = scdf.sort_values("End").groupby(by).End
        starts = g.last() - abs(start)
    else:
        g = scdf.sort_values("Start").groupby(by).Start
        starts = g.first() + abs(start)

    if (strand == "-" and end < 0) or (strand == "+" and end >= 0):
        g = scdf.sort_values("End").groupby(by).End
        ends = g.last() - abs(end)
    else:
        g = scdf.sort_values("Start").groupby(by).Start
        ends = g.first() + abs(end)

    counts = g.counts()
    scdf.insert("__min__", scdf.shape[1], np.repeat(starts, counts))
    scdf.insert("__max__", scdf.shape[1], np.repeat(ends, counts))

    r = scdf[~((scdf.Start > scdf.__max__) | (scdf.End < scdf.__min__))]
    r.Start = np.maximum(r.Start, r.__min__)
    r.End = np.minimum(r.End, r.__max__)
    r = r.drop(["__min__", "__max__"], axis=1)

    return scdf
