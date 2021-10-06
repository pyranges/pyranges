#!/usr/bin/env python

import numpy as np

def _subseq(scdf, **kwargs):

    if scdf.empty:
        return None

    scdf = scdf.copy().reset_index(drop=True)

    by = kwargs.get("by", scdf.index.get_level_values())
    start = kwargs["start"]
    end = kwargs["end"]

    strand_should_be_used_and_is_neg = "Strand" in scdf and kwargs.get("strand", None) and scdf.Strand.head(1).iloc[0] == "-"
    if strand_should_be_used_and_is_neg:
        g = scdf.sort_values("End").groupby(by).End
        starts = g.last() - start
        ends = starts + start - end
        counts = g.count()
    else:
        g = scdf.sort_values("Start").groupby(by).Start
        starts = g.first() + start
        ends = g.first() - start + end
        counts = g.count()

    scdf.insert("__min__", scdf.shape[1], np.repeat(starts, counts))
    scdf.insert("__max__", scdf.shape[1], np.repeat(ends, counts))

    r = scdf[~((scdf.Start > scdf.__max__) | (scdf.End < scdf.__min__))]
    r.Start = np.maximum(r.Start, r.__min__)
    r.End = np.minimum(r.End, r.__max__)
    r = r.drop(["__min__", "__max__"], axis=1)

    return scdf
