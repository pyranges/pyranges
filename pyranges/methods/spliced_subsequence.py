#!/usr/bin/env python

import pandas as pd
import numpy as np

def _spliced_subseq(scdf, **kwargs):
    if scdf.empty:
        return None

    scdf = scdf.copy()

    by = kwargs.get("by") if kwargs.get("by") else "__index__"
    start = kwargs.get("start") if kwargs.get("start") else 0
    end = kwargs.get("end") if kwargs.get("end") else scdf.End.max()
    length = end - start

    strand = kwargs.get("strand") #, None)
    # at this point, strand is False if 1. spliced_subsequence was called with strand=False or
    #                                   2. it was called with strand=None and self is not stranded
    # or it can be '+' or '-': it was input as True to spliced_subsequence and passed
    #  to pyrange_apply_single as True, which updates it to '-' or '+' before calling _spliced_subseq
    if strand:
        assert "Strand" in scdf, "Cannot have strand=True on unstranded pyranges!"
    #use_strand = True if strand or (strand is None and "Strand" in scdf) else False

    if strand == "-": #use_strand and 
        scdf = scdf.sort_values(["Start", "End"], ascending=False)

    scdf.insert(scdf.shape[1], "__length__", scdf.End - scdf.Start)
    scdf.insert(scdf.shape[1], "__index__", np.arange(len(scdf)))
    g = scdf.groupby(by)

    minstart_idx = g.__index__.first()
    cumsum_lengths = g.__length__.cumsum()

    cs_start = cumsum_lengths.shift(1, fill_value=0).values
    #scdf.insert(scdf.shape[1], "_cs_start", cs_start)    
    
    cs_start[minstart_idx] = 0
    #scdf.insert(scdf.shape[1], "_cs_start2", cs_start)    
    
    cs_end = cumsum_lengths.values
    #scdf.insert(scdf.shape[1], "__cs_end", cs_end)    

    if strand == "-": # and use_strand:
        start_adjustments = start - cs_start
        adjust_start = start_adjustments > 0
        scdf.loc[adjust_start, "End"] -= start_adjustments[adjust_start].astype(int)

        end_adjustments = cs_end - end
        adjust_end = end_adjustments > 0
        scdf.loc[adjust_end, "Start"] += end_adjustments[adjust_end]
    else:
        start_adjustments = start - cs_start
        adjust_start = start_adjustments > 0
        scdf.loc[adjust_start, "Start"] += start_adjustments[adjust_start].astype(int)

        end_adjustments = cs_end - end
        adjust_end = end_adjustments > 0
        scdf.loc[adjust_end, "End"] -= end_adjustments[adjust_end] 

    scdf = scdf[(scdf.Start < scdf.End) & (scdf.Start >= 0) & (scdf.End > 0)]

    return scdf.drop(["__index__", "__length__"], axis=1).sort_index()

