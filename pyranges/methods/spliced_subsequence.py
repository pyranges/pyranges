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

    strand = kwargs.get("strand") #, None)
    # at this point, strand is False if 1. spliced_subsequence was called with strand=False or
    #                                   2. it was called with strand=None and self is not stranded
    # or it can be '+' or '-': it was input as True to spliced_subsequence and passed
    #  to pyrange_apply_single as True, which updates it to '-' or '+' before calling _spliced_subseq
    if strand:
        assert "Strand" in scdf, "Cannot have strand=True on unstranded pyranges!"

    if strand == "-": 

        scdf = scdf.sort_values(["Start", "End"], ascending=False)

    scdf.insert(scdf.shape[1], "__length__", scdf.End - scdf.Start)
    scdf.insert(scdf.shape[1], "__index__", np.arange(len(scdf)))
    g = scdf.groupby(by)

    minstart_idx = g.__index__.first()
    cumsum_lengths = g.__length__.cumsum()
    
    if start < 0 or (end is not None and end < 0):
        exon_count = g.__index__.count()
        # transc_len has one row per group (transcript) with total sum of exon length
        transc_len= cumsum_lengths.iloc[ g.__index__.last() ].values
        
        if start < 0:
            actual_start=transc_len + start 
            start = np.repeat(actual_start, exon_count) 
        if end is not None and end < 0:
            actual_end = transc_len + end
            end = np.repeat(actual_end, exon_count) 
        
    cs_start = cumsum_lengths.shift(1, fill_value=0).values
    cs_start[minstart_idx] = 0
    
    cs_end = cumsum_lengths.values

    ## NOTE
    #  here below, start is a scalar if originally provided > 0, or a Series if < 0
    #  same for end
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

