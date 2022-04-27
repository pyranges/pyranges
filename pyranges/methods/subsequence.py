#!/usr/bin/env python

import pandas as pd
import numpy as np

def _subseq(scdf, **kwargs):

    if scdf.empty:
        return None

    scdf = scdf.copy()
    orig_order=scdf.index.copy()
    scdf.insert(0, "__i__",  orig_order)
    
    by = kwargs.get("by") if kwargs.get("by") else "__i__"
    if not type(by) is list:
        by=[by]
    start = kwargs.get("start") if kwargs.get("start") else 0
    end = kwargs.get("end") if kwargs.get("end") else scdf.End.max()

    # no strand is treated the same way as positive strand
    strand = scdf.Strand.iloc[0] if "Strand" in scdf else "+"    

    # creating j which holds the boundaries per group
    # j contains one row per group; columns: Start  End (+ by columns); indexed by __i__    
    agg_dict={"__i__": "first", "Start": "min", "End": "max"}
    for b in by:
        agg_dict[b]="first"

    
    if kwargs.get("by"):
        j=scdf.groupby(by)[ ["Start", "End", "__i__"] + by ].agg(agg_dict).set_index("__i__")
    else:
        j=scdf.groupby(by)[ ["Start", "End", "__i__"]  ].agg(agg_dict).set_index("__i__")
        j.insert(0, "__i__",  j.index)
        j.index.name = None
    
    # below we add columns starts, ends to j
    # start and ends define the desired start and end, with one row per group (transcript).
    # they may be out of bounds of the interval, though (this is dealt with later)    
    
    if (strand == "+" and start >= 0) or (strand == "-" and start < 0):
        j['starts'] = j.Start + abs(start)
    else:
        j['starts'] = j.End - abs(start)

    if (strand == "+" and end >= 0) or (strand == "-" and end < 0):
        j['ends'] = j.Start + abs(end)
    else:
        j['ends'] = j.End - abs(end)
        
    if strand == "-":
        j.rename(columns={'ends':'__min__', 'starts':'__max__'}, inplace=True)
    else:
        j.rename(columns={'starts':'__min__', 'ends':'__max__'}, inplace=True)

    # I'm maintaing the original row order    
    scdf=scdf.merge(j[by+['__min__', '__max__']], on=by).set_index('__i__').loc[orig_order]
    
    # instead of simply using starts and ends as computed above, we're dealing here with potential out of bounds:
    r = scdf[~((scdf.Start >= scdf.__max__) | (scdf.End <= scdf.__min__))].copy()
    r.loc[:, "Start"] = np.maximum(r.Start, r.__min__)
    r.loc[:, "End"] = np.minimum(r.End, r.__max__)

    r = r.drop(["__min__", "__max__"], axis=1)

    return r


def _old_corrected_subseq2(scdf, **kwargs):
    """ not used: _subseq above adopts a min/max strategy rather than double sorting as here, and it is faster"""
    
    # from easyterm import write
    if scdf.empty:
        return None

    scdf = scdf.copy()
    #scdf.insert(0, "__index__", np.arange(len(scdf)))
    orig_order=scdf.index.copy()
    scdf.insert(0, "__i__",  orig_order)
    
    #by = kwargs.get("by") if kwargs.get("by") else "__index__"
    by = kwargs.get("by") if kwargs.get("by") else "__i__"
    if not type(by) is list:
        by=[by]
    start = kwargs.get("start") if kwargs.get("start") else 0
    end = kwargs.get("end") if kwargs.get("end") else scdf.End.max()

    # no strand is treated the same way as positive strand
    strand = scdf.Strand.iloc[0] if "Strand" in scdf else "+"

    ### not sure this is optimal: isn't it faster to just ask for min (start) and max (end) per group? --> it is, switching to max/min strategy
    
    gs = scdf.sort_values("Start").groupby(by)
    #gs = gs[["Start", "__index__"]].agg({"__index__": ["count", "first"], "Start": "first"})
    #gs = gs[["Start", "__i__"]].agg({"__i__": "first", "Start": "first"}).rename(
    #    columns={"__i___first": "__i__", "Start_first": "Start"}).set_index("__i__")
    gs = gs[["Start", "__i__"]+by].first().set_index('__i__') #keeping "by" columns
    #agg({"__i__": "first", "Start": "first"})    .rename(
    #    columns={"__i___first": "__i__", "Start_first": "Start"}).set_index("__i__")
    #write(gs, how='magenta')
          
    #gs.columns = ['_'.join(col) if type(col) is tuple else col for col in gs.columns.values]
    #gs = gs.rename(columns={"__index___count": "counts", "__index___first": "__index__", "Start_first": "Start"}).set_index("__index__")
    
    ge = scdf.sort_values("End").groupby(by)[["End", "__i__"]].agg({"__i__": "first", "End": "last"}).set_index("__i__") #dropping by
    #write(ge, how='red')
    j = gs.join(ge).loc[gs.index]
    # j contains one row per group; columns: counts  Start  End
    #  counts is the number of exons; Start and End are the boundaries of the whole group
    
    if (strand == "+" and start >= 0) or (strand == "-" and start < 0):
        j['starts'] = j.Start + abs(start)
    else:
        j['starts'] = j.End - abs(start)

    if (strand == "+" and end >= 0) or (strand == "-" and end < 0):
        j['ends'] = j.Start + abs(end)
    else:
        j['ends'] = j.End - abs(end)
 
    #write(j, how='red')       
        
    # start and ends define the desired start and end, with one row per group (transcript).
    # they may be out of bounds of the interval, though (this is dealt with later)

    #### Bug below (solved): here it assumed that exons of the same transcripts are in consecutive rows, can't do that!
    # Here below: repeat the desired start and end to have one row per interval
    #starts = np.repeat(starts, j.counts).reset_index(drop=True)
    #ends = np.repeat(ends, j.counts).reset_index(drop=True)  
    
    if strand == "-":
        j.rename(columns={'ends':'__min__', 'starts':'__max__'}, inplace=True)
        #ends, starts = starts, ends
    else:
        j.rename(columns={'starts':'__min__', 'ends':'__max__'}, inplace=True)
        
    scdf=scdf.merge(j[by+['__min__', '__max__']], on=by).set_index('__i__').loc[orig_order]
    #write(scdf, how='reverse')    
    #write(scdf.index)
    
    # instead of simply using starts and ends as computed above, we're dealing here with potential out of bounds:
    r = scdf[~((scdf.Start >= scdf.__max__) | (scdf.End <= scdf.__min__))].copy()
    r.loc[:, "Start"] = np.maximum(r.Start, r.__min__)
    r.loc[:, "End"] = np.minimum(r.End, r.__max__)

    r = r.drop(["__min__", "__max__"], axis=1)

    return r



