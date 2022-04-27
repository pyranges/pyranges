#!/usr/bin/env python

import pandas as pd
import numpy as np

def _bounds(scdf, **kwargs):

    if scdf.empty:
        return None

    col_order=[c for c in scdf.columns]
    
    # scdf = scdf.copy()
    # orig_order=scdf.index.copy()
    #scdf.insert(0, "__i__",  orig_order)
    
    by = kwargs.get("group_by") 
    if not type(by) is list:
        by=[by]
    
    agg_dict = kwargs.get("agg") if kwargs.get("agg") else {}
    agg_dict.update( {"Start":"min", "End":"max", "Chromosome":"first"} )
    if "Strand" in scdf.columns:
        agg_dict["Strand"]="first"
            
    res=scdf.groupby(by).agg(agg_dict).reset_index()
    res=res.reindex(columns=[c for c in col_order if c in res.columns])
    
    return res
        
    # # start = kwargs.get("start") if kwargs.get("start") else 0
    # # end = kwargs.get("end") if kwargs.get("end") else scdf.End.max()

    # # no strand is treated the same way as positive strand
    # #strand = scdf.Strand.iloc[0] if "Strand" in scdf else "+"    

    # # creating j which holds the boundaries per group
    # # j contains one row per group; columns: Start  End (+ by columns); indexed by __i__    
    # agg_dict={"__i__": "first", "Start": "min", "End": "max"}
    # for b in by:
    #     agg_dict[b]="first"
    # j=scdf.groupby(by)[ ["Start", "End", "__i__"] + by ].agg(agg_dict).set_index("__i__")
    
    # # below we add columns starts, ends to j
    # # start and ends define the desired start and end, with one row per group (transcript).
    # # they may be out of bounds of the interval, though (this is dealt with later)    
    
    # if (strand == "+" and start >= 0) or (strand == "-" and start < 0):
    #     j['starts'] = j.Start + abs(start)
    # else:
    #     j['starts'] = j.End - abs(start)

    # if (strand == "+" and end >= 0) or (strand == "-" and end < 0):
    #     j['ends'] = j.Start + abs(end)
    # else:
    #     j['ends'] = j.End - abs(end)
        
    # if strand == "-":
    #     j.rename(columns={'ends':'__min__', 'starts':'__max__'}, inplace=True)
    # else:
    #     j.rename(columns={'starts':'__min__', 'ends':'__max__'}, inplace=True)

    # # I'm maintaing the original row order
    # scdf=scdf.merge(j[by+['__min__', '__max__']], on=by).set_index('__i__').loc[orig_order]
    
    # # instead of simply using starts and ends as computed above, we're dealing here with potential out of bounds:
    # r = scdf[~((scdf.Start >= scdf.__max__) | (scdf.End <= scdf.__min__))].copy()
    # r.loc[:, "Start"] = np.maximum(r.Start, r.__min__)
    # r.loc[:, "End"] = np.minimum(r.End, r.__max__)

    # r = r.drop(["__min__", "__max__"], axis=1)

    # return r
