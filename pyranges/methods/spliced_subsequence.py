#!/usr/bin/env python

import pandas as pd
import numpy as np
#from easyterm import write

def _spliced_subseq(scdf, **kwargs):
    if scdf.empty:
        return None

    scdf = scdf.copy()

    by = kwargs.get("by") if kwargs.get("by") else "__i__"
    #else "__index__"
    if not type(by) is list: 
        by=[by]              
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
        #scdf = scdf.sort_values( by + ["Start", "End"], ascending=False, ignore_index=True) #debug
        scdf = scdf.sort_values(["Start", "End"], ascending=False)         
    #else:   #debug
        #scdf = scdf.sort_values( by + ["Start", "End"], ascending=True, ignore_index=True) #debug        
        
    scdf.insert(scdf.shape[1], "__length__", scdf.End - scdf.Start)
    #scdf.insert(scdf.shape[1], "__index__", np.arange(len(scdf)))
    scdf.insert(scdf.shape[1], "__i__",  scdf.index)
    
    g = scdf.groupby(by)
    scdf.insert(scdf.shape[1], "__cumsum__",  g.__length__.cumsum())

    #minstart_idx = g.__index__.first()
    minstart_idx = g.__i__.first()
    later_dropped=["__i__", "__length__", "__cumsum__"]
    
    #cumsum_lengths = g.__length__.cumsum()  ## remove
    # write('scdf', how='reverse')
    # write(scdf, how='reverse')
    
    # write('minstart_idx', how='green')
    # write(    g.__i__.first(),  how='green')
        
    if start < 0 or (end is not None and end < 0):
        # exon_count is indexed by transcript:
        exon_count = g.__i__.count()
        
        #len_per_transc is total sum of exon length per transcript 
        len_per_transc=scdf.loc[ g.__i__.last(), by+['__cumsum__'] ].rename(
            columns={'__cumsum__':'__totexonlen__'})
        
        # exp_itransc_len has same rows of scdf with total sum of exon length
        # had to add bits to keep the order of rows right, or merge would destroy it
        exp_len_per_transc=(scdf.loc[:,by+['__i__']].merge(len_per_transc, on=by).set_index('__i__').loc[scdf.index]  )
        #write(('exp_len_per_transc\n', exp_len_per_transc), how='magenta') 
        
        if start < 0:
            #actual_start=transc_len + start
            ## BUG below (solved): assumed that exons of the same transcripts are in consecutive rows, can't do that!
            #start = np.repeat(actual_start, exon_count) 
            #start=exp_len_per_transc['__totexonlen__'].values + start
            start=exp_len_per_transc['__totexonlen__'] + start
            
        if end is not None and end < 0:
            #actual_end = transc_len + end
            ## BUG below (solved): assumed that exons of the same transcripts are in consecutive rows, can't do that!
            #end = np.repeat(actual_end, exon_count)
            #end = exp_len_per_transc['__totexonlen__'].values + end
            end = exp_len_per_transc['__totexonlen__'] + end            
            
        # write( ('start\n', start), how='reverse,green') 
        # write( ('end\n', end), how='reverse,magenta')               
                       

    ### Bug below! assumes exons of the same transcripts are in consecutive rows, can't do that!    
    #cs_start = scdf['__cumsum__'].shift(1, fill_value=0) #.values
    cs_start = g.__cumsum__.shift(1, fill_value=0) 
    cs_start.loc[minstart_idx] = 0

    #cs_start[minstart_idx] = 0
    #print ('minstart_idx', minstart_idx)
    #write( 'cs_start', how='green')
    #write( cs_start, how='green')
    
    #cs_start[g.index.first()]=0
    #print ('cs_start\n', cs_start),
    #print ('not val cs_start\n',     scdf['__cumsum__'].shift(1, fill_value=0) )
    
    cs_end = scdf['__cumsum__']
    #write( 'cs_end', how='magenta')
    #write( cs_end, how='magenta')
        
    ## NOTE
    #  here below, start is a scalar if originally provided > 0, or a Series if < 0
    #  same for end
    if strand == "-": # and use_strand:        
        start_adjustments = start - cs_start  
        #write( 'start_adjustments', how='green')
        #write( start_adjustments, how='green')        
        
        adjust_start = start_adjustments > 0
        scdf.loc[adjust_start, "End"] -= start_adjustments[adjust_start].astype(int)

        end_adjustments = cs_end - end
        #write( 'end_adjustments', how='magenta')
        #write( end_adjustments, how='magenta')
        adjust_end = end_adjustments > 0
        scdf.loc[adjust_end, "Start"] += end_adjustments[adjust_end]
    else:
        start_adjustments = start - cs_start
        #write( 'start_adjustments', how='green')
        #write( start_adjustments, how='green')        
        
        adjust_start = start_adjustments > 0
        scdf.loc[adjust_start, "Start"] += start_adjustments[adjust_start].astype(int)

        end_adjustments = cs_end - end
        #write( 'end_adjustments', how='magenta')
        #write( end_adjustments, how='magenta')
        
        adjust_end = end_adjustments > 0
        scdf.loc[adjust_end, "End"] -= end_adjustments[adjust_end] 

    scdf = scdf[(scdf.Start < scdf.End) & (scdf.Start >= 0) & (scdf.End > 0)]

    #x= scdf.drop(later_dropped, axis=1).sort_index()
    #write(scdf.to_csv(sep='\t'), how='reverse,red')
    
    return scdf.drop(later_dropped, axis=1).sort_index()

