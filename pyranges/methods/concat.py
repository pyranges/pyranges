from typing import Optional, List

import pandas as pd

import pyranges as pr


def concat(pyranges: List["pr.PyRanges"], strand: Optional[bool] = None) -> "pr.PyRanges":
    non_empty_pyranges = [gr for gr in pyranges if not gr.empty]
    if not non_empty_pyranges:
        return pr.empty()

    consider_strand = all(gr.stranded for gr in non_empty_pyranges) if strand is None else strand

    if consider_strand and not all(gr.stranded for gr in non_empty_pyranges):
        raise ValueError("Cannot do stranded concat, not all pyranges contain strand info.")

    dfs_to_concat = [gr.df.drop(columns="Strand", errors='ignore') if not consider_strand else gr.df
                     for gr in non_empty_pyranges]

    return pr.PyRanges(pd.concat(dfs_to_concat))