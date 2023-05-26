import pandas as pd

import pyranges as pr


def concat(pyranges, strand=None):
    if not pyranges:
        return None

    pyranges = [gr for gr in pyranges if not gr.empty]

    strand_info = [gr.stranded for gr in pyranges]

    if strand is None:
        strand = all(strand_info)

    if strand:
        assert all([gr.stranded for gr in pyranges]), "Cannot do stranded concat, not all pyranges contain strand info."

    return pr.PyRanges(pd.concat([gr.df for gr in pyranges]))
