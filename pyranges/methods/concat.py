import pyranges as pr
import pandas as pd

from collections import defaultdict


def concat(pyranges, strand=False):

    grs_per_chromosome = defaultdict(list)

    if strand:
        assert all([
            gr.stranded for gr in pyranges
        ]), "Cannot do stranded concat, not all pyranges contain strand info."

    for gr in pyranges:
        for k, df in gr.dfs.items():
            grs_per_chromosome[k].append(df)

    new_pyrange = {}

    for k, v in grs_per_chromosome.items():
        new_pyrange[k] = pd.concat(v)

    return pr.PyRanges(new_pyrange)
