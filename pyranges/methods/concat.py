import pyranges as pr
import pandas as pd

from collections import defaultdict


def concat(pyranges, strand=None):

    if not pyranges:
        return None

    # from pydbg import dbg
    pyranges = [pr for pr in pyranges if not pr.empty]
    # dbg(pyranges)
    # dbg([p.df.dtypes for p in pyranges])
    grs_per_chromosome = defaultdict(list)

    strand_info = [gr.stranded for gr in pyranges]

    if strand is None:
        strand = all(strand_info)

    if strand:
        assert all([
            gr.stranded for gr in pyranges
        ]), "Cannot do stranded concat, not all pyranges contain strand info."

        for gr in pyranges:
            for k, df in gr.dfs.items():
                # dbg(df)
                grs_per_chromosome[k].append(df)
    else:
        for gr in pyranges:
            for chromosome in gr.chromosomes:
                df = gr[chromosome].df
                grs_per_chromosome[chromosome].append(df)

    new_pyrange = {}

    for k, v in grs_per_chromosome.items():
        new_pyrange[k] = pd.concat(v, sort=False)

    res = pr.multithreaded.process_results(new_pyrange.values(),
                                           new_pyrange.keys())

    if any(strand_info) and not all(strand_info):
        new_res = {}
        for k, v in res.items():
            new_res[k] = v.assign(Strand=v.Strand.fillna("."))
        res = pr.PyRanges(new_res)
        res.Strand = res.Strand
    else:
        res = pr.PyRanges(res)

    return  res
