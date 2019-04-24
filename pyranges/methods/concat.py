import pyranges as pr
import pandas as pd

from collections import defaultdict


def concat(pyranges, strand=False):

    # from pydbg import dbg
    pyranges = [pr for pr in pyranges if not pr.empty]
    # dbg(pyranges)
    # dbg([p.df.dtypes for p in pyranges])
    grs_per_chromosome = defaultdict(list)

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
                # dbg(gr)
                # dbg(gr[chromosome])
                df = gr[chromosome].df
                # dbg(df.dtypes)
                grs_per_chromosome[chromosome].append(df)

    new_pyrange = {}

    for k, v in grs_per_chromosome.items():
        # dbg([_v.dtypes for _v in v])
        new_pyrange[k] = pd.concat(v, sort=False)
        # dbg(new_pyrange[k].dtypes)

    res = pr.multithreaded.process_results(new_pyrange.values(),
                                           new_pyrange.keys())

    # dbg([r.dtypes for r in res.values()])

    return pr.PyRanges(res)
