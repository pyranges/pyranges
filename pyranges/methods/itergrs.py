import pandas as pd

from collections import defaultdict


def itergrs(prs, strand=None, keys=False):

    if strand is None:
        strand = all([gr.stranded for gr in prs])

    if strand is False and any([gr.stranded for gr in prs]):
        prs = [gr.unstrand() for gr in prs]

    grs_per_chromosome = defaultdict(list)
    set_keys = set()
    for gr in prs:
        set_keys.update(gr.dfs.keys())

    empty_dfs = [pd.DataFrame(columns=gr.columns) for gr in prs]
    for gr, empty in zip(prs, empty_dfs):
        for k in set_keys:
            df = gr.dfs.get(k, empty)
            grs_per_chromosome[k].append(df)

    if not keys:
        return iter(grs_per_chromosome.values())
    else:
        return iter(grs_per_chromosome.items())
