import pandas as pd
from ncls import NCLS  # type: ignore


def create_ncls(df):
    return NCLS(df.Start.values, df.End.values, df.index.values)


def find_overlaps(df, start, end):
    n = create_ncls(df)

    idxes = []
    for r in n.find_overlap(start, end):
        idxes.append(r[2])

    return idxes


def get_slice(self, val):
    # 100:999

    d = {}

    for k, df in self.items():
        start = val.start or 0
        stop = val.stop or max(df.End.max(), start)
        idxes = find_overlaps(df, start, stop)
        d[k] = df.reindex(idxes)

    return d


def get_string(self, val):
    if val in self.chromosomes:
        if self.stranded:
            return {k: self.dfs[k] for k in self.keys() if k[0] == val}
        else:
            return {val: self.dfs[val]}

    elif val in "+ -".split():
        return {k: v for k, v in self.items() if k[1] == val}
    else:
        return {}


def get_tuple(self, val):
    if len(val) == 2:
        dfs = get_double(self, val)
    elif len(val) == 3:
        dfs = get_triple(self, val)

    return dfs


def get_double(self, val):
    if len(val) == 2 and val[0] in self.chromosomes and isinstance(val[1], slice):
        chromosome, loc = val
        start = loc.start or 0
        if self.stranded:
            dfs = {k: df for k, df in self.items() if k[0] == chromosome}
            max_end = max([df.End.max() for df in dfs.values()])
        else:
            dfs = {val[0]: self.dfs[val[0]]}
            max_end = list(dfs.values())[0].End.max()

        # in case 1:None
        stop = loc.stop or max(max_end, start)

        dfs2 = {}
        for k, df in dfs.items():
            idxes = find_overlaps(df, start, stop)
            if idxes:
                dfs2[k] = df.loc[idxes]

        return dfs2

    # "+", 5:10
    if len(val) == 2 and val[0] in "+ -".split() and isinstance(val[1], slice):
        strand, loc = val
        start = loc.start or 0

        dfs = {k: df for k, df in self.items() if k[1] == strand}
        max_end = max([df.End.max() for df in dfs.values()])

        stop = loc.stop or max(max_end, start)

        dfs2 = {}
        for k, df in dfs.items():
            idxes = find_overlaps(df, start, stop)
            if idxes:
                dfs2[k] = df.loc[idxes]

        return dfs2

    # "chr1", "+"
    if len(val) == 2 and val[1] in "+ -".split():
        chromosome, strand = val

        if (chromosome, strand) in self.dfs:
            return {(chromosome, strand): self.dfs[chromosome, strand]}
        else:
            return {}


def get_triple(self, val):
    # "chr1", "+", 5:10
    chromosome, strand, loc = val
    start = loc.start or 0

    if strand not in "+ -".split():
        raise Exception("Strand '{}' invalid.".format(val))

    r = self[chromosome, strand].values()
    if len(r):
        df = r[0]
    else:
        df = pd.DataFrame(columns="Chromosome Start End".split())
        return df

    max_end = df.End.max()

    stop = loc.stop or max(max_end, start)

    idxes = find_overlaps(df, start, stop)
    return {(chromosome, strand): df.reindex(idxes)}


def get_booldict(self, df):
    _overlapping = set(self.dfs.keys()).intersection(set(df.keys()))

    new_dfs = {}
    for k in _overlapping:
        new_dfs[k] = self.dfs[k][df[k]]

    return new_dfs
