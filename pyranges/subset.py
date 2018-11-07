from pyranges.ncls import find_overlaps

def get_slice(self, val):

    # 100:999

    start = val.start or 0
    stop = val.stop or max(self.df.End.max(), start)

    idxes = []
    for it in self.__ncls__.values():
        idxes.extend([r[2] for r in it.find_overlap(start, stop)])

    return self.df.loc[idxes]


def get_string(self, val):

    if val in self.chromosomes:
        if self.stranded:
            return {k: self.dfs[k] for k in self.keys if k[0] == val}
        else:
            return self.dfs[val]

    elif val in "+ -".split():
        return self.df.loc[self.df.Strand == val]
    else:
        raise Exception("Chromosome or Strand '{}' not found in PyRanges.".format(val))



def get_tuple(self, val):

    if len(val) == 2:
        df = get_double(self, val)
    elif len(val) == 3:
        df = get_triple(self, val)

    return df


def get_double(self, val):

    print(val)
    # "chr1", 5:10
    if len(val) == 2 and val[0] in self.chromosomes and isinstance(val[1], slice):
        chromosome, loc = val
        start = loc.start or 0
        df = self.dfs[chromosome].values
        stop = loc.stop or max(df.End.max(), start)

        idxes = [r[2] for r in find_overlaps(df, start, stop)]
        # idxes = [r[2] for r in self.__ncls__[chromosome, "+"].find_overlap(start, stop)] + \
        #         [r[2] for r in self.__ncls__[chromosome, "-"].find_overlap(start, stop)]

        return df

    # "+", 5:10
    if len(val) == 2 and val[0] in "+ -".split() and isinstance(val[1], slice):
        strand, loc = val
        start = loc.start or 0
        stop = loc.stop or max(self.df.loc[self.df.Chromosome == chromosome].End.max(), start)

        dfs = []
        for df in self.values():

            idxes = [r[2] for r in find_overlaps(df, start, stop)]
            dfs.append(df.reindex(idxes))

        return dfs

    # "chr1", "+"
    if len(val) == 2 and val[1] in "+ -".split():

        chromosome, strand = val

        return self.dfs[chromosome, strand]


def get_triple(self, val):

    # "chr1", "+", 5:10
    chromosome, strand, loc = val
    start = loc.start or 0
    stop = loc.stop or max(self.df.loc[(self.df.Chromosome == chromosome) & (self.df.Strand == strand)].End.max(), start)
    if strand not in "+ -".split():
        raise Exception("Strand '{}' invalid.".format(val))

    if (chromosome, strand) in self.__ncls__:
        idxes = [r[2] for r in self.__ncls__[chromosome, strand].find_overlap(start, stop)]
        return self.df.loc[idxes]
    else:
        raise Exception("Chromosome and Strand pair {}, {} not found in PyRanges.".format(chromosome, strand))
