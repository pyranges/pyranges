cimport cython

import pandas as pd

from time import time
import datetime

from collections import defaultdict

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef c_overlap(self, other, strandedness, invert):

    df = self.df
    other_strand = {"+": "-", "-": "+"}
    cdef int length, i
    cdef int start, end, idx

    cdef long [::1] starts
    cdef long [::1] ends
    cdef long [::1] indexes

    other_strand = {"+": "-", "-": "+"}
    same_strand = {"+": "+", "-": "-"}

    if strandedness == "same":
        strand_dict = same_strand
    elif strandedness == "opposite":
        strand_dict = other_strand

    if "Strand" in df and "Strand" in other.df and strandedness:

        indexes_of_overlapping = []
        for (chromosome, strand), cdf in df.groupby("Chromosome Strand".split()):

            starts = cdf.Start.values
            ends = cdf.End.values
            indexes = cdf.index.values

            strand = strand_dict[strand]
            it = other.__ncls__[chromosome, strand]

            indexes_of_overlapping.extend(it.has_overlaps(starts, ends, indexes))

    else:

        indexes_of_overlapping = set()
        for chromosome, cdf in df.groupby("Chromosome"):

            starts = cdf.Start.values
            ends = cdf.End.values
            indexes = cdf.index.values

            it = other.__ncls__[chromosome, "+"]
            it2 = other.__ncls__[chromosome, "-"]

            idx1 = it.has_overlaps(starts, ends, indexes)
            idx2 = it2.has_overlaps(starts, ends, indexes)

            indexes_of_overlapping.update(set(idx1).union(idx2))

    if not invert:
        return df.loc[df.index.isin(indexes_of_overlapping)]
    else:
        return df[~df.index.isin(indexes_of_overlapping)]



@cython.boundscheck(False)
@cython.wraparound(False)
cpdef both_indexes(self, other, strandedness):

    assert strandedness in ["same", "opposite", False, None]

    df = self.df

    other_strand = {"+": "-", "-": "+"}
    same_strand = {"+": "+", "-": "-"}

    if strandedness == "same":
        strand_dict = same_strand
    elif strandedness == "opposite":
        strand_dict = other_strand

    self_d = defaultdict(list)
    other_d = defaultdict(list)

    if "Strand" in df and "Strand" in other.df and strandedness:

        for (chromosome, strand), cdf in df.groupby("Chromosome Strand".split()):

            starts = cdf.Start.values
            ends = cdf.End.values
            indexes = cdf.index.values

            strand = strand_dict[strand]
            it = other.__ncls__[chromosome, strand]

            _self_indexes, _other_indexes = it.all_overlaps_both(starts, ends, indexes)
            self_d[chromosome, strand].extend(_self_indexes)
            other_d[chromosome, strand].extend(_other_indexes)

    else:

        for chromosome, cdf in df.groupby("Chromosome"):

            it = other.__ncls__[chromosome, "+"]
            it2 = other.__ncls__[chromosome, "-"]

            starts = cdf.Start.values
            ends = cdf.End.values
            indexes = cdf.index.values

            for _it in [it, it2]:
                _self_indexes, _other_indexes = _it.all_overlaps_both(starts, ends, indexes)
                self_d[chromosome].extend(_self_indexes)
                other_d[chromosome].extend(_other_indexes)

    return self_d, other_d


if __name__ == "__main__":

    print("hello")
