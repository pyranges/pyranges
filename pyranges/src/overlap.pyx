cimport cython

import pandas as pd

# from libcpp.vector cimport vector

# try:
#     dummy = profile
# except:
#     profile = lambda x: x


# @profile

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

    indexes_of_overlapping = []

    if "Strand" in df and "Strand" in other.df and strandedness == "same":

        for (chromosome, strand), cdf in df.groupby("Chromosome Strand".split()):

            starts = cdf.Start.values
            ends = cdf.End.values
            indexes = cdf.index.values

            opposite_strand = other_strand[strand]
            it = other.__ncls__[chromosome, strand]

            indexes_of_overlapping.extend(it.has_overlaps(starts, ends, indexes))


    elif "Strand" in df and "Strand" in other.df and strandedness == "opposite":

        for (chromosome, strand), cdf in df.groupby("Chromosome Strand".split()):

            starts = cdf.Start.values
            ends = cdf.End.values
            indexes = cdf.index.values

            opposite_strand = other_strand[strand]

            it = other.__ncls__[chromosome, opposite_strand]
            indexes_of_overlapping.extend(it.has_overlaps(starts, ends, indexes))

    else:

        for chromosome, cdf in df.groupby("Chromosome"):

            starts = cdf.Start.values
            ends = cdf.End.values
            indexes = cdf.index.values

            it = other.__ncls__[chromosome, "+"]
            it2 = other.__ncls__[chromosome, "-"]

            idx1 = it.has_overlaps(starts, ends, indexes)
            idx2 = it2.has_overlaps(starts, ends, indexes)

            indexes_of_overlapping.extend(list(set(idx1).union(idx2)))


    if not invert:
        return df.loc[indexes_of_overlapping]
    else:
        return df[~df.index.isin(indexes_of_overlapping)]


if __name__ == "__main__":

    print("hello")
