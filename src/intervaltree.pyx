
import numpy as np

cimport cython

from quicksect import IntervalTree


try:
    dummy = profile
except:
    profile = lambda x: x


@profile
@cython.boundscheck(False)
@cython.wraparound(False)
def create_intervaltree(long [::1] indexes, long [::1] starts, long [::1] ends):

    intervalt = IntervalTree()

    cdef int i = 0
    cdef int start
    cdef int end
    cdef int idx

    # print(indexes)
    # print(range(len(indexes)))
    # print(i)

    for i in range(len(indexes)):
        start = starts[i]
        end = ends[i]
        idx = indexes[i]
        i += 1

        intervalt.add(start, end, (start, end, idx))

    return intervalt
