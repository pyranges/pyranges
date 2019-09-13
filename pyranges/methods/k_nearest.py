import numpy as np

from sorted_nearest import nearest_previous_nonoverlapping_all, nearest_next_nonoverlapping_all

try:
    dummy = profile
except NameError:
    profile = lambda x: x

np.random.seed(0)


def nearest_previous_idx(d1, d2, k):

    d1 = d1.sort_values("Start")
    d2 = d2.sort_values("End")

    d1s = d1.Start
    d2e = d2.End

    # ix where d1s could be inserted into d2e to preserve order of d2e
    ix = np.searchsorted(d2e, d1s, side="right") - 1

    valid = ix >= 0
    ix = ix[valid]

    d1s = d1s.iloc[valid]

    lidx, ridx_pos, dist = nearest_previous_nonoverlapping_all(
        d1s.values, d2e.values, d1s.index.values, ix, k)

    ridx = d2e.iloc[ridx_pos].index

    return lidx, ridx, dist


def nearest_next_idx(d1, d2, k):

    d1 = d1.sort_values("End")
    d2 = d2.sort_values("Start")

    d1e = d1.End
    d2s = d2.Start

    # ix where d1e could be inserted into d2s to preserve order of d2s
    ix = np.searchsorted(d2s, d1e, side="left")

    print(ix)
    valid = ix < len(d2s)
    ix = ix[valid]

    d1e = d1e.iloc[valid]

    print(ix)
    print(d1)
    print(d1e)
    print(d2)

    lidx, ridx_pos, dist = nearest_next_nonoverlapping_all(
        d1e.values, d2s.values, d1e.index.values, ix, k)

    ridx = d2s.iloc[ridx_pos].index

    return lidx, ridx, dist


def nearest_previous(d1, d2, kwargs):

    k = kwargs.get("k", 1)
    suffix = kwargs.get("suffix", "_b")

    lidx, ridx, dist = nearest_previous_idx(d1, d2, k)

    d1 = d1.reindex(lidx)
    d2 = d2.reindex(ridx)
    d1.index = range(len(d1))
    d2.index = range(len(d1))
    df = d1.join(d2, rsuffix=suffix)
    df.insert(df.shape[1], "Distance", dist)
    return df


def nearest_next(d1, d2, kwargs):

    k = kwargs.get("k", 1)
    suffix = kwargs.get("suffix", "_b")

    lidx, ridx, dist = nearest_next_idx(d1, d2, k)

    d1 = d1.reindex(lidx)
    d2 = d2.reindex(ridx)
    d1.index = range(len(d1))
    d2.index = range(len(d1))
    df = d1.join(d2, rsuffix=suffix)
    df.insert(df.shape[1], "Distance", dist)
    return df


@profile
def _nearest(d1, d2, kwargs):

    if d1.empty or d2.empty:
        return None

    # how = kwargs["how"]
    # suffix = kwargs["suffix"]
    # k = kwargs["k"]

    # return nearest_previous(d1, d2, kwargs)
    return nearest_next(d1, d2, kwargs)
    # lidx, ridx, dist = nearest_previous_idx(d1, d2, k)

    # if how == "upstream" or how == "downstream":
    #     strand = d1.Strand.iloc[0]
    #     __nearest = nearest_method[how, strand]
    # else:
    #     __nearest = nearest_idx

    # x = __nearest(d1, d2)

    # x = x[x.Dist >= 0]  # negative dist denotes "had no nearest"
    # x.Dist += 1

    # x = x.drop_duplicates()

    # print("x" * 10)
    # print(x)
    # if overlap:
    #     overlapping_ix = overlap_for_nearest(d1, d2)
    #     print("overlap " * 10)
    #     print(overlapping_ix)
    #     x = pd.concat([overlapping_ix, x])

    # d1 = d1.loc[x.D1X]
    # d1.index = range(len(d1))
    # d2 = d2.loc[x.D2X]
    # d2.index = range(len(d2))

    # j = d1.join(d2, rsuffix=suffix)

    # j = _insert_distance(j, x.Dist.values, suffix)

    # return j
