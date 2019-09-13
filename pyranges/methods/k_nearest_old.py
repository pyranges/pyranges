import numpy as np
import pandas as pd

from .join import _both_indexes

try:
    dummy = profile
except NameError:
    profile = lambda x: x

np.random.seed(0)


def sort_one_by_one(d, col1, col2):
    """
    Equivalent to pd.sort_values(by=[col1, col2]), but faster.
    """

    d = d.sort_values(by=[col2])
    return d.sort_values(by=[col1], kind='mergesort')


def _insert_distance(ocdf, dist, suffix):

    if "Distance" not in ocdf:
        distance_column_name = "Distance"
    elif "Distance" + suffix not in ocdf:
        distance_column_name = "Distance" + suffix
    else:
        i = 1
        while "Distance" + str(i) in ocdf:
            i += 1
        distance_column_name = "Distance" + str(i)

    ocdf.insert(ocdf.shape[1], distance_column_name,
                pd.Series(dist, index=ocdf.index).fillna(-1).astype(int))

    return ocdf


def nearest_next_idx(d1, d2, k=1):
    """Return k indexes from d1, d2 and their distance.

    d1.End <= d2.Start, i.e. next in forward direction.

    dist negative means invalid, i.e. there were less than k nearest intervals."""

    d1e = d1.End.sort_values()
    d2s = d2.Start.sort_values()

    ix = np.searchsorted(d2s, d1e, side="left")

    if k != 1:
        new_ix = np.ones(len(ix) * (k), dtype=d1e.dtype) * -1
        new_ix[::k] = ix
        for _k in range(1, k):
            _k_m1 = _k - 1
            new_ix[_k::k] = new_ix[_k_m1::k] + 1

        ix = new_ix

    ix[ix >= len(d2s)] = len(d2s) - 1

    d1_idx = d1e.index
    if k != 1:
        r = np.repeat(np.arange(len(d1e)), k)
        d1e = d1e.iloc[r]
        d1_idx = d1e.index

    d2_idx = d2s.iloc[ix].index

    d2s = d2s[d2_idx].values
    dist = d2s - d1e

    return pd.DataFrame({"D1X": d1_idx, "D2X": d2_idx, "Dist": dist})


def nearest_next(d1, d2, kwargs):

    k = kwargs["k"]
    suffix = kwargs["suffix"]

    x = nearest_next_idx(d1, d2, k)
    x = x[x.Dist >= 0]
    d1 = d1.loc[x.D1X]
    d1.index = range(len(d1))
    d2 = d2.loc[x.D2X]
    d2.insert(d2.shape[1], "Distance", x.Dist)
    d2.index = range(len(d2))
    return d1.join(d2, rsuffix=suffix)


def nearest_previous_idx(d1, d2, k=1):

    d1s = d1.Start.sort_values()
    d2e = d2.End.sort_values()

    ix = np.searchsorted(d2e, d1s, side="right") - 1
    ix[ix < 0] = 0

    d2_idx = d2e.iloc[ix].index

    if k != 1:
        new_ix = np.ones(len(ix) * (k), dtype=d1s.dtype) * -1
        new_ix[::k] = ix
        for _k in range(1, k):
            _k_m1 = _k - 1
            new_ix[_k::k] = new_ix[_k_m1::k] - 1

        ix = new_ix

    ix[ix < 0] = 0

    d1_idx = d1s.index
    if k != 1:
        r = np.repeat(np.arange(len(d1s)), k)
        d1s = d1s.iloc[r]
        d1_idx = d1s.index

    d2_idx = d2e.iloc[ix].index

    dist = d1s - d2e[d2_idx].values

    return pd.DataFrame({"D1X": d1_idx, "D2X": d2_idx, "Dist": dist})


def nearest_previous(d1, d2, kwargs):

    suffix = kwargs["suffix"]

    x = nearest_previous_idx(d1, d2, k)
    x = x[x.Dist >= 0]
    d1 = d1.loc[x.D1X]
    d1.index = range(len(d1))
    d2 = d2.loc[x.D2X]
    d2.insert(d2.shape[1], "Distance", x.Dist)
    d2.index = range(len(d2))
    return d1.join(d2, rsuffix=suffix)


def nearest_idx(d1, d2, k=1):

    n = nearest_next_idx(d1, d2, k)
    p = nearest_previous_idx(d1, d2, k)

    df = pd.concat([n, p])
    df = df[df.Dist >= 0]
    df = sort_one_by_one(df, "D1X", "Dist")
    # df = df.groupby("D1X", sort=False).head(k)

    return df


nearest_method = {
    ("upstream", "+"): nearest_previous_idx,
    ("upstream", "-"): nearest_next_idx,
    ("downstream", "+"): nearest_next_idx,
    ("downstream", "-"): nearest_previous_idx
}


@profile
def overlap_for_nearest(scdf, ocdf):

    six, oix = _both_indexes(scdf, ocdf)

    return pd.DataFrame(data={"D1X": six, "D2X": oix, "Dist": 0}, index=six)


@profile
def _nearest(d1, d2, kwargs):

    if d1.empty or d2.empty:
        return None

    how = kwargs["how"]
    # k = kwargs["k"]
    suffix = kwargs["suffix"]

    if how == "upstream" or how == "downstream":
        strand = d1.Strand.iloc[0]
        __nearest = nearest_method[how, strand]
    else:
        __nearest = nearest_idx

    x = __nearest(d1, d2)

    x = x[x.Dist >= 0]  # negative dist denotes "had no nearest"
    x.Dist += 1

    x = x.drop_duplicates()

    print("x" * 10)
    print(x)
    if overlap:
        overlapping_ix = overlap_for_nearest(d1, d2)
        print("overlap " * 10)
        print(overlapping_ix)
        x = pd.concat([overlapping_ix, x])

    d1 = d1.loc[x.D1X]
    d1.index = range(len(d1))
    d2 = d2.loc[x.D2X]
    d2.index = range(len(d2))

    j = d1.join(d2, rsuffix=suffix)

    j = _insert_distance(j, x.Dist.values, suffix)

    return j


# def _k_nearest(d1, d2, kwargs):

# def _nearest(d1, d2, kwargs):

#     how = kwargs["how"]
#     k = kwargs["k"]
#     suffix = kwargs["suffix"]
#     overlap = kwargs["overlap"]

#     overlapping, nonoverlapping = d1.overlap_and_nonoverlap_tuple(d2)
