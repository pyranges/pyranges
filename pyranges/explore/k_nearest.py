import numpy as np
import pandas as pd

from pyranges.methods.join import _both_dfs

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


large = True

if large:
    size_chr1 = 248956422
    n = int(1e6)
    s1 = np.random.randint(size_chr1, size=n)
    e1 = s1 + np.random.randint(10000, size=n)
    s2 = np.random.randint(size_chr1, size=n)
    e2 = s2 + np.random.randint(10000, size=n)
else:
    n = 10
    n2 = n * 2
    s1 = np.random.randint(100, size=n2)
    e1 = s1 + 10
    s2 = np.random.randint(100, size=n)
    e2 = s2 + 10




d1 = pd.DataFrame({"Start": s1, "End": e1})
d2 = pd.DataFrame({"Start": s2, "End": e2})




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

def nearest_next(d1, d2, k, suffix):

    x = nearest_next_idx(d1, d2, k)
    d1 = d1.loc[x.D1X]
    d1.index = range(len(d1))
    d2 = d2.loc[x.D2X]
    d2 = _insert_distance(d2, x.Dist, suffix)
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


def nearest_idx(d1, d2, k=1):

    n = nearest_next_idx(d1, d2, k)
    p = nearest_previous_idx(d1, d2, k)

    df = pd.concat([n, p])
    df = df[df.Dist >= 0]
    df = sort_one_by_one(df, "D1X", "Dist")
    df = df.groupby("D1X", sort=False).head(k)

    return df


n = nearest_idx(d1, d2, k=1)

# print(n)
# print(d1.loc[n.D1X].head())
# print(d2.loc[n.D2X].head())
# print(n.Dist.head())


# def _overlapping_for_nearest(scdf, ocdf, suffix):

#     nearest_df = pd.DataFrame(columns="Chromosome Start End Strand".split())

#     scdf2, ocdf2 = _both_dfs(scdf, ocdf, how="first")

#     if not ocdf2.empty:
#         original_idx = scdf.index

#         idxs = scdf2.index
#         original_idx = scdf.index.copy(deep=True)
#         missing_idxs = ~original_idx.isin(idxs)
#         missing_overlap = scdf.index[missing_idxs]

#         df_to_find_nearest_in = scdf.reindex(missing_overlap)

#         odf = ocdf.reindex(ocdf2.index)
#         odf.index = idxs
#         sdf = scdf.reindex(idxs)

#         nearest_df = sdf.join(odf, rsuffix=suffix)
#         nearest_df = _insert_distance(nearest_df, 0, suffix)
#     else:
#         df_to_find_nearest_in = scdf

#     return nearest_df, df_to_find_nearest_in


# def _nearest(scdf, ocdf, kwargs):

#     if scdf.empty or ocdf.empty:
#         return None

#     overlap = kwargs["overlap"]
#     how = kwargs["how"]
#     suffix = kwargs["suffix"]
#     k = kwargs["k"]

#     if how == "upstream":
#         strand = scdf.Strand.iloc[0]
#         how = {"+": "previous", "-": "next"}[strand]
#     elif how == "downstream":
#         strand = scdf.Strand.iloc[0]
#         how = {"+": "next", "-": "previous"}[strand]

#     if overlap:
#         nearest_df, df_to_find_nearest_in = _overlapping_for_nearest(
#             scdf, ocdf, suffix)
#     else:
#         df_to_find_nearest_in = scdf

#     if not df_to_find_nearest_in.empty:

#         if how == "next":
#             df = nearest_next(scdf, df_to_find_nearest_in, k, suffix)

#             print(df)
#             raise

# print(_nearest(d1, d2, {"overlap": False, "suffix": "hooo", "k": 1, "how": "next"}))

# d1x, d2x, dist = nearest_next(d1, d2, k=2)
# d1x, d2x, dist = nearest_previous(d1, d2, k=2)


# print("d1")
# print(d1)
# print("d2")
# print(d2)

# print("d1 right order")
# print(d1.loc[d1x])
# print("d2 right order")
# print(d2.loc[d2x])
# print(dist)





# print(len(d1))
# print(d1)
# print(len(res))
# print(res)




# a1 = np.sort(a1)
# # CPU times: user 9.78 s, sys: 648 ms, total: 10.4 s
# # Wall time: 951 ms
# a2_s = np.sort(a2)
# # CPU times: user 9.82 s, sys: 544 ms, total: 10.4 s
# # Wall time: 956 ms


# d1s = sort_one_by_one(d1, "Start", "End")
# # CPU times: user 48.2 s, sys: 3.88 s, total: 52.1 s
# # Wall time: 4.22 s
# # time d1.sort_values(["Start", "End"])
# # CPU times: user 1min, sys: 3.92 s, total: 1min 4s
# # Wall time: 25.3 s
# d2s = sort_one_by_one(d2, "Start", "End")

# r = np.searchsorted(d1s.Start, d2s.End)
