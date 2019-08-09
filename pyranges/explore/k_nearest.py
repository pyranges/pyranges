import numpy as np
import pandas as pd

np.random.seed(0)

def sort_one_by_one(d, col1, col2):
    """
    Equivalent to pd.sort_values(by=[col1, col2]), but faster.
    """

    d = d.sort_values(by=[col2])
    return d.sort_values(by=[col1], kind='mergesort')

large = False

if large:
    size_chr1 = 248956422
    n = int(1e7)
    s1 = np.random.randint(size_chr1, size=n)
    e1 = s1 + np.random.randint(10000, size=n)
    s2 = np.random.randint(size_chr1, size=n)
    e2 = s2 + np.random.randint(10000, size=n)
else:
    n = 10
    s1 = np.random.randint(100, size=n)
    e1 = s1 + 10
    s2 = np.random.randint(100, size=n*2)
    e2 = s2 + 10


d1 = pd.DataFrame({"Start": s1, "End": e1}).sort_values("Start End".split())
d2 = pd.DataFrame({"Start": s2, "End": e2}).sort_values("Start End".split())


def nearest_next(d1, d2, k=1):
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
        print("k " * 5, k)
        r = np.repeat(np.arange(len(d1e)), k)
        print(r)
        d1e = d1e.iloc[r]

    d2_idx = d2s.iloc[ix].index

    print("len idx" * 5)
    print(len(d1_idx))
    print(len(d2_idx))
    print("len dfs" * 5)
    print(len(d1e))
    print(len(d2s))
    print("ds " * 5)
    print(d1e)
    print(d2s)
    print("se " * 5)
    print(d1e[d1_idx])
    print(d2s[d2_idx])
    dist = d2s[d2_idx].values - d1e[d1_idx].values

    return d1_idx, d2_idx, dist

d1x, d2x, dist = nearest_next(d1, d2, k=2)

print(d1.loc[d1x])
print(d2.loc[d2x])
print(dist)
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
