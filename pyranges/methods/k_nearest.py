import numpy as np
import pandas as pd

"""Problems with overlap:

If we find overlapping then how many others do we need? Just fetch all?

Number different lengths rather.

Simplest, then reduce number afterwards :)"""

from sorted_nearest import nearest_previous_nonoverlapping, nearest_next_nonoverlapping

try:
    dummy = profile
except NameError:
    profile = lambda x: x



def nearest_previous_idx(d1, d2, k, ties=None):

    d1 = d1.sort_values("Start")
    d2 = d2.sort_values("End")

    d1s = d1.Start
    d2e = d2.End

    # ix where d1s could be inserted into d2e to preserve order of d2e
    ix = np.searchsorted(d2e, d1s, side="right") - 1

    valid = ix >= 0
    ix = ix[valid]

    d1s = d1s.iloc[valid]

    lidx, ridx_pos, dist = nearest_previous_nonoverlapping(
        d1s.values, d2e.values, d1s.index.values, ix, k, ties)

    ridx = d2e.iloc[ridx_pos].index

    return lidx, ridx, dist


def nearest_next_idx(d1, d2, k, ties=None):

    d1 = d1.sort_values("End")
    d2 = d2.sort_values("Start")

    d1e = d1.End
    d2s = d2.Start

    # ix where d1e could be inserted into d2s to preserve order of d2s
    ix = np.searchsorted(d2s, d1e, side="left")

    #- print(ix)
    valid = ix < len(d2s)
    ix = ix[valid]

    d1e = d1e.iloc[valid]

    lidx, ridx_pos, dist = nearest_next_nonoverlapping(
        d1e.values, d2s.values, d1e.index.values, ix, k, ties)

    ridx = d2s.iloc[ridx_pos].index

    return lidx, ridx, dist


def nearest(d1, d2, kwargs):

    suffix = kwargs.get("suffix", "_b")
    ties = kwargs.get("ties", None)

    # print(f"{ties=}")
    # print("ties", ties)

    plidx, pridx, pdist = nearest_previous_idx(d1, d2, d1.__k__.values, ties)
    nlidx, nridx, ndist = nearest_next_idx(d1, d2, d1.__k__.values, ties)

    pk = d1.__k__.reindex(plidx)
    nk = d1.__k__.reindex(nlidx)

    p = pd.DataFrame({"LX": plidx, "RX": pridx, "D": pdist, "PN": "P", "k": pk})
    n = pd.DataFrame({"LX": nlidx, "RX": nridx, "D": ndist, "PN": "N", "k": nk})

    df = pd.concat([p, n]).sort_values(["LX", "D"])

    k_dfs = []
    k_grpby = df.groupby("k", sort=False)
    if not ties == "different":
        for k, kdf in k_grpby:
            k_dfs.append(kdf.groupby("LX", sort=False).head(k))

        xdf = pd.concat(k_dfs)
    else:
        xdf = df

    # row_indexer = xdf.PN == "P"
    # xdf.loc[row_indexer, "D"] = -xdf[row_indexer].D

    d1 = d1.reindex(xdf.LX)
    d2 = d2.reindex(xdf.RX)
    d1.index = range(len(d1))
    d2.index = range(len(d1))
    d2 = d2.drop("Chromosome", 1)
    df = d1.join(d2, rsuffix=suffix)
    df.insert(df.shape[1], "Distance", xdf.D.values)

    # to_drop = [c for c in df.columns if "__k__" in c]

    # df = df.drop(to_drop, axis=1)

    return df



def nearest_previous(d1, d2, kwargs):

    suffix = kwargs.get("suffix", "_b")
    ties = kwargs.get("ties", None)

    lidx, ridx, dist = nearest_previous_idx(d1, d2, d1.__k__.values, ties)

    d1 = d1.reindex(lidx)
    d2 = d2.reindex(ridx)
    d1.index = range(len(d1))
    d2.index = range(len(d1))

    d2 = d2.drop("Chromosome", 1)
    df = d1.join(d2, rsuffix=suffix)
    df.insert(df.shape[1], "Distance", -dist)

    # to_drop = [c for c in df.columns if "__k__" in c]
    # df = df.drop(to_drop, axis=1)

    return df


def nearest_next(d1, d2, kwargs):

    suffix = kwargs.get("suffix", "_b")
    ties = kwargs.get("ties", None)

    lidx, ridx, dist = nearest_next_idx(d1, d2, d1.__k__.values, ties)

    d1 = d1.reindex(lidx)
    d2 = d2.reindex(ridx)
    d2 = d2.drop("Chromosome", 1)
    d1.index = range(len(d1))
    d2.index = range(len(d1))
    df = d1.join(d2, rsuffix=suffix)
    df.insert(df.shape[1], "Distance", dist)

    # to_drop = [c for c in df.columns if "__k__" in c]
    # df = df.drop(to_drop, axis=1)

    return df


@profile
def _nearest(d1, d2, kwargs):

    if d1.empty or d2.empty:
        return None

    how = kwargs["how"]

    if how in ["upstream", "downstream"] and kwargs["stranded"]:
        strand = d1.Strand.iloc[0]
        __nearest = {("+", "upstream"): nearest_previous,
                     ("-", "upstream"): nearest_next,
                     ("+", "downstream"): nearest_next,
                     ("-", "downstream"): nearest_previous}[strand, how]
    elif how in ["upstream", "downstream"] and not kwargs["stranded"]:
        __nearest = {"upstream": nearest_previous, "downstream": nearest_next}[how]
    else:
        __nearest = nearest

    df = __nearest(d1, d2, kwargs)

    return df


if __name__ == "__main__":

    import pyranges as pr
    import numpy as np
    np.random.seed(0)
    chrM = pr.data.chromsizes()
    # chrM = chrM[chrM.Chromosome == "chrM"]
    size = int(5e5)
    print(np.log10(size))
    half_size = int(size / 2)
    strand = True

    # print("size")
    # print(chrM)

    gr = pr.random(size, chromsizes=chrM, strand=strand).sort()
    gr2 = pr.random(size, chromsizes=chrM, strand=strand).sort()
    gr.ID = np.arange(len(gr))
    gr2.ID = np.arange(len(gr2))

    # print(gr.nearest(gr2))
    from time import time
    start = time()
    ks = np.array([1, 2] * half_size, dtype=int)
    result = gr.k_nearest(gr2, k=ks, strandedness=None, overlap=True, ties="different")
    end = time()

    print(end - start)
    print(result)
    # print(result.sort("ID").msp())
    ids = result[result.Distance == 0].ID
    print(result[result.ID.isin(ids)])
