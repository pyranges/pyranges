import pandas as pd
import numpy as np

import pyranges as pr
from pyranges.multithreaded import pyrange_apply

from pyranges.methods.statistics import _relative_distance


def simes(df, groupby, pcol):

    if isinstance(groupby, str):
        groupby = [groupby]

    sorter = groupby + [pcol]

    sdf = df[sorter].sort_values(sorter)
    g = sdf.groupby(groupby)

    ranks = g.cumcount().values + 1
    size = g.size().values
    size = np.repeat(size, size)
    multiplied = (sdf[pcol].values * size)

    simes = multiplied / ranks

    sdf.insert(sdf.shape[1], "Simes", simes)

    simes = sdf.groupby(groupby).Simes.min().reset_index()

    return simes


def rowbased_spearman(x, y):

    x = np.asarray(x)
    y = np.asarray(y)

    rx = rowbased_rankdata(x)
    ry = rowbased_rankdata(y)

    return rowbased_pearson(rx, ry)

def rowbased_pearson(x, y):

    # Thanks to https://github.com/dengemann

    def ss(a, axis):
        return np.sum(a * a, axis=axis)

    x = np.asarray(x)
    y = np.asarray(y)

    mx = x.mean(axis=-1)
    my = y.mean(axis=-1)

    xm, ym = x - mx[..., None], y - my[..., None]

    r_num = np.add.reduce(xm * ym, axis=-1)
    r_den = np.sqrt(ss(xm, axis=-1) * ss(ym, axis=-1))

    with np.errstate(divide='ignore', invalid="ignore"):

        r = r_num / r_den

    return r


def rowbased_rankdata(data):

    """Row-based rankdata using method=mean"""

    dc = np.asarray(data).copy()
    sorter = np.apply_along_axis(np.argsort, 1, data)

    inv = np.empty(data.shape, np.intp)

    ranks = np.tile(np.arange(data.shape[1]), (len(data), 1))

    np.put_along_axis(inv, sorter, ranks, axis=1)

    dc = np.take_along_axis(dc, sorter, 1)

    res = np.apply_along_axis(lambda r: r[1:] != r[:-1], 1, dc)

    obs = np.column_stack([np.ones(len(res), dtype=bool), res])

    dense = np.take_along_axis(np.apply_along_axis(np.cumsum, 1, obs), inv, 1)

    len_r = obs.shape[1]

    nonzero = np.count_nonzero(obs, axis=1)
    obs = pd.DataFrame(obs)
    nonzero = pd.Series(nonzero)
    dense = pd.DataFrame(dense)

    ranks = []
    for _nonzero, nzdf in obs.groupby(nonzero, sort=False):

        nz = np.apply_along_axis(lambda r: np.nonzero(r)[0], 1, nzdf)

        _count = np.column_stack([nz, np.ones(len(nz)) * len_r])
        _dense = dense.reindex(nzdf.index).values

        _result = 0.5 * (np.take_along_axis(_count, _dense, 1) + np.take_along_axis(_count, _dense - 1, 1) + 1)

        result = pd.DataFrame(_result, index=nzdf.index)
        ranks.append(result)

    final = pd.concat(ranks).sort_index(kind="mergesort")

    return final


def fdr(p_vals):

    from scipy.stats import rankdata
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1

    return fdr


def fisher_exact(n1, d1, n2, d2, **kwargs):
    try:
        from fisher import pvalue_npy
    except:
        import sys
        print("fisher needs to be installed to use fisher exact. pip install fisher or conda install -c bioconda fisher.")
        sys.exit(-1)

    pseudocount = kwargs.get("pseudocount", 0)
    fe_type = kwargs.get("alternative", "twosided")

    n1 = np.array(n1, dtype=np.uint)
    n2 = np.array(n2, dtype=np.uint)
    d1 = np.array(d1, dtype=np.uint)
    d2 = np.array(d2, dtype=np.uint)

    left, right, twosided = pvalue_npy(n1, d1, n2, d2)

    if fe_type == "twosided":
        p_vals = twosided
    elif fe_type == "left":
        p_vals = left
    elif fe_type == "right":
        p_vals = right
    else:
        raise Exception("alternative must be twosided, left or right")


    OR = ((n1 + pseudocount) / (d2 + pseudocount)) / ((n2 + pseudocount) / (d1 + pseudocount))

    df = pd.DataFrame({"OR": OR, "P": p_vals})

    return df


def chromsizes_as_int(chromsizes):
        if isinstance(chromsizes, int):
            pass
        elif isinstance(chromsizes, dict):
            chromsizes = sum(chromsizes.values())
        elif isinstance(chromsizes, (pd.DataFrame, pr.PyRanges)):
            chromsizes = chromsizes.End.sum()

        return chromsizes


class StatisticsMethods():

    pr = None

    def __init__(self, pr):

        self.pr = pr


    # def __tetrachoric(self, other, chromsizes, **kwargs):

    #     self = self.pr

    #     chromsizes = chromsizes_as_int(chromsizes)

    #     kwargs["new_pos"] = "intersection"
    #     strand = True if kwargs.get("strandedness") else False

    #     ss = self.merge(strand=strand)
    #     so = other.merge(strand=strand)
    #     a = ss.intersect(so, **kwargs).length
    #     b = ss.subtract(so, **kwargs).length
    #     c = so.subtract(ss, **kwargs).length

    #     m = pr.concat([ss, so]).merge(strand=strand).length

    #     d = chromsizes - m

    #     from math import cos, sqrt

    #     _tetrachoric = cos(180/(1 + sqrt((b * c) / (a * d))))

    #     return _tetrachoric


    def forbes(self, other, chromsizes, **kwargs):

        chromsizes = chromsizes_as_int(chromsizes)

        self = self.pr

        kwargs["sparse"] = {"self": True, "other": True}
        kwargs = pr.pyranges.fill_kwargs(kwargs)
        strand = True if kwargs.get("strandedness") else False
        kwargs["new_pos"] = "intersection"

        reference_length = self.merge(strand=strand).length
        query_length = other.merge(strand=strand).length

        intersection_sum = sum(
            v.sum()
            for v in self.set_intersect(other, **kwargs).lengths(as_dict=True).values())

        forbes = chromsizes * intersection_sum / (reference_length * query_length)

        return forbes


    def jaccard(self, other, **kwargs):

        self = self.pr

        kwargs["sparse"] = {"self": True, "other": True}
        kwargs = pr.pyranges.fill_kwargs(kwargs)
        strand = True if kwargs.get("strandedness") else False

        kwargs["new_pos"] = "intersection"
        intersection_sum = sum(
            v.sum()
            for v in self.set_intersect(other, **kwargs).lengths(as_dict=True).values())

        union_sum = 0
        for gr in [self, other]:
            union_sum += sum(
                v.sum() for v in gr.merge(strand=strand).lengths(as_dict=True).values())

        denominator = (union_sum - intersection_sum)
        if denominator == 0:
            return 1
        else:
            jc = intersection_sum / denominator

        return jc

    def relative_distance(self, other, **kwargs):

        self = self.pr

        kwargs["sparse"] = {"self": True, "other": True}
        kwargs = pr.pyranges.fill_kwargs(kwargs)

        result = pyrange_apply(_relative_distance, self, other, **kwargs)  # pylint: disable=E1132

        result = pd.Series(np.concatenate(list(result.values())))

        not_nan = ~np.isnan(result)
        result.loc[not_nan] = np.floor(result[not_nan] * 100) / 100
        vc = result.value_counts(dropna=False).to_frame().reset_index()
        vc.columns = "reldist count".split()
        vc.insert(vc.shape[1], "total", len(result))
        vc.insert(vc.shape[1], "fraction", vc["count"] / len(result))
        vc = vc.sort_values("reldist", ascending=True)
        vc = vc.reset_index(drop=True)

        return vc

from math import sqrt
def _mcc(tp, fp, tn, fn):

    # https://stackoverflow.com/a/56875660/992687
    x = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    return ((tp * tn) - (fp * fn)) / sqrt(x)


def mcc(grs, genome=None, labels=None, strand=False, verbose=False):
    import sys
    from itertools import combinations_with_replacement, chain

    if labels is None:
        _labels = list(range(len(grs)))
        _labels = combinations_with_replacement(_labels, r=2)
    else:
        assert len(labels) == len(grs)
        _labels = combinations_with_replacement(labels, r=2)

    # remove all non-loc columns before computation
    grs = [gr.merge(strand=strand) for gr in grs]

    if genome is not None:
        try:
            genome_length = int(genome)
        except (TypeError, ValueError):
            genome_length = int(genome.End.sum())

        if verbose:
            # check that genome definition does not have many more
            # chromosomes than datafiles
            gr_cs = set(chain(*[gr.chromosomes for gr in grs]))

            g_cs = set(genome.chromosomes)
            surplus = g_cs - gr_cs
            if len(surplus):
                print("The following chromosomes are in the genome, but not the PyRanges:", ", ".join(surplus), file=sys.stderr)


        if strand:
            def make_stranded(df):
                df = df.copy()
                df2 = df.copy()
                df.insert(df.shape[1], "Strand", "+")
                df2.insert(df2.shape[1], "Strand", "-")
                return pd.concat([df, df2])

            genome = genome.apply(make_stranded)

        genome_length = genome.length
    else:
        if len(grs) == 2:
            print("If you do not have a genome and the number of compared pyranges is two, mcc will not give any true negatives.", file=sys.stderr)
        genome_length = pr.concat(grs).merge().length


    strandedness = "same" if strand else None

    rowdicts = []
    for (lt, lf), (t, f) in zip(_labels, combinations_with_replacement(grs, r=2)):
        if verbose:
            print(lt, lf, file=sys.stderr)

        if lt == lf:

            if not strand:
                tp = t.length
                fn = 0
                tn = genome_length - tp
                fp = 0
                rowdicts.append({"T": lt, "F": lf, "TP": tp, "FP": fp, "TN": tn, "FN": fn, "MCC": 1})
            else:
                for strand in "+ -".split():
                    tp = t[strand].length
                    fn = 0
                    tn = genome_length - tp
                    fp = 0
                    rowdicts.append({"T": lt, "F": lf, "Strand": strand, "TP": tp, "FP": fp, "TN": tn, "FN": fn, "MCC": 1})
            continue

        else:
            j = t.join(f, strandedness=strandedness)
            tp_gr = j.new_position("intersection").merge(strand=strand)
            if strand:
                for strand in "+ -".split():
                    tp = tp_gr[strand].length
                    fp = f[strand].length - tp
                    fn = t[strand].length - tp
                    tn = genome_length - (tp + fp + fn)
                    mcc = _mcc(tp, fp, tn, fn)
                    rowdicts.append({"T": lt, "F": lf, "Strand": strand, "TP": tp, "FP": fp, "TN": tn, "FN": fn, "MCC": mcc})
                    rowdicts.append({"T": lf, "F": lt, "Strand": strand, "TP": tp, "FP": fn, "TN": tn, "FN": fp, "MCC": mcc})
            else:
                tp = tp_gr.length
                fp = f.length - tp
                fn = t.length - tp
                tn = genome_length - (tp + fp + fn)
                mcc = _mcc(tp, fp, tn, fn)

                rowdicts.append({"T": lt, "F": lf, "TP": tp, "FP": fp, "TN": tn, "FN": fn, "MCC": mcc})
                rowdicts.append({"T": lf, "F": lt, "TP": tp, "FP": fn, "TN": tn, "FN": fp, "MCC": mcc})

    df = pd.DataFrame.from_dict(rowdicts).sort_values(["T", "F"])

    return df
