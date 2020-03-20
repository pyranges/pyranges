"""Fast implementation of  bb"""

import pandas as pd
import numpy as np

import pyranges as pr
from pyranges.multithreaded import pyrange_apply

from pyranges.methods.statistics import _relative_distance

__all__ = ["simes", "fisher_exact", "StatisticsMethods", "fdr", "rowbased_rankdata", "rowbased_pearson", "rowbased_spearman", "mcc"]

def simes(df, groupby, pcol, keep_position=False):

    """Apply Simes method for giving dependent events a p-value.

    Parameters
    ----------
    df : pandas.DataFrame

        Data to analyse with Simes.

    groupby : str or list of str

        Features equal in these columns will be merged with Simes.

    pcol : str

        Name of column with p-values.

    keep_position : bool, default False

        Keep columns "Chromosome", "Start", "End" and "Strand" if they exist.

    See Also
    --------

    pr.stats.fdr : correct for multiple testing

    Examples
    --------

    >>> s = '''Chromosome Start End Strand Gene PValue
    ... 1 10 20 + P53 0.0001
    ... 1 20 20 + P53 0.0002
    ... 1 30 20 + P53 0.0003
    ... 2 60 65 - FOX 0.05
    ... 2 70 75 - FOX 0.0000001
    ... 2 80 90 - FOX 0.0000021'''

    >>> gr = pr.from_string(s)
    >>> gr
    +--------------+-----------+-----------+--------------+------------+-------------+
    |   Chromosome |     Start |       End | Strand       | Gene       |      PValue |
    |   (category) |   (int32) |   (int32) | (category)   | (object)   |   (float64) |
    |--------------+-----------+-----------+--------------+------------+-------------|
    |            1 |        10 |        20 | +            | P53        |     0.0001  |
    |            1 |        20 |        20 | +            | P53        |     0.0002  |
    |            1 |        30 |        20 | +            | P53        |     0.0003  |
    |            2 |        60 |        65 | -            | FOX        |     0.05    |
    |            2 |        70 |        75 | -            | FOX        |     1e-07   |
    |            2 |        80 |        90 | -            | FOX        |     2.1e-06 |
    +--------------+-----------+-----------+--------------+------------+-------------+
    Stranded PyRanges object has 6 rows and 6 columns from 2 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.

    >>> simes = pr.stats.simes(gr.df, "Gene", "PValue")
    >>> simes
      Gene         Simes
    0  FOX  3.000000e-07
    1  P53  3.000000e-04

    >>> gr.apply(lambda df:
    ... pr.stats.simes(df, "Gene", "PValue", keep_position=True))
    +--------------+-----------+-----------+-------------+------------+------------+
    |   Chromosome |     Start |       End |       Simes | Strand     | Gene       |
    |     (object) |   (int32) |   (int32) |   (float64) | (object)   | (object)   |
    |--------------+-----------+-----------+-------------+------------+------------|
    |            1 |        10 |        20 |      0.0001 | +          | P53        |
    |            2 |        60 |        90 |      1e-07  | -          | FOX        |
    +--------------+-----------+-----------+-------------+------------+------------+
    Stranded PyRanges object has 2 rows and 6 columns from 2 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    if isinstance(groupby, str):
        groupby = [groupby]

    positions = []
    if "Strand" in df:
        stranded = True

    if keep_position:
        positions += ["Chromosome", "Start", "End"]
        if stranded:
            positions += ["Strand"]

    sorter = groupby + [pcol]

    sdf = df[positions + sorter].sort_values(sorter)
    g = sdf.groupby(positions + groupby)

    ranks = g.cumcount().values + 1
    size = g.size().values
    size = np.repeat(size, size)
    multiplied = (sdf[pcol].values * size)

    simes = multiplied / ranks

    sdf.insert(sdf.shape[1], "Simes", simes)

    if keep_position:

        grpby_dict = {"Chromosome": "first", "Start": "min", "End": "max", "Simes": "min"}

        if stranded:
            grpby_dict["Strand"] = "first"

        simes = sdf.groupby(groupby).agg(grpby_dict).reset_index()
        columns = list(simes.columns)
        columns.append(columns[0])
        del columns[0]
        simes = simes[columns]
    else:
        simes = sdf.groupby(groupby).Simes.min().reset_index()

    return simes


def rowbased_spearman(x, y):

    """Fast row-based Spearman's correlation.

    Parameters
    ----------
    x : matrix-like

        2D numerical matrix. Same size as y.

    y : matrix-like

        2D numerical matrix. Same size as x.

    Returns
    -------
    numpy.array

        Array with same length as input, where values are P-values.

    See Also
    --------

    pyranges.statistics.rowbased_pearson : fast row-based Pearson's correlation.
    pr.stats.fdr : correct for multiple testing

    Examples
    --------

    >>> np.random.seed(0)
    >>> x = np.random.randint(10, size=(10, 10))
    >>> y = np.random.randint(10, size=(10, 10))

    Perform Spearman's correlation pairwise on each row in 10x10 matrixes:

    >>> pr.stats.rowbased_spearman(x, y)
    array([ 0.07523548, -0.24838724,  0.03703774,  0.24194052,  0.04778621,
           -0.23913505,  0.12923138,  0.26840486,  0.13292204, -0.29846295])
    """

    x = np.asarray(x)
    y = np.asarray(y)

    rx = rowbased_rankdata(x)
    ry = rowbased_rankdata(y)

    return rowbased_pearson(rx, ry)


def rowbased_pearson(x, y):

    """Fast row-based Pearson's correlation.

    Parameters
    ----------
    x : matrix-like

        2D numerical matrix. Same size as y.

    y : matrix-like

        2D numerical matrix. Same size as x.

    Returns
    -------
    numpy.array

        Array with same length as input, where values are P-values.

    See Also
    --------

    pyranges.statistics.rowbased_spearman : fast row-based Spearman's correlation.
    pr.stats.fdr : correct for multiple testing

    Examples
    --------

    >>> np.random.seed(0)
    >>> x = np.random.randint(10, size=(10, 10))
    >>> y = np.random.randint(10, size=(10, 10))

    Perform Pearson's correlation pairwise on each row in 10x10 matrixes:

    >>> pr.stats.rowbased_pearson(x, y)
    array([ 0.20349603, -0.01667236, -0.01448763, -0.00442322,  0.06527234,
           -0.36710862,  0.14978726,  0.32360286,  0.17209191, -0.08902829])
    """

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

    """Rank order of entries in each row.

    Same as SciPy rankdata with method=mean.

    Parameters
    ----------
    data : matrix-like

        The data to find order of.

    Returns
    -------
    Pandas.DataFrame

        DataFrame where values are order of data.

    Examples
    --------

    >>> np.random.seed(0)
    >>> x = np.random.randint(10, size=(3, 10))
    >>> x
    array([[5, 0, 3, 3, 7, 9, 3, 5, 2, 4],
           [7, 6, 8, 8, 1, 6, 7, 7, 8, 1],
           [5, 9, 8, 9, 4, 3, 0, 3, 5, 0]])
    >>> pr.stats.rowbased_rankdata(x)
         0    1    2    3    4     5    6    7    8    9
    0  7.5  1.0  4.0  4.0  9.0  10.0  4.0  7.5  2.0  6.0
    1  6.0  3.5  9.0  9.0  1.5   3.5  6.0  6.0  9.0  1.5
    2  6.5  9.5  8.0  9.5  5.0   3.5  1.5  3.5  6.5  1.5
    """

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

    """Adjust p-values with Benjamini-Hochberg.

    Parameters
    ----------
    data : array-like


    Returns
    -------
    Pandas.DataFrame

        DataFrame where values are order of data.

    Examples
    --------

    >>> np.random.seed(0)
    >>> x = np.random.random(10) / 100

    >>> gr = pr.random(10)
    >>> gr.PValue = x
    >>> gr
    +--------------+-----------+-----------+--------------+----------------------+
    | Chromosome   | Start     | End       | Strand       | PValue               |
    | (category)   | (int32)   | (int32)   | (category)   | (float64)            |
    |--------------+-----------+-----------+--------------+----------------------|
    | chr1         | 176601938 | 176602038 | +            | 0.005488135039273248 |
    | chr1         | 155082851 | 155082951 | -            | 0.007151893663724195 |
    | chr2         | 211134424 | 211134524 | -            | 0.006027633760716439 |
    | chr9         | 78826761  | 78826861  | -            | 0.005448831829968969 |
    | ...          | ...       | ...       | ...          | ...                  |
    | chr22        | 16728001  | 16728101  | +            | 0.003834415188257777 |
    | chr19        | 17333425  | 17333525  | +            | 0.009636627605010294 |
    | chr17        | 8085927   | 8086027   | -            | 0.008917730007820798 |
    | chr16        | 52216522  | 52216622  | +            | 0.004375872112626925 |
    +--------------+-----------+-----------+--------------+----------------------+
    Stranded PyRanges object has 10 rows and 5 columns from 9 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.

    >>> gr.FDR = pr.stats.fdr(gr.PValue)
    >>> gr.print(formatting={"PValue": "{:.4f}", "FDR": "{:.4}"})
    +--------------+-----------+-----------+--------------+-------------+-------------+
    | Chromosome   | Start     | End       | Strand       | PValue      | FDR         |
    | (category)   | (int32)   | (int32)   | (category)   | (float64)   | (float64)   |
    |--------------+-----------+-----------+--------------+-------------+-------------|
    | chr1         | 176601938 | 176602038 | +            | 0.0055      | 0.01098     |
    | chr1         | 155082851 | 155082951 | -            | 0.0072      | 0.00894     |
    | chr2         | 211134424 | 211134524 | -            | 0.0060      | 0.01005     |
    | chr9         | 78826761  | 78826861  | -            | 0.0054      | 0.01362     |
    | ...          | ...       | ...       | ...          | ...         | ...         |
    | chr22        | 16728001  | 16728101  | +            | 0.0038      | 0.03834     |
    | chr19        | 17333425  | 17333525  | +            | 0.0096      | 0.009637    |
    | chr17        | 8085927   | 8086027   | -            | 0.0089      | 0.009909    |
    | chr16        | 52216522  | 52216622  | +            | 0.0044      | 0.01459     |
    +--------------+-----------+-----------+--------------+-------------+-------------+
    Stranded PyRanges object has 10 rows and 6 columns from 9 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.
    """

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
