"""Statistics useful for genomics."""

import pandas as pd
import numpy as np

import pyranges as pr
from pyranges.multithreaded import pyrange_apply

from pyranges.methods.statistics import _relative_distance

from collections import defaultdict

__all__ = ["simes", "fisher_exact", "StatisticsMethods", "fdr", "rowbased_rankdata", "rowbased_pearson", "rowbased_spearman", "mcc"]

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
    | chr16        | 52216522  | 52216622  | +            | 0.004375872112626925 |
    | chr17        | 8085927   | 8086027   | -            | 0.008917730007820798 |
    | chr19        | 17333425  | 17333525  | +            | 0.009636627605010294 |
    | chr22        | 16728001  | 16728101  | +            | 0.003834415188257777 |
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
    | chr16        | 52216522  | 52216622  | +            | 0.0044      | 0.01459     |
    | chr17        | 8085927   | 8086027   | -            | 0.0089      | 0.009909    |
    | chr19        | 17333425  | 17333525  | +            | 0.0096      | 0.009637    |
    | chr22        | 16728001  | 16728101  | +            | 0.0038      | 0.03834     |
    +--------------+-----------+-----------+--------------+-------------+-------------+
    Stranded PyRanges object has 10 rows and 6 columns from 9 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    from scipy.stats import rankdata
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1

    return fdr


def fisher_exact(tp, fp, fn, tn, pseudocount=0):

    """Fisher's exact for contingency tables.

    Computes the hypotheses two-sided, less and greater at the same time.

    The odds-ratio is

    Parameters
    ----------
    tp : array-like of int

        Top left square of contingency table (true positives).

    fp : array-like of int

        Top right square of contingency table (false positives).

    fn : array-like of int

        Bottom left square of contingency table (false negatives).

    tn : array-like of int

        Bottom right square of contingency table (true negatives).

    pseudocount : float, default 0

        Values > 0 allow Odds Ratio to always be a finite number.

    Notes
    -----

    The odds-ratio is computed thusly:

    ``((tp + pseudocount) / (fp + pseudocount)) / ((fn + pseudocount) / (tn + pseudocount))``

    Returns
    -------
    pandas.DataFrame

        DataFrame with columns OR and P, PLeft and PRight.

    See Also
    --------

    pr.stats.fdr : correct for multiple testing

    Examples
    --------

    >>> d = {"TP": [12, 0], "FP": [5, 12], "TN": [29, 10], "FN": [2, 2]}
    >>> df = pd.DataFrame(d)
    >>> df
       TP  FP  TN  FN
    0  12   5  29   2
    1   0  12  10   2

    >>> pr.stats.fisher_exact(df.TP, df.FP, df.TN, df.FN)
             OR         P     PLeft    PRight
    0  0.165517  0.080269  0.044555  0.994525
    1  0.000000  0.000067  0.000034  1.000000
    """


    try:
        from fisher import pvalue_npy
    except:
        import sys
        print("fisher needs to be installed to use fisher exact. pip install fisher or conda install -c bioconda fisher.")
        sys.exit(-1)

    tp = np.array(tp, dtype=np.uint)
    fp = np.array(fp, dtype=np.uint)
    fn = np.array(fn, dtype=np.uint)
    tn = np.array(tn, dtype=np.uint)

    left, right, twosided = pvalue_npy(tp, fp, fn, tn)

    OR = ((tp + pseudocount) / (fp + pseudocount)) / ((fn + pseudocount) / (tn + pseudocount))

    df = pd.DataFrame({"OR": OR, "P": twosided, "PLeft": left, "PRight": right})

    return df


def mcc(grs, genome=None, labels=None, strand=False, verbose=False):

    """Compute Matthew's correlation coefficient for PyRanges overlaps.

    Parameters
    ----------
    grs : list of PyRanges

        PyRanges to compare.

    genome : DataFrame or dict, default None

        Should contain chromosome sizes. By default, end position of the
        rightmost intervals are used as proxies for the chromosome size, but
        it is recommended to use a genome.

    labels : list of str, default None

        Names to give the PyRanges in the output.

    strand : bool, default False

        Whether to compute correlations per strand.

    verbose : bool, default False

        Warn if some chromosomes are in the genome, but not in the PyRanges.

    Examples
    --------
    >>> np.random.seed(0)
    >>> chromsizes = {"chrM": 16000}
    >>> grs = [pr.random(chromsizes=chromsizes) for _ in range(3)]
    >>> labels = ["a", "b", "c"]
    >>> mcc = pr.stats.mcc(grs, labels=labels, genome=chromsizes)
    >>> mcc
       T  F     TP  FP  TN  FN       MCC
    0  a  a  15920   0  80   0  1.000000
    1  a  b  15875  65  15  45  0.213109
    3  a  c  15896  72   8  24  0.155496
    2  b  a  15875  45  15  65  0.213109
    5  b  b  15940   0  60   0  1.000000
    6  b  c  15916  52   8  24  0.180354
    4  c  a  15896  24   8  72  0.155496
    7  c  b  15916  24   8  52  0.180354
    8  c  c  15968   0  32   0  1.000000

    To create a symmetric matrix (useful for heatmaps of correlations):

    >>> mcc.set_index(["T", "F"]).MCC.unstack()
    F         a         b         c
    T
    a  1.000000  0.213109  0.155496
    b  0.213109  1.000000  0.180354
    c  0.155496  0.180354  1.000000"""

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
        if isinstance(genome, (pd.DataFrame, pr.PyRanges)):
            genome_length = int(genome.End.sum())
        else:
            genome_length = sum(genome.values())

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

    else:
        d = defaultdict(int)
        for gr in grs:
            for k, v in gr:
                d[k] = max(d[k], v.End.max())

        genome_length = sum(d.values())

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



def chromsizes_as_int(chromsizes):
        if isinstance(chromsizes, int):
            pass
        elif isinstance(chromsizes, dict):
            chromsizes = sum(chromsizes.values())
        elif isinstance(chromsizes, (pd.DataFrame, pr.PyRanges)):
            chromsizes = chromsizes.End.sum()

        return chromsizes


class StatisticsMethods():

    """Namespace for statistical comparsion-operations.

    Accessed with gr.stats."""

    pr = None

    def __init__(self, pr):

        self.pr = pr


    def forbes(self, other, chromsizes, strandedness=None):

        """Compute Forbes coefficient.

        Ratio which represents observed versus expected co-occurence.

        Described in ``Forbes SA (1907): On the local distribution of certain Illinois fishes: an essay in statistical ecology.``

        Parameters
        ----------
        other : PyRanges

            Intervals to compare with.

        chromsizes : int, dict, DataFrame or PyRanges

            Integer representing genome length or mapping from chromosomes
            to its length.

        strandedness : {None, "same", "opposite", False}, default None, i.e. "auto"

            Whether to compute without regards to strand or on same or opposite.

        Returns
        -------
        float

            Ratio of observed versus expected co-occurence.

        See Also
        --------

        pyranges.statistics.jaccard : compute the jaccard coefficient

        Examples
        --------

        >>> gr, gr2 = pr.data.chipseq(), pr.data.chipseq_background()
        >>> chromsizes = pr.data.chromsizes()
        >>> gr.stats.forbes(gr2, chromsizes=chromsizes)
        1.7168314674978278"""

        chromsizes = chromsizes_as_int(chromsizes)

        self = self.pr

        kwargs = {}
        kwargs["sparse"] = {"self": True, "other": True}
        kwargs = pr.pyranges.fill_kwargs(kwargs)
        strand = True if kwargs.get("strandedness") else False

        reference_length = self.merge(strand=strand).length
        query_length = other.merge(strand=strand).length

        intersection_sum = sum(
            v.sum()
            for v in self.set_intersect(
                    other, strandedness=strandedness).lengths(as_dict=True).values())

        forbes = chromsizes * intersection_sum / (reference_length * query_length)

        return forbes


    def jaccard(self, other, **kwargs):

        """Compute Jaccards coefficient.

        Ratio of the intersection and union of two sets.

        Parameters
        ----------
        other : PyRanges

            Intervals to compare with.

        chromsizes : int, dict, DataFrame or PyRanges

            Integer representing genome length or mapping from chromosomes
            to its length.

        strandedness : {None, "same", "opposite", False}, default None, i.e. "auto"

            Whether to compute without regards to strand or on same or opposite.

        Returns
        -------
        float

            Ratio of the intersection and union of two sets.

        See Also
        --------

        pyranges.statistics.forbes : compute the forbes coefficient

        Examples
        --------

        >>> gr, gr2 = pr.data.chipseq(), pr.data.chipseq_background()
        >>> chromsizes = pr.data.chromsizes()
        >>> gr.stats.jaccard(gr2, chromsizes=chromsizes)
        6.657941988519211e-05"""

        self = self.pr

        kwargs["sparse"] = {"self": True, "other": True}
        kwargs = pr.pyranges.fill_kwargs(kwargs)
        strand = True if kwargs.get("strandedness") else False

        intersection_sum = sum(
            v.sum()
            for v in self.set_intersect(other).lengths(as_dict=True).values())

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

    def relative_distance(self, other):

        """Compute spatial correllation between two sets.

        Metric which describes relative distance between each interval in one
        set and two closest intervals in another.

        Parameters
        ----------
        other : PyRanges

            Intervals to compare with.

        chromsizes : int, dict, DataFrame or PyRanges

            Integer representing genome length or mapping from chromosomes
            to its length.

        strandedness : {None, "same", "opposite", False}, default None, i.e. "auto"

            Whether to compute without regards to strand or on same or opposite.

        Returns
        -------
        pandas.DataFrame

            DataFrame containing the frequency of each relative distance.

        See Also
        --------

        pyranges.statistics.jaccard : compute the jaccard coefficient
        pyranges.statistics.forbes : compute the forbes coefficient

        Examples
        --------

        >>> gr, gr2 = pr.data.chipseq(), pr.data.chipseq_background()
        >>> chromsizes = pr.data.chromsizes()
        >>> gr.stats.relative_distance(gr2)
            reldist  count  total  fraction
        0      0.00    264   9956  0.026517
        1      0.01    226   9956  0.022700
        2      0.02    206   9956  0.020691
        3      0.03    235   9956  0.023604
        4      0.04    194   9956  0.019486
        5      0.05    241   9956  0.024207
        6      0.06    201   9956  0.020189
        7      0.07    191   9956  0.019184
        8      0.08    192   9956  0.019285
        9      0.09    191   9956  0.019184
        10     0.10    186   9956  0.018682
        11     0.11    203   9956  0.020390
        12     0.12    218   9956  0.021896
        13     0.13    209   9956  0.020992
        14     0.14    201   9956  0.020189
        15     0.15    178   9956  0.017879
        16     0.16    202   9956  0.020289
        17     0.17    197   9956  0.019787
        18     0.18    208   9956  0.020892
        19     0.19    202   9956  0.020289
        20     0.20    191   9956  0.019184
        21     0.21    188   9956  0.018883
        22     0.22    213   9956  0.021394
        23     0.23    192   9956  0.019285
        24     0.24    199   9956  0.019988
        25     0.25    181   9956  0.018180
        26     0.26    172   9956  0.017276
        27     0.27    191   9956  0.019184
        28     0.28    190   9956  0.019084
        29     0.29    192   9956  0.019285
        30     0.30    201   9956  0.020189
        31     0.31    212   9956  0.021294
        32     0.32    213   9956  0.021394
        33     0.33    177   9956  0.017778
        34     0.34    197   9956  0.019787
        35     0.35    163   9956  0.016372
        36     0.36    191   9956  0.019184
        37     0.37    198   9956  0.019888
        38     0.38    160   9956  0.016071
        39     0.39    188   9956  0.018883
        40     0.40    200   9956  0.020088
        41     0.41    188   9956  0.018883
        42     0.42    230   9956  0.023102
        43     0.43    197   9956  0.019787
        44     0.44    224   9956  0.022499
        45     0.45    184   9956  0.018481
        46     0.46    198   9956  0.019888
        47     0.47    187   9956  0.018783
        48     0.48    200   9956  0.020088
        49     0.49    194   9956  0.019486
        """

        self = self.pr

        kwargs = {}
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
