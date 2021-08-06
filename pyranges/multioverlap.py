import pyranges as pr
import numpy as np


def count_overlaps(grs, features=None, strandedness=None, how=None,  nb_cpu=1):

    """Count overlaps in multiple pyranges.

    Parameters
    ----------
    grs : dict of PyRanges

        The PyRanges to use as queries.

    features : PyRanges, default None

        The PyRanges to use as subject in the query. If None, the PyRanges themselves are used as a query.

    strandedness : {None, "same", "opposite", False}, default None, i.e. auto

        Whether to compare PyRanges on the same strand, the opposite or ignore strand
        information. The default, None, means use "same" if both PyRanges are stranded,
        otherwise ignore the strand information.

     how : {None, "all", "containment"}, default None, i.e. all

        What intervals to report. By default reports all overlapping intervals. "containment"
        reports intervals where the overlapping is contained within it.

    nb_cpu : int, default 1

        How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
        Will only lead to speedups on large datasets.

    Examples
    --------

    >>> a = '''Chromosome Start End
    ... chr1    6    12
    ... chr1    10    20
    ... chr1    22    27
    ... chr1    24    30'''

    >>> b = '''Chromosome Start End
    ... chr1    12    32
    ... chr1    14    30'''

    >>> c = '''Chromosome Start End
    ... chr1    8    15
    ... chr1    10    14
    ... chr1    32    34'''

    >>> grs = {n: pr.from_string(s) for n, s in zip(["a", "b", "c"], [a, b, c])}
    >>> for k, v in grs.items():
    ...     print("Name: " + k)
    ...     print(v)
    Name: a
    +--------------+-----------+-----------+
    | Chromosome   |     Start |       End |
    | (category)   |   (int32) |   (int32) |
    |--------------+-----------+-----------|
    | chr1         |         6 |        12 |
    | chr1         |        10 |        20 |
    | chr1         |        22 |        27 |
    | chr1         |        24 |        30 |
    +--------------+-----------+-----------+
    Unstranded PyRanges object has 4 rows and 3 columns from 1 chromosomes.
    For printing, the PyRanges was sorted on Chromosome.
    Name: b
    +--------------+-----------+-----------+
    | Chromosome   |     Start |       End |
    | (category)   |   (int32) |   (int32) |
    |--------------+-----------+-----------|
    | chr1         |        12 |        32 |
    | chr1         |        14 |        30 |
    +--------------+-----------+-----------+
    Unstranded PyRanges object has 2 rows and 3 columns from 1 chromosomes.
    For printing, the PyRanges was sorted on Chromosome.
    Name: c
    +--------------+-----------+-----------+
    | Chromosome   |     Start |       End |
    | (category)   |   (int32) |   (int32) |
    |--------------+-----------+-----------|
    | chr1         |         8 |        15 |
    | chr1         |        10 |        14 |
    | chr1         |        32 |        34 |
    +--------------+-----------+-----------+
    Unstranded PyRanges object has 3 rows and 3 columns from 1 chromosomes.
    For printing, the PyRanges was sorted on Chromosome.

    >>> pr.count_overlaps(grs)
    +--------------+-----------+-----------+-----------+-----------+-----------+
    | Chromosome   | Start     | End       | a         | b         | c         |
    | (object)     | (int32)   | (int32)   | (int32)   | (int32)   | (int32)   |
    |--------------+-----------+-----------+-----------+-----------+-----------|
    | chr1         | 6         | 8         | 1         | 0         | 0         |
    | chr1         | 8         | 10        | 1         | 0         | 1         |
    | chr1         | 10        | 12        | 2         | 0         | 2         |
    | chr1         | 12        | 14        | 1         | 1         | 2         |
    | ...          | ...       | ...       | ...       | ...       | ...       |
    | chr1         | 24        | 27        | 2         | 2         | 0         |
    | chr1         | 27        | 30        | 1         | 2         | 0         |
    | chr1         | 30        | 32        | 0         | 1         | 0         |
    | chr1         | 32        | 34        | 0         | 0         | 1         |
    +--------------+-----------+-----------+-----------+-----------+-----------+
    Unstranded PyRanges object has 12 rows and 6 columns from 1 chromosomes.
    For printing, the PyRanges was sorted on Chromosome.

    >>> gr = pr.PyRanges(chromosomes=["chr1"] * 4, starts=[0, 10, 20, 30], ends=[10, 20, 30, 40])
    >>> gr
    +--------------+-----------+-----------+
    | Chromosome   |     Start |       End |
    | (category)   |   (int32) |   (int32) |
    |--------------+-----------+-----------|
    | chr1         |         0 |        10 |
    | chr1         |        10 |        20 |
    | chr1         |        20 |        30 |
    | chr1         |        30 |        40 |
    +--------------+-----------+-----------+
    Unstranded PyRanges object has 4 rows and 3 columns from 1 chromosomes.
    For printing, the PyRanges was sorted on Chromosome.

    >>> pr.count_overlaps(grs, gr)
    +--------------+-----------+-----------+-----------+-----------+-----------+
    | Chromosome   |     Start |       End |         a |         b |         c |
    | (category)   |   (int32) |   (int32) |   (int32) |   (int32) |   (int32) |
    |--------------+-----------+-----------+-----------+-----------+-----------|
    | chr1         |         0 |        10 |         1 |         0 |         1 |
    | chr1         |        10 |        20 |         2 |         2 |         2 |
    | chr1         |        20 |        30 |         2 |         2 |         0 |
    | chr1         |        30 |        40 |         0 |         1 |         1 |
    +--------------+-----------+-----------+-----------+-----------+-----------+
    Unstranded PyRanges object has 4 rows and 6 columns from 1 chromosomes.
    For printing, the PyRanges was sorted on Chromosome.
    """

    kwargs = {"as_pyranges": False, "nb_cpu": nb_cpu, "strandedness": strandedness, "how": how, "nb_cpu": nb_cpu}
    names = list(grs.keys())

    if features is None:
        features = pr.concat(grs.values()).split(between=True)
    else:
        features = features.copy()

    from pyranges.methods.intersection import _count_overlaps

    for name, gr in grs.items():

        gr = gr.drop()

        kwargs["name"] = name
        res = features.apply_pair(gr, _count_overlaps, **kwargs)

    def to_int(df):
        df.loc[:, names] = df[names].astype(np.int32)
        return df

    features = features.apply(to_int)

    return features
