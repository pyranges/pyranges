from __future__ import print_function


try:
    import mkl
    mkl.set_num_threads(1)
except ImportError:
    pass

import pandas as pd
import numpy as np

import pkg_resources

from pyranges.pyranges import PyRanges
from pyranges.readers import read_gtf, read_bam, read_bed, read_gff3
from pyranges import data
from pyranges.methods.concat import concat

from pyrle import PyRles, Rle

from pyranges.version import __version__

get_example_path = data.get_example_path

read_gff = read_gtf

def from_dict(d, int64=False):

    return PyRanges(pd.DataFrame(d), int64=int64)

def from_string(s, int64=False):

    from io import StringIO
    df = pd.read_csv(StringIO(s), sep=r"\s+", index_col=None)

    return PyRanges(df, int64=int64)


import pyranges.genomicfeatures.genomicfeatures as gf

random = gf.random

from pyranges.methods.itergrs import itergrs
iter = itergrs

from pyranges.methods.multioverlap import count_overlaps#, interval_split

from pyranges import statistics
stats = statistics

# __all__ = [read_gff3, read_bed, read_bam, read_gtf, concat, iter, count_overlaps]

def to_bigwig(gr, path, chromosome_sizes):

    """Write df to bigwig.

    Must contain the columns Chromosome, Start, End and Score. All others are ignored.

    Parameters
    ----------
    path : str

        Where to write bigwig.

    chromosome_sizes : PyRanges or dict

        If dict: map of chromosome names to chromosome length.

    Examples
    --------

    Extended example with how to prepare your data for writing bigwigs:

    >>> d =  {'Chromosome': ['chr1', 'chr1', 'chr1'], 'Start': [1, 4, 6],
    ...       'End': [7, 8, 10], 'Strand': ['+', '-', '-'],
    ...       'Value': [10, 20, 30]}
    >>> import pyranges as pr
    >>> gr = pr.from_dict(d)
    >>> hg19 = pr.data.chromsizes()
    >>> print(hg19)
    +--------------+-----------+-----------+
    | Chromosome   | Start     | End       |
    | (category)   | (int32)   | (int32)   |
    |--------------+-----------+-----------|
    | chr1         | 0         | 249250621 |
    | chr2         | 0         | 243199373 |
    | chr3         | 0         | 198022430 |
    | chr4         | 0         | 191154276 |
    | ...          | ...       | ...       |
    | chrY         | 0         | 59373566  |
    | chrX         | 0         | 155270560 |
    | chrM         | 0         | 16571     |
    | chr22        | 0         | 51304566  |
    +--------------+-----------+-----------+
    Unstranded PyRanges object has 25 rows and 3 columns from 25 chromosomes.
    For printing, the PyRanges was sorted on Chromosome.

    >>> to_bigwig(gr, "outpath.bw", hg19)
    Traceback (most recent call last):
    ...
    AssertionError: Can only write one strand at a time. Use an unstranded PyRanges or subset on strand first.

    >>> to_bigwig(gr["-"], "outpath.bw", hg19)
    Traceback (most recent call last):
    ...
    AssertionError: Intervals must not overlap.

    >>> gr
    +--------------+-----------+-----------+--------------+-----------+
    | Chromosome   |     Start |       End | Strand       |     Value |
    | (category)   |   (int32) |   (int32) | (category)   |   (int64) |
    |--------------+-----------+-----------+--------------+-----------|
    | chr1         |         1 |         7 | +            |        10 |
    | chr1         |         4 |         8 | -            |        20 |
    | chr1         |         6 |        10 | -            |        30 |
    +--------------+-----------+-----------+--------------+-----------+
    Stranded PyRanges object has 3 rows and 5 columns from 1 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.

    >>> value = gr.to_rle(rpm=False, value_col="Value")
    >>> value
    chr1 +
    --
    +--------+-----+------+
    | Runs   | 1   | 6    |
    |--------+-----+------|
    | Values | 0.0 | 10.0 |
    +--------+-----+------+
    Rle of length 7 containing 2 elements (avg. length 3.5)
    <BLANKLINE>
    chr1 -
    --
    +--------+-----+------+------+------+
    | Runs   | 4   | 2    | 2    | 2    |
    |--------+-----+------+------+------|
    | Values | 0.0 | 20.0 | 50.0 | 30.0 |
    +--------+-----+------+------+------+
    Rle of length 10 containing 4 elements (avg. length 2.5)
    PyRles object with 2 chromosomes/strand pairs.

    >>> raw = gr.to_rle(rpm=False)
    >>> raw
    chr1 +
    --
    +--------+-----+-----+
    | Runs   | 1   | 6   |
    |--------+-----+-----|
    | Values | 0.0 | 1.0 |
    +--------+-----+-----+
    Rle of length 7 containing 2 elements (avg. length 3.5)
    <BLANKLINE>
    chr1 -
    --
    +--------+-----+-----+-----+-----+
    | Runs   | 4   | 2   | 2   | 2   |
    |--------+-----+-----+-----+-----|
    | Values | 0.0 | 1.0 | 2.0 | 1.0 |
    +--------+-----+-----+-----+-----+
    Rle of length 10 containing 4 elements (avg. length 2.5)
    PyRles object with 2 chromosomes/strand pairs.

    >>> result = (value / raw).apply_values(np.log10)
    >>> result
    chr1 +
    --
    +--------+-----+-----+
    | Runs   | 1   | 6   |
    |--------+-----+-----|
    | Values | nan | 1.0 |
    +--------+-----+-----+
    Rle of length 7 containing 2 elements (avg. length 3.5)
    <BLANKLINE>
    chr1 -
    --
    +--------+-----+--------------------+--------------------+--------------------+
    | Runs   | 4   | 2                  | 2                  | 2                  |
    |--------+-----+--------------------+--------------------+--------------------|
    | Values | nan | 1.3010300397872925 | 1.3979400396347046 | 1.4771212339401245 |
    +--------+-----+--------------------+--------------------+--------------------+
    Rle of length 10 containing 4 elements (avg. length 2.5)
    PyRles object with 2 chromosomes/strand pairs.

    >>> out = result.numbers_only().to_ranges()
    >>> out
    +--------------+-----------+-----------+-------------+------------+
    | Chromosome   |     Start |       End |       Score | Strand     |
    | (object)     |   (int32) |   (int32) |   (float64) | (object)   |
    |--------------+-----------+-----------+-------------+------------|
    | chr1         |         1 |         7 |     1       | +          |
    | chr1         |         4 |         6 |     1.30103 | -          |
    | chr1         |         6 |         8 |     1.39794 | -          |
    | chr1         |         8 |        10 |     1.47712 | -          |
    +--------------+-----------+-----------+-------------+------------+
    Stranded PyRanges object has 4 rows and 5 columns from 1 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.

    >>> to_bigwig(out["-"], "deleteme_reverse.bw", hg19)
    >>> to_bigwig(out["+"], "deleteme_forward.bw", hg19)
    """

    try:
        import pyBigWig
    except ModuleNotFoundError:
        print("pybigwig must be installed to create bigwigs. Use `conda install -c bioconda pybigwig` or `pip install pybigwig` to install it.")
        import sys
        sys.exit(1)

    assert len(gr.strands) <= 1, "Can only write one strand at a time. Use an unstranded PyRanges or subset on strand first."
    assert np.sum(gr.lengths()) == gr.merge().length, "Intervals must not overlap."

    df = gr.df

    unique_chromosomes = list(df.Chromosome.drop_duplicates())

    if not isinstance(chromosome_sizes, dict):
        size_df = chromosome_sizes.df
        chromosome_sizes = {k: v for k, v in zip(size_df.Chromosome, size_df.End)}

    header = [(c, int(chromosome_sizes[c])) for c in unique_chromosomes]

    bw = pyBigWig.open(path, "w")
    bw.addHeader(header)

    chromosomes = df.Chromosome.tolist()
    starts = df.Start.tolist()
    ends = df.End.tolist()
    values = df.Score.tolist()

    bw.addEntries(chromosomes, starts, ends=ends, values=values)
