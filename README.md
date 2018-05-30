# pyranges

[![Build Status](https://travis-ci.org/endrebak/pyranges.svg?branch=master)](https://travis-ci.org/endrebak/pyranges) [![Documentation Status](https://readthedocs.org/projects/pyranges/badge/?version=latest)](https://pyranges.readthedocs.io/?badge=latest) [![PyPI version](https://badge.fury.io/py/pyranges.svg)](https://badge.fury.io/py/pyranges)

GenomicRanges for Python.

This library tries to be a thin, but extremely useful wrapper around genomic data contained in pandas dataframes. This allows for all the wonderful functionality of bedtools/bedops and/or GenomicRanges, while being able to use the the enormous universe of Python datascience libraries to manipulate and do computations on the data.

PyRanges also contains a run-length encoding library for extremely efficient arithmetic computation of scores associated with genomic intervals.

### Paper

Being written here: https://github.com/endrebak/pyranges-paper

Timings will be added to docs when ready.

### Install

```bash
pip install pyranges
```

### TODO

Necessary for biorxiv-paper:

* Benchmarks comparing with Bioconductor GenomicRanges

For the future:

* Test code with hypothesis (https://hypothesis.readthedocs.io/en/latest/index.html)
* Allow writing to different UCSC genome-browser compatible formats such as
  bigwig, bedgraph, histograms, colored bed etc
* Add visualization capabilites?
* Enable annotation with featurefetch?
* Add dtypes to PyRanges-header?
* Find sequences of ranges with pyfaidx

```
import pyranges as pr

>>> cs = pr.load_dataset("chipseq")

>>> cs

+--------------|-----------|-----------|--------|---------|----------+
| Chromosome   | Start     | End       | Name   | Score   | Strand   |
|--------------|-----------|-----------|--------|---------|----------|
| chr8         | 28510032  | 28510057  | U0     | 0       | -        |
| chr7         | 107153363 | 107153388 | U0     | 0       | -        |
| chr5         | 135821802 | 135821827 | U0     | 0       | -        |
| ...          | ...       | ...       | ...    | ...     | ...      |
| chr6         | 89296757  | 89296782  | U0     | 0       | -        |
| chr1         | 194245558 | 194245583 | U0     | 0       | +        |
| chr8         | 57916061  | 57916086  | U0     | 0       | +        |
+--------------|-----------|-----------|--------|---------|----------+
PyRanges object has 10000 sequences from 24 chromosomes.

>>> bg = pr.load_dataset("chipseq_background")

>>> cs.nearest(bg, suffix="_IP")

+--------------|----------|----------|--------|---------|----------|-----------------|------------|----------|-----------|------------|-------------|------------+
| Chromosome   | Start    | End      | Name   | Score   | Strand   | Chromosome_IP   | Start_IP   | End_IP   | Name_IP   | Score_IP   | Strand_IP   | Distance   |
|--------------|----------|----------|--------|---------|----------|-----------------|------------|----------|-----------|------------|-------------|------------|
| chr1         | 1325303  | 1325328  | U0     | 0       | -        | chr1            | 1041102    | 1041127  | U0        | 0          | +           | 284176     |
| chr1         | 1541598  | 1541623  | U0     | 0       | +        | chr1            | 1770383    | 1770408  | U0        | 0          | -           | 228760     |
| chr1         | 1599121  | 1599146  | U0     | 0       | +        | chr1            | 1770383    | 1770408  | U0        | 0          | -           | 171237     |
| ...          | ...      | ...      | ...    | ...     | ...      | ...             | ...        | ...      | ...       | ...        | ...         | ...        |
| chrY         | 21910706 | 21910731 | U0     | 0       | -        | chrY            | 20557165   | 20557190 | U0        | 0          | +           | 1353516    |
| chrY         | 22054002 | 22054027 | U0     | 0       | -        | chrY            | 20557165   | 20557190 | U0        | 0          | +           | 1496812    |
| chrY         | 22210637 | 22210662 | U0     | 0       | -        | chrY            | 20557165   | 20557190 | U0        | 0          | +           | 1653447    |
+--------------|----------|----------|--------|---------|----------|-----------------|------------|----------|-----------|------------|-------------|------------+
PyRanges object has 10000 sequences from 24 chromosomes.

>>> cs.set_intersection(bg, strandedness="opposite")

+--------------|-----------|-----------|----------+
| Chromosome   |     Start |       End | Strand   |
|--------------|-----------|-----------|----------|
| chr1         | 226987603 | 226987617 | +        |
| chr8         |  38747236 |  38747251 | -        |
+--------------|-----------|-----------|----------+
PyRanges object has 2 sequences from 2 chromosomes.

>>> cv = cs.coverage(stranded=True)
>>> cv

chr1 +
+--------|-----------|------|---------|------|-----------|---------|------|-----------|------|-----------|------+
| Runs   |   1541598 |   25 |   57498 |   25 |   1904886 |  ...    |   25 |   2952580 |   25 |   1156833 |   25 |
|--------|-----------|------|---------|------|-----------|---------|------|-----------|------|-----------|------|
| Values |         0 |    1 |       0 |    1 |         0 | ...     |    1 |         0 |    1 |         0 |    1 |
+--------|-----------|------|---------|------|-----------|---------|------|-----------|------|-----------|------+
Rle of length 247134924 containing 894 elements
...
chrY -
+--------|-----------|------|----------|------|----------|---------|------|----------|------|----------|------+
| Runs   |   7046809 |   25 |   358542 |   25 |   296582 |  ...    |   25 |   143271 |   25 |   156610 |   25 |
|--------|-----------|------|----------|------|----------|---------|------|----------|------|----------|------|
| Values |         0 |    1 |        0 |    1 |        0 | ...     |    1 |        0 |    1 |        0 |    1 |
+--------|-----------|------|----------|------|----------|---------|------|----------|------|----------|------+
Rle of length 22210662 containing 32 elements
PyRles object with 48 chromosomes/strand pairs.

>>> cv + 10.42

chr1 +
+--------|-----------|-------|---------|-------|-----------|---------|-------|-----------|-------|-----------|-------+
| Runs   |   1541598 |    25 |   57498 |    25 |   1904886 |  ...    |    25 |   2952580 |    25 |   1156833 |    25 |
|--------|-----------|-------|---------|-------|-----------|---------|-------|-----------|-------|-----------|-------|
| Values |     10.42 | 11.42 |   10.42 | 11.42 |     10.42 | ...     | 11.42 |     10.42 | 11.42 |     10.42 | 11.42 |
+--------|-----------|-------|---------|-------|-----------|---------|-------|-----------|-------|-----------|-------+
Rle of length 247134924 containing 894 elements
...
chrY -
+--------|-----------|-------|----------|-------|----------|---------|-------|----------|-------|----------|-------+
| Runs   |   7046809 |    25 |   358542 |    25 |   296582 |  ...    |    25 |   143271 |    25 |   156610 |    25 |
|--------|-----------|-------|----------|-------|----------|---------|-------|----------|-------|----------|-------|
| Values |     10.42 | 11.42 |    10.42 | 11.42 |    10.42 | ...     | 11.42 |    10.42 | 11.42 |    10.42 | 11.42 |
+--------|-----------|-------|----------|-------|----------|---------|-------|----------|-------|----------|-------+
Rle of length 22210662 containing 32 elements
PyRles object with 48 chromosomes/strand pairs.

>>> bg_cv = bg.coverage()

>>> cv - bg_cv
chr1
+--------|----------|------|----------|------|---------|---------|------|----------|------|----------|------+
| Runs   |   887771 |   25 |   106864 |   25 |   46417 |  ...    |   25 |   730068 |   25 |   259250 |   25 |
|--------|----------|------|----------|------|---------|---------|------|----------|------|----------|------|
| Values |        0 |   -1 |        0 |   -1 |       0 | ...     |    1 |        0 |   -1 |        0 |    1 |
+--------|----------|------|----------|------|---------|---------|------|----------|------|----------|------+
Rle of length 247134924 containing 3242 elements
...
chrY
+--------|-----------|------|----------|------|----------|---------|------|----------|------|------------|------+
| Runs   |   7046809 |   25 |   147506 |   25 |   211011 |  ...    |   25 |   156610 |   25 |   35191552 |   25 |
|--------|-----------|------|----------|------|----------|---------|------|----------|------|------------|------|
| Values |         0 |    1 |        0 |    1 |        0 | ...     |    1 |        0 |    1 |          0 |   -1 |
+--------|-----------|------|----------|------|----------|---------|------|----------|------|------------|------+
Rle of length 57402239 containing 60 elements
Unstranded PyRles object with 25 chromosomes.
```

# See also

* https://github.com/endrebak/pyrle
* https://github.com/vsbuffalo/BioRanges/tree/master/BioRanges
* https://github.com/daler/pybedtools
* http://bedtools.readthedocs.io/en/latest/
* https://github.com/phaverty/RLEVectors.jl
* https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
* https://bedops.readthedocs.io/en/latest/
