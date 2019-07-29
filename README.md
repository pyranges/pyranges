# pyranges

[![Coverage Status](https://img.shields.io/coveralls/github/biocore-ntnu/pyranges.svg)](https://coveralls.io/github/biocore-ntnu/pyranges?branch=master) [![Codacy Badge](https://api.codacy.com/project/badge/Grade/b61a53346d764a8d8f0ab2a6afd7b100)](https://www.codacy.com/app/endrebak/pyranges?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=biocore-ntnu/pyranges&amp;utm_campaign=Badge_Grade) [![Build Status](https://travis-ci.org/biocore-ntnu/pyranges.svg?branch=master)](https://travis-ci.org/biocore-ntnu/pyranges) [![hypothesis tested](graphs/hypothesis-tested-brightgreen.svg)](http://hypothesis.readthedocs.io/) [![PyPI version](https://badge.fury.io/py/pyranges.svg)](https://badge.fury.io/py/pyranges) [![MIT](https://img.shields.io/pypi/l/pyranges.svg?color=green)](https://opensource.org/licenses/MIT) ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/pyranges.svg) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/pyranges/README.html)


## Introduction

GenomicRanges and genomic Rle-objects for Python.

*"Finally ... This was what Python badly needed for years."* - Heng Li

## Release

PyRanges is in a beta state. We are extremely responsive to bug-reports, so if you
have problems or come across unexpected behavior, please create an issue.

## Asking for help

Feel free to ask questions of the type "how do I do X with pyranges?" on public forums like Stack Overflow, Bioinformatics Stack Exchange or Biostars. You can use endrebak85#gmail.com to point me to the question(s).

## Quick example

```python
import pyranges as pr
# load example datasets
exons, cpg = pr.data.exons(), pr.data.cpg()

# subsetting pyranges is easy
exons["chrY", "-",  15591259:27197945]
# +--------------|-----------|-----------|----------------------------------------|-----------|--------------+
# | Chromosome   | Start     | End       | Name                                   | Score     | Strand       |
# | (category)   | (int32)   | (int32)   | (object)                               | (int64)   | (category)   |
# |--------------|-----------|-----------|----------------------------------------|-----------|--------------|
# | chrY         | 15591393  | 15592550  | NR_047610_exon_27_0_chrY_15591394_r    | 0         | -            |
# | chrY         | 15591393  | 15592550  | NR_047607_exon_29_0_chrY_15591394_r    | 0         | -            |
# | chrY         | 15591393  | 15592550  | NM_001258269_exon_29_0_chrY_15591394_r | 0         | -            |
# | chrY         | 15591393  | 15592550  | NR_047599_exon_28_0_chrY_15591394_r    | 0         | -            |
# | ...          | ...       | ...       | ...                                    | ...       | ...          |
# | chrY         | 25336491  | 25336631  | NM_004081_exon_22_0_chrY_25336492_r    | 0         | -            |
# | chrY         | 26952215  | 26952307  | NM_020364_exon_16_0_chrY_26952216_r    | 0         | -            |
# | chrY         | 27197822  | 27197945  | NM_004678_exon_7_0_chrY_27197823_r     | 0         | -            |
# | chrY         | 27197822  | 27197945  | NM_001002760_exon_7_0_chrY_27197823_r  | 0         | -            |
# +--------------|-----------|-----------|----------------------------------------|-----------|--------------+
# Stranded PyRanges object has 22 rows and 6 columns from 1 chromosomes.

# you can use your pandas-skills with pyranges
exons[~exons.Name.str.startswith("NR")] # all rows where the name column does not start with "NR"
# +--------------|-----------|-----------|----------------------------------------|-----------|--------------+
# | Chromosome   | Start     | End       | Name                                   | Score     | Strand       |
# | (category)   | (int32)   | (int32)   | (object)                               | (int64)   | (category)   |
# |--------------|-----------|-----------|----------------------------------------|-----------|--------------|
# | chrX         | 135574120 | 135574598 | NM_001727_exon_2_0_chrX_135574121_f    | 0         | +            |
# | chrX         | 47868945  | 47869126  | NM_205856_exon_4_0_chrX_47868946_f     | 0         | +            |
# | chrX         | 77294333  | 77294480  | NM_000052_exon_17_0_chrX_77294334_f    | 0         | +            |
# | chrX         | 91090459  | 91091043  | NM_001168360_exon_0_0_chrX_91090460_f  | 0         | +            |
# | ...          | ...       | ...       | ...                                    | ...       | ...          |
# | chrY         | 15481135  | 15481229  | NM_182659_exon_16_0_chrY_15481136_r    | 0         | -            |
# | chrY         | 25325872  | 25325936  | NM_004081_exon_18_0_chrY_25325873_r    | 0         | -            |
# | chrY         | 15560896  | 15560946  | NM_001258258_exon_25_0_chrY_15560897_r | 0         | -            |
# | chrY         | 15467254  | 15467278  | NM_001258270_exon_13_0_chrY_15467255_r | 0         | -            |
# +--------------|-----------|-----------|----------------------------------------|-----------|--------------+
# Stranded PyRanges object has 847 rows and 6 columns from 2 chromosomes.

# the API allows for easy and terse chaining
(cpg # use the cpg dataset
  .join(exons, suffix="_xn") # join with exons, use suffix _xn for duplicate cols
  .subset(lambda df: df.CpG > 30) # keep only rows with a CpG score over 30
  .sort(nb_cpu=2) # sort on Chromosome, Start and End
                  # note that virtually all pyranges-methods take a nb_cpu argument
                  # to use multiple cores, you need to install ray with pip install ray
  .pc() # print, while keeping the chain (c) going. Try .sp(), msp(), rp(), spc(), mspc() also :)
  ["chrX"] # keep only chromosome X
  .assign("CpGDecile", lambda df: df.CpG / 10) # Insert new column
  .unstrand()) # remove the strand info
# +--------------|-----------|-----------|-----------|------------|-----------|----------------------------------------|-----------+
# | Chromosome   | Start     | End       | CpG       | Start_xn   | End_xn    | Name                                   | Score     |
# | (category)   | (int32)   | (int32)   | (int64)   | (int32)    | (int32)   | (object)                               | (int64)   |
# |--------------|-----------|-----------|-----------|------------|-----------|----------------------------------------|-----------|
# | chrX         | 584563    | 585326    | 66        | 585078     | 585337    | NM_000451_exon_0_0_chrX_585079_f       | 0         |
# | chrX         | 1510501   | 1511838   | 173       | 1510791    | 1511039   | NM_001636_exon_3_0_chrX_1510792_r      | 0         |
# | chrX         | 2846195   | 2847511   | 92        | 2847272    | 2847416   | NM_001669_exon_9_0_chrX_2847273_r      | 0         |
# | chrX         | 13587648  | 13588221  | 49        | 13587693   | 13588054  | NM_001167890_exon_0_0_chrX_13587694_f  | 0         |
# | ...          | ...       | ...       | ...       | ...        | ...       | ...                                    | ...       |
# | chrY         | 15591259  | 15591720  | 33        | 15591393   | 15592550  | NR_047607_exon_29_0_chrY_15591394_r    | 0         |
# | chrY         | 15591259  | 15591720  | 33        | 15591393   | 15592550  | NM_001258269_exon_29_0_chrY_15591394_r | 0         |
# | chrY         | 15591259  | 15591720  | 33        | 15591393   | 15592550  | NR_047599_exon_28_0_chrY_15591394_r    | 0         |
# | chrY         | 16941822  | 16942188  | 32        | 16941609   | 16942399  | NM_014893_exon_4_0_chrY_16941610_f     | 0         |
# +--------------|-----------|-----------|-----------|------------|-----------|----------------------------------------|-----------+
# Unstranded PyRanges object has 65 rows and 8 columns from 2 chromosomes.
# +--------------|-----------|-----------|-----------|------------|-----------|-----------------------------------------|-----------|-------------+
# | Chromosome   | Start     | End       | CpG       | Start_xn   | End_xn    | Name                                    | Score     | CpGDecile   |
# | (category)   | (int32)   | (int32)   | (int64)   | (int32)    | (int32)   | (object)                                | (int64)   | (float64)   |
# |--------------|-----------|-----------|-----------|------------|-----------|-----------------------------------------|-----------|-------------|
# | chrX         | 584563    | 585326    | 66        | 585078     | 585337    | NM_000451_exon_0_0_chrX_585079_f        | 0         | 6.6         |
# | chrX         | 1510501   | 1511838   | 173       | 1510791    | 1511039   | NM_001636_exon_3_0_chrX_1510792_r       | 0         | 17.3        |
# | chrX         | 2846195   | 2847511   | 92        | 2847272    | 2847416   | NM_001669_exon_9_0_chrX_2847273_r       | 0         | 9.2         |
# | chrX         | 13587648  | 13588221  | 49        | 13587693   | 13588054  | NM_001167890_exon_0_0_chrX_13587694_f   | 0         | 4.9         |
# | ...          | ...       | ...       | ...       | ...        | ...       | ...                                     | ...       | ...         |
# | chrX         | 153068787 | 153070353 | 134       | 153067622  | 153070355 | NM_032512_exon_0_0_chrX_153067623_r     | 0         | 13.4        |
# | chrX         | 153284685 | 153285655 | 94        | 153284647  | 153284779 | NM_001025243_exon_10_0_chrX_153284648_r | 0         | 9.4         |
# | chrX         | 153598874 | 153600604 | 164       | 153599240  | 153599729 | NM_001456_exon_45_0_chrX_153599241_r    | 0         | 16.4        |
# | chrX         | 153990840 | 153991831 | 105       | 153991030  | 153991256 | NM_001363_exon_0_0_chrX_153991031_f     | 0         | 10.5        |
# +--------------|-----------|-----------|-----------|------------|-----------|-----------------------------------------|-----------|-------------+
# Unstranded PyRanges object has 58 rows and 9 columns from 1 chromosomes.

cpg
# +--------------|-----------|-----------|-----------+
# | Chromosome   | Start     | End       | CpG       |
# | (category)   | (int32)   | (int32)   | (int64)   |
# |--------------|-----------|-----------|-----------|
# | chrX         | 64181     | 64793     | 62        |
# | chrX         | 69133     | 70029     | 100       |
# | chrX         | 148685    | 149461    | 85        |
# | chrX         | 166504    | 167721    | 96        |
# | ...          | ...       | ...       | ...       |
# | chrY         | 28555535  | 28555932  | 32        |
# | chrY         | 28773315  | 28773544  | 25        |
# | chrY         | 59213794  | 59214183  | 36        |
# | chrY         | 59349266  | 59349574  | 29        |
# +--------------|-----------|-----------|-----------+
# Unstranded PyRanges object has 1,077 rows and 4 columns from 2 chromosomes.

cpg_rle = cpg.to_rle(value_col="CpG") # ignore value_col empty for regular coverage
cpg_rle
# chrX
# ----
# +--------|---------|-------|--------|-------|---------|---------|-------|----------|-------|----------|-------+
# | Runs   | 64181   | 612   | 4340   | 896   | 78656   |  ...    | 607   | 268069   | 389   | 135083   | 308   |
# |--------|---------|-------|--------|-------|---------|---------|-------|----------|-------|----------|-------|
# | Values | 0.0     | 62.0  | 0.0    | 100.0 | 0.0     | ...     | 44.0  | 0.0      | 36.0  | 0.0      | 29.0  |
# +--------|---------|-------|--------|-------|---------|---------|-------|----------|-------|----------|-------+
# Rle of length 155246568 containing 1792 elements
#
# chrY
# ----
# +--------|---------|-------|--------|-------|---------|---------|-------|------------|-------|----------|-------+
# | Runs   | 14181   | 612   | 4340   | 896   | 78656   |  ...    | 229   | 30440250   | 389   | 135083   | 308   |
# |--------|---------|-------|--------|-------|---------|---------|-------|------------|-------|----------|-------|
# | Values | 0.0     | 62.0  | 0.0    | 100.0 | 0.0     | ...     | 25.0  | 0.0        | 36.0  | 0.0      | 29.0  |
# +--------|---------|-------|--------|-------|---------|---------|-------|------------|-------|----------|-------+
# Rle of length 59349574 containing 362 elements
# Unstranded PyRles object with 2 chromosomes.

cpg_rle * exons.to_rle() # only keep cpg_rle scores where exons has coverage
# chrX
# ----
# +--------|----------|-------|----------|-------|---------|---------|----------|-------|----------|-------|-----------+
# | Runs   | 585078   | 248   | 925465   | 248   | 42875   |  ...    | 314461   | 489   | 391301   | 226   | 1255312   |
# |--------|----------|-------|----------|-------|---------|---------|----------|-------|----------|-------|-----------|
# | Values | 0.0      | 66.0  | 0.0      | 173.0 | 0.0     | ...     | 0.0      | 164.0 | 0.0      | 105.0 | 0.0       |
# +--------|----------|-------|----------|-------|---------|---------|----------|-------|----------|-------|-----------+
# Rle of length 155246568 containing 131 elements
#
# chrY
# ----
# +--------|----------|-------|---------|-------|-----------|---------|-----------|-------|------------|-------|------------+
# | Runs   | 155399   | 137   | 89131   | 585   | 1117968   |  ...    | 1350102   | 366   | 10037778   | 150   | 32369458   |
# |--------|----------|-------|---------|-------|-----------|---------|-----------|-------|------------|-------|------------|
# | Values | 0.0      | 19.0  | 0.0     | 310.0 | 0.0       | ...     | 0.0       | 32.0  | 0.0        | 21.0  | 0.0        |
# +--------|----------|-------|---------|-------|-----------|---------|-----------|-------|------------|-------|------------+
# Rle of length 59349574 containing 15 elements
# Unstranded PyRles object with 2 chromosomes.

cpg_rle["chrX"] * 5
# +--------|---------|-------|--------|-------|---------|---------|-------|----------|-------|----------|-------+
# | Runs   | 64181   | 612   | 4340   | 896   | 78656   |  ...    | 607   | 268069   | 389   | 135083   | 308   |
# |--------|---------|-------|--------|-------|---------|---------|-------|----------|-------|----------|-------|
# | Values | 0.0     | 310.0 | 0.0    | 500.0 | 0.0     | ...     | 220.0 | 0.0      | 180.0 | 0.0      | 145.0 |
# +--------|---------|-------|--------|-------|---------|---------|-------|----------|-------|----------|-------+
# Rle of length 155246568 containing 1792 elements

cpg_rle["chrX"][64100:100000]
# +--------|------|-------|--------|-------|---------+
# | Runs   | 81   | 612   | 4340   | 896   | 29971   |
# |--------|------|-------|--------|-------|---------|
# | Values | 0.0  | 62.0  | 0.0    | 100.0 | 0.0     |
# +--------|------|-------|--------|-------|---------+
# Rle of length 35900 containing 5 elements
```

## Features

- fast (also in single-core mode)
- supports multiple cores
- memory-efficient
- featureful
- pythonic/pandastic
- supports chaining with a terse syntax
- uses Pandas DataFrames, so the whole Python data science stack works on PyRanges.

## Documentation

<https://biocore-ntnu.github.io/pyranges/>

(Might be slightly out of date; watch the CHANGELOG too)

## Install

```bash
pip install pyranges
```

## Paper/Cite

https://www.biorxiv.org/content/10.1101/609396v1

## TODO

For the future:

*   groupby
*   settings
*   write docstrings, autogenerate API-docs
*   K-nearest
*   write bam

Potentially:

*   support mixes of stranded and unstranded data

PyRanges should always be the fastest general-purpose genomics library for
Python. So I will happily change the multithreading library and overlap
datastructures sometime in the future, if rigorous tests show that the proposed
alternatives are indeed faster. (As the multithreading requires about 30 lines
of code and the overlap queries about 15, this will not be hard.)

## Performance

<img src="./graphs/main_paper_annotation_binary.png" />

Comprehensive set of graphs for many types of functions on different datasets are here:

[Time](https://github.com/endrebak/pyranges-paper/blob/master/supplementary_paper/time.md)

[Memory](https://github.com/endrebak/pyranges-paper/blob/master/supplementary_paper/memory.md)

The exact code tested is found [here](https://github.com/endrebak/pyranges-paper/tree/master/supplementaries).
