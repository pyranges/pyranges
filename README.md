# pyranges

[![Coverage Status](https://img.shields.io/coveralls/github/biocore-ntnu/pyranges.svg)](https://coveralls.io/github/biocore-ntnu/pyranges?branch=master) [![Codacy Badge](https://api.codacy.com/project/badge/Grade/b61a53346d764a8d8f0ab2a6afd7b100)](https://www.codacy.com/app/endrebak/pyranges?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=biocore-ntnu/pyranges&amp;utm_campaign=Badge_Grade) [![Build Status](https://travis-ci.org/biocore-ntnu/pyranges.svg?branch=master)](https://travis-ci.org/biocore-ntnu/pyranges) [![hypothesis tested](graphs/hypothesis-tested-brightgreen.svg)](http://hypothesis.readthedocs.io/) [![PyPI version](https://badge.fury.io/py/pyranges.svg)](https://badge.fury.io/py/pyranges) [![MIT](https://img.shields.io/pypi/l/pyranges.svg?color=green)](https://opensource.org/licenses/MIT)

## Release

PyRanges is in a beta state.

## Introduction

GenomicRanges and genomic Rle-objects for Python.

## Changelog

```
# 0.0.22 (14.04.19)
Additions:
  - pr.PyRanges() returns empty PyRange # before you needed pr.PyRanges({})
  - pyranges are now callable. Examples: gr("df.Score > 0") and gr("df.A.astype(str) + mysuffix")
  - can subset PyRanges with a dict of boolean vectors
  - pr.data.exons(), pr.data.cpg()
  - gr.unstrand() removes strand information from a PyRanges
  - throw exception if trying to drop Strand from df without setting drop_strand=True
  - adding a Strand column to the PyRanges makes it stranded

Changes:
  - write dtype as category, not int8/int16/...

Fixes:
  - remove empty dfs in the dict given to the PyRanges constructor

# 0.0.21 (14.04.19)
Additions:
  - gr.cluster(): assign ID to each cluster found by merge
  - gr.columns: return the columns in the pyranges
  - gr.drop: drop columns based on regex or list
  - gr[["Score", "Name"]]: select subset of columns
Fixes:
  - gr.stranded errored if chromosomes were ints
  - gr.join errored if other had duplicate indexes
```

## Quick example

```
import pyranges as pr
# load example datasets
exons, cpg = pr.data.exons(), pr.data.cpg()

# subsetting pyranges is easy
exons["chrY", "-",  15591259:27197945]
# +--------------|-----------|-----------|----------------------------------------|-----------|----------+
# | Chromosome   | Start     | End       | Name                                   | Score     | Strand   |
# | (int8)       | (int32)   | (int32)   | (object)                               | (int64)   | (int8)   |
# |--------------|-----------|-----------|----------------------------------------|-----------|----------|
# | chrY         | 15591393  | 15592550  | NR_047610_exon_27_0_chrY_15591394_r    | 0         | -        |
# | chrY         | 15591393  | 15592550  | NR_047607_exon_29_0_chrY_15591394_r    | 0         | -        |
# | chrY         | 15591393  | 15592550  | NM_001258269_exon_29_0_chrY_15591394_r | 0         | -        |
# | ...          | ...       | ...       | ...                                    | ...       | ...      |
# | chrY         | 26952215  | 26952307  | NM_020364_exon_16_0_chrY_26952216_r    | 0         | -        |
# | chrY         | 27197822  | 27197945  | NM_004678_exon_7_0_chrY_27197823_r     | 0         | -        |
# | chrY         | 27197822  | 27197945  | NM_001002760_exon_7_0_chrY_27197823_r  | 0         | -        |
# +--------------|-----------|-----------|----------------------------------------|-----------|----------+
# PyRanges object has 22 sequences from 1 chromosomes.

# plenty of other ways to subset:
# or exons["chrY", 50:20000], or exons["-"] or exons["chrX", "+"] or ...

# chaining operations a delight

# use the cpg dataset
(cpg
  # to join with the exons dataset
  # and give colnames in exons that exist in cpg the suffix "_xn"
  .join(exons, suffix="_xn")
  # keep only rows with a CpG score over 30
  ('df.CpG > 30')
  .sort()
  # drop rows where the name starts with NR
  ('~df.Name.str.startswith("NR")')
  # remove the strand info
  .unstrand())
```

## Documentation

<https://biocore-ntnu.github.io/pyranges/>

(Might be slightly out of date; watch the CHANGELOG too)

## Install

```bash
pip install pyranges
```

## Paper

Being written here: <https://github.com/endrebak/pyranges-paper>

## TODO

For the future:

*   K-nearest
*   write bam

## Performance

<img src="./graphs/main_paper_annotation_binary.png" />

Comprehensive set of graphs for many types of functions on different datasets are here:

[Time](https://github.com/endrebak/pyranges-paper/blob/master/supplementary_paper/time.md))

[Memory](https://github.com/endrebak/pyranges-paper/blob/master/supplementary_paper/memory.md)

The exact code tested is found [here](https://github.com/endrebak/pyranges-paper/tree/master/supplementaries).
