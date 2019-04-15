# pyranges

[![Coverage Status](https://img.shields.io/coveralls/github/biocore-ntnu/pyranges.svg)](https://coveralls.io/github/biocore-ntnu/pyranges?branch=master) [![Codacy Badge](https://api.codacy.com/project/badge/Grade/b61a53346d764a8d8f0ab2a6afd7b100)](https://www.codacy.com/app/endrebak/pyranges?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=biocore-ntnu/pyranges&amp;utm_campaign=Badge_Grade) [![Build Status](https://travis-ci.org/biocore-ntnu/pyranges.svg?branch=master)](https://travis-ci.org/biocore-ntnu/pyranges) [![hypothesis tested](graphs/hypothesis-tested-brightgreen.svg)](http://hypothesis.readthedocs.io/) [![PyPI version](https://badge.fury.io/py/pyranges.svg)](https://badge.fury.io/py/pyranges) [![MIT](https://img.shields.io/pypi/l/pyranges.svg?color=green)](https://opensource.org/licenses/MIT)

## Release

PyRanges is in a beta state.

## Introduction

GenomicRanges and genomic Rle-objects for Python.

## Quick example

```python
import pyranges as pr
# load example datasets
exons, cpg = pr.data.exons(), pr.data.cpg()

# subsetting pyranges is easy
exons["chrY", "-",  15591259:27197945]
# +--------------|-----------|-----------|----------------------------------------|-----------|--------------+
# | Chromosome   | Start     | End       | Name                                   | Score     | Strand       |
# | (category)   | (int64)   | (int64)   | (object)                               | (int64)   | (category)   |
# |--------------|-----------|-----------|----------------------------------------|-----------|--------------|
# | chrY         | 15591393  | 15592550  | NR_047610_exon_27_0_chrY_15591394_r    | 0         | -            |
# | chrY         | 15591393  | 15592550  | NR_047607_exon_29_0_chrY_15591394_r    | 0         | -            |
# | chrY         | 15591393  | 15592550  | NM_001258269_exon_29_0_chrY_15591394_r | 0         | -            |
# | ...          | ...       | ...       | ...                                    | ...       | ...          |
# | chrY         | 26952215  | 26952307  | NM_020364_exon_16_0_chrY_26952216_r    | 0         | -            |
# | chrY         | 27197822  | 27197945  | NM_004678_exon_7_0_chrY_27197823_r     | 0         | -            |
# | chrY         | 27197822  | 27197945  | NM_001002760_exon_7_0_chrY_27197823_r  | 0         | -            |
# +--------------|-----------|-----------|----------------------------------------|-----------|--------------+
# PyRanges object has 22 sequences from 1 chromosomes.

# the API allows for easy and terse chaining
(cpg # use the cpg dataset
  .join(exons, suffix="_xn") # join with exons, use suffix _xn for duplicate cols
  (lambda df: df.CpG > 30) # keep only rows with a CpG score over 30
  .sort()
  (lambda df: ~df.Name.str.startswith("NR")) # drop rows where the name starts with NR
  .unstrand()) # remove the strand info
# +--------------|-----------|-----------|-----------|------------|-----------|----------------------------------------|-----------+
# | Chromosome   | Start     | End       | CpG       | Start_xn   | End_xn    | Name                                   | Score     |
# | (category)   | (int64)   | (int64)   | (int64)   | (int64)    | (int64)   | (object)                               | (int64)   |
# |--------------|-----------|-----------|-----------|------------|-----------|----------------------------------------|-----------|
# | chrX         | 584563    | 585326    | 66        | 585078     | 585337    | NM_000451_exon_0_0_chrX_585079_f       | 0         |
# | chrX         | 13587648  | 13588221  | 49        | 13587693   | 13588054  | NM_001167890_exon_0_0_chrX_13587694_f  | 0         |
# | chrX         | 17755053  | 17756648  | 117       | 17755568   | 17755800  | NM_001037535_exon_0_0_chrX_17755569_f  | 0         |
# | ...          | ...       | ...       | ...       | ...        | ...       | ...                                    | ...       |
# | chrY         | 16941822  | 16942188  | 32        | 16941609   | 16942399  | NM_014893_exon_4_0_chrY_16941610_f     | 0         |
# | chrY         | 241398    | 245968    | 310       | 244667     | 245252    | NM_013239_exon_0_0_chrY_244668_r       | 0         |
# | chrY         | 15591259  | 15591720  | 33        | 15591393   | 15592550  | NM_001258269_exon_29_0_chrY_15591394_r | 0         |
# +--------------|-----------|-----------|-----------|------------|-----------|----------------------------------------|-----------+
# PyRanges object has 57 sequences from 2 chromosomes.
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
