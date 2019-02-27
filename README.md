# pyranges

[![Build Status](https://travis-ci.org/biocore-ntnu/pyranges.svg?branch=master)](https://travis-ci.org/biocore-ntnu/pyranges) [![hypothesis tested](graphs/hypothesis-tested-brightgreen.svg)](http://hypothesis.readthedocs.io/) [![PyPI version](https://badge.fury.io/py/pyranges.svg)](https://badge.fury.io/py/pyranges)


## Release

PyRanges is slated for an early 2019 release, late February at the latest.

#### Introduction

GenomicRanges for Python.

This library is a thin, transparent wrapper around genomic data contained in
pandas dataframes. This allows for all the wonderful functionality of
bedtools/bedops and/or GenomicRanges, while being able to use the the enormous
universe of Python datascience libraries to manipulate and do computations on
the data.

PyRanges also contains a run-length encoding library for extremely efficient
arithmetic computation of scores associated with genomic intervals.


### Documentation

https://biocore-ntnu.github.io/pyranges/


### Install

```bash
pip install pyranges
```

### Paper

Being written here: https://github.com/endrebak/pyranges-paper

### TODO

For the future:

* Allow writing to different UCSC genome-browser compatible formats such as
  bigwig, bedgraph, histograms, colored bed etc
* Add visualization capabilites?
* Find sequences of ranges with pyfaidx
* K-nearest


#### Performance


<img src="./graphs/main_paper_annotation_binary.png" />

Comprehensive set of graphs for many types of functions on different datasets are here:

[Time](https://github.com/endrebak/pyranges-paper/blob/master/supplementary_paper/time.md)

[Memory](https://github.com/endrebak/pyranges-paper/blob/master/supplementary_paper/memory.md)

The exact code tested is found [here](https://github.com/endrebak/pyranges-paper/tree/master/supplementaries).

### See also

* https://github.com/endrebak/pyrle
* https://github.com/vsbuffalo/BioRanges/tree/master/BioRanges
* https://github.com/daler/pybedtools
* http://bedtools.readthedocs.io/en/latest/
* https://github.com/phaverty/RLEVectors.jl
* https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
* https://bedops.readthedocs.io/en/latest/
