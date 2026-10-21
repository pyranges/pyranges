# pyranges0

## Deprecation note!

This repo currently hosts the "version 0" of pyranges. A redesigned, faster version 1 is available at https://github.com/pyranges/pyranges1
Version 1 will soon replace version 0. We encourage existing v0 users to migrate to v1 (see guide at https://pyranges1.readthedocs.io/en/latest/migration_guide.html), and new users to directly adopt v1. Read about v1 in our preprint at https://doi.org/10.64898/2025.12.11.693639

## Recent changes

Version 0.1.6 supports pandas 3 (pandas 1.5 and 2.x are still supported). With pandas 3, `to_rle` and coverage need pyrle >= 0.0.45. See [CHANGELOG.txt](CHANGELOG.txt) for all changes.

Install coverage support with `pip install 'pyranges[coverage]'`.

## Introduction

PyRanges is a Python library specifically designed for efficient and intuitive manipulation of genomics data,
particularly genomic intervals (like genes, genomic features, or reads).
The library is optimized for fast querying and manipulation of genomic annotations.

*"Finally ... This was what Python badly needed for years."* - Heng Li

## Documentation

The pyranges0 documentation, including installation instructions, API, tutorial, and how-to pages, is available at <https://pyranges0.readthedocs.io/en/latest/>

## Features

  - fast
  - memory-efficient
  - featureful
  - pythonic/pandastic
  - supports chaining with a terse syntax
  - uses Pandas DataFrames, so the whole Python data science stack works on PyRanges

## Paper/Cite

Stovner EB, Sætrom P (2020) PyRanges: efficient comparison of genomic intervals in Python. *Bioinformatics 36(3):918-919*  http://dx.doi.org/10.1093/bioinformatics/btz615

## Supporting pyranges

  - most importantly, cite pyranges if you use it. It is the main metric funding sources care about.
  - use pyranges in Stack Overflow/biostars questions and answers
  - star the repo (possibly important for github visibility and as a proxy for project popularity)
  - if you are a business using pyranges, please give to one of the charities listed at https://www.givewell.org/

## Asking for help

If you encounter bugs, or the documentation is not enough and you cannot accomplish a specific task of interest, or if you'd like new features implemented, open an issue on GitHub: https://github.com/pyranges/pyranges0/issues

## Contributing to pyranges

Pyranges accepts code contributions in the form of pull requests. For details, visit [https://pyranges0.readthedocs.io/en/latest/developer_guide.html](https://pyranges0.readthedocs.io/en/latest/developer_guide.html)

