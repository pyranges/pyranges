import pytest

from pyranges.pyranges import PyRanges

import pandas as pd

from io import StringIO


def test_instantiation_with_list_like(chip_10_plus_one):

    # mix series, lists and array
    seqnames = chip_10_plus_one.Chromosome.values
    starts = chip_10_plus_one.Start.values
    ends = chip_10_plus_one.End
    strands = chip_10_plus_one.Strand.tolist()

    pr = PyRanges(seqnames=seqnames, starts=starts, ends=ends, strands=strands)


def test_instantiation_with_str_for_strands_seqnames(chip_10_plus_one):

    # mix series, lists and array
    seqnames = "chr1"
    starts = chip_10_plus_one.Start.values
    ends = chip_10_plus_one.End
    strands = "+"

    pr = PyRanges(seqnames=seqnames, starts=starts, ends=ends, strands=strands)
