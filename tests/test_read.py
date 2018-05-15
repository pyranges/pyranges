import pytest

from pyranges.pyranges import PyRanges
import pyranges as pr

import pandas as pd

from io import StringIO


def test_read_bam():

    pr.read_bam("tests/test_data/test_sorted.bam")


def test_read_gtf():

    gr = pr.read_gtf("tests/test_data/ensembl.gtf")

    assert list(gr.df.columns[:4]) == "Chromosome Start End Strand".split()


def test_read_bed():

    pr.read_bed("pyranges/example_data/chipseq.bed")
