import pytest

import pandas as pd

from pyranges import PyRanges


@pytest.fixture
def names():
    return "Chromosome  Start  End  Name Score Strand".split()

@pytest.fixture
def chip_10(names):

    df = pd.read_table("tests/chip_10.bed", header=None, names=names)

    gr = PyRanges(df)

    assert gr.stranded

    return gr

@pytest.fixture
def f1(names):

    df = pd.read_table("tests/f1.bed", sep="\t", header=None, names="Chromosome  Start  End  Name Score Strand".split())

    return PyRanges(df)


@pytest.fixture
def f2(names):

    df = pd.read_table("tests/f2.bed", sep="\t", header=None, names=names)

    return PyRanges(df)
