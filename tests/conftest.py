import pytest

from tests.helpers import assert_df_equal

from pyranges import PyRanges
import pyranges as pr

import pandas as pd

from io import StringIO


@pytest.fixture
def cs():

    cs = pr.load_dataset("chipseq")
    return cs



@pytest.fixture
def bg():

    bg = pr.load_dataset("chipseq_background")
    return bg


@pytest.fixture
def exons():

    df = pd.read_table("tests/exon.txt", sep="\t", header=None)
    df.columns = "Chromosome Start End".split() + list(df.columns[3:])

    return PyRanges(df)


@pytest.fixture
def names():
    return "Chromosome  Start  End  Name Score Strand".split()


@pytest.fixture
def introns():

    df = pd.read_table("tests/intron.txt", sep="\t", header=None)

    print(df.head())
    print(df.shape)
    df.columns = "Chromosome Start End".split() + list(df.columns[3:])
    print(df.columns)

    return PyRanges(df)


@pytest.fixture
def f1():

    df = pd.read_table("tests/f1.bed", sep="\t", header=None, names="Chromosome  Start  End  Name Score Strand".split())

    return PyRanges(df)


@pytest.fixture
def f2():

    df = pd.read_table("tests/f2.bed", sep="\t", header=None, names=names)

    return PyRanges(df)


def assert_df_equal(df1, df2):

    df1 = df1.reset_index(drop=True)
    df2 = df2.reset_index(drop=True)

    return df1.equals(df2)


@pytest.fixture
def chip_10(names):

    df = pd.read_table("tests/chip_10.bed", header=None, names=names)

    gr = PyRanges(df)

    assert gr.stranded

    return gr


@pytest.fixture
def chip_10_plus_one(names):

    df = pd.read_table("tests/chip_10_plus_one.bed", header=None, names=names)

    gr = PyRanges(df)

    assert gr.stranded

    return gr


@pytest.fixture
def input_10(names):

    df = pd.read_table("tests/input_10.bed", header=None, names=names)

    gr = PyRanges(df)

    assert gr.stranded

    return gr



@pytest.fixture
def chip_10_no_strand(names):

    df = pd.read_table("tests/chip_10.bed", header=None, names=names)
    df = df.drop("Strand", 1)

    return PyRanges(df)


@pytest.fixture
def input_10_no_strand(names):

    df = pd.read_table("tests/input_10.bed", header=None, names=names)
    df = df.drop("Strand", 1)

    return PyRanges(df)
