import pytest

from tests.helpers import assert_df_equal

from pyranges import PyRanges
import pyranges as pr

import pandas as pd
import numpy as np

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
def gtf():

    return pr.read_gtf("tests/test_data/ensembl.gtf")



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
def f1(names):

    df = pd.read_table("tests/f1.bed", sep="\t", header=None, names="Chromosome  Start  End  Name Score Strand".split())

    return PyRanges(df)


@pytest.fixture
def f2(names):

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

##### HYPOTHESIS

import hypothesis.strategies as st
from hypothesis import given, settings, reproduce_failure, unlimited, HealthCheck, seed
from hypothesis.extra.pandas import data_frames, columns, range_indexes, column, indexes

from hypothesis.extra.numpy import arrays


positions = st.integers(min_value=0, max_value=int(1e7))
lengths = st.integers(min_value=1, max_value=int(1e7))
small_lengths = st.integers(min_value=1, max_value=int(1e4))
strands = st.sampled_from("+ -".split())
names = st.text("abcdefghijklmnopqrstuvxyz", min_size=1)
scores = st.integers(min_value=0, max_value=256)

better_df_minsize = 1
better_dfs_min = data_frames(index=indexes(dtype=np.int64, min_size=better_df_minsize, unique=True),
                             columns=[column("Chromosome", cs),
                                      column("Start", elements=lengths),
                                      column("End", elements=small_lengths),
                                      # column("Name", elements=names),
                                      # column("Score", elements=scores),
                                      column("Strand", strands)])

@st.composite
def dfs_min(draw):
    df = draw(better_dfs_min)
    df.loc[:, "End"] += df.Start

    df.insert(3, "Name", "a")
    df.insert(4, "Score", 0)
    gr = PyRanges(df)

    # do not sort like this, use pyranges sort
    # np.random.seed(draw(st.integers(min_value=0, max_value=int(1e6))))
    # gr.df = df.reindex(np.random.permutation(df.index.values))

    return gr
