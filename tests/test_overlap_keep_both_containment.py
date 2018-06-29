
import pytest

from tests.helpers import assert_df_equal
from pyranges.pyranges import PyRanges
import pyranges as pr

import pandas as pd

from io import StringIO



def test_advanced_overlap_containment(cs, bg):

    result = cs.join(bg, how="containment")

    assert result.df.empty


@pytest.fixture
def contained():

    c = """chr1	2	4	U0	0	+
chr8	7	8	U0	0	-
chr15	11	13	U0	0	+"""

    df = pd.read_table(StringIO(c), header=None, names="Chromosome Start End Name Score Strand".split(), sep="\s+")

    return PyRanges(df)




@pytest.fixture
def container():

    c = """chr1	1 9	U0	0	+
chr8	6 9	U0	0	+
chr15	12	13	U0	0	+"""

    df = pd.read_table(StringIO(c), header=None, names="Chromosome Start End Name Score Strand".split(), sep="\s+")

    return PyRanges(df)


@pytest.fixture
def expected_result_unstranded():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b
0 chr1 2 4 U0 0 + 1 9 U0 0 +
1 chr8 7 8 U0 0 - 6 9 U0 0 +"""

    df = pd.read_table(StringIO(c), header=0, sep="\s+")

    return PyRanges(df)


def test_simple_containment(contained, container, expected_result_unstranded):

    result = contained.join(container, how="containment")
    print(result.df.to_csv(sep=" "))

    assert_df_equal(result.df, expected_result_unstranded.df)

@pytest.fixture
def expected_result_opposite_strand():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b
chr8 7 8 U0 0 - 6 9 U0 0 +"""

    df = pd.read_table(StringIO(c), header=0, sep="\s+")

    return PyRanges(df)



def test_simple_containment_opposite_strand(contained, container, expected_result_opposite_strand):

    result = contained.join(container, strandedness="opposite", how="containment")
    print(result.df.to_csv(sep=" "))

    assert_df_equal(result.df, expected_result_opposite_strand.df)



@pytest.fixture
def expected_result_same_strand():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b
chr1 2 4 U0 0 + 1 9 U0 0 +"""

    df = pd.read_table(StringIO(c), header=0, sep="\s+")

    return PyRanges(df)



def test_simple_containment_same_strand(contained, container, expected_result_same_strand):

    result = contained.join(container, strandedness="same", how="containment")
    print(result.df.to_csv(sep=" "))

    assert_df_equal(result.df, expected_result_same_strand.df)
