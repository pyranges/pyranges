import pytest

from tests.helpers import assert_df_equal
from pyranges.pyranges import PyRanges
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
def expected_result_regular_overlap():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b
chr8	38747226	38747251	U0	0	-	38747236	38747261	U0	0	+
chr1	226987592	226987617	U0	0	+	226987603	226987628	U0	0	-
chr15	26105515	26105540	U0	0	+	26105493	26105518	U0	0	+"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+"))


def test_advanced_overlap(cs, bg, expected_result_regular_overlap):

    result = cs.join(bg)

    print("result")
    print(result)
    print("expected_result_regular_overlap")
    print(expected_result_regular_overlap)

    assert assert_df_equal(expected_result_regular_overlap.df, result.df)

@pytest.fixture
def expected_result_regular_overlap_intersection():

    c = """Chromosome Start End Start_a End_a Name_a Score_a Strand Start_b End_b Name_b Score_b Strand_b
chr1 226987603 226987617 226987592 226987617 U0 0 + 226987603 226987628 U0 0 -
chr8 38747236 38747251 38747226 38747251 U0 0 - 38747236 38747261 U0 0 +
chr15 26105515 26105518 26105515 26105540 U0 0 + 26105493 26105518 U0 0 +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+"))

def test_advanced_overlap_new_pos_intersect(cs, bg, expected_result_regular_overlap_intersection):

    result = cs.join(bg, new_pos="intersection")

    assert assert_df_equal(result.df, expected_result_regular_overlap_intersection.df)
    # assert expected_result_regular_overlap_intersection.df.equals(result.df)



@pytest.fixture
def expected_result_overlap_same_strand():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b
chr15	26105515	26105540	U0	0	+	26105493	26105518	U0	0	+"""

    df = pd.read_table(StringIO(c), header=0, sep="\s+")
    return PyRanges(df)

def test_advanced_overlap_stranded(cs, bg, expected_result_overlap_same_strand):

    result = cs.join(bg, strandedness="same")

    print("result")
    print(result)
    print("expected result")
    print(expected_result_overlap_same_strand)

    assert assert_df_equal(expected_result_overlap_same_strand.df, result.df)


@pytest.fixture
def expected_result_overlap_opposite_strand():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b
chr8	38747226	38747251	U0	0	-	38747236	38747261	U0	0	+
chr1	226987592	226987617	U0	0	+	226987603	226987628	U0	0	-"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+"))


def test_advanced_overlap_opposite_stranded(cs, bg, expected_result_overlap_opposite_strand):

    result = cs.join(bg, strandedness="opposite")

    print(result.df)
    print(expected_result_overlap_opposite_strand)

    assert assert_df_equal(expected_result_overlap_opposite_strand.df, result.df)



@pytest.fixture
def expected_result_overlap_opposite_strand_suffixes():

    c = """Chromosome Start End Name Score Strand Starthiya Endhiya Namehiya Scorehiya Strandhiya
chr8	38747226	38747251	U0	0	-	38747236	38747261	U0	0	+
chr1	226987592	226987617	U0	0	+	226987603	226987628	U0	0	-"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+"))


def test_advanced_overlap_opposite_stranded_suffixes(cs, bg, expected_result_overlap_opposite_strand_suffixes):

    result = cs.join(bg, strandedness="opposite", suffixes=["", "hiya"])

    print(result.df)
    print(expected_result_overlap_opposite_strand_suffixes)

    assert assert_df_equal(expected_result_overlap_opposite_strand_suffixes.df, result.df)


@pytest.fixture
def expected_result_regular_overlap_union():

    c = """Chromosome Start End Start_a End_a Name_a Score_a Strand Start_b End_b Name_b Score_b Strand_b
0 chr1 226987592 226987628 226987592 226987617 U0 0 + 226987603 226987628 U0 0 -
1 chr8 38747226 38747261 38747226 38747251 U0 0 - 38747236 38747261 U0 0 +
2 chr15 26105493 26105540 26105515 26105540 U0 0 + 26105493 26105518 U0 0 +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+"))

def test_advanced_overlap_new_pos_union(cs, bg, expected_result_regular_overlap_union):

    result = cs.join(bg, new_pos="union")

    print("result")
    print(result.df.to_csv(sep=" "))
    print("expected_result_regular_overlap")
    print(expected_result_regular_overlap_union)

    assert assert_df_equal(expected_result_regular_overlap_union.df, result.df)
