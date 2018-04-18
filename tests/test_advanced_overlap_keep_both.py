import pytest

from pyranges.pyranges import GRanges
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

    c = """Chromosome Start End Name Score Strand Start_chipseq_background End_chipseq_background Name_chipseq_background Score_chipseq_background Strand_chipseq_background
chr8	38747226	38747251	U0	0	-	38747236	38747261	U0	0	+
chr1	226987592	226987617	U0	0	+	226987603	226987628	U0	0	-
chr15	26105515	26105540	U0	0	+	26105493	26105518	U0	0	+"""

    return GRanges(pd.read_table(StringIO(c), sep="\s+"))


def test_advanced_overlap(cs, bg, expected_result_regular_overlap):

    result = cs.overlap_join(bg)

    print("result")
    print(result)
    print("expected_result_regular_overlap")
    print(expected_result_regular_overlap)

    assert expected_result_regular_overlap.df.equals(result.df)

@pytest.fixture
def expected_result_regular_overlap_intersection():

    c = """Chromosome Start End Start_chipseq End_chipseq Name_chipseq Score_chipseq Strand_chipseq Start_chipseq_background End_chipseq_background Name_chipseq_background Score_chipseq_background Strand_chipseq_background
chr1 226987603 226987617 226987592 226987617 U0 0 + 226987603 226987628 U0 0 -
chr8 38747236 38747251 38747226 38747251 U0 0 - 38747236 38747261 U0 0 +
chr15 26105515 26105518 26105515 26105540 U0 0 + 26105493 26105518 U0 0 +"""

    return GRanges(pd.read_table(StringIO(c), sep="\s+"))

def test_advanced_overlap_new_pos_intersect(cs, bg, expected_result_regular_overlap_intersection):

    result = cs.overlap_join(bg, new_pos="intersect")

    print("result")
    print(result.df.to_csv(sep=" "))
    print("expected_result_regular_overlap")
    print(expected_result_regular_overlap_intersection)

    assert expected_result_regular_overlap_intersection.df.equals(result.df)



@pytest.fixture
def expected_result_overlap_same_strand():

    c = """Chromosome Start End Name Score Strand Start_chipseq_background End_chipseq_background Name_chipseq_background Score_chipseq_background Strand_chipseq_background
chr15	26105515	26105540	U0	0	+	26105493	26105518	U0	0	+"""

    df = pd.read_table(StringIO(c), header=0, sep="\s+")
    return GRanges(df)

def test_advanced_overlap_stranded(cs, bg, expected_result_overlap_same_strand):

    result = cs.overlap_join(bg, strandedness="same")

    print("result")
    print(result)
    print("expected result")
    print(expected_result_overlap_same_strand)

    assert expected_result_overlap_same_strand.df.equals(result.df)


@pytest.fixture
def expected_result_overlap_opposite_strand():

    c = """Chromosome Start End Name Score Strand Start_chipseq_background End_chipseq_background Name_chipseq_background Score_chipseq_background Strand_chipseq_background
chr8	38747226	38747251	U0	0	-	38747236	38747261	U0	0	+
chr1	226987592	226987617	U0	0	+	226987603	226987628	U0	0	-"""

    return GRanges(pd.read_table(StringIO(c), sep="\s+"))


def test_advanced_overlap_opposite_stranded(cs, bg, expected_result_overlap_opposite_strand):

    result = cs.overlap_join(bg, strandedness="opposite")

    print(result.df)
    print(expected_result_overlap_opposite_strand)

    assert expected_result_overlap_opposite_strand.df.equals(result.df)


@pytest.fixture
def expected_result_regular_overlap_union():

    c = """Chromosome Start End Start_chipseq End_chipseq Name_chipseq Score_chipseq Strand_chipseq Start_chipseq_background End_chipseq_background Name_chipseq_background Score_chipseq_background Strand_chipseq_background
0 chr1 226987592 226987628 226987592 226987617 U0 0 + 226987603 226987628 U0 0 -
1 chr8 38747226 38747261 38747226 38747251 U0 0 - 38747236 38747261 U0 0 +
2 chr15 26105493 26105540 26105515 26105540 U0 0 + 26105493 26105518 U0 0 +"""

    return GRanges(pd.read_table(StringIO(c), sep="\s+"))

def test_advanced_overlap_new_pos_union(cs, bg, expected_result_regular_overlap_union):

    result = cs.overlap_join(bg, new_pos="union")

    print("result")
    print(result.df.to_csv(sep=" "))
    print("expected_result_regular_overlap")
    print(expected_result_regular_overlap_union)

    assert expected_result_regular_overlap_union.df.equals(result.df)
