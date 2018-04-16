
import pytest

from pyranges.pyranges import GRanges

import pandas as pd

from io import StringIO


@pytest.fixture
def simple_gr1():

    c = """Chromosome Start End Strand Score
chr1 3 6 + 5
chr1 5 7 - 7
chr1 8 9 + 1"""

    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return GRanges(df)



@pytest.fixture
def simple_gr2():

    c = """Chromosome Start End Strand Score
chr1 1 2 + 1
chr1 6 7 - 2"""
    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return GRanges(df)


@pytest.fixture
def expected_result_subtract_simple_granges():

    c = """Chromosome Start End Strand Score
chr1	3	6	+ 5
chr1	8	9	+ 1"""
    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return GRanges(df)



def test_overlap_invert_simple_granges(simple_gr1, simple_gr2, expected_result_subtract_simple_granges):

    result = simple_gr1.overlap(simple_gr2, invert=True, strandedness=False)

    print(result)

    assert result.df.equals(expected_result_subtract_simple_granges.df)


@pytest.fixture
def expected_result_overlap_same_strand_simple_granges():

    c = """Chromosome Start End Strand Score
chr1 5 7 - 7"""
    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return GRanges(df)


def test_overlap_same_strand_simple_granges(simple_gr1, simple_gr2, expected_result_overlap_same_strand_simple_granges):

    result = simple_gr1.overlap(simple_gr2, strandedness="same")

    print(result)

    assert result.df.equals(expected_result_overlap_same_strand_simple_granges.df)



def test_overlap_opposite_strand_simple_granges(simple_gr1, simple_gr2):

    result = simple_gr1.overlap(simple_gr2, strandedness="opposite")

    print(result)

    assert result.df.empty


def test_default_overlap_simple_granges(simple_gr1, simple_gr2, expected_result_overlap_same_strand_simple_granges):

    print(simple_gr1)
    print(simple_gr2)

    result = simple_gr1.overlap(simple_gr2)

    print(result.df.dtypes)
    print(expected_result_overlap_same_strand_simple_granges.df.dtypes)

    print(result.df.index)
    print(expected_result_overlap_same_strand_simple_granges.df.index)

    assert result.df.equals(expected_result_overlap_same_strand_simple_granges.df)



# @pytest.fixture
# def expected_result_intersect_simple_granges():

#     c = """Chromosome Start End
# chr1 3 6
# chr1 5 7"""

#     df = pd.read_table(StringIO(c), sep="\s+", header=0)
#     return GRanges(df)


# def test_intersect_simple_granges(simple_gr1, simple_gr2, expected_result_intersect_simple_granges):

#     result = simple_gr1 | simple_gr2
#     print(result)

#     assert result.df.equals(expected_result_intersect_simple_granges.df)



# def test_str_rep_grange(gr):

#     result = str(gr)

#     print(result)

#     expected_result = """+--------------+----------+----------+--------------------+
# | Chromosome   | Start    | End      | Score              |
# +--------------+----------+----------+--------------------|
# | chr1         | 11323785 | 11617177 | 0.86217008797654   |
# | chr1         | 12645605 | 13926923 | 0.9348914858096831 |
# | chr1         | 14750216 | 15119039 | 0.9459459459459459 |
# | ...          | ...      | ...      | ...                |
# | chr1         | 46974417 | 47547260 | 0.91701244813278   |
# | chr1         | 51234626 | 51418366 | 0.83453237410072   |
# | chr1         | 56953387 | 57050584 | 0.91304347826087   |
# +--------------+----------+----------+--------------------+
# GRanges object with 15 sequences from 1 chromosomes."""

#     assert result == expected_result


# def test_overlap_introns_exons(introns, exons):

#     igr = introns | exons

#     # bedtools agrees
#     assert igr.df.empty



# def test_overlap_introns_exons(introns, exons):

#     igr = introns | exons

#     # bedtools agrees
#     assert igr.df.empty


# @pytest.fixture
# def expected_result_f1_f2():

#     c = """Chromosome  Start  End  Name Score Strand
# 1      1    2  .  0  +"""

#     df = pd.read_table(StringIO(c), sep="\s+")

#     return df
