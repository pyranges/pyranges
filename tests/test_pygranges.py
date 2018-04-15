import pytest

from pyranges.pyranges import GRanges

import pandas as pd

from io import StringIO


@pytest.fixture
def simple_gr1():

    c = """Chromosome Start End
    chr1 3 6
    chr1 5 7
    chr1 8 9"""

    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return GRanges(df)



@pytest.fixture
def simple_gr2():

    c = """Chromosome Start End
chr1 1 2
chr1 6 6"""
    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return GRanges(df)


@pytest.fixture
def gr():
    c = """Chromosome Start End Score
chr1	11323785	11617177	0.86217008797654
chr1	12645605	13926923	0.934891485809683
chr1	14750216	15119039	0.945945945945946
chr1	18102157	19080189	0.895174708818636
chr1	29491029	30934636	0.892526250772082
chr1	33716472	35395979	0.911901081916538
chr1	36712462	37685238	0.95655951346655
chr1	37838094	38031209	0.944206008583691
chr1	38272060	39078902	0.940932642487047
chr1	42197905	42388954	0.759162303664921
chr1	42583551	42913605	0.801857585139319
chr1	45281207	45416320	0.885714285714286
chr1	46974417	47547260	0.91701244813278
chr1	51234626	51418366	0.83453237410072
chr1	56953387	57050584	0.91304347826087"""

    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return GRanges(df)



@pytest.fixture
def expected_result_subtract_simple_granges():

    c = """Chromosome Start End
chr1 8 9"""
    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return GRanges(df)


def test_subtract_simple_granges(simple_gr1, simple_gr2, expected_result_subtract_simple_granges):

    print(simple_gr1[9:])
    result = simple_gr1 - simple_gr2

    print(result)

    assert result.df.equals(expected_result_subtract_simple_granges.df)



@pytest.fixture
def expected_result_intersect_simple_granges():

    c = """Chromosome Start End
chr1 3 6
chr1 5 7"""

    df = pd.read_table(StringIO(c), sep="\s+", header=0)
    return GRanges(df)


def test_intersect_simple_granges(simple_gr1, simple_gr2, expected_result_intersect_simple_granges):

    result = simple_gr1 | simple_gr2
    print(result)

    assert result.df.equals(expected_result_intersect_simple_granges.df)



def test_str_rep_grange(gr):

    result = str(gr)

    print(result)

    expected_result = """+-----+--------------+----------+----------+--------------------+
|     | Chromosome   | Start    | End      | Score              |
|-----+--------------+----------+----------+--------------------|
| 0   | chr1         | 11323785 | 11617177 | 0.86217008797654   |
| 1   | chr1         | 12645605 | 13926923 | 0.9348914858096831 |
| 2   | chr1         | 14750216 | 15119039 | 0.9459459459459459 |
| ... | ...          | ...      | ...      | ...                |
| 12  | chr1         | 46974417 | 47547260 | 0.91701244813278   |
| 13  | chr1         | 51234626 | 51418366 | 0.83453237410072   |
| 14  | chr1         | 56953387 | 57050584 | 0.91304347826087   |
+-----+--------------+----------+----------+--------------------+
GRanges object with 15 sequences from 1 chromosomes."""

    assert result == expected_result


def test_intersection_introns_exons(introns, exons):

    igr = introns | exons

    # bedtools agrees
    assert igr.df.empty



def test_intersection_introns_exons(introns, exons):

    igr = introns | exons

    # bedtools agrees
    assert igr.df.empty


@pytest.fixture
def expected_result_f1_f2():

    c = """Chromosome  Start  End  Name Score Strand
1      1    2  .  0  +"""

    df = pd.read_table(StringIO(c), sep="\s+")

    return df
