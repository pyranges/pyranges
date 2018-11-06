
import pytest
from tests.helpers import assert_dfs_equal

from pyranges.pyranges import PyRanges
import pyranges as pr

import pandas as pd

from io import StringIO

@pytest.fixture
def expected_result_simple_cluster():

    c = """chr1	3	6	+
chr1	8	9	+
chr1	5	7	-"""

    return pr.PyRanges(pd.read_table(StringIO(c), header=None, names="Chromosome Start End Strand".split()))


def test_simple_cluster(f1, expected_result_simple_cluster):

    result = f1.cluster(strand=True)

    print(result)

    print(expected_result_simple_cluster)

    assert_dfs_equal(result, expected_result_simple_cluster)


@pytest.fixture
def expected_result_advanced():

    c = """chrY	7046809	7046834
chrY	7194340	7194365
chrY	7405376	7405401
chrY	7463444	7463469
chrY	7701983	7702008
chrY	7761026	7761051
chrY	8010951	8010976
chrY	8316773	8316798
chrY	11942770	11942795
chrY	12930373	12930398
chrY	13216614	13216639
chrY	13517892	13517917
chrY	14774053	14774078
chrY	15224235	15224260
chrY	15548022	15548047
chrY	16045242	16045267
chrY	16495497	16495522
chrY	21559181	21559206
chrY	21707662	21707687
chrY	21751211	21751236
chrY	21910706	21910731
chrY	22054002	22054027
chrY	22210637	22210662"""


    return pr.PyRanges(pd.read_table(StringIO(c), header=None, names="Chromosome Start End".split()))


def test_advanced_cluster(cs, expected_result_advanced):

    chrY = cs["chrY"]

    result = chrY.cluster()

    print(result.values)

    assert_dfs_equal(result, expected_result_advanced)

# chrY	7046809	7046834
# chrY	7194340	7194365
# chrY	7405376	7405401
# chrY	7463444	7463469
# chrY	7701983	7702008
# chrY	7761026	7761051
# chrY	8010951	8010976
# chrY	8316773	8316798
# chrY	11942770	11942795
# chrY	12930373	12930398
# chrY	13216614	13216639
# chrY	13517892	13517917
# chrY	14774053	14774078
# chrY	15224235	15224260
# chrY	15548022	15548047
# chrY	16045242	16045267
# chrY	16495497	16495522
# chrY	21559181	21559206
# chrY	21707662	21707687
# chrY	21751211	21751236
# chrY	21910706	21910731
# chrY	22054002	22054027
# chrY	22210637	22210662
