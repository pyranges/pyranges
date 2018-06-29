

import pytest
from tests.helpers import assert_df_equal

from pyranges.pyranges import PyRanges
import pyranges as pr

import pandas as pd

from io import StringIO


@pytest.fixture()
def expected_result_counterexample1():

    c = """chr1 0 1 + 5726225 5726228 + 5726225
chr1 0 1 + 5726225 5726228 + 5726225
chr1 0 1 + 5726225 5726228 + 5726225
chr1 0 1 + 5726225 5726228 + 5726225
chr1 0 1 + 5726225 5726228 + 5726225
chr1 0 1 + 5726225 5726228 + 5726225
chr1 0 1 + 5726225 5726228 + 5726225
chr1 0 1 + 5726225 5726228 + 5726225
chr1 0 5726225 + 5726225 5726228 + 1
chr1 3538885 3832293 + 5726225 5726228 + 1893933
chr1 4426346 9655531 + 0 1 + 4426346"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))


@pytest.fixture()
def hyp1():

    c = """chr1        0        1      +
chr1        0        1      +
chr1        0        1      +
chr1        0        1      +
chr1        0        1      +
chr1        0        1      +
chr1        0        1      +
chr1        0        1      +
chr1        0  5726225      +
chr1  3538885  3832293      +
chr1  4426346  9655531      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp2():

    c = """chr1        0        1      +
chr1        0        1      +
chr1        0        1      +
chr1        0        1      +
chr1        0        1      +
chr1        0        1      +
chr1        0        1      +
chr1        0        1      +
chr1        0        1      +
chr1        0        1      +
chr1  5726225  5726228      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))

def test_hypothesis_counterexample(hyp1, hyp2, expected_result_counterexample1):

     print(hyp1)
     print(hyp2)

     result = hyp1.nearest(hyp2, overlap=False)

     print(result)
     # 0
     # print(PyRanges(result.df.tail()))
     assert_df_equal(result.df, expected_result_counterexample1.df)
