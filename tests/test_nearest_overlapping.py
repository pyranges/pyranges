

import pytest
from tests.helpers import assert_df_equal

from pyranges.pyranges import PyRanges
import pyranges as pr

import pandas as pd

from io import StringIO


@pytest.fixture()
def expected_result_counterexample1():

    c = """chr1	0	1	+	0	1	+	0
chr1	0	1	+	0	1	+	0
chr1	0	1	+	0	1	+	0
chr1	0	1	+	0	1	+	0
chr1	0	1	+	0	1	+	0
chr1	0	1	+	0	1	+	0
chr1	0	1	+	0	1	+	0
chr1	0	1	+	0	1	+	0
chr1	0	5726225	+	0	1	+	0
chr1	3538885	3832293	+	5726225	5726228	+	1893933
chr1	4426346	9655531	+	5726225	5726228	+	0"""

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

     result = hyp1.nearest(hyp2, overlap=True)

     print(result)
     # assert 0
     # print(PyRanges(result.df.tail()))
     assert_df_equal(result.df, expected_result_counterexample1.df)


def test_hypothesis_counterexample_nearest_next(hyp1, hyp2, expected_result_counterexample1):

     print(hyp1)
     print(hyp2)

     result = hyp1.nearest(hyp2, overlap=True, how="next")

     print(result)
     # print(PyRanges(result.df.tail()))
     assert_df_equal(result.df, expected_result_counterexample1.df)


@pytest.fixture()
def expected_result_nearest_previous():

    c = """chr1 0 1 + 0 1 + 0
chr1 0 1 + 0 1 + 0
chr1 0 1 + 0 1 + 0
chr1 0 1 + 0 1 + 0
chr1 0 1 + 0 1 + 0
chr1 0 1 + 0 1 + 0
chr1 0 1 + 0 1 + 0
chr1 0 1 + 0 1 + 0
chr1 0 5726225 + 0 1 + 0
chr1 3538885 3832293 + 0 1 + 3538885
chr1 4426346 9655531 + 5726225 5726228 + 0"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))


def test_hypothesis_counterexample_nearest_previous(hyp1, hyp2, expected_result_nearest_previous):

     print(hyp1)
     print(hyp2)

     result = hyp1.nearest(hyp2, overlap=True, how="previous")

     print(result)
     # print(PyRanges(result.df.tail()))
     assert_df_equal(result.df, expected_result_nearest_previous.df)
