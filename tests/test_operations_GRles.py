import pytest

import pandas as pd
import numpy as np

from io import StringIO

from pyrle import PyRles

import pyranges as pr

@pytest.fixture
def chip():

    c = """Chromosome Start End Strand
chr1 5 7 +
chr1 3 10 -
chr2 1 3 +"""

    return PyRles(pd.read_table(StringIO(c), sep="\s+"))


@pytest.fixture
def background():

    c = """Chromosome Start End Strand
chr1 1 4 +
chr1 2 5 -
chr2 0 1 +"""

    return PyRles(pd.read_table(StringIO(c), sep="\s+"))


@pytest.fixture
def expected_result():

    c = """Chromosome  Start  End  Score
chr1      0    1    0.0
chr1      1    2   -1.0
chr1      2    3   -2.0
chr1      3    4   -1.0
chr1      4    5    0.0
chr1      5    7    2.0
chr1      7   10    1.0
chr2      0    1   -1.0
chr2      1    3    1.0"""

    grles =  PyRles(pd.read_table(StringIO(c), sep="\s+"), value_col="Score")

    return grles


def test_subtraction(chip, background, expected_result):

    result = chip - background

    print("result\n", result, "\n")
    print("expected_result\n", expected_result)
    assert result == expected_result




def test_subtraction_multicpu(chip, background, expected_result):

    result = chip.sub(background, n_jobs=5)

    print("result\n", result, "\n")
    print("expected_result\n", expected_result)
    assert result == expected_result
