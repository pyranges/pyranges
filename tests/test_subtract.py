import pytest

import pandas as pd
import numpy as np

from io import StringIO

from pyrle import Rle
from pyrle.methods import _to_ranges, coverage

from pyranges import PyRanges

@pytest.fixture
def simple_rle():

    r = pd.Series([1, 2, 3, 4], dtype=np.int16)
    r1 = Rle(r, r)

    return r1

@pytest.fixture
def simple_rle2():

    r = pd.Series([1, 2, 3, 4], dtype=np.int16)
    r2 = Rle(r * 2, r * 2)

    return r2


@pytest.fixture
def shorter_rle():

    runs = pd.Series([1, 2, 3])
    values = pd.Series([1, 0, 1])
    r = Rle(runs, values)

    return r

@pytest.fixture
def long_rle():

    r = pd.Series([6, 5, 4, 3, 2, 1], dtype=np.int16)
    r2 = Rle(r, r)

    return r2

# numeric-Rle of length 20 with 4 runs
#   Lengths: 2 4 6 8
#   Values : 2 4 6 8
# numeric-Rle of length 20 with 7 runs
#   Lengths:  1  1  1  3  4  2  8
#   Values : -1  0 -2 -1 -2 -6 -8

def test_subtract_simple(simple_rle, simple_rle2):

    result = simple_rle - simple_rle2

    print(result.runs)
    print(result.values)

    expected_runs = np.array([1, 1, 1, 3, 4, 2, 8])
    expected_values = np.array([-1,  0, -2, -1, -2, -6,  -8])

    assert all(np.equal(result.runs, expected_runs))
    assert all(np.equal(result.values, expected_values))



def test_subtract_simple2(simple_rle, simple_rle2):

# numeric-Rle of length 20 with 7 runs
#   Lengths: 1 1 1 3 4 2 8
#   Values : 1 0 2 1 2 6 8

    result2 = simple_rle2 - simple_rle

    print(result2.runs)
    print(result2.values)
    expected_runs2 = [1, 1, 1, 3, 4, 2, 8]
    expected_values2 = [1, 0, 2, 1, 2, 6, 8]
    assert all(np.equal(result2.runs, expected_runs2))
    assert all(np.equal(result2.values, expected_values2))


def test_subtract_advanced(shorter_rle, long_rle):

# numeric-Rle of length 21 with 7 runs
#   Lengths:  1  2  8  4  3  2  1
#   Values : -5 -6 -5 -4 -3 -2 -1
    expected_runs = [1, 2, 8, 4, 3, 2, 1]
    expected_values = [-5, -6, -5, -4, -3, -2, -1]
    result = shorter_rle - long_rle
    print(result.runs)
    print(result.values)

    assert all(np.equal(result.runs, expected_runs))
    assert all(np.equal(result.values, expected_values))

# numeric-Rle of length 21 with 7 runs
#   Lengths: 1 2 8 4 3 2 1
#   Values : 5 6 5 4 3 2 1
def test_subtract_advanced(shorter_rle, long_rle):

    expected_runs = [1, 2, 8, 4, 3, 2, 1]
    expected_values2 = [5, 6, 5, 4, 3, 2, 1]
    result2 = long_rle - shorter_rle

    print("result2.runs", result2.runs)
    print("result2.values", result2.values)

    assert all(np.equal(result2.runs, expected_runs))
    assert all(np.equal(result2.values, expected_values2))


@pytest.fixture
def chip():

    c = """Chromosome Start End Strand
chr2 1 3 +"""

    return coverage(PyRanges(pd.read_table(StringIO(c), sep="\s+")))



@pytest.fixture
def background():

    c = """Chromosome Start End Strand
chr2 0 1 +"""

    return coverage(PyRanges(pd.read_table(StringIO(c), sep="\s+")))


def test_subtract_result_same_start(chip, background):

    print(chip)
    print(background)

    result = chip - background

    print(result)

    assert result == Rle([1, 2], [-1, 1])




@pytest.fixture
def chip_chr1():

    c = """Chromosome Start End Strand
chr1 5 7 +
chr1 3 10 -"""

    return coverage(PyRanges(pd.read_table(StringIO(c), sep="\s+")))


@pytest.fixture
def background_chr1():

    c = """Chromosome Start End Strand
chr1 1 4 +
chr1 2 5 -"""

    return coverage(PyRanges(pd.read_table(StringIO(c), sep="\s+")))



def test_subtract_result_chr1(chip_chr1, background_chr1):

    result = chip_chr1 - background_chr1

    print(result)

    assert len(result.values) == 7
