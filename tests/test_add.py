import pytest

import pandas as pd
import numpy as np

from pyrle import Rle


@pytest.fixture
def rle():

    r = pd.Series([1], dtype=np.int16)
    r1 = Rle(r, r)

    return r1

@pytest.fixture
def rle2():

    r = pd.Series([2, 1], dtype=np.int16)
    v = pd.Series([0, 1], dtype=np.int16)
    r2 = Rle(r, v)

    return r2

def test_hypothesis_counterexample(rle, rle2):

    result = rle + rle2

    expected_runs = [1, 1, 1]
    expected_values = [1, 0, 1]

    print(result.runs)
    print(result.values)

    assert list(result.runs) == expected_runs
    assert list(result.values) == expected_values

    result2 = rle2 + rle

    print(result2.runs)
    print(result2.values)

    assert list(result2.runs) == expected_runs
    assert list(result2.values) == expected_values

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

# numeric-Rle of length 21 with 8 runs
#   Lengths: 1 2 3 5 4 3 2 1
#   Values : 7 6 7 5 4 3 2 1

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


def test_add_simple(simple_rle, simple_rle2):

    result = simple_rle + simple_rle2

    expected_runs = np.array([1, 1, 1, 3, 4, 2, 8])
    expected_values = np.array([3, 4, 6, 7, 10, 6, 8])

    print(result.runs)
    print(result.values)

    assert all(np.equal(result.runs, expected_runs))
    assert all(np.equal(result.values, expected_values))

    result2 = simple_rle2 + simple_rle

    print(result2.runs)
    print(result2.values)

    print(result2)

    assert all(np.equal(result2.runs, expected_runs))
    assert all(np.equal(result2.values, expected_values))


def test_add_advanced(shorter_rle, long_rle):

    result = shorter_rle + long_rle
    print("result runs", result.runs)
    print("result values", result.values)

    result2 = long_rle + shorter_rle

    print("result runs", result2.runs)
    print("result values", result2.values)

    expected_runs = np.array([1, 2, 3, 5, 4, 3, 2, 1])
    expected_values = np.array([7.0, 6.0, 7.0, 5.0, 4.0, 3.0, 2.0, 1.0])

    assert all(np.equal(result.runs, expected_runs))
    assert all(np.equal(result.values, expected_values))

    assert all(np.equal(result2.runs, expected_runs))
    assert all(np.equal(result2.values, expected_values))
