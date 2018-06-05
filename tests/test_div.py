
import pytest

import pandas as pd
import numpy as np

from pyrle import Rle

@pytest.fixture
def hypothesis_rle():

    return Rle([66561, 253], [1, np.inf])

@pytest.fixture
def hypothesis_rle2():

    return Rle([66561], [1])


@pytest.fixture
def expected_result_hypothesis():

    return Rle([66561, 253], [1, np.inf])

def test_hypothesis(hypothesis_rle, hypothesis_rle2, expected_result_hypothesis):

    result = hypothesis_rle / hypothesis_rle2

    print("result\n", result)
    print("expected_result\n", expected_result_hypothesis)

    assert result == expected_result_hypothesis

@pytest.fixture
def hypothesis2_rle():

    return Rle([9793796], [627117])

@pytest.fixture
def hypothesis2_rle2():

    return Rle([1553679], [742893])


@pytest.fixture
def expected_result_hypothesis2():

    return Rle([1553679, 8240117], [0.844155, np.inf])

def test_hypothesis2(hypothesis2_rle, hypothesis2_rle2, expected_result_hypothesis2):

    print(hypothesis2_rle, "\n", hypothesis2_rle2)

    result = hypothesis2_rle / hypothesis2_rle2

    print("result\n", result)
    print("expected_result\n", expected_result_hypothesis2)

    assert result == expected_result_hypothesis2


@pytest.fixture
def weird_rle():

    return Rle([10, 20, 30, 40], [1, 2, 3, 4])



@pytest.fixture
def weird_rle2():

    return Rle([1], [1])


@pytest.fixture
def expected_result_weird():

    return Rle([1, 99], [1, 0])

def test_weird(weird_rle, weird_rle2, expected_result_weird):

    result = weird_rle2 / weird_rle

    print(result)

    assert result == expected_result_weird


@pytest.fixture
def expected_result_weird_opposite():

    return Rle([1, 99], [1, np.inf])


def test_weird_opposite(weird_rle, weird_rle2, expected_result_weird_opposite):

    result = weird_rle / weird_rle2

    print(result)

    assert result == expected_result_weird_opposite


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

@pytest.fixture
def long_rle_neg_vals():

    r = pd.Series([6, 5, 4, 3, 2, 1], dtype=np.int16)
    r2 = Rle(r, r * -1)

    return r2

@pytest.fixture
def expected_result_div_simple():

    return Rle([1, 1, 1, 3, 4, 10], [0.5, 1, 0.5, 0.75, 0.66666667, 0.0])


def test_div_simple(simple_rle, simple_rle2, expected_result_div_simple):

    result = simple_rle / simple_rle2

    print(result)
    print(expected_result_div_simple)


    print(result.values)
    print(expected_result_div_simple.values)

    assert result == expected_result_div_simple



@pytest.fixture
def expected_result_div_simple_opposite():

    return Rle([1, 1, 1, 3, 4, 10], [2.0, 1.0, 2.0, 1.333333, 1.5, np.inf])


def test_div_simple_opposite(simple_rle, simple_rle2, expected_result_div_simple_opposite):

    result = simple_rle2 / simple_rle

    print(result.values)
    print(result.values)

    print(result.runs)

    assert result == expected_result_div_simple_opposite


@pytest.fixture
def expected_result_div_advanced():

    return Rle([1, 2, 3, 15], [0.16666667, 0., 0.16666667, 0.])


def test_div_advanced(shorter_rle, long_rle, expected_result_div_advanced):

    result = shorter_rle / long_rle
    print("result runs", result.runs)
    print("result values", result.values)

    assert result == expected_result_div_advanced




@pytest.fixture
def expected_result_div_advanced_opposite():

    return Rle([1, 2, 3, 15], [6, np.inf, 6, np.inf])

def test_div_advanced_opposite(shorter_rle, long_rle, expected_result_div_advanced_opposite):

    result = long_rle / shorter_rle
    print("result runs", result.runs)
    print("result values", result.values)

    assert result == expected_result_div_advanced_opposite



@pytest.fixture
def expected_result_div_advanced_opposite():

    return Rle([1, 2, 3, 15], [6, np.inf, 6, np.inf])



@pytest.fixture
def expected_result_div_advanced_neg_vals():

    return Rle([1, 2, 3, 15], [-0.16666667, 0., -0.16666667, 0.])


def test_div_advanced_neg_vals(shorter_rle, long_rle_neg_vals, expected_result_div_advanced_neg_vals):

    result = shorter_rle / long_rle_neg_vals
    print("result runs", result.runs)
    print("result values", result.values)

    assert result == expected_result_div_advanced_neg_vals


@pytest.fixture
def expected_result_div_advanced_opposite_neg_vals():

    return Rle([1, 2, 3, 15], [-6, -np.inf, -6, -np.inf])


def test_div_advanced_opposite_neg_vals(shorter_rle, long_rle_neg_vals, expected_result_div_advanced_opposite_neg_vals):

    result = long_rle_neg_vals / shorter_rle
    print("result runs", result.runs)
    print("result values", result.values)

    assert result == expected_result_div_advanced_opposite_neg_vals
