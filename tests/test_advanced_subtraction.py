from io import StringIO

from tests.helpers import assert_df_equal
from pyranges.pyranges import PyRanges

import pandas as pd

import pytest

@pytest.fixture
def expected_result_unstranded(names):

    c = """chr1    9916    9988    HWI-ST216_313:3:1203:10227:6568 1       -
chr1    9939    9988    HWI-ST216_313:3:2301:15791:16298        1       +
chr1    9951    9988    HWI-ST216_313:3:2205:20086:33508        1       -
chr1    9953    9988    HWI-ST216_313:3:1305:6975:102491        1       +
chr1    9978    9988    HWI-ST216_313:3:1204:5599:113305        1       -"""

    df = pd.read_table(StringIO(c), sep="\s+", names=names, header=None)
    return PyRanges(df)


@pytest.fixture
def expected_result_self_unstranded(names):

    c = """chr1    9916    9988    HWI-ST216_313:3:1203:10227:6568 1
chr1    9939    9988    HWI-ST216_313:3:2301:15791:16298        1
chr1    9951    9988    HWI-ST216_313:3:2205:20086:33508        1
chr1    9953    9988    HWI-ST216_313:3:1305:6975:102491        1
chr1    9978    9988    HWI-ST216_313:3:1204:5599:113305        1"""

    df = pd.read_table(StringIO(c), sep="\s+", names=names[:-1], header=None)
    return PyRanges(df)


@pytest.fixture
def expected_result_no_strand_plus_one(names):

    c = """chr1    9916    9988    HWI-ST216_313:3:1203:10227:6568 1       -
chr1    9939    9988    HWI-ST216_313:3:2301:15791:16298        1       +
chr1    9951    9988    HWI-ST216_313:3:2205:20086:33508        1       -
chr1    9953    9988    HWI-ST216_313:3:1305:6975:102491        1       +
chr1    9978    9988    HWI-ST216_313:3:1204:5599:113305        1       -
chr1    10348   10445   HWI-ST216_313:3:1207:4315:142177        1       +"""

    df = pd.read_table(StringIO(c), sep="\s+", names=names, header=None)
    return PyRanges(df)



@pytest.fixture
def expected_result_same_strand(names):

    c = """chr1    9916    9988    HWI-ST216_313:3:1203:10227:6568 1       -
chr1    9939    10073   HWI-ST216_313:3:2301:15791:16298        1       +
chr1    9951    9988    HWI-ST216_313:3:2205:20086:33508        1       -
chr1    9953    10073   HWI-ST216_313:3:1305:6975:102491        1       +
chr1    9978    9988    HWI-ST216_313:3:1204:5599:113305        1       -
chr1    10024   10073   HWI-ST216_313:3:2201:5209:155139        1       +
chr1    10348   10440   HWI-ST216_313:3:1302:4516:156396        1       -
chr1    10272   10280   HWI-ST216_313:3:1207:4315:142177        1       +"""

    df = pd.read_table(StringIO(c), sep="\s+", names=names, header=None)
    return PyRanges(df)

@pytest.fixture
def expected_result_same_strand_plus_one(names):

    c = """chr1    9916    9988    HWI-ST216_313:3:1203:10227:6568 1       -
chr1    9939    10073   HWI-ST216_313:3:2301:15791:16298        1       +
chr1    9951    9988    HWI-ST216_313:3:2205:20086:33508        1       -
chr1    9953    10073   HWI-ST216_313:3:1305:6975:102491        1       +
chr1    9978    9988    HWI-ST216_313:3:1204:5599:113305        1       -
chr1    10024   10073   HWI-ST216_313:3:2201:5209:155139        1       +
chr1    10348   10440   HWI-ST216_313:3:1302:4516:156396        1       -
chr1    10272   10280   HWI-ST216_313:3:1207:4315:142177        1       +
chr1    10348   10445   HWI-ST216_313:3:1207:4315:142177        1       +"""

    df = pd.read_table(StringIO(c), sep="\s+", names=names, header=None)
    return PyRanges(df)

@pytest.fixture
def expected_result_opposite_strand(names):

    c = """chr1    9916    10073   HWI-ST216_313:3:1203:10227:6568 1       -
chr1    9939    9988    HWI-ST216_313:3:2301:15791:16298        1       +
chr1    9951    10073   HWI-ST216_313:3:2205:20086:33508        1       -
chr1    9953    9988    HWI-ST216_313:3:1305:6975:102491        1       +
chr1    9978    10073   HWI-ST216_313:3:1204:5599:113305        1       -
chr1    10001   10073   HWI-ST216_313:3:1102:14019:151362       1       -
chr1    10272   10280   HWI-ST216_313:3:2207:7406:122346        1       -
chr1    10272   10280   HWI-ST216_313:3:1302:4516:156396        1       -
chr1    10348   10445   HWI-ST216_313:3:1207:4315:142177        1       +"""


    df = pd.read_table(StringIO(c), sep="\s+", names=names, header=None)
    return PyRanges(df)

@pytest.fixture
def expected_result_opposite_strand_plus_one(names):

    c = """chr1    9916    10073   HWI-ST216_313:3:1203:10227:6568 1       -
chr1    9939    9988    HWI-ST216_313:3:2301:15791:16298        1       +
chr1    9951    10073   HWI-ST216_313:3:2205:20086:33508        1       -
chr1    9953    9988    HWI-ST216_313:3:1305:6975:102491        1       +
chr1    9978    10073   HWI-ST216_313:3:1204:5599:113305        1       -
chr1    10001   10073   HWI-ST216_313:3:1102:14019:151362       1       -
chr1    10272   10280   HWI-ST216_313:3:2207:7406:122346        1       -
chr1    10272   10280   HWI-ST216_313:3:1302:4516:156396        1       -
chr1    10348   10445   HWI-ST216_313:3:1207:4315:142177        1       +
chr1    110246  110445  HWI-ST216_313:3:1207:4315:142177        1       +"""


    df = pd.read_table(StringIO(c), sep="\s+", names=names, header=None)
    return PyRanges(df)



def test_advanced_subtraction_unstranded(chip_10, input_10, expected_result_unstranded):

    result = chip_10.subtraction(input_10)

    print(result)

    assert_df_equal(result.df, expected_result_unstranded.df)


def test_advanced_subtraction_same_strand(chip_10, input_10, expected_result_same_strand):

    result = chip_10.subtraction(input_10, strandedness="same")

    assert_df_equal(result.df, expected_result_same_strand.df)


def test_advanced_subtraction_opposite_strand(chip_10, input_10, expected_result_opposite_strand):

    result = chip_10.subtraction(input_10, strandedness="opposite")

    assert_df_equal(result.df, expected_result_opposite_strand.df)

def test_advanced_subtraction_no_strand_plus_one(chip_10_plus_one, input_10, expected_result_no_strand_plus_one):

    result = chip_10_plus_one.subtraction(input_10)

    assert_df_equal(result.df, expected_result_no_strand_plus_one.df)

def test_advanced_subtraction_same_strand_plus_one(chip_10_plus_one, input_10, expected_result_same_strand_plus_one):

    result = chip_10_plus_one.subtraction(input_10, strandedness="same")

    assert_df_equal(result.df, expected_result_same_strand_plus_one.df)

def test_advanced_subtraction_opposite_strand_plus_one(chip_10_plus_one, input_10, expected_result_opposite_strand_plus_one):

    result = chip_10_plus_one.subtraction(input_10, strandedness="opposite")

    assert_df_equal(result.df, expected_result_opposite_strand_plus_one.df)


def test_advanced_subtraction_unstranded_self_no_strand(chip_10_no_strand, input_10, expected_result_self_unstranded):

    result = chip_10_no_strand.subtraction(input_10)

    assert_df_equal(result.df, expected_result_self_unstranded.df)


def test_advanced_subtraction_unstranded_other_no_strand(chip_10, input_10_no_strand, expected_result_unstranded):

    result = chip_10.subtraction(input_10_no_strand)

    assert_df_equal(result.df, expected_result_unstranded.df)
