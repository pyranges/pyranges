
import pytest
from tests.helpers import assert_df_equal

from pyranges.pyranges import PyRanges
import pyranges as pr

import pandas as pd

from io import StringIO


@pytest.fixture
def expected_result_unstranded():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b Distance
0 chr1 3 6 h 0 + 6.0 7.0 f 0.0 - 0
1 chr1 5 7 h 0 - 6.0 7.0 f 0.0 - 0"""

    return PyRanges(pd.read_table(StringIO(c), sep=" "))


def test_nearest_bed_unstranded(f1, f2, expected_result_unstranded):

    result = f1.nearest(f2, strandedness=False, suffix="_b", how="next_nonoverlapping")

    print(result.df.to_csv(sep=" "))

    assert result.df.empty


@pytest.fixture
def expected_result_opposite_stranded():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b Distance
0 chr1 9939 10138 HWI-ST216_313:3:2301:15791:16298 1 + 10149 10348 HWI-ST216:427:D29R1ACXX:2:2306:7654:33038 1 - 11
1 chr1 9953 10152 HWI-ST216_313:3:1305:6975:102491 1 + 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 9806
2 chr1 10024 10223 HWI-ST216_313:3:2201:5209:155139 1 + 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 9735
3 chr1 10246 10445 HWI-ST216_313:3:1207:4315:142177 1 + 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 9513
4 chr1 9916 10115 HWI-ST216_313:3:1203:10227:6568 1 - 10280 10479 HWI-ST216:427:D29R1ACXX:2:2103:15965:69480 1 + 165
5 chr1 9951 10150 HWI-ST216_313:3:2205:20086:33508 1 - 10280 10479 HWI-ST216:427:D29R1ACXX:2:2103:15965:69480 1 + 130
6 chr1 9978 10177 HWI-ST216_313:3:1204:5599:113305 1 - 10280 10479 HWI-ST216:427:D29R1ACXX:2:2103:15965:69480 1 + 103
7 chr1 10001 10200 HWI-ST216_313:3:1102:14019:151362 1 - 10280 10479 HWI-ST216:427:D29R1ACXX:2:2103:15965:69480 1 + 80
8 chr1 10127 10326 HWI-ST216_313:3:2207:7406:122346 1 - 16056 16255 HWI-ST216:427:D29R1ACXX:2:1102:7604:12113 1 + 5730
9 chr1 10241 10440 HWI-ST216_313:3:1302:4516:156396 1 - 16056 16255 HWI-ST216:427:D29R1ACXX:2:1102:7604:12113 1 + 5616"""

    return PyRanges(pd.read_table(StringIO(c), sep=" "))


def test_nearest_bed_opposite_stranded(chip_10_plus_one, input_10, expected_result_opposite_stranded):

    result = chip_10_plus_one.nearest(input_10, strandedness="opposite", suffix="_b", how="next_nonoverlapping")

    print("result", result.df.to_csv(sep=" "))

    assert assert_df_equal(result.df, expected_result_opposite_stranded.df)


@pytest.fixture
def expected_result_same_stranded():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b Distance
0 chr1 9916 10115 HWI-ST216_313:3:1203:10227:6568 1 - 10149 10348 HWI-ST216:427:D29R1ACXX:2:2306:7654:33038 1 - 34
1 chr1 9939 10138 HWI-ST216_313:3:2301:15791:16298 1 + 10280 10479 HWI-ST216:427:D29R1ACXX:2:2103:15965:69480 1 + 142
2 chr1 9951 10150 HWI-ST216_313:3:2205:20086:33508 1 - 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 9808
3 chr1 9953 10152 HWI-ST216_313:3:1305:6975:102491 1 + 10280 10479 HWI-ST216:427:D29R1ACXX:2:2103:15965:69480 1 + 128
4 chr1 9978 10177 HWI-ST216_313:3:1204:5599:113305 1 - 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 9781
5 chr1 10001 10200 HWI-ST216_313:3:1102:14019:151362 1 - 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 9758
6 chr1 10024 10223 HWI-ST216_313:3:2201:5209:155139 1 + 10280 10479 HWI-ST216:427:D29R1ACXX:2:2103:15965:69480 1 + 57
7 chr1 10127 10326 HWI-ST216_313:3:2207:7406:122346 1 - 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 9632
8 chr1 10241 10440 HWI-ST216_313:3:1302:4516:156396 1 - 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 9518
9 chr1 10246 10445 HWI-ST216_313:3:1207:4315:142177 1 + 16056 16255 HWI-ST216:427:D29R1ACXX:2:1102:7604:12113 1 + 5611"""

    return PyRanges(pd.read_table(StringIO(c), sep=" "))


def test_nearest_bed_same_stranded(chip_10_plus_one, input_10, expected_result_same_stranded):

    result = chip_10_plus_one.nearest(input_10, strandedness="same", suffix="_b", how="next_nonoverlapping")

    print("result", result.df.to_csv(sep=" "))

    assert assert_df_equal(result.df, expected_result_same_stranded.df)
