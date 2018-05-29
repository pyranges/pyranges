

import pytest
from tests.helpers import assert_df_equal

from pyranges.pyranges import PyRanges
import pyranges as pr

import pandas as pd

from io import StringIO


@pytest.fixture
def expected_result_unstranded():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b Distance
0 chr1 3 6 h 0 + 1 2 f 0 + 1
1 chr1 5 7 h 0 - 1 2 f 0 + 3
2 chr1 8 9 h 0 + 6 7 f 0 - 1"""

    return PyRanges(pd.read_table(StringIO(c), sep=" "))


def test_nearest_bed_unstranded(f1, f2, expected_result_unstranded):

    result = f1.nearest(f2, strandedness=False, suffix="_b", how="previous_nonoverlapping")

    print(result.df.to_csv(sep=" "))

    assert assert_df_equal(result.df, expected_result_unstranded.df)


@pytest.fixture
def expected_result_opposite_stranded():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b Distance
0 chr1 10246 10445 HWI-ST216_313:3:1207:4315:142177 1 + 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 59
1 chr1 110246 110445 HWI-ST216_313:3:1207:4315:142177 1 + 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 90089"""

    return PyRanges(pd.read_table(StringIO(c), sep=" "))


def test_nearest_bed_opposite_stranded(chip_10_plus_one, input_10, expected_result_opposite_stranded):

    result = chip_10_plus_one.nearest(input_10, strandedness="opposite", suffix="_b", how="previous_nonoverlapping")

    print("result", result.df.to_csv(sep=" "))

    assert assert_df_equal(result.df, expected_result_opposite_stranded.df)


@pytest.fixture
def expected_result_same_stranded():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b Distance
0 chr1 10241 10440 HWI-ST216_313:3:1302:4516:156396 1 - 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 54
1 chr1 110246 110445 HWI-ST216_313:3:1207:4315:142177 1 + 16109 16308 HWI-ST216:427:D29R1ACXX:2:2110:12286:25379 1 + 93938"""

    return PyRanges(pd.read_table(StringIO(c), sep=" "))


def test_nearest_bed_same_stranded(chip_10_plus_one, input_10, expected_result_same_stranded):

    result = chip_10_plus_one.nearest(input_10, strandedness="same", suffix="_b", how="previous_nonoverlapping")

    print("result", result.df.to_csv(sep=" "))

    assert assert_df_equal(result.df, expected_result_same_stranded.df)
