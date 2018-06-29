
import pytest
from tests.helpers import assert_df_equal

from pyranges.pyranges import PyRanges
import pyranges as pr

import pandas as pd

from io import StringIO

@pytest.fixture
def expected_result_unstranded():

 c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b Distance
chr1 9916 10115 HWI-ST216_313:3:1203:10227:6568 1 - 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
chr1 9939 10138 HWI-ST216_313:3:2301:15791:16298 1 + 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
chr1 9951 10150 HWI-ST216_313:3:2205:20086:33508 1 - 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
chr1 9953 10152 HWI-ST216_313:3:1305:6975:102491 1 + 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
chr1 9978 10177 HWI-ST216_313:3:1204:5599:113305 1 - 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
chr1 10001 10200 HWI-ST216_313:3:1102:14019:151362 1 - 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
chr1 10024 10223 HWI-ST216_313:3:2201:5209:155139 1 + 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
chr1 10127 10326 HWI-ST216_313:3:2207:7406:122346 1 - 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
chr1 10241 10440 HWI-ST216_313:3:1302:4516:156396 1 - 10073 10272 HWI-ST216:427:D29R1ACXX:2:2302:14161:85418 1 + 0
chr1 10246 10445 HWI-ST216_313:3:1207:4315:142177 1 + 10073 10272 HWI-ST216:427:D29R1ACXX:2:2302:14161:85418 1 + 0
chr1 110246 110445 HWI-ST216_313:3:1207:4315:142177 1 + 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 90090"""

 return PyRanges(pd.read_table(StringIO(c), sep=" "))


def test_nearest_bed_unstranded(chip_10_plus_one, input_10, expected_result_unstranded):

    result = chip_10_plus_one.nearest(input_10, strandedness=False, suffix="_b")

    # print(result.df.to_csv(sep=" "))

    assert_df_equal(result.df, expected_result_unstranded.df)



@pytest.fixture
def expected_result_same_stranded():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b Distance
chr1 9939 10138 HWI-ST216_313:3:2301:15791:16298 1 + 10073 10272 HWI-ST216:427:D29R1ACXX:2:2302:14161:85418 1 + 0
chr1 9953 10152 HWI-ST216_313:3:1305:6975:102491 1 + 10073 10272 HWI-ST216:427:D29R1ACXX:2:2302:14161:85418 1 + 0
chr1 10024 10223 HWI-ST216_313:3:2201:5209:155139 1 + 10073 10272 HWI-ST216:427:D29R1ACXX:2:2302:14161:85418 1 + 0
chr1 10246 10445 HWI-ST216_313:3:1207:4315:142177 1 + 10073 10272 HWI-ST216:427:D29R1ACXX:2:2302:14161:85418 1 + 0
chr1 110246 110445 HWI-ST216_313:3:1207:4315:142177 1 + 16109 16308 HWI-ST216:427:D29R1ACXX:2:2110:12286:25379 1 + 93939
chr1 9916 10115 HWI-ST216_313:3:1203:10227:6568 1 - 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
chr1 9951 10150 HWI-ST216_313:3:2205:20086:33508 1 - 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
chr1 9978 10177 HWI-ST216_313:3:1204:5599:113305 1 - 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
chr1 10001 10200 HWI-ST216_313:3:1102:14019:151362 1 - 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
chr1 10127 10326 HWI-ST216_313:3:2207:7406:122346 1 - 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
chr1 10241 10440 HWI-ST216_313:3:1302:4516:156396 1 - 10079 10278 HWI-ST216:427:D29R1ACXX:2:1314:10333:38924 1 - 0"""

    return PyRanges(pd.read_table(StringIO(c), sep=" "))


def test_nearest_bed_same_stranded(chip_10_plus_one, input_10, expected_result_same_stranded):

    result = chip_10_plus_one.nearest(input_10, strandedness="same", suffix="_b")

    print(result.df.to_csv(sep=" "))

    assert_df_equal(result.df, expected_result_same_stranded.df)


@pytest.fixture
def expected_result_opposite_stranded():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b Distance
chr1 9939 10138 HWI-ST216_313:3:2301:15791:16298 1 + 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
chr1 9953 10152 HWI-ST216_313:3:1305:6975:102491 1 + 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
chr1 10024 10223 HWI-ST216_313:3:2201:5209:155139 1 + 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
chr1 10246 10445 HWI-ST216_313:3:1207:4315:142177 1 + 10079 10278 HWI-ST216:427:D29R1ACXX:2:1314:10333:38924 1 - 0
chr1 110246 110445 HWI-ST216_313:3:1207:4315:142177 1 + 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 90090
chr1 9916 10115 HWI-ST216_313:3:1203:10227:6568 1 - 10073 10272 HWI-ST216:427:D29R1ACXX:2:2302:14161:85418 1 + 0
chr1 9951 10150 HWI-ST216_313:3:2205:20086:33508 1 - 10073 10272 HWI-ST216:427:D29R1ACXX:2:2302:14161:85418 1 + 0
chr1 9978 10177 HWI-ST216_313:3:1204:5599:113305 1 - 10073 10272 HWI-ST216:427:D29R1ACXX:2:2302:14161:85418 1 + 0
chr1 10001 10200 HWI-ST216_313:3:1102:14019:151362 1 - 10073 10272 HWI-ST216:427:D29R1ACXX:2:2302:14161:85418 1 + 0
chr1 10127 10326 HWI-ST216_313:3:2207:7406:122346 1 - 10073 10272 HWI-ST216:427:D29R1ACXX:2:2302:14161:85418 1 + 0
chr1 10241 10440 HWI-ST216_313:3:1302:4516:156396 1 - 10073 10272 HWI-ST216:427:D29R1ACXX:2:2302:14161:85418 1 + 0"""

    return PyRanges(pd.read_table(StringIO(c), sep=" "))


def test_nearest_bed_opposite_stranded(chip_10_plus_one, input_10, expected_result_opposite_stranded):

    result = chip_10_plus_one.nearest(input_10, strandedness="opposite", suffix="_b", overlap=True)

    print(result.df.to_csv(sep=" "))

    assert_df_equal(result.df, expected_result_opposite_stranded.df)


@pytest.fixture
def expected_result_other_no_strand():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Distance
0 chr1 9916 10115 HWI-ST216_313:3:1203:10227:6568 1 - 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 0
1 chr1 9939 10138 HWI-ST216_313:3:2301:15791:16298 1 + 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 0
2 chr1 9951 10150 HWI-ST216_313:3:2205:20086:33508 1 - 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 0
3 chr1 9953 10152 HWI-ST216_313:3:1305:6975:102491 1 + 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 0
4 chr1 9978 10177 HWI-ST216_313:3:1204:5599:113305 1 - 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 0
5 chr1 10001 10200 HWI-ST216_313:3:1102:14019:151362 1 - 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 0
6 chr1 10024 10223 HWI-ST216_313:3:2201:5209:155139 1 + 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 0
7 chr1 10127 10326 HWI-ST216_313:3:2207:7406:122346 1 - 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 0
8 chr1 10241 10440 HWI-ST216_313:3:1302:4516:156396 1 - 10073 10272 HWI-ST216:427:D29R1ACXX:2:2302:14161:85418 1 0
9 chr1 10246 10445 HWI-ST216_313:3:1207:4315:142177 1 + 10073 10272 HWI-ST216:427:D29R1ACXX:2:2302:14161:85418 1 0
10 chr1 110246 110445 HWI-ST216_313:3:1207:4315:142177 1 + 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 90090"""

    return PyRanges(pd.read_table(StringIO(c), sep=" "))

def test_nearest_other_unstranded(chip_10_plus_one, input_10_no_strand, expected_result_other_no_strand):

    result = chip_10_plus_one.nearest(input_10_no_strand, strandedness=False, suffix="_b")

    print(result.df.to_csv(sep=" "))

    assert_df_equal(result.df, expected_result_other_no_strand.df)


@pytest.fixture
def expected_result_self_no_strand():

    c = """Chromosome Start End Name Score Start_b End_b Name_b Score_b Strand Distance
0 chr1 9916 10115 HWI-ST216_313:3:1203:10227:6568 1 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
1 chr1 9939 10138 HWI-ST216_313:3:2301:15791:16298 1 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
2 chr1 9951 10150 HWI-ST216_313:3:2205:20086:33508 1 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
3 chr1 9953 10152 HWI-ST216_313:3:1305:6975:102491 1 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
4 chr1 9978 10177 HWI-ST216_313:3:1204:5599:113305 1 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
5 chr1 10001 10200 HWI-ST216_313:3:1102:14019:151362 1 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
6 chr1 10024 10223 HWI-ST216_313:3:2201:5209:155139 1 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
7 chr1 10127 10326 HWI-ST216_313:3:2207:7406:122346 1 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 0
8 chr1 10241 10440 HWI-ST216_313:3:1302:4516:156396 1 10073 10272 HWI-ST216:427:D29R1ACXX:2:2302:14161:85418 1 + 0
9 chr1 10246 10445 HWI-ST216_313:3:1207:4315:142177 1 10073 10272 HWI-ST216:427:D29R1ACXX:2:2302:14161:85418 1 + 0"""

    return PyRanges(pd.read_table(StringIO(c), sep=" "))


def test_nearest_self_unstranded(chip_10_no_strand, input_10, expected_result_self_no_strand):

    result = chip_10_no_strand.nearest(input_10, strandedness=False, suffix="_b")

    print(result.df.to_csv(sep=" "))

    assert_df_equal(result.df, expected_result_self_no_strand.df)



def test_nearest_self_unstranded_strandedness(chip_10_no_strand, input_10, expected_result_self_no_strand):

    with pytest.raises(AssertionError):
        result = chip_10_no_strand.nearest(input_10, strandedness="opposite", suffix="_b")


# @pytest.fixture()
# def hyp1():

#     c = """chr1	0	1	+
# chr1	0	1	+
# chr1	4426346	9655531	+
# chr1	0	1	+
# chr1	0	5726225	+
# chr1	0	1	+
# chr1	0	1	+
# chr1	0	1	+
# chr1	3538885	3832293	+
# chr1	0	1	+
# chr1	0	1	+"""

#     return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))

# @pytest.fixture()
# def hyp2():

#     c = """chr1	0	1	+
# chr1	0	1	+
# chr1	0	1	+
# chr1	0	1	+
# chr1	0	1	+
# chr1	0	1	+
# chr1	0	1	+
# chr1	0	1	+
# chr1	0	1	+
# chr1	5767167	5923191	+
# chr1	0	1	+"""

#     return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


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

     result = hyp1.nearest(hyp2)

     print(result)
     # 0
     # print(PyRanges(result.df.tail()))
     assert_df_equal(result.df, expected_result_counterexample1.df)




@pytest.fixture()
def hyp3():

    c = """chr1      3    4      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp4():

    c = """chr1      0    3      +
chr1      3    4      +
chr1      1    2      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))

@pytest.fixture()
def expected_result_counterexample2():

    c = """chr1	3	4	+	3	4	+	0"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))


def test_hypothesis_counterexample2(hyp3, hyp4, expected_result_counterexample2):

     print(hyp3)
     print(hyp4)

     result = hyp3.nearest(hyp4)

     print(result)
     # print(PyRanges(result.df.tail()))
     assert_df_equal(result.df, expected_result_counterexample2.df)




@pytest.fixture()
def hyp5():

    c = """chr1      0    1      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp6():

    c = """chr1      1    2      +
chr1      2    3      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))

@pytest.fixture()
def expected_result_counterexample3():

    c = """chr1	0	1	+	1	2	+	1"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))


def test_hypothesis_counterexample3(hyp5, hyp6, expected_result_counterexample3):

     print(hyp5)
     print(hyp6)

     result = hyp5.nearest(hyp6)

     print(result)
     # print(PyRanges(result.df.tail()))
     assert_df_equal(result.df, expected_result_counterexample3.df)



@pytest.fixture()
def hyp7():

    c = """chr1      3    4      +
chr1      0    1      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp8():

    c = """chr1      1    2      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))

@pytest.fixture()
def expected_result_counterexample4():

    c = """chr1  0  1  +  1   2   +         1
chr1  3  4  +  1   2   +         2"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))


def test_hypothesis_counterexample4(hyp7, hyp8, expected_result_counterexample4):

     print(hyp7)
     print(hyp8)

     result = hyp7.nearest(hyp8)

     print(result)
     # print(PyRanges(result.df.tail()))
     assert_df_equal(result.df, expected_result_counterexample4.df)



@pytest.fixture()
def hyp9():

    c = """chr2      0    1      +
chr2      2    3      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp10():

    c = """chr2      1    2      +
chr2      1    3      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))

@pytest.fixture()
def expected_result_counterexample5():

    c = """chr2  0  1  +  1   2   +         1
chr2  2  3  +  1   3   +         0"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))


def test_hypothesis_counterexample5(hyp9, hyp10, expected_result_counterexample5):

     print(hyp9)
     print(hyp10)

     result = hyp9.nearest(hyp10)

     print(result)
     # print(PyRanges(result.df.tail()))
     assert_df_equal(result.df, expected_result_counterexample5.df)
