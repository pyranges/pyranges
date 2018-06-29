
import pytest
from tests.helpers import assert_df_equal

from pyranges.pyranges import PyRanges
import pyranges as pr

import pandas as pd

from io import StringIO


@pytest.fixture
def expected_result_unstranded():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b Distance
0 chr1 3 6 h 0 + 6 7 f 0 - 1"""

    return PyRanges(pd.read_table(StringIO(c), sep=" "))


def test_nearest_bed_unstranded(f1, f2, expected_result_unstranded):

    result = f1.nearest(f2, strandedness=False, suffix="_b", how="next", overlap=False)


    assert_df_equal(result.df, expected_result_unstranded.df)


@pytest.fixture
def expected_result_opposite_stranded():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b Distance
0 chr1 9939 10138 HWI-ST216_313:3:2301:15791:16298 1 + 10149 10348 HWI-ST216:427:D29R1ACXX:2:2306:7654:33038 1 - 12
1 chr1 9953 10152 HWI-ST216_313:3:1305:6975:102491 1 + 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 9807
2 chr1 10024 10223 HWI-ST216_313:3:2201:5209:155139 1 + 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 9736
3 chr1 10246 10445 HWI-ST216_313:3:1207:4315:142177 1 + 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 9514
4 chr1 9916 10115 HWI-ST216_313:3:1203:10227:6568 1 - 10280 10479 HWI-ST216:427:D29R1ACXX:2:2103:15965:69480 1 + 166
5 chr1 9951 10150 HWI-ST216_313:3:2205:20086:33508 1 - 10280 10479 HWI-ST216:427:D29R1ACXX:2:2103:15965:69480 1 + 131
6 chr1 9978 10177 HWI-ST216_313:3:1204:5599:113305 1 - 10280 10479 HWI-ST216:427:D29R1ACXX:2:2103:15965:69480 1 + 104
7 chr1 10001 10200 HWI-ST216_313:3:1102:14019:151362 1 - 10280 10479 HWI-ST216:427:D29R1ACXX:2:2103:15965:69480 1 + 81
8 chr1 10127 10326 HWI-ST216_313:3:2207:7406:122346 1 - 16056 16255 HWI-ST216:427:D29R1ACXX:2:1102:7604:12113 1 + 5731
9 chr1 10241 10440 HWI-ST216_313:3:1302:4516:156396 1 - 16056 16255 HWI-ST216:427:D29R1ACXX:2:1102:7604:12113 1 + 5617"""

    gr = PyRanges(pd.read_table(StringIO(c), sep=" "))

    return gr


def test_nearest_bed_opposite_stranded(chip_10_plus_one, input_10, expected_result_opposite_stranded):

    result = chip_10_plus_one.nearest(input_10, strandedness="opposite", suffix="_b", how="next", overlap=False)

    print("result", result.df.to_csv(sep=" "))

    assert_df_equal(result.df, expected_result_opposite_stranded.df)


@pytest.fixture
def expected_result_same_stranded():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b Distance
0 chr1 9916 10115 HWI-ST216_313:3:1203:10227:6568 1 - 10149 10348 HWI-ST216:427:D29R1ACXX:2:2306:7654:33038 1 - 35
1 chr1 9939 10138 HWI-ST216_313:3:2301:15791:16298 1 + 10280 10479 HWI-ST216:427:D29R1ACXX:2:2103:15965:69480 1 + 143
2 chr1 9951 10150 HWI-ST216_313:3:2205:20086:33508 1 - 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 9809
3 chr1 9953 10152 HWI-ST216_313:3:1305:6975:102491 1 + 10280 10479 HWI-ST216:427:D29R1ACXX:2:2103:15965:69480 1 + 129
4 chr1 9978 10177 HWI-ST216_313:3:1204:5599:113305 1 - 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 9782
5 chr1 10001 10200 HWI-ST216_313:3:1102:14019:151362 1 - 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 9759
6 chr1 10024 10223 HWI-ST216_313:3:2201:5209:155139 1 + 10280 10479 HWI-ST216:427:D29R1ACXX:2:2103:15965:69480 1 + 58
7 chr1 10127 10326 HWI-ST216_313:3:2207:7406:122346 1 - 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 9633
8 chr1 10241 10440 HWI-ST216_313:3:1302:4516:156396 1 - 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 9519
9 chr1 10246 10445 HWI-ST216_313:3:1207:4315:142177 1 + 16056 16255 HWI-ST216:427:D29R1ACXX:2:1102:7604:12113 1 + 5612"""

    return PyRanges(pd.read_table(StringIO(c), sep=" "))


def test_nearest_bed_same_stranded(chip_10_plus_one, input_10, expected_result_same_stranded):

    result = chip_10_plus_one.nearest(input_10, strandedness="same", suffix="_b", how="next", overlap=False)

    print("result", result.df.to_csv(sep=" "))

    assert_df_equal(result.df, expected_result_same_stranded.df)



@pytest.fixture()
def expected_result_counterexample1():

    c = """chr1  1  2  +  2   3   +         1"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))


@pytest.fixture()
def hyp1():

    c = """chr1	0	3	+
chr1	1	2	+"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp2():

    c = """chr1	2	3	+"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))

def test_hypothesis_counterexample(hyp1, hyp2, expected_result_counterexample1):

     print(hyp1)
     print(hyp2)

     result = hyp1.nearest(hyp2, how="next", overlap=False)

     print(result)
     # assert 0
     # print(PyRanges(result.df.tail()))
     assert_df_equal(result.df, expected_result_counterexample1.df)



@pytest.fixture()
def expected_result_counterexample2():

    c = """chr1  1  2  +  2   3   +         1"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))


@pytest.fixture()
def hyp3():

    c = """chr1	0	4	+
chr1	1	2	+"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp4():

    c = """chr1	2	3	+"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))

def test_hypothesis_counterexample2(hyp3, hyp4, expected_result_counterexample2):

     print(hyp3)
     print(hyp4)

     result = hyp3.nearest(hyp4, how="next", overlap=False)

     print(result)
     # assert 0
     # print(PyRanges(result.df.tail()))
     assert_df_equal(result.df, expected_result_counterexample2.df)


@pytest.fixture()
def expected_result_counterexample3():

    c = """chr1 0 1 + 2 3 + 2
chr1 0 2 + 2 3 + 1
chr1 1 2 + 2 3 + 1"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))



@pytest.fixture()
def hyp5():

    c = """chr1	0	4	+
chr1	0	1	+
chr1	0	2	+
chr1	1	2	+"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp6():

    c = """chr1	2	3	+"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))

def test_hypothesis_counterexample3(hyp5, hyp6, expected_result_counterexample3):

     print(hyp5)
     print(hyp6)

     result = hyp5.nearest(hyp6, how="next", overlap=False)

     print(result)
     # assert 0
     # print(PyRanges(result.df.tail()))
     assert_df_equal(result.df, expected_result_counterexample3.df)



@pytest.fixture()
def hyp7():

    c = """chr1	0	1	+"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp8():

    c = """chr1      0    1      +
chr1      1    2      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def expected_result_counterexample4():

    c = """chr1  0  1  +  1   2   +         1"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))

def test_hypothesis_counterexample4(hyp7, hyp8, expected_result_counterexample4):

     print(hyp7)
     print(hyp8)

     result = hyp7.nearest(hyp8, how="next", overlap=False)

     print(result.df)
     # assert 0
     # print(PyRanges(result.df.tail()))

     assert_df_equal(result.df, expected_result_counterexample4.df)



@pytest.fixture()
def minus_chip():

    c = """chr1	9916	10115	HWI-ST216_313:3:1203:10227:6568	1	-
chr1	9951	10150	HWI-ST216_313:3:2205:20086:33508	1	-
chr1	9978	10177	HWI-ST216_313:3:1204:5599:113305	1	-
chr1	10001	10200	HWI-ST216_313:3:1102:14019:151362	1	-
chr1	10127	10326	HWI-ST216_313:3:2207:7406:122346	1	-
chr1	10241	10440	HWI-ST216_313:3:1302:4516:156396	1	-"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Name Score Strand".split()))


@pytest.fixture()
def minus_input():

    c = """chr1	9988	10187	HWI-ST216:427:D29R1ACXX:2:1205:6095:16532	1	-
chr1	10079	10278	HWI-ST216:427:D29R1ACXX:2:1314:10333:38924	1	-
chr1	10082	10281	HWI-ST216:427:D29R1ACXX:2:2213:10285:5151	1	-
chr1	10149	10348	HWI-ST216:427:D29R1ACXX:2:2306:7654:33038	1	-
chr1	19958	20157	HWI-ST216:427:D29R1ACXX:2:1313:6283:67310	1	-"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Name Score Strand".split()))


@pytest.fixture()
def expected_result_minus():

    c = """chr1 9916 10115 HWI-ST216_313:3:1203:10227:6568 1 - 10149 10348 HWI-ST216:427:D29R1ACXX:2:2306:7654:33038 1 - 35
chr1 9951 10150 HWI-ST216_313:3:2205:20086:33508 1 - 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 9809
chr1 9978 10177 HWI-ST216_313:3:1204:5599:113305 1 - 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 9782
chr1 10001 10200 HWI-ST216_313:3:1102:14019:151362 1 - 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 9759
chr1 10127 10326 HWI-ST216_313:3:2207:7406:122346 1 - 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 9633
chr1 10241 10440 HWI-ST216_313:3:1302:4516:156396 1 - 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 9519"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b Distance".split()))

def test_minus(minus_chip, minus_input, expected_result_minus):

     print(minus_chip)
     print(minus_input)

     result = minus_chip.nearest(minus_input, how="next", overlap=False)

     print(result.df)
     # assert 0
     # print(PyRanges(result.df.tail()))

     assert_df_equal(result.df, expected_result_minus.df)


# @pytest.fixture()
# def plus_input():

#     c = """chr1	10073	10272	HWI-ST216:427:D29R1ACXX:2:2302:14161:85418	1	+
# chr1	10280	10479	HWI-ST216:427:D29R1ACXX:2:2103:15965:69480	1	+
# chr1	16056	16255	HWI-ST216:427:D29R1ACXX:2:1102:7604:12113	1	+
# chr1	16064	16263	HWI-ST216:427:D29R1ACXX:2:2105:18202:77798	1	+
# chr1	16109	16308	HWI-ST216:427:D29R1ACXX:2:2110:12286:25379	1	+"""

#     return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Name Score Strand".split()))



# def test_opposite(minus_chip, plus_input, expected_result_minus):

#      print(minus_chip)
#      print(plus_input)

#      result = minus_chip.nearest(plus_input, how="next", overlap=False)

#      print(result.df)
#      # assert 0
#      # print(PyRanges(result.df.tail()))

#      assert 0
     # assert_df_equal(result.df, expected_result_minus.df)



@pytest.fixture()
def hyp9():

    c = """chr1    2       3       +
chr1    1       4       +
chr1    1       4       +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp10():

    c = """chr1    3       4       +
chr1    4       5       +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def expected_result_counterexample4():

    c = """chr1  1  4  +  4   5   +         1
chr1  1  4  +  4   5   +         1
chr1  2  3  +  3   4   +         1"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))

def test_hypothesis_counterexample4(hyp9, hyp10, expected_result_counterexample4):

     print(hyp9)
     print(hyp10)

     result = hyp9.nearest(hyp10, how="next", overlap=False)

     print(result.df)
     # assert 0
     # print(PyRanges(result.df.tail()))

     assert_df_equal(result.df, expected_result_counterexample4.df)



@pytest.fixture()
def hyp11():

    c = """chr1    2       3       +
chr1    1       4       +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp12():

    c = """chr1    3       4       +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def expected_result_counterexample5():

    c = """chr1  2  3  +  3   4   +         1"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))

def test_hypothesis_counterexample5(hyp11, hyp12, expected_result_counterexample5):

     print(hyp11)
     print(hyp12)

     result = hyp11.nearest(hyp12, how="next", overlap=False)

     print(result.df)
     # assert 0
     # print(PyRanges(result.df.tail()))

     assert_df_equal(result.df, expected_result_counterexample5.df)



@pytest.fixture()
def hyp13():

    c = """chr1    1       2       +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp14():

    c = """chr1    1       2       +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))

def test_hypothesis_counterexample6(hyp13, hyp14):

     print(hyp13)
     print(hyp14)

     result = hyp13.nearest(hyp14, how="next", overlap=False)

     print(result.df)
     # assert 0
     # print(PyRanges(result.df.tail()))

     assert result.df.empty



@pytest.fixture()
def expected_result_counterexample6():

    c = """chr1 0 1 + 5726225 5726228 + 5726225
chr1 0 1 + 5726225 5726228 + 5726225
chr1 0 1 + 5726225 5726228 + 5726225
chr1 0 1 + 5726225 5726228 + 5726225
chr1 0 1 + 5726225 5726228 + 5726225
chr1 0 1 + 5726225 5726228 + 5726225
chr1 0 1 + 5726225 5726228 + 5726225
chr1 0 1 + 5726225 5726228 + 5726225
chr1 0 5726225 + 5726225 5726228 + 1
chr1 3538885 3832293 + 5726225 5726228 + 1893933"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))


@pytest.fixture()
def hyp15():

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
def hyp16():

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

def test_hypothesis_counterexample6(hyp15, hyp16, expected_result_counterexample6):

     print(hyp15)
     print(hyp16)

     result = hyp15.nearest(hyp16, how="next", overlap=False)

     print(result)
     # assert 0
     # print(PyRanges(result.df.tail()))
     assert_df_equal(result.df, expected_result_counterexample6.df)
