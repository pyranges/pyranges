

import pytest
from tests.helpers import assert_df_equal

from pyranges.pyranges import PyRanges
import pyranges as pr

import pandas as pd

from io import StringIO


@pytest.fixture
def expected_result_unstranded():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b Distance
0 chr1 3 6 h 0 + 1 2 f 0 + 2
1 chr1 5 7 h 0 - 1 2 f 0 + 4
2 chr1 8 9 h 0 + 6 7 f 0 - 2"""

    return PyRanges(pd.read_table(StringIO(c), sep=" "))


def test_nearest_bed_unstranded(f1, f2, expected_result_unstranded):

    print(f1)
    print(f2)

    result = f1.nearest(f2, strandedness=False, suffix="_b", how="previous", overlap=False)

    print(result.df.to_csv(sep=" "))

    assert_df_equal(result.df, expected_result_unstranded.df)


@pytest.fixture
def expected_result_opposite_stranded():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b Distance
0 chr1 10246 10445 HWI-ST216_313:3:1207:4315:142177 1 + 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 60
1 chr1 110246 110445 HWI-ST216_313:3:1207:4315:142177 1 + 19958 20157 HWI-ST216:427:D29R1ACXX:2:1313:6283:67310 1 - 90090"""

    return PyRanges(pd.read_table(StringIO(c), sep=" "))


def test_nearest_bed_opposite_stranded(chip_10_plus_one, input_10, expected_result_opposite_stranded):

    result = chip_10_plus_one.nearest(input_10, strandedness="opposite", suffix="_b", how="previous", overlap=False)

    # print("result", result.df.to_csv(sep=" "))

    assert_df_equal(result.df, expected_result_opposite_stranded.df)


@pytest.fixture
def expected_result_same_stranded():

    c = """Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b Distance
0 chr1 10241 10440 HWI-ST216_313:3:1302:4516:156396 1 - 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 55
1 chr1 110246 110445 HWI-ST216_313:3:1207:4315:142177 1 + 16109 16308 HWI-ST216:427:D29R1ACXX:2:2110:12286:25379 1 + 93939"""

    return PyRanges(pd.read_table(StringIO(c), sep=" "))


def test_nearest_bed_same_stranded(chip_10_plus_one, input_10, expected_result_same_stranded):

    print(chip_10_plus_one)
    print(input_10)

    result = chip_10_plus_one.nearest(input_10, strandedness="same", suffix="_b", how="previous", overlap=False)

    assert_df_equal(result.df, expected_result_same_stranded.df)




@pytest.fixture()
def expected_result_counterexample1():

    c = """chr1  2  3  +  1   2   +         1"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))


@pytest.fixture()
def hyp1():

    c = """chr1	2	3	+"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp2():

    c = """chr1	1	2	+
chr1	0	3	+"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))

def test_hypothesis_counterexample(hyp1, hyp2, expected_result_counterexample1):

     print(hyp1)
     print(hyp2)

     result = hyp1.nearest(hyp2, how="previous", overlap=False)

     print(result)
     # 0
     # print(PyRanges(result.df.tail()))
     assert_df_equal(result.df, expected_result_counterexample1.df)



@pytest.fixture()
def expected_result_counterexample2():

    c = """chr1  2  3  +  1   2   +         1"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))


@pytest.fixture()
def hyp3():

    c = """chr1      2    3      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp4():

    c = """chr1      1    2      +
chr1      0    4      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))

def test_hypothesis_counterexample2(hyp3, hyp4, expected_result_counterexample2):

     print(hyp3)
     print(hyp4)

     result = hyp3.nearest(hyp4, how="previous", overlap=False)

     print(result)
     # 0
     # print(PyRanges(result.df.tail()))
     assert_df_equal(result.df, expected_result_counterexample2.df)



@pytest.fixture()
def hyp5():

    c = """chr1      0    1      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp6():

    c = """chr1      0    1      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))

def test_hypothesis_counterexample3(hyp5, hyp6):

     print(hyp5)
     print(hyp6)

     result = hyp5.nearest(hyp6, how="previous", overlap=False)

     print(result)
     # 0
     # print(PyRanges(result.df.tail()))
     result.df.empty



@pytest.fixture()
def hyp7():

    c = """chr1	1	2	+
chr1	1	2	+"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp8():

    c = """chr1      0    1      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def expected_result_counterexample4():

    c = """chr1      1    2      +        0      1        +         1
chr1      1    2      +        0      1        +         1"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))


def test_hypothesis_counterexample4(hyp7, hyp8, expected_result_counterexample4):

     print(hyp7)
     print(hyp8)

     result = hyp7.nearest(hyp8, how="previous", overlap=False)

     print(result.df)
     # 0
     # print(PyRanges(result.df.tail()))
     assert_df_equal(result.df, expected_result_counterexample4.df)



@pytest.fixture()
def hyp9():

    c = """chr1	0	1	+
chr1	0	1	+"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp10():

    c = """chr1      0    1      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))

def test_hypothesis_counterexample5(hyp9, hyp10):

     print(hyp9)
     print(hyp10)

     result = hyp9.nearest(hyp10, how="previous", overlap=False)

     print(result.df)
     # 0
     # print(PyRanges(result.df.tail()))
     result.df.empty



@pytest.fixture()
def hyp11():

    c = """chr1	1	2	+
chr1	1	2	+"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp12():

    c = """chr1      0    1      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def expected_result_counterexample6():

    c = """chr1      1    2      +        0      1        +         1
chr1      1    2      +        0      1        +         1"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))

def test_hypothesis_counterexample6(hyp11, hyp12, expected_result_counterexample6):

     print(hyp11)
     print(hyp12)

     result = hyp11.nearest(hyp12, how="previous", overlap=False)

     print(result.df)
     # 0
     # print(PyRanges(result.df.tail()))

     assert_df_equal(result.df, expected_result_counterexample6.df)


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

    c = """0 chr1 10241 10440 HWI-ST216_313:3:1302:4516:156396 1 - 9988 10187 HWI-ST216:427:D29R1ACXX:2:1205:6095:16532 1 - 55"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Name Score Strand Start_b End_b Name_b Score_b Strand_b Distance".split()))

def test_minus(minus_chip, minus_input, expected_result_minus):

     print(minus_chip)
     print(minus_input)

     result = minus_chip.nearest(minus_input, how="previous", overlap=False)

     print(result.df)
     # 0
     # print(PyRanges(result.df.tail()))

     assert_df_equal(result.df, expected_result_minus.df)



@pytest.fixture()
def hyp13():

    c = """chr1	1	2	+
chr1	1	2	+
chr1	1	2	+
chr1	1	2	+
chr1	1	2	+
chr1	1	2	+
chr1	1	2	+"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp14():

    c = """chr1      0    1      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def expected_result_counterexample6():

    c = """chr1      1    2      +        0      1        +         1
chr1      1    2      +        0      1        +         1
chr1      1    2      +        0      1        +         1
chr1      1    2      +        0      1        +         1
chr1      1    2      +        0      1        +         1
chr1      1    2      +        0      1        +         1
chr1      1    2      +        0      1        +         1"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))

def test_hypothesis_counterexample6(hyp13, hyp14, expected_result_counterexample6):

     print(hyp13)
     print(hyp14)

     result = hyp13.nearest(hyp14, how="previous", overlap=False)

     print(result.df)
     # 0
     # print(PyRanges(result.df.tail()))

     assert_df_equal(result.df, expected_result_counterexample6.df)



@pytest.fixture()
def hyp15():

    c = """chr1	0	3	+
chr1	1	2	+"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp16():

    c = """chr1	2	3	+"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


# @pytest.fixture()
# def expected_result_counterexample7():

#     c = """chr1  1  2  +  2   3   +         1"""

#     return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
#                                   names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))

def test_hypothesis_counterexample7(hyp15, hyp16):

     print(hyp15)
     print(hyp16)

     result = hyp15.nearest(hyp16, how="previous", overlap=False)

     print(result.df)
     # 0
     # print(PyRanges(result.df.tail()))

     result.df.empty # assert_df_equal(result.df, expected_result_counterexample7.df)



@pytest.fixture()
def hyp17():

    c = """chr1      1    2      +
chr1      2    3      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp18():

    c = """chr1      0    1      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def expected_result_counterexample8():

    c = """chr1  1  2  +  0   1   +         1
chr1  2  3  +  0   1   +         2"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))

def test_hypothesis_counterexample8(hyp17, hyp18, expected_result_counterexample8):

     print(hyp17)
     print(hyp18)

     result = hyp17.nearest(hyp18, how="previous", overlap=False)

     print(result.df)
     # 0
     # print(PyRanges(result.df.tail()))

     assert_df_equal(result.df, expected_result_counterexample8.df)



@pytest.fixture()
def hyp19():

    c = """chr1      1    2      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp20():

    c = """chr1      0    1      +
chr1      0    2      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def expected_result_counterexample9():

    c = """chr1  1  2  +  0   1   +         1"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))

def test_hypothesis_counterexample9(hyp19, hyp20, expected_result_counterexample9):

     print(hyp19)
     print(hyp20)

     result = hyp19.nearest(hyp20, how="previous", overlap=False)

     print(result.df)
     # 0
     # print(PyRanges(result.df.tail()))

     assert_df_equal(result.df, expected_result_counterexample9.df)



@pytest.fixture()
def hyp21():

    c = """chr1	1	2	+
chr1	2	3	+"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp22():

    c = """chr1      0    2      +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def expected_result_counterexample10():

    c = """chr1  2  3  +  0   2   +         1"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))

def test_hypothesis_counterexample10(hyp21, hyp22, expected_result_counterexample10):

     print(hyp21)
     print(hyp22)

     result = hyp21.nearest(hyp22, how="previous", overlap=False)

     print(result.df)
     # 0
     # print(PyRanges(result.df.tail()))

     assert_df_equal(result.df, expected_result_counterexample10.df)



@pytest.fixture()
def hyp23():

    c = """chr1    1       2       +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp24():

    c = """chr1    1       2       +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


def test_hypothesis_counterexample11(hyp23, hyp24):

     print(hyp23)
     print(hyp24)

     result = hyp23.nearest(hyp24, how="previous", overlap=False)

     print(result.df)
     # 0
     # print(PyRanges(result.df.tail()))

     result.df.empty



@pytest.fixture()
def hyp25():

    c = """chr1    3       4       +
chr1    4       5       +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def hyp26():

    c = """chr1    0       4       +
chr1    1       3       +
chr1    2       3       +"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None, names="Chromosome Start End Strand".split()))


@pytest.fixture()
def expected_result_counterexample12():

    c = """chr1  3  4  +  2   3   +         1
chr1  4  5  +  0   4   +         1"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))


def test_hypothesis_counterexample12(hyp25, hyp26, expected_result_counterexample12):

     print(hyp25)
     print(hyp26)

     result = hyp25.nearest(hyp26, how="previous", overlap=False)

     print(result.df)
     # 0
     # print(PyRanges(result.df.tail()))

     assert_df_equal(result.df, expected_result_counterexample12.df)



@pytest.fixture()
def expected_result_counterexample13():

    c = """chr1 3538885 3832293 + 0 1 + 3538885
chr1 4426346 9655531 + 0 1 + 4426346"""

    return PyRanges(pd.read_table(StringIO(c), sep="\s+", header=None,
                                  names="Chromosome Start End Strand Start_b End_b Strand_b Distance".split()))


@pytest.fixture()
def hyp27():

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
def hyp28():

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

def test_hypothesis_counterexample13(hyp27, hyp28, expected_result_counterexample13):

     print(hyp27)
     print(hyp2)

     result = hyp27.nearest(hyp28, how="previous", overlap=False)

     print(result)

     assert_df_equal(result.df, expected_result_counterexample13.df)
