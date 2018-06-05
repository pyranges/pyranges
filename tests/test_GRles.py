
import pytest

import pandas as pd
import numpy as np

from pyrle import Rle
from io import StringIO

from pyrle import PyRles
from pyrle.methods import coverage, _to_ranges, to_ranges

import pyranges as pr


@pytest.fixture
def chipseq_dataset():

    gr = pr.load_dataset("chipseq")

    return gr


@pytest.fixture
def teensy():

    c = """Chromosome Start End Name Score Strand
chr2 13611 13636 U0 2 -
chr2 32620 32645 U0 4 -
chr2 33241 33266 U0 1 +
chr2 788268 788293 U0 0 -
chr2 1150665 1150690 U0 -1 -"""

    return pr.PyRanges(pd.read_table(StringIO(c), sep="\s+"))


# def test_roundtrip_to_ranges_grles(chipseq_dataset):

#     cv = PyRles(chipseq_dataset)

#     chipseq_dataset.df.to_csv("input_dataset.txt", sep=" ")
#     gr = to_ranges(cv)

#     print(chipseq_dataset)
#     print(gr)
#     gr.df.to_csv("result.txt", sep=" ")

#     assert list(gr.df.Start) == list(chipseq_dataset.df.Start)
#     assert list(gr.df.Chromosome) == list(chipseq_dataset.df.Chromosome)
#     assert list(gr.df.End) == list(chipseq_dataset.df.End)


def test_roundtrip_to_ranges_single_rle_teensy(teensy):

    cv = coverage(teensy, value_col="Score")

    starts, ends, values = _to_ranges(cv)

    df = pd.concat([pd.Series(a) for a in [starts, ends, values]], axis=1)
    df.columns = "Start End Score".split()
    df.insert(0, "Chromosome", "chr2")

    gr = pr.PyRanges(df)

    # print(", ".join(str(i) for i in starts))
    # print(", ".join(str(f) for f in ends))
    assert list(starts) == [0, 13611, 13636, 32620, 32645, 33241, 33266, 1150665]
    assert list(ends) == [13611, 13636, 32620, 32645, 33241, 33266, 1150665, 1150690]


@pytest.fixture
def teensy_duplicated():

    "Contains duplicates"

    c = """Chromosome     Start       End Name  Score Strand
180       chr2  42058716  42058741   U0      0      +
181       chr2  42130511  42130536   U0      0      +
182       chr2  42593165  42593190   U0      0      -
183       chr2  42593165  42593190   U0      0      -
184       chr2  42635413  42635438   U0      0      -
185       chr2  43357333  43357358   U0      0      +
186       chr2  43854685  43854710   U0      0      +"""

    df = pd.read_table(StringIO(c), sep="\s+")

    return pr.PyRanges(df)



@pytest.fixture
def expected_result_teensy_duplicated():

    "Contains duplicates"

    c = """Chromosome     Start       End Name  Score Strand
180       chr2  42058716  42058741   U0      1      +
181       chr2  42130511  42130536   U0      1      +
183       chr2  42593165  42593190   U0      2      -
184       chr2  42635413  42635438   U0      1      -
185       chr2  43357333  43357358   U0      1      +
186       chr2  43854685  43854710   U0      1      +"""

    df = pd.read_table(StringIO(c), sep="\s+")

    return pr.PyRanges(df)


def test_roundtrip_to_ranges_single_rle_teensy_duplicated(teensy_duplicated, expected_result_teensy_duplicated):

    gr = teensy_duplicated
    cv = coverage(teensy_duplicated)
    # print(cv.values)

    starts, ends, values = _to_ranges(cv)

    # print(gr)
    # print(pr.PyRanges(gr.df.drop_duplicates()))
    # print("len(starts)", len(starts))
    # print("starts")
    # print(starts[:5])
    # print(starts[-5:])
    # print("ends")
    # print(ends[:5])
    # print(ends[-5:])

    assert list(starts) == [0, 42058716, 42058741, 42130511, 42130536, 42593165, 42593190, 42635413, 42635438, 43357333, 43357358, 43854685]
    assert list(ends) == [42058716, 42058741, 42130511, 42130536, 42593165, 42593190, 42635413, 42635438, 43357333, 43357358, 43854685, 43854710]


@pytest.fixture
def expected_result_single_chromosome(chipseq_dataset):

    return pr.PyRanges(chipseq_dataset["chr2"].df.drop_duplicates())


@pytest.fixture
def overlapping_gr():

    c = """Chromosome      Start        End Name  Score Strand
chr2  1  7   U0      0      -
chr2  4  10   U0      0      +"""

    df = pd.read_table(StringIO(c), sep="\s+")

    return pr.PyRanges(df)


def test_roundtrip_to_ranges_single_rle_overlapping(overlapping_gr):

    gr = overlapping_gr
    print(gr)
    cv = coverage(gr)
    print(cv)

    starts, ends, values = _to_ranges(cv)

    df = pd.concat([pd.Series(a) for a in [starts, ends, cv.values]], axis=1)
    df.columns = "Start End Score".split()
    df.insert(0, "Chromosome", "chr2")
    print(pr.PyRanges(df))
    # print(gr)
    # print(pr.PyRanges(gr.df.drop_duplicates()))
    # print("len(starts)", len(starts))
    # print("starts")
    # print(starts[:5])
    # print("ends")
    # print(ends[:5])


    assert list(starts) == [0, 1, 4, 7]
    assert list(ends) == [1, 4, 7, 10]

    # assert (starts == expected_result_problematic_gr.df.Start).all()
    # assert (ends == expected_result_problematic_gr.df.End).all()
    # assert (scores == expected_result_problematic_gr.df.Score).all()

# gr = PyRanges('toyChr',IRanges(cumsum(c(0,runLength(toyData)[-nrun(toyData)])),
#                               width=runLength(toyData)),
#              toyData = runValue(toyData))




# def test_roundtrip_to_ranges_single_chromosome(chipseq_dataset, expected_result_single_chromosome):

#     chr2 = chipseq_dataset["chr2"]
#     print(chr2)
#     print(expected_result_single_chromosome)

#     cv = coverage(chr2)
#     print(cv.values)

#     starts, ends, scores = _to_ranges(cv)

#     print(chr2.df.loc[list(range(600, 640))])
#     print(starts[183])
#     print(ends[183])

#     # gr = pr.PyRanges(pd.concat([chr2.df.head(5), chr2.df.tail(5)]))
#     # print(gr)
#     print("len(starts)", len(starts))
#     print("starts")
#     print(starts[:5])
#     print(starts[-5:])
#     print("ends")
#     print(ends[:5])
#     print(ends[-5:])
#     print("scores")
#     print(scores[:5])
#     print(scores[-5:])

#     assert list(starts) == list(expected_result_single_chromosome.df.Start)
#     assert (ends == expected_result_single_chromosome.df.End).all()
#     assert (scores == expected_result_single_chromosome.df.Score).all()



def test_create_grles(chipseq_dataset):

    grles = PyRles(chipseq_dataset)

    print(grles)



def test_create_stranded_grles(chipseq_dataset):


    grles = PyRles(chipseq_dataset, stranded=True)

    print(grles)


def test_create_stranded_grles_multicpu(chipseq_dataset):


    grles = PyRles(chipseq_dataset, stranded=True, n_jobs=5)

    print(grles)


# def test_grange_to_grle_and_back():







@pytest.fixture
def expected_result_coverage():

    runs = """9739215      25  463205      25 3069430      25    9143      25  993038 25  142071      25  260968      25   71512      25   18072      25 103292 25""".split()
    runs = [int(s) for s in runs]
    values = [0., 1., 0., 1., 0., 1., 0., 1., 0., 1., 0., 1., 0., 1., 0., 1., 0., 1., 0., 1.]

    return runs, values


@pytest.fixture
def df():

    c = """chr21	9739215	9739240	U0	0	-
chr21	10202445	10202470	U0	0	+
chr21	13271900	13271925	U0	0	-
chr21	13281068	13281093	U0	0	-
chr21	14274131	14274156	U0	0	-
chr21	14416227	14416252	U0	0	+
chr21	14677220	14677245	U0	0	-
chr21	14748757	14748782	U0	0	+
chr21	14766854	14766879	U0	0	+
chr21	14870171	14870196	U0	0	-"""

    return pd.read_table(StringIO(c), header=None, names="Chromosome Start End Name Score Strand".split())


@pytest.fixture
def grle1():

    r1 = Rle([1, 5, 10], [1, 2, 3])
    r2 = Rle([4, 7, 9], [0.01, 0.02, 0.03])

    d = {"chr1": r1, "chr2": r2}
    print(d)
    g = PyRles(d)
    print(g)
    return g

@pytest.fixture
def grle1_stranded(grle1):

    return PyRles({("chr1", "+"): grle1["chr1"], ("chr1", "-"): grle1["chr1"], ("chr2", "-"): grle1["chr2"]})


@pytest.fixture()
def grle2():

    r1 = Rle([1, 2, 3], [1, 2, 3])
    r2 = Rle([5, 4, 2], [-0.1, -0.2, -0.3])

    d = {"chr1": r1, "chr2": r2}

    return PyRles(d)

@pytest.fixture
def expected_result():

    r1 = Rle([1, 2, 3, 10], [2, 4, 5, 3])
    r2 = Rle([4, 1, 4, 2, 9], [-0.09, -0.08, -0.18, -0.28, 0.03])

    d = {"chr1": r1, "chr2": r2}

    return PyRles(d)


@pytest.fixture
def d1():

    r1 = Rle([1, 5, 10], [1, 2, 3])
    r2 = Rle([4, 7, 9], [0.01, 0.02, 0.03])
    print(r1)

    d = {"chr1": r1, "chr2": r2}

    return d


def test_create_PyRles(d1):

    # No errors
    result = PyRles(d1)



def test_add_PyRles(grle1, grle2, expected_result):

    result = grle1.add(grle2)
    print(result)
    print("----")
    print(expected_result)

    assert result == expected_result


def test_create_stranded_grles(grle1_stranded):

    result = to_ranges(grle1_stranded)




# df = read.table(text="Chromosome Start End Name Score Strand
# 625 chr2 175474407 175474432 U0 0 -
# 626 chr2 175474427 175474452 U0 0 +", sep=" ")

# ir = IRanges(start=df$Start, end=df$End)

# cv = coverage(ir)


# ir2 = IRanges(cumsum(c(0,runLength(cv)[-nrun(cv)])),
#               width=runLength(cv))
