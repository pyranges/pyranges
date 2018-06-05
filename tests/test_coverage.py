import pytest

import pandas as pd
import numpy as np

from pyrle import Rle
from io import StringIO

from pyrle.methods import coverage

import pyranges as pr

# expected result for whole chr21 according to Rle
#   [1] 9739215      25  463205      25 3069430      25    9143      25  993038
#  [10]      25  142071      25  260968      25   71512      25   18072      25
#  [19]  103292      25  211969      25  275065      25  337961      25   28264
#  [28]      25   85316      25 1327122      25  942292      25  515321      25
#  [37]   51031      25  505279      25  106831      25  525295      25  428835
#  [46]      25 1075916      25   59077      25  397488      25  146032      25
#  [55]  241044      25   40723      25    6325      25  790758      25  636326
#  [64]      25  265715      25  707648      25  442479      25  623307      25
#  [73]   39301      25   56624      25   37674      25    6412      25   75632
#  [82]      25    7020      25   91867      25  516098      25  455342      25
#  [91]  207316      25   89376      25  220415      25   63873      25  178563
# [100]      25  208833      25  231384      25  233483      25  210876      25
# [109]   13625      25  698897      25    9427      25  199410      25 1334739
# [118]      25  181942      25   17516      25  186634      25   18445      25
# [127]   56362      25    6113      25    8690      25  314362      25  427565
# [136]      25  194922      25   25119      25   96491      25  814824      25
# [145]  329690      25   77988      25  715491      25   23244      25  305720
# [154]      25   45296      25  332040      25  174008      25  163489      25
# [163]   51692      25  622758      25  582426      25  167975      25   99400
# [172]      25   14499      25  481575      25  224239      25  109893      25
# [181]  229084      25  481166      25   76961      25  104924      25  262629
# [190]      25  123925      25   72451      25  423954      25  622114      25
# [199]  870208      25  291275      25   58970      25  189900      25  972143
# [208]      25  532150      25  157577      25  360979      25  122030      25
# [217]  365189      25 1376353      25  251038      25  338889      25


# def test_chr2_to_coverage():

#     gr = pr.load_dataset("chipseq")
#     gr_chr2 = gr["chr2"]

#     # print(gr_chr2)
#     # gr_chr2.df.to_csv("gr_chr2", sep=" ", index=False)

#     coverage(gr_chr2)

#     assert 0

@pytest.fixture
def simple():
  c = """Chromosome Start End Score
chr2      0    1   -1.0
chr2      1    3    1.0"""

  return pd.read_table(StringIO(c), sep="\s+")


def test_coverage_simple(simple):

    result = coverage(simple, value_col="Score")

    print(result)

    assert result == Rle([1, 2], [-1, 1])




@pytest.fixture()
def expected_result_coverage():

    runs = """9739215      25  463205      25 3069430      25    9143      25  993038 25  142071      25  260968      25   71512      25   18072      25 103292 25""".split()
    runs = [int(s) for s in runs]
    values = [0., 1., 0., 1., 0., 1., 0., 1., 0., 1., 0., 1., 0., 1., 0., 1., 0., 1., 0., 1.]

    return runs, values


@pytest.fixture()
def df():

    c = """chr21	9739215	9739240	U0	0	-
chr21	10202445	10202470	U0	0	+
chr21	14416227	14416252	U0	0	+
chr21	14677220	14677245	U0	0	-
chr21	14748757	14748782	U0	0	+
chr21	13271900	13271925	U0	0	-
chr21	13281068	13281093	U0	0	-
chr21	14274131	14274156	U0	0	-
chr21	14870171	14870196	U0	0	-
chr21	14766854	14766879	U0	0	+"""

    return pd.read_table(StringIO(c), header=None, names="Chromosome Start End Name Score Strand".split())


def test_coverage(df, expected_result_coverage):

    expected_runs, expected_values = expected_result_coverage

    # result = coverage(df)
    result = coverage(df)

    print(result.runs)
    print(result.values)


    print(len(result.runs))
    print(len(result.values))

    print(expected_runs)
    print(expected_values)

    assert list(result.runs) == expected_runs
    assert list(result.values) == expected_values



@pytest.fixture
def supersimple_bed():

    c = """Start End
2 3"""

    return pd.read_table(StringIO(c), sep="\s+", header=0)


# def test_coverage(supersimple_bed):

#     result = coverage(supersimple_bed)

#     assert list(result.runs) == [2, 1]
#     assert list(result.values) == [0, 1]


# @pytest.fixture
# def empty_bed():

#     c = """Start End
# 2 2
# 4 4"""

#     return pd.read_table(StringIO(c), sep="\s+", header=0)


# @pytest.fixture()
# def expected_result_supersimple_bed2():

#     runs = np.array([2, 1, 1, 1], dtype=np.int)
#     values = np.array([0, 1, 0, 1], dtype=np.float)

#     return Rle(runs, values)


# def test_supersimple_bed2(empty_bed, expected_result_supersimple_bed2):

#     result = coverage(empty_bed)
#     print(result.runs)
#     print(result.values)

#     assert list(result.runs) == [2, 2]


@pytest.fixture
def simple_bed():

    c = """Start End Value
3 6 2.4
4 7 0.9
5 6 3.33"""

    return pd.read_table(StringIO(c), sep="\s+", header=0)


@pytest.fixture()
def expected_result_simple_bed():

    runs = np.array([3, 1, 1, 1, 1], dtype=np.int)
    values = np.array([0, 1, 2, 3, 1], dtype=np.float)

    return Rle(runs, values)

@pytest.fixture()
def expected_result_simple_bed_values():

    runs = np.array([3, 1, 1, 1, 1], dtype=np.int)
    values = np.array([0., 2.4, 3.3, 6.63, 0.9], dtype=np.float)

    return Rle(runs, values)


def test_simple_bed(simple_bed, expected_result_simple_bed):

    result = coverage(simple_bed)
    print(result.runs, expected_result_simple_bed.runs)
    print(result.values, expected_result_simple_bed.values)
    assert list(result.runs) == list(expected_result_simple_bed.runs)
    assert list(result.values) == list(expected_result_simple_bed.values)


def test_simple_bed_with_scores(simple_bed, expected_result_simple_bed_values):

    result = coverage(simple_bed, value_col="Value")
    print(result.runs, expected_result_simple_bed_values.runs)
    print(result.values, expected_result_simple_bed_values.values)
    assert list(result.runs) == list(expected_result_simple_bed_values.runs)
    assert np.allclose(result.values, expected_result_simple_bed_values.values)

# df = read.table("test.bed", header=FALSE)
# bg_df = read.table("control.bed", header=FALSE)
# > gr = GRanges(seqnames = df$V1,
# +              ranges = IRanges(start = df$V2, end = df$V3), strand = df$V6)
# >
# > gr
# GRanges object with 10000 ranges and 0 metadata columns:
#           seqnames                 ranges strand
#              <Rle>              <IRanges>  <Rle>
#       [1]     chr8 [ 28510032,  28510057]      -
#       [2]     chr7 [107153363, 107153388]      -
#       [3]     chr5 [135821802, 135821827]      -
#       [4]    chr14 [ 19418999,  19419024]      -
#       [5]    chr12 [106679761, 106679786]      -
#       ...      ...                    ...    ...
#    [9996]     chr4 [153155301, 153155326]      +
#    [9997]     chr9 [120803448, 120803473]      +
#    [9998]     chr6 [ 89296757,  89296782]      -
#    [9999]     chr1 [194245558, 194245583]      +
#   [10000]     chr8 [ 57916061,  57916086]      +
#   -------
#   seqinfo: 24 sequences from an unspecified genome; no seqlengths
# > bg_df
# > bg_gr = GRanges(seqnames = bg_df$V1,
# +              ranges = IRanges(start = bg_df$V2, end = bg_df$V3), strand = bg_df$V6)
# >
# > cbg = coverage(bg_gr)
# > c = coverage(gr)
# > c + cbg
# RleList of length 25
# $chr1
# integer-Rle of length 247134924 with 3242 runs
#   Lengths: 887770     26 106863     26  46416 ... 730067     26 259249     26
#   Values :      0      1      0      1      0 ...      0      1      0      1

# $chr10
# integer-Rle of length 135184308 with 1830 runs
#   Lengths: 122196     26 183504     26  44691 ...  43354     26  80887     26
#   Values :      0      1      0      1      0 ...      0      1      0      1

# $chr11
# integer-Rle of length 134407943 with 1884 runs
#   Lengths: 183570     26 114611     26  19735 ...   2576     26  97466     26
#   Values :      0      1      0      1      0 ...      0      1      0      1

# $chr12
# integer-Rle of length 132153042 with 1916 runs
#   Lengths: 124014     26 464843     26  32405 ...  16454     26   2826     26
#   Values :      0      1      0      1      0 ...      0      1      0      1

# $chr13
# integer-Rle of length 114064253 with 1236 runs
#   Lengths: 18236898       26   157375       26 ...       26    29952       26
#   Values :        0        1        0        1 ...        1        0        1

# ...
# <20 more elements>
# There were 26 warnings (use warnings() to see them)
# > result = c + cbg
# There were 26 warnings (use warnings() to see them)
# > result["chrY"]
# RleList of length 1
# $chrY
# integer-Rle of length 153847044 with 1096 runs
#   Lengths: 2820035      26   34804      26 ...  113694      26   16734      26
#   Values :       0       1       0       2 ...       0       1       0       1

# > runValues(result["chrY"])
# Error: could not find function "runValues"
# > runValue(result["chrY"])
# IntegerList of length 1
# [["chrY"]] 0 1 0 2 0 2 0 2 0 1 0 1 0 1 0 1 ... 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
# > runValue(result[["chrY"]])
#    [1] 0 1 0 2 0 2 0 2 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0
#   [38] 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
#   [75] 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 2 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0
#  [112] 1 0 1 0 1 0 1 0 1 0 1 0 1 0 2 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
#  [149] 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 2 0 1 0
#  [186] 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
#  [223] 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 2 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0
#  [260] 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 2 0 1 0 1 0 1 0 2 0 1 0 1 0 1 0 1
#  [297] 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0
#  [334] 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
#  [371] 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0
#  [408] 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
#  [445] 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0
#  [482] 1 0 1 0 1 0 1 0 1 0 2 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
#  [519] 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0
#  [556] 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
#  [593] 0 1 0 1 0 1 0 2 0 1 0 1 0 1 0 2 0 1 0 2 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0
#  [630] 2 0 1 0 1 0 1 0 1 0 1 0 2 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 2 0 1 0 2 0 1 0 1
#  [667] 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0
#  [704] 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 2 0 1 0 1 0 1 0 2 0 1 0 1 0 1 0 1 0 1 0 2
#  [741] 0 1 0 1 0 1 0 1 0 2 0 1 0 1 0 1 0 1 0 2 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0
#  [778] 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 2 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
#  [815] 0 2 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0
#  [852] 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 2 0 1 0 1 0 2 0 1 0 1
#  [889] 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0
#  [926] 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
#  [963] 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0
# [1000] 1 0 1 0 2 0 1 0 1 0 1 0 1 0 2 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
# [1037] 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 2 0 1 0
# [1074] 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
# > sum(runLength(c))
#      chr1     chr10     chr11     chr12     chr13     chr14     chr15     chr16
# 247134924 135184308 134310451 130956155 114064253 106021297  99484340  87620182
#     chr17     chr18     chr19      chr2     chr20     chr21     chr22      chr3
#  76607176  76034254  63775899 242351150  62164654  45986669  48886512 199362640
#      chr4      chr5      chr6      chr7      chr8      chr9      chrX      chrY
# 190843795 180513820 169893543 158684288 146247528 140101127 153874131  22210662
# > sum(runLength(c["chrY"]))
#     chrY
# 22210662
# > sum(runLength(cbg["chrY"]))
#     chrY
# 57402239
# > sum(runLength(cbg["chrY"])) - sum(runLength(cbg["chrY"]))
# chrY
#    0
# > sum(runLength(cbg["chrY"])) - sum(runLength(c["chrY"]))
#     chrY
# 35191577
# > Rlec["chrY"]
