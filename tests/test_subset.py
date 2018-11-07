import pytest

from pyranges.pyranges import PyRanges

from tests.helpers import assert_dfs_equal, string_to_pyrange
import pandas as pd

from io import StringIO

@pytest.fixture
def subset_chr_strand():

    s = """Chromosome  Start  End Name  Score Strand
chr1      5    7    h      0      -"""

    return string_to_pyrange(s)


def test_subset_chr_strand(f1, subset_chr_strand):

    subset = f1["chr1", "-"]

    assert_dfs_equal(subset, subset_chr_strand)



@pytest.fixture
def subset_chr_slice():

     s = """Chromosome     Start       End Name  Score Strand
chrY  13517892  13517917   U0      0      -"""
     return string_to_pyrange(s)


def test_subset_chr_slice(cs, subset_chr_slice):

    subset = cs["chrY", 13517892:13517917]

    assert_dfs_equal(subset_chr_slice, subset)

@pytest.fixture
def subset_chr():

     s = """Chromosome     Start       End Name  Score Strand
2575       chrY  12930373  12930398   U0      0      +
2757       chrY  15548022  15548047   U0      0      +
3731       chrY   7194340   7194365   U0      0      +
6890       chrY  21559181  21559206   U0      0      +
7639       chrY  11942770  11942795   U0      0      +
8432       chrY   8316773   8316798   U0      0      +
8657       chrY   7463444   7463469   U0      0      +
242        chrY  13216614  13216639   U0      0      -
329        chrY  21751211  21751236   U0      0      -
359        chrY   7701983   7702008   U0      0      -
419        chrY  21910706  21910731   U0      0      -
913        chrY  22054002  22054027   U0      0      -
1209       chrY  16045242  16045267   U0      0      -
1704       chrY  21707662  21707687   U0      0      -
2002       chrY   7761026   7761051   U0      0      -
2044       chrY  22210637  22210662   U0      0      -
2085       chrY  14774053  14774078   U0      0      -
2151       chrY  16495497  16495522   U0      0      -
3023       chrY   7046809   7046834   U0      0      -
3131       chrY  15224235  15224260   U0      0      -
3816       chrY  13517892  13517917   U0      0      -
3897       chrY   8010951   8010976   U0      0      -
9570       chrY   7405376   7405401   U0      0      -"""
     return string_to_pyrange(s)


def test_subset_chr(cs, subset_chr):

    subset = cs["chrY"]

    print(subset.as_df())

    assert_dfs_equal(subset, subset_chr)

@pytest.fixture
def subset_strand():

    s = """Chromosome  Start  End Name  Score Strand
0       chr1      3    6    h      0      +
2       chr1      8    9    h      0      +"""

    return string_to_pyrange(s)


def test_subset_strand(f1, subset_strand):

    subset = f1["+"]
    print(subset.as_df())
    assert_dfs_equal(subset, subset_strand)

@pytest.fixture
def subset_strand_slice():

    s = """Chromosome Start End Name Score Strand
0 chr1 3 6 h 0 +"""

    return string_to_pyrange(s)

def test_subset_strand_slice(f1, subset_strand_slice):

    subset = f1["+", 1:5]
    print(subset.as_df())
    assert_dfs_equal(subset, subset_strand_slice)




def test_subset_chromosome_strand_slice(f1, subset_strand_slice):

    subset = f1["chr1", "+", 1:5]
    print(subset.as_df())
    assert_dfs_equal(subset, subset_strand_slice)


@pytest.fixture
def subset_slice():

    s = """Chromosome  Start  End Name  Score Strand
chr1      8    9    h      0      +
chr1      5    7    h      0      -"""

    return string_to_pyrange(s)


def test_subset_slice(f1, subset_slice):

    subset = f1[6:10]
    print(subset.as_df())
    assert_dfs_equal(subset, subset_slice)
