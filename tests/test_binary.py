

import pytest

from hypothesis import given, settings, reproduce_failure, unlimited, HealthCheck, seed
from hypothesis.extra.pandas import data_frames, columns, range_indexes, column, indexes
from hypothesis.extra.numpy import arrays
import hypothesis.strategies as st

from itertools import product
import tempfile
import subprocess
from io import StringIO

from pyrle import Rle
import pyranges as pr

import pandas as pd
import numpy as np


from tests.helpers import assert_df_equal
from tests.hypothesis_helper import dfs_min



import numpy as np

from os import environ

if environ.get("TRAVIS"):
    max_examples = 100
    deadline = None
else:
    max_examples = 100
    deadline = None


strandedness = [False, "same", "opposite"]
no_opposite = [False, "same"]


def run_bedtools(command, gr, gr2, strandedness, nearest_overlap=False):

    bedtools_strand = {False: "", "same": "-s", "opposite": "-S"}[strandedness]
    bedtools_overlap = {True: "", False: "-io"}[nearest_overlap]

    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = "{}/f1.bed".format(temp_dir)
        f2 = "{}/f2.bed".format(temp_dir)
        gr.df.to_csv(f1, sep="\t", header=False, index=False)
        gr2.df.to_csv(f2, sep="\t", header=False, index=False)

        cmd = command.format(f1=f1, f2=f2, strand=bedtools_strand, overlap=bedtools_overlap)
        print(cmd)

        result = subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()

    return result


def read_bedtools_result_set_op(bedtools_result, strandedness):

    if strandedness:
        usecols = [0, 1, 2, 5]
        names = "Chromosome Start End Strand".split()
    else:
        usecols = [0, 1, 2]
        names = "Chromosome Start End".split()

    return pd.read_table(StringIO(bedtools_result), header=None, usecols=usecols, names=names)


def compare_results(bedtools_df, result):

    if not bedtools_df.empty:
        assert_df_equal(result.df, bedtools_df)
    else:
        assert bedtools_df.empty == result.df.empty


def compare_results_nearest(bedtools_df, result):

    print("1 " * 50)
    print(bedtools_df)
    if not bedtools_df.empty:
        bedtools_df = bedtools_df[bedtools_df.Distance != -1]

    print("2 " * 50)
    print(result)

    result = result.df

    if not len(result) == 0:
        result_df = result["Chromosome Start End Strand Distance".split()]
        assert_df_equal(result_df, bedtools_df)
    else:
        assert bedtools_df.empty


@pytest.mark.bedtools
@pytest.mark.parametrize("strandedness", no_opposite)
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min(), gr2=dfs_min())
def test_set_intersect(gr, gr2, strandedness):

    set_intersect_command = "bedtools intersect {strand} -a <(sort -k1,1 -k2,2n {f1} | bedtools merge {strand} -c 4,5,6 -o first -i -) -b <(sort -k1,1 -k2,2n {f2} | bedtools merge {strand} -c 4,5,6 -o first -i -)"
    bedtools_result = run_bedtools(set_intersect_command, gr, gr2, strandedness)

    bedtools_df = read_bedtools_result_set_op(bedtools_result, strandedness)

    result = gr.set_intersect(gr2, strandedness=strandedness)

    compare_results(bedtools_df, result)


@pytest.mark.bedtools
@pytest.mark.parametrize("strandedness", no_opposite)
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min(), gr2=dfs_min())
def test_set_union(gr, gr2, strandedness):

    set_union_command = "cat {f1} {f2} | bedtools sort | bedtools merge {strand} -c 4,5,6 -o first -i -"    # set_union_command = "bedtools merge {strand} -c 4,5,6 -o first -i {f1}"
    bedtools_result = run_bedtools(set_union_command, gr, gr2, strandedness)

    bedtools_df = read_bedtools_result_set_op(bedtools_result, strandedness)

    result = gr.set_union(gr2, strandedness=strandedness)

    compare_results(bedtools_df, result)


@pytest.mark.bedtools
@pytest.mark.parametrize("strandedness", strandedness)
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min(), gr2=dfs_min())
def test_overlap(gr, gr2, strandedness):

    overlap_command = "bedtools intersect -wa {strand} -a {f1} -b {f2}"

    bedtools_result = run_bedtools(overlap_command, gr, gr2, strandedness)

    bedtools_df = pd.read_table(StringIO(bedtools_result), header=None, names="Chromosome Start End Name Score Strand".split())

    result = gr.overlap(gr2, strandedness=strandedness)

    compare_results(bedtools_df, result)


@pytest.mark.bedtools
@pytest.mark.parametrize("strandedness", strandedness)
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min(), gr2=dfs_min())
def test_intersect(gr, gr2, strandedness):

    intersect_command = "bedtools intersect {strand} -a {f1} -b {f2}"

    bedtools_result = run_bedtools(intersect_command, gr, gr2, strandedness)

    bedtools_df = pd.read_table(StringIO(bedtools_result), header=None, names="Chromosome Start End Name Score Strand".split())

    result = gr.intersect(gr2, strandedness=strandedness)

    compare_results(bedtools_df, result)



@pytest.mark.bedtools
@pytest.mark.parametrize("strandedness", ["same", "opposite", False]) #
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min(), gr2=dfs_min())
def test_subtraction(gr, gr2, strandedness):

    subtract_command = "bedtools subtract {strand} -a {f1} -b {f2}"

    bedtools_result = run_bedtools(subtract_command, gr, gr2, strandedness)

    bedtools_df = pd.read_table(StringIO(bedtools_result), header=None, names="Chromosome Start End Name Score Strand".split())

    result = gr.subtract(gr2, strandedness=strandedness)

    compare_results(bedtools_df, result)


                    # "bedtools closest {} -t first -io -d -a <(sort -k1,1 -k2,2n {}) -b <(sort -k1,1 -k2,2n {})"]
nearest_hows = [None, "previous", "next"]
overlaps = [True, False]
strandedness = [False, "same", "opposite"]


@pytest.mark.bedtools
@pytest.mark.parametrize("nearest_how,overlap,strandedness", product(nearest_hows, overlaps, strandedness))
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min(), gr2=dfs_min())
def test_nearest(gr, gr2, nearest_how, overlap, strandedness):

    nearest_command = "bedtools closest {strand} {overlap} -t first -d -a <(sort -k1,1 -k2,2n {f1}) -b <(sort -k1,1 -k2,2n {f2})"

    bedtools_result = run_bedtools(nearest_command, gr, gr2, strandedness, overlap)

    bedtools_df = pd.read_table(StringIO(bedtools_result), header=None, names="Chromosome Start End Strand Distance".split(), usecols=[0, 1, 2, 5, 12])

    result = gr.nearest(gr2, strandedness=strandedness, overlap=overlap)

    compare_results_nearest(bedtools_df, result)



@pytest.mark.bedtools
@pytest.mark.parametrize("strandedness", no_opposite)
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min(), gr2=dfs_min())
def test_jaccard(gr, gr2, strandedness):

    # jaccard_command = "bedtools jaccard {strand}  -a <(sort -k1,1 -k2,2n {f1}) -b <(sort -k1,1 -k2,2n {f2})"

    # bedtools_result = run_bedtools(jaccard_command, gr, gr2, strandedness)

    # bedtools_jaccard = float(bedtools_result.split("\n")[1].split()[2])

    result = gr.jaccard(gr2, strandedness=strandedness)

    # there is a bug in bedtools, so cannot use as oracle
    assert 0 <= result <= 1
    # assert abs(result - bedtools_jaccard) < 0.001

# jaccard example, which should be 1, but bedtools considers 0.5:
# +--------------+-----------+-----------+------------+-----------+----------+
# | Chromosome   |     Start |       End | Name       |     Score | Strand   |
# | (int8)       |   (int32) |   (int32) | (object)   |   (int64) | (int8)   |
# |--------------+-----------+-----------+------------+-----------+----------|
# | chr1         |         1 |         2 | a          |         0 | +        |
# +--------------+-----------+-----------+------------+-----------+----------+
# PyRanges object has 1 sequences from 1 chromosomes.,
# =+--------------+-----------+-----------+------------+-----------+----------+
# | Chromosome   |     Start |       End | Name       |     Score | Strand   |
# | (int8)       |   (int32) |   (int32) | (object)   |   (int64) | (int8)   |
# |--------------+-----------+-----------+------------+-----------+----------|
# | chr1         |         1 |         2 | a          |         0 | +        |
# | chr1         |         1 |         2 | a          |         0 | -        |
# +--------------+-----------+-----------+------------+-----------+----------+
# PyRanges object has 2 sequences from 1 chromosomes., strandedness='same')
#
# bedtools jaccard -s  -a <(sort -k1,1 -k2,2n /tmp/tmpuxs7jtw0/f1.bed) -b <(sort -k1,1 -k2,2n /tmp/tmpuxs7jtw0/f2.bed)
# intersection	union-intersection	jaccard	n_intersections
# 1	2	0.5	1




@pytest.mark.bedtools
@pytest.mark.parametrize("strandedness", strandedness)
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min(), gr2=dfs_min())
def test_join(gr, gr2, strandedness):

    join_command = "bedtools intersect {strand} -wo -a {f1} -b {f2}"

    bedtools_result = run_bedtools(join_command, gr, gr2, strandedness)

    bedtools_df = pd.read_table(StringIO(bedtools_result), header=None,
                                names="Chromosome Start End Name Score Strand Chromosome_b Start_b End_b Name_b Score_b Strand_b Overlap".split(),
                                dtype={"Chromosome": "category", "Strand": "category"}).drop("Chromosome_b Overlap".split(), axis=1)

    result = gr.join(gr2, strandedness=strandedness)

    if result.df.empty:
        assert bedtools_df.empty
    else:
        assert_df_equal(result.df, bedtools_df)
