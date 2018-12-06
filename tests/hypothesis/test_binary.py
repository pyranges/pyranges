

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
from tests.hypothesis.hypothesis_helper import dfs_min



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


def run_bedtools(command, gr, gr2, strandedness):

    bedtools_strand = {False: "", "same": "-s", "opposite": "-S"}[strandedness]

    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = "{}/f1.bed".format(temp_dir)
        f2 = "{}/f2.bed".format(temp_dir)
        gr.df.to_csv(f1, sep="\t", header=False, index=False)
        gr2.df.to_csv(f2, sep="\t", header=False, index=False)

        cmd = command.format(f1=f1, f2=f2, strand=bedtools_strand)

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
@reproduce_failure('3.59.0', b'AXicY2RAA4woFBAAAABfAAQ=')
def test_set_union(gr, gr2, strandedness):

    set_union_command = "cat {f1} {f2} | bedtools sort | bedtools merge -c 4,5,6 -o first -i -"    # set_union_command = "bedtools merge {strand} -c 4,5,6 -o first -i {f1}"
    bedtools_result = run_bedtools(set_union_command, gr, gr2, strandedness)

    bedtools_df = read_bedtools_result_set_op(bedtools_result, strandedness)

    result = gr.set_union(gr2, strandedness=strandedness)
    print("r " * 30)
    print(result)
    # result_df

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
@pytest.mark.parametrize("strandedness", no_opposite)
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

    bedtools_strand = {False: "", "same": "-s", "opposite": "-S"}[strandedness]
    subtraction_command = "bedtools subtract {} -a {} -b {}"

    result_df = None
    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = "{}/f1.bed".format(temp_dir)
        f2 = "{}/f2.bed".format(temp_dir)
        gr.df.to_csv(f1, sep="\t", header=False, index=False)
        gr2.df.to_csv(f2, sep="\t", header=False, index=False)

        cmd = subtraction_command.format(bedtools_strand, f1, f2)

        result = subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()

        bedtools_df = pd.read_table(StringIO(result), header=None, squeeze=True, names="Chromosome Start End Name Score Strand".split())

    result = gr.subtraction(gr2, strandedness=strandedness)
    print("result\n", result)
    print("bedtools_df\n", PyRanges(bedtools_df))

    if not bedtools_df.empty or not result.df.empty:
        assert_df_equal(result.df, bedtools_df)
    else:
        assert bedtools_df.empty == result.df.empty


                    # "bedtools closest {} -t first -io -d -a <(sort -k1,1 -k2,2n {}) -b <(sort -k1,1 -k2,2n {})"]
nearest_hows = [None, "previous", "next"]
overlaps = [True, False]
strandedness = [False, "same", "opposite"]


@pytest.mark.bedtools
@pytest.mark.parametrize("nearest_how,overlap,strandedness", product(nearest_hows, overlaps, strandedness))
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min(), gr2=dfs_min())
def test_nearest(gr, gr2, nearest_how, overlap, strandedness):

    bedtools_strand = {False: "", "same": "-s", "opposite": "-S"}[strandedness]
    bedtools_overlap = {True: "", False: "-io"}[overlap]

    nearest_command = "bedtools closest {} {} -t first -d -a <(sort -k1,1 -k2,2n {}) -b <(sort -k1,1 -k2,2n {})"
    # sort_values = "Distance Start End".split()
    result_df = None
    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = "{}/f1.bed".format(temp_dir)
        f2 = "{}/f2.bed".format(temp_dir)
        gr.as_df().to_csv(f1, sep="\t", header=False, index=False)
        gr2.as_df().to_csv(f2, sep="\t", header=False, index=False)

        cmd = nearest_command.format(bedtools_overlap, bedtools_strand, f1, f2)
        # print(cmd)
        result = subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()

        bedtools_df = pd.read_table(StringIO(result), header=None, squeeze=True, names="Chromosome Start End Name Score Strand Chromosome_b Start_b End_b Name_b Score_b Strand_b Distance".split())
        # raise
        if not bedtools_df.empty:
            bedtools_df = bedtools_df[bedtools_df.Distance != -1]["Chromosome Start End Strand Distance".split()]

    result = gr.nearest(gr2, overlap=overlap, strandedness=strandedness).as_df()
    # print(result)
    if not len(result) == 0:
        result_df = result["Chromosome Start End Strand Distance".split()]

        # resetting index, because bedtools sorts result on start, end
        # while we want order in original file
        # print("before" * 50, result_df)
        # result_df.index = range(len(result_df))
        # print("after" * 50, result_df)
        assert_df_equal(result_df, bedtools_df)
    else:
        result_df = pd.DataFrame(columns="Chromosome Start End Strand Distance".split())
        assert bedtools_df.empty

    # if not result.df.empty:
    #     pyranges_distances = result_df.Distance.tolist()
    # else:
    #     pyranges_distances = []

    # print("bedtools", bedtools_distances)



@pytest.mark.bedtools
@pytest.mark.parametrize("strandedness", no_opposite)
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min(), gr2=dfs_min())
def test_jaccard(gr, gr2, strandedness):

    print(gr.df)
    print(gr2.df)
    jaccard_command = "bedtools jaccard {}  -a <(sort -k1,1 -k2,2n {}) -b <(sort -k1,1 -k2,2n {})"

    print(bedtools_strand)

    result_df = None
    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = "{}/f1.bed".format(temp_dir)
        f2 = "{}/f2.bed".format(temp_dir)
        gr.df.to_csv(f1, sep="\t", header=False, index=False)
        gr2.df.to_csv(f2, sep="\t", header=False, index=False)

        cmd = jaccard_command.format(bedtools_strand, f1, f2)
        print(f1)
        print(f2)
        print(cmd)
        result = subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()

        bedtools_jaccard = float(result.split("\n")[1].split()[2])

    print("result bedtools", bedtools_jaccard)
    result = gr.jaccard(gr2, strandedness=strandedness)
    print("pyranges result", result)

    assert 0 <= result <=1
    # sicne bedtools is buggy for this func?
    # assert np.isclose(result, bedtools_jaccard)


join_command = "bedtools intersect {} -wo -a {} -b {}"

@pytest.mark.bedtools
@pytest.mark.parametrize("strandedness", strandedness)
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min(), gr2=dfs_min())
def test_join(gr, gr2, strandedness):

    print(gr.df)
    print(gr2.df)

    bedtools_strand = {False: "", "same": "-s", "opposite": "-S"}[strandedness]

    result_df = None
    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = "{}/f1.bed".format(temp_dir)
        f2 = "{}/f2.bed".format(temp_dir)
        gr.df.to_csv(f1, sep="\t", header=False, index=False)
        gr2.df.to_csv(f2, sep="\t", header=False, index=False)

        cmd = join_command.format(bedtools_strand, f1, f2)
        print(cmd)
        result = subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()

        bedtools_df = pd.read_table(StringIO(result), header=None, squeeze=True, names="Chromosome Start End Name Score Strand Chromosome_b Start_b End_b Name_b Score_b Strand_b Overlap".split(), dtype={"Chromosome": "category", "Strand": "category"} )
        bedtools_df = bedtools_df.drop(["Overlap", "Chromosome_b"], 1)

    print("gr\n", gr)
    print("gr2\n", gr2)

    result = gr.join(gr2, strandedness=strandedness)

    print("result\n", result)
    print("bedtools\n", PyRanges(bedtools_df))

    if result.df.empty:
        assert bedtools_df.empty
    else:
        assert_df_equal(result.df, bedtools_df)
