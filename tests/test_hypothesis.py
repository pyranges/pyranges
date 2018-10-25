import pytest
from hypothesis import given, settings, reproduce_failure, unlimited, HealthCheck, seed
from hypothesis.extra.pandas import data_frames, columns, range_indexes, column, indexes
from hypothesis.extra.numpy import arrays
import hypothesis.strategies as st

from itertools import product
import tempfile
import subprocess
from io import StringIO

from pyranges import PyRanges

# using assert df equal, because we want to consider output from bedtools and
# pyranges equal even if they have different sort order
from tests.helpers import assert_df_equal

import numpy as np

from os import environ

if environ.get("TRAVIS"):
    max_examples = 100
    deadline = None
else:
    max_examples = 100
    deadline = None


def mysort(tp):

    if tp[1] == tp[2]:
        tp = (tp[0], tp[1], tp[2] + 1, tp[3])

    key = [-1, tp[1], tp[2], int(1e10)]

    return [x for _, x in sorted(zip(key, tp))]

chromosomes = st.sampled_from(["chr{}".format(str(e)) for e in list(range(1, 23)) + "X Y M".split()])
chromosomes_small = st.sampled_from(["chr1"])
cs = st.one_of(chromosomes, chromosomes_small)

positions = st.integers(min_value=0, max_value=int(1e7))
lengths = st.integers(min_value=1, max_value=int(1e7))
small_lengths = st.integers(min_value=1, max_value=int(1e4))
strands = st.sampled_from("+ -".split())
names = st.text("abcdefghijklmnopqrstuvxyz", min_size=1)
scores = st.integers(min_value=0, max_value=256)


# dfs = data_frames(columns=columns("Chromosome Start End Strand".split(),
#                                   dtype=int), rows=st.tuples(chromosomes, positions, positions,
#                                                              strands).map(mysort))

df_minsize = 1
nonempty_dfs = data_frames(index=indexes(dtype=np.int64, min_size=df_minsize, unique=True),
                           columns=columns("Chromosome Start End Strand".split(), dtype=int),
                           rows=st.tuples(chromosomes, positions, positions, strands).map(mysort))


better_df_minsize = 1
better_dfs = data_frames(index=indexes(dtype=np.int64, min_size=better_df_minsize, unique=True),
                         columns=[column("Chromosome", chromosomes),
                                  column("Start", elements=positions),
                                  column("End", elements=lengths),
                                  column("Strand", strands)])

better_dfs_min = data_frames(index=indexes(dtype=np.int64, min_size=better_df_minsize, unique=True),
                             columns=[column("Chromosome", cs),
                                      column("Start", elements=lengths),
                                      column("End", elements=small_lengths),
                                      # column("Name", elements=names),
                                      # column("Score", elements=scores),
                                      column("Strand", strands)])

better_dfs_min_single_chromosome = data_frames(index=indexes(dtype=np.int64, min_size=better_df_minsize, unique=True),
                                               columns=[column("Chromosome", chromosomes_small),
                                                        column("Start", elements=lengths),
                                                        column("End", elements=small_lengths),
                                                        column("Name", elements=names),
                                                        column("Score", elements=scores),
                                                        column("Strand", strands)])

better_dfs_min_nostrand = data_frames(index=indexes(dtype=np.int64, min_size=better_df_minsize, unique=True),
                             columns=[column("Chromosome", chromosomes_small),
                                      column("Start", elements=small_lengths),
                                      column("End", elements=small_lengths)])


large_df_minsize = 5
large_dfs = data_frames(index=indexes(dtype=np.int64, min_size=large_df_minsize, unique=True),
                           columns=columns("Chromosome Start End Name Score Strand".split(), dtype=int),
                           rows=st.tuples(chromosomes, positions, positions, names, small_lengths, strands).map(mysort))

runlengths = data_frames(index=indexes(dtype=np.int64, min_size=df_minsize, unique=True),
                         columns=[column("Runs", st.integers(min_value=1, max_value=int(1e7))),
                                  # must have a min/max on floats because R S4vectors translates too big ones into inf.
                                  # which is unequal to eg -1.79769e+308 so the tests fail
                                  column("Values", st.integers(min_value=-int(1e7), max_value=int(1e7)))])

# runlengths_same_length_floats = data_frames(index=range_indexes(min_size=df_minsize),
#                                             columns=[column("Runs", st.integers(min_value=1, max_value=int(1e7))),
#                                                      column("Values", st.floats(min_value=0.001, max_value=int(1e7))),
#                                                      column("Values2", st.floats(min_value=0.001, max_value=int(1e7)))])


runlengths_same_length_integers = data_frames(index=indexes(dtype=np.int64, min_size=df_minsize, unique=True),
                                              columns=[column("Runs", st.integers(min_value=1, max_value=int(1e4))),
                                              column("Values", st.integers(min_value=1, max_value=int(1e4))),
                                              column("Values2", st.integers(min_value=1, max_value=int(1e4)))])


# runs = arrays(st.integers(min_value=1, max_value=int(1e5)), shape=)
# values = arrays(dtype=np.float)
@st.composite
def dfs_min(draw):
    df = draw(better_dfs_min)
    df.loc[:, "End"] += df.Start

    df.insert(3, "Name", "a")
    df.insert(4, "Score", 0)
    gr = PyRanges(df)

    # do not sort like this, use pyranges sort
    # np.random.seed(draw(st.integers(min_value=0, max_value=int(1e6))))
    # gr.df = df.reindex(np.random.permutation(df.index.values))

    return gr


@st.composite
def dfs_min_single_chromosome(draw):
    df = draw(better_dfs_min_single_chromosome)
    df.loc[:, "End"] += df.Start
    df.insert(3, "Name", "a")
    df.insert(4, "Score", 0)

    return df


import pyranges as pr
import pandas as pd
from pyrle import Rle
import numpy as np


strandedness = [False, "same", "opposite"]
methods = ["join", "intersection", "overlap"]

# strandedness = strandedness[1:2]
# methods = methods[1:2]


coverage_cmd = "Rscript --vanilla tests/compute_coverage.R {} {}"

@pytest.mark.r
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(df=dfs_min_single_chromosome())
def test_coverage(df):

    print("---" * 10)
    p = pr.PyRanges(df)
    print("pyranges\n", p)

    c = p.coverage()["chr1"]

    result_df = None
    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = "{}/f1.txt".format(temp_dir)
        outfile = "{}/result.txt".format(temp_dir)
        R_df = df
        R_df.End = R_df.End - 1
        R_df.to_csv(f1, sep="\t", index=False)

        cmd = coverage_cmd.format(f1, outfile) + " 2>/dev/null"

        subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()

        result = pd.read_table(outfile)[["Runs.value", "Values.value"]]
        result.columns = "Runs Values".split()
        result = pd.concat([pd.DataFrame(index=[0], data={"Runs": 1, "Values": 0}), result], ignore_index=True)
        s4vectors_result = Rle(result.Runs, result.Values)

    print("pyranges result\n", c)
    print("s4vectors result\n", s4vectors_result)
    print(str(c == s4vectors_result) + " " * 10, c == s4vectors_result)

    assert np.all(np.equal(c.runs, s4vectors_result.runs))
    assert np.all(np.equal(c.values, s4vectors_result.values))


rle_operations = "+ - / *".split()
# rle_operations = ["-"]

rle_operation_cmd = "Rscript --vanilla tests/compute_Rle.R {} {} '{}' {}"

@pytest.mark.r
@given(runlengths=runlengths, runlengths2=runlengths)
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@pytest.mark.parametrize("operation", rle_operations)
def test_rle(runlengths, runlengths2, operation):

    # It is only compared against bioc with integers because float equality is hard,
    # for both libraries, sometimes end up with slightly different runlengths
    # [257.492105948544, 257.492075430654] gives
    # pyranges result
    #  +--------+---------+---------+---------+
    # | Runs   |       1 |       1 |       1 |
    # |--------+---------+---------+---------|
    # | Values | 257.492 | 257.492 | 257.492 |
    # +--------+---------+---------+---------+
    # Rle of length 3 containing 3 elements
    # s4vectors result
    #  +--------+---------+---------+
    # | Runs   |       2 |       1 |
    # |--------+---------+---------|
    # | Values | 257.492 | 257.492 |
    # +--------+---------+---------+
    # Rle of length 3 containing 2 elements

    pyop = {"+": "__add__", "-": "__sub__", "*": "__mul__", "/": "__truediv__"}[operation]

    print("runlengths", runlengths)
    print("runlengths2", runlengths2)

    r = Rle(runlengths.Runs, runlengths.Values)
    r2 = Rle(runlengths2.Runs, runlengths2.Values)

    print("r\n", r)
    print("r2\n", r2)

    m = getattr(r, pyop)
    result_pyranges = m(r2)

    print("pyranges result\n", result_pyranges)
    # f1 = "f1.txt"
    # f2 = "f2.txt"
    # runlengths.to_csv(f1, sep="\t", index=False)
    # runlengths2.to_csv(f2, sep="\t", index=False)
    result_df = None
    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = "{}/f1.txt".format(temp_dir)
        f2 = "{}/f2.txt".format(temp_dir)
        outfile = "{}/result.txt".format(temp_dir)
        runlengths.to_csv(f1, sep="\t", index=False)
        runlengths2.to_csv(f2, sep="\t", index=False)

        cmd = rle_operation_cmd.format(f1, f2, operation, outfile) + " 2>/dev/null"

        subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()

        result = pd.read_table(outfile)
        s4vectors_result = Rle(result.Runs, result.Values)

    print("pyranges result\n", result_pyranges)
    print("s4vectors result\n", s4vectors_result)

    assert np.allclose(result_pyranges.runs, s4vectors_result.runs, equal_nan=False)
    assert np.allclose(result_pyranges.values, s4vectors_result.values, equal_nan=True)

set_intersection_command = "bedtools intersect -a <(sort -k1,1 -k2,2n {} | bedtools merge -i -) -b <(sort -k1,1 -k2,2n {} | bedtools merge -i -)"


def name_score(df, df2):
    df.insert(3, "Name", "a")
    df.insert(4, "Score", 0)

    df2.insert(3, "Name", "b")
    df2.insert(4, "Score", 0)

    return df, df2

merge_command = "bedtools merge -o first -c 6 {} -i <(sort -k1,1 -k2,2n {})"
@pytest.mark.bedtools
@pytest.mark.parametrize("strand", [True, False])
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min())
def test_cluster(gr, strand):

    bedtools_strand = {True: "-s", False: ""}[strand]

    print(gr)

    result_df = None
    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = "{}/f1.bed".format(temp_dir)
        gr.df.to_csv(f1, sep="\t", header=False, index=False)

        cmd = merge_command.format(bedtools_strand, f1)
        print(cmd)

        result = subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()

        if not strand:
            bedtools_df = pd.read_table(StringIO(result), header=None, squeeze=True, usecols=[0, 1, 2], names="Chromosome Start End".split(), dtype={"Chromosome": "category"})
        else:
            bedtools_df = pd.read_table(StringIO(result), header=None, squeeze=True, names="Chromosome Start End Strand".split(), dtype={"Chromosome": "category"})

    print("bedtools_df\n", bedtools_df)
    result = gr.cluster(strand=strand)
    print("result\n", result.df)

    if not bedtools_df.empty:
        assert_df_equal(result.df, bedtools_df)
    else:
        assert bedtools_df.empty == result.df.empty


@pytest.mark.bedtools
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min(), gr2=dfs_min())
def test_set_intersection(gr, gr2):

    print(gr)
    print(gr2)

    result_df = None
    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = "{}/f1.bed".format(temp_dir)
        f2 = "{}/f2.bed".format(temp_dir)
        gr.df.to_csv(f1, sep="\t", header=False, index=False)
        gr2.df.to_csv(f2, sep="\t", header=False, index=False)

        cmd = set_intersection_command.format(f1, f2)

        result = subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()

        bedtools_df = pd.read_table(StringIO(result), header=None, squeeze=True, names="Chromosome Start End".split())

    print("bedtools_df\n", bedtools_df)
    result = gr.set_intersection(gr2)
    print("result\n", result.df)

    if not bedtools_df.empty:
        assert_df_equal(result.df, bedtools_df)
    else:
        assert bedtools_df.empty == result.df.empty


overlap_command = "bedtools intersect -u {} -a {} -b {}"

@pytest.mark.bedtools
@pytest.mark.parametrize("strandedness", strandedness)
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min(), gr2=dfs_min())
def test_overlap(gr, gr2, strandedness):

    bedtools_strand = {False: "", "same": "-s", "opposite": "-S"}[strandedness]

    print(gr.df)
    print(gr2.df)

    result_df = None
    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = "{}/f1.bed".format(temp_dir)
        f2 = "{}/f2.bed".format(temp_dir)
        gr.df.to_csv(f1, sep="\t", header=False, index=False)
        gr2.df.to_csv(f2, sep="\t", header=False, index=False)

        cmd = overlap_command.format(bedtools_strand, f1, f2)

        result = subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()

        bedtools_df = pd.read_table(StringIO(result), header=None, squeeze=True, names="Chromosome Start End Name Score Strand".split())

    result = gr.overlap(gr2, strandedness=strandedness)
    # print("result\n", result.df)
    # print("bedtools_df\n", bedtools_df)

    if not bedtools_df.empty:
        assert_df_equal(result.df, bedtools_df)
    else:
        assert bedtools_df.empty == result.df.empty

intersection_command = "bedtools intersect {} -a {} -b {}"

@pytest.mark.bedtools
@pytest.mark.parametrize("strandedness", strandedness)
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min(), gr2=dfs_min())
def test_intersection(gr, gr2, strandedness):

    bedtools_strand = {False: "", "same": "-s", "opposite": "-S"}[strandedness]

    print(gr.df)
    print(gr2.df)

    result_df = None
    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = "{}/f1.bed".format(temp_dir)
        f2 = "{}/f2.bed".format(temp_dir)
        gr.df.to_csv(f1, sep="\t", header=False, index=False)
        gr2.df.to_csv(f2, sep="\t", header=False, index=False)

        cmd = intersection_command.format(bedtools_strand, f1, f2)

        result = subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()

        bedtools_df = pd.read_table(StringIO(result), header=None, squeeze=True, names="Chromosome Start End Name Score Strand".split())

    result = gr.intersection(gr2, strandedness=strandedness)
    print("result\n", result.df)
    print("bedtools_df\n", bedtools_df)

    if not bedtools_df.empty:
        assert_df_equal(result.df, bedtools_df)
    else:
        assert bedtools_df.empty == result.df.empty


set_union_command = "bedtools merge -i <(sort -k1,1 -k2,2n --merge <(sort -k1,1 -k2,2n {} | bedtools merge -i -) <(sort -k1,1 -k2,2n {} | bedtools merge -i -))"

@pytest.mark.bedtools
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min(), gr2=dfs_min())
def test_set_union(gr, gr2):

    result_df = None
    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = "{}/f1.bed".format(temp_dir)
        f2 = "{}/f2.bed".format(temp_dir)
        gr.df.to_csv(f1, sep="\t", header=False, index=False)
        gr2.df.to_csv(f2, sep="\t", header=False, index=False)

        cmd = set_union_command.format(f1, f2)

        result = subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()

        bedtools_df = pd.read_table(StringIO(result), header=None, squeeze=True, names="Chromosome Start End".split())

    result = gr.set_union(gr2)
    print("result\n", result)
    print("bedtools_df\n", PyRanges(bedtools_df))

    if not bedtools_df.empty:
        assert_df_equal(result.df, bedtools_df)
    else:
        assert bedtools_df.empty == result.df.empty

subtraction_command = "bedtools subtract {} -a {} -b {}"

@pytest.mark.bedtools
@pytest.mark.parametrize("strandedness", ["same"]) # , False
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min(), gr2=dfs_min())
def test_subtraction(gr, gr2, strandedness):

    # print("gr\n", gr)
    # print("gr2\n", gr2)
    bedtools_strand = {False: "", "same": "-s", "opposite": "-S"}[strandedness]

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

    if not bedtools_df.empty:
        assert_df_equal(result.df, bedtools_df)
    else:
        assert bedtools_df.empty == result.df.empty


how = [False, "nonoverlapping", "previous_nonoverlapping", "next_nonoverlapping", "next", "previous"]

rle_commute_how = ["__add__", "__mul__"]

@pytest.mark.parametrize("how", rle_commute_how)
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min(), gr2=dfs_min())
def test_commutative_rles(gr, gr2, how):

    cv = gr.coverage(stranded=True)
    cv2 = gr2.coverage(stranded=True)

    method = getattr(cv, how)
    method2 = getattr(cv2, how)

    result = method(cv2)
    result2 = method2(cv)

    assert result == result2, "\n".join([str(e) for e in [cv, cv2, result, result2, "---" * 10]])

@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(df=runlengths_same_length_integers)
def test_inverse_div_mul_rles(df):

    """Testing with small integers, since small value floating points might lead to
mul then div not being equal to identity function because of float equality."""

    print(df)
    runlength = df.Runs.sum()

    cv = Rle(df.Runs.values, df.Values.values)

    newruns = np.random.permutation(df.Runs.values)
    print("newruns", newruns)
    cv2 = Rle(newruns, df.Values2.values)

    print("cv\n", cv)
    print("cv2\n", cv2)

    assert runlength == np.sum(cv.runs) and runlength == np.sum(cv2.runs)

    result = cv / cv2

    result2 = result * cv2

    print("result\n", result)
    print("result2\n", result2)

    assert np.all(np.equal(result2.runs, cv.runs))
    assert np.allclose(result2.values, cv.values)


@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(df=runlengths_same_length_integers)
def test_inverse_add_sub_rles(df):

    """Testing with small integers, since small value floating points might lead to
mul then div not being equal to identity function because of float equality."""

    cv = Rle(df.Runs.values, df.Values.values)

    cv2 = Rle(np.random.permutation(df.Runs.values), df.Values2.values)

    result = cv + cv2

    result2 = result - cv2

    assert np.all(np.equal(result2.runs, cv.runs))
    assert np.allclose(result2.values, cv.values)


nearest_command = "bedtools closest {} {} -t first -d -a <(sort -k1,1 -k2,2n {}) -b <(sort -k1,1 -k2,2n {})"
                    # "bedtools closest {} -t first -io -d -a <(sort -k1,1 -k2,2n {}) -b <(sort -k1,1 -k2,2n {})"]
nearest_hows = [None, "previous", "next"]
overlaps = [True, False]
strandedness = [False, "same", "opposite"]


@pytest.mark.bedtools
@pytest.mark.parametrize("nearest_how,overlap,strandedness", product(nearest_hows, overlaps, strandedness))
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min(), gr2=dfs_min())
def test_nearest(gr, gr2, nearest_how, overlap, strandedness):

    print(gr.df)
    print(gr2.df)

    bedtools_strand = {False: "", "same": "-s", "opposite": "-S"}[strandedness]
    bedtools_overlap = {True: "", False: "-io"}[overlap]

    result_df = None
    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = "{}/f1.bed".format(temp_dir)
        f2 = "{}/f2.bed".format(temp_dir)
        gr.df.to_csv(f1, sep="\t", header=False, index=False)
        gr2.df.to_csv(f2, sep="\t", header=False, index=False)

        cmd = nearest_command.format(bedtools_overlap, bedtools_strand, f1, f2)
        # print(cmd)
        result = subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()

        bedtools_df = pd.read_table(StringIO(result), header=None, squeeze=True, names="Chromosome Start End Name Score Strand Chromosome_b Start_b End_b Name_b Score_b Strand_b Distance".split())
        print("bedtools_df", bedtools_df)
        # raise
        if not bedtools_df.empty:
            bedtools_df = bedtools_df[bedtools_df.Distance != -1]["Chromosome Start End Strand Distance".split()].sort_values("Distance")

    result = gr.nearest(gr2, overlap=overlap, strandedness=strandedness)
    if not result.df.empty:

        print("pyranges_df", result.df)
        result_df = result.df.sort_values("Distance")["Chromosome Start End Strand Distance".split()]
        # print("bedtools_df", bedtools_df)
        # print("pyranges", pyranges_distances)

        assert_df_equal(result_df, bedtools_df)
    else:
        result_df = pd.DataFrame(columns="Chromosome Start End Strand Distance".split())
        assert bedtools_df.empty

    # if not result.df.empty:
    #     pyranges_distances = result_df.Distance.tolist()
    # else:
    #     pyranges_distances = []

    # print("bedtools", bedtools_distances)



jaccard_strandedness = [False, "same"]
jaccard_command = "bedtools jaccard {}  -a <(sort -k1,1 -k2,2n {}) -b <(sort -k1,1 -k2,2n {})"
@pytest.mark.bedtools
@pytest.mark.parametrize("strandedness", jaccard_strandedness)
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min(), gr2=dfs_min())
def test_jaccard(gr, gr2, strandedness):

    print(gr.df)
    print(gr2.df)

    bedtools_strand = {False: "", "same": "-s", "opposite": "-S"}[strandedness]
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
join_strandedness = [False, "opposite", "same"]

@pytest.mark.bedtools
@pytest.mark.parametrize("strandedness", join_strandedness)
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
