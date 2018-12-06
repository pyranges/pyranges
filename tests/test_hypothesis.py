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
use_strand = st.booleans()

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

better_dfs_min = data_frames(index=indexes(dtype=np.int64, min_size=better_df_minsize, unique=True, elements=lengths),
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
                                                        # column("Name", elements=names),
                                                        # column("Score", elements=scores),
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
    # strand = draw(use_strand)
    df.loc[:, "End"] += df.Start

    df.insert(3, "Name", "a")
    df.insert(4, "Score", 0)

    # if not strand:
    #     df = df.drop("Strand", axis=1)

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
@pytest.mark.parametrize("strandedness", ["same", "opposite", False]) #
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(gr=dfs_min(), gr2=dfs_min())
def test_subtraction(gr, gr2, strandedness):

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

    if not bedtools_df.empty or not result.df.empty:
        assert_df_equal(result.df, bedtools_df)
    else:
        assert bedtools_df.empty == result.df.empty


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

    bedtools_strand = {False: "", "same": "-s", "opposite": "-S"}[strandedness]
    bedtools_overlap = {True: "", False: "-io"}[overlap]

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
