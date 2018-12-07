import pytest

from hypothesis import given, settings, reproduce_failure, unlimited, HealthCheck, seed
from hypothesis.extra.pandas import data_frames, columns, range_indexes, column, indexes
from hypothesis.extra.numpy import arrays
import hypothesis.strategies as st

from pyranges import PyRanges

import numpy as np

lengths = st.integers(min_value=1, max_value=int(1e7))
small_lengths = st.integers(min_value=1, max_value=int(1e4))

strands = st.sampled_from("+ -".split())
names = st.text("abcdefghijklmnopqrstuvxyz", min_size=1)
scores = st.integers(min_value=0, max_value=256)

chromosomes = st.sampled_from(["chr{}".format(str(e)) for e in list(range(1, 23)) + "X Y M".split()])
chromosomes_small = st.sampled_from(["chr1"])
cs = st.one_of(chromosomes, chromosomes_small)

runlengths = data_frames(index=indexes(dtype=np.int64, min_size=1, unique=True),
                         columns=[column("Runs", st.integers(min_value=1, max_value=int(1e7))),
                                  # must have a min/max on floats because R S4vectors translates too big ones into inf.
                                  # which is unequal to eg -1.79769e+308 so the tests fail
                                  column("Values", st.integers(min_value=-int(1e7), max_value=int(1e7)))])




better_dfs_min = data_frames(index=indexes(dtype=np.int64, min_size=1, unique=True, elements=lengths),
                             columns=[column("Chromosome", cs),
                                      column("Start", elements=lengths),
                                      column("End", elements=small_lengths),
                                      # column("Name", elements=names),
                                      # column("Score", elements=scores),
                                      column("Strand", strands)])


better_dfs_min_single_chromosome = data_frames(index=indexes(dtype=np.int64, min_size=1, unique=True, elements=lengths),
                                               columns=[column("Chromosome", chromosomes_small),
                                                        column("Start", elements=lengths),
                                                        column("End", elements=small_lengths),
                                                        # column("Name", elements=names),
                                                        # column("Score", elements=scores),
                                      column("Strand", strands)])


runlengths_same_length_integers = data_frames(index=indexes(dtype=np.int64, min_size=1, unique=True),
                                              columns=[column("Runs", st.integers(min_value=1, max_value=int(1e4))),
                                              column("Values", st.integers(min_value=1, max_value=int(1e4))),
                                              column("Values2", st.integers(min_value=1, max_value=int(1e4)))])


@st.composite
def dfs_min(draw):
    df = draw(better_dfs_min)
    # strand = draw(use_strand)
    df.loc[:, "End"] += df.Start

    df.insert(3, "Name", "a")
    df.insert(4, "Score", 0)

    # if not strand:
    #     df = df.drop("Strand", axis=1)

    gr = PyRanges(df, extended=True)
    # gr = PyRanges(df)

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
