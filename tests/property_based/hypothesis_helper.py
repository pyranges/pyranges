from os import environ

import hypothesis.strategies as st
import numpy as np
import pandas as pd
from hypothesis.extra.pandas import column, data_frames, indexes

import pyranges as pr
from pyranges import PyRanges

max_examples = 15
slow_max_examples = 5
deadline = None

lengths = st.integers(min_value=1, max_value=int(1e7))
small_lengths = st.integers(min_value=1, max_value=int(1e4))

strands = st.sampled_from("+ -".split())
single_strand = st.sampled_from(["+"])
names = st.text("abcdefghijklmnopqrstuvxyz", min_size=1)
scores = st.integers(min_value=0, max_value=256)

datatype = st.sampled_from([pd.Series, np.array, list])

feature_data = st.sampled_from(["ensembl_gtf", "gencode_gtf", "ucsc_bed"])

chromosomes = st.sampled_from(["chr{}".format(str(e)) for e in list(range(1, 23)) + "X Y M".split()])
chromosomes_small = st.sampled_from(["chr1"])
cs = st.one_of(chromosomes, chromosomes_small)

runlengths = data_frames(
    index=indexes(dtype=np.int64, min_size=1, unique=True),
    columns=[
        column("Runs", st.integers(min_value=1, max_value=int(1e7))),
        # must have a min/max on floats because R S4vectors translates too big ones into inf.
        # which is unequal to eg -1.79769e+308 so the tests fail
        column("Values", st.integers(min_value=-int(1e7), max_value=int(1e7))),
    ],
)

better_dfs_no_min = data_frames(
    index=indexes(dtype=np.int64, min_size=0, unique=True, elements=lengths),
    columns=[
        column("Chromosome", cs),
        column("Start", elements=lengths),
        column("End", elements=small_lengths),
        # column("Name", elements=names),
        # column("Score", elements=scores),
        column("Strand", strands),
    ],
)

better_dfs_min = data_frames(
    index=indexes(dtype=np.int64, min_size=1, unique=True, elements=lengths),
    columns=[
        column("Chromosome", cs),
        column("Start", elements=lengths),
        column("End", elements=small_lengths),
        # column("Name", elements=names),
        # column("Score", elements=scores),
        column("Strand", strands),
    ],
)

better_dfs_min_2 = data_frames(
    index=indexes(dtype=np.int64, min_size=2, unique=True, elements=lengths),
    columns=[
        column("Chromosome", chromosomes_small),
        column("Start", elements=lengths),
        column("End", elements=small_lengths),
        # column("Name", elements=names),
        # column("Score", elements=scores),
        column("Strand", single_strand),
    ],
)

better_dfs_min_single_chromosome = data_frames(
    index=indexes(dtype=np.int64, min_size=1, unique=True, elements=lengths),
    columns=[
        column("Chromosome", chromosomes_small),
        column("Start", elements=lengths),
        column("End", elements=small_lengths),
        # column("Name", elements=names),
        # column("Score", elements=scores),
        column("Strand", strands),
    ],
)

runlengths_same_length_integers = data_frames(
    index=indexes(dtype=np.int64, min_size=1, unique=True),
    columns=[
        column("Runs", st.integers(min_value=1, max_value=int(1e4))),
        column("Values", st.integers(min_value=1, max_value=int(1e4))),
        column("Values2", st.integers(min_value=1, max_value=int(1e4))),
    ],
)


@st.composite
def dfs_min2(draw):  # nosec
    df = draw(better_dfs_min_2)
    # strand = draw(use_strand)
    df.loc[:, "End"] += df.Start

    df.insert(3, "Name", "a")
    df.insert(4, "Score", 0)

    # if not strand:
    #     df = df.drop("Strand", axis=1)

    gr = PyRanges(df, int64=True)
    # gr = PyRanges(df)

    # do not sort like this, use pyranges sort
    # np.random.seed(draw(st.integers(min_value=0, max_value=int(1e6))))
    # gr.df = df.reindex(np.random.permutation(df.index.values))

    return gr


@st.composite
def dfs_min(draw):  # nosec
    df = draw(better_dfs_min)
    # strand = draw(use_strand)
    df.loc[:, "End"] += df.Start

    df.insert(3, "Name", "a")
    df.insert(4, "Score", 0)

    # df.Start = df.Start.astype(np.int32)
    # df.End = df.End.astype(np.int32)
    # print(df.dtypes)
    # stranded = draw(st.booleans())
    # if not strand:
    #     df = df.drop("Strand", axis=1)

    gr = PyRanges(df, int64=True)
    # print(gr)
    # raise
    # gr = PyRanges(df)

    # do not sort like this, use pyranges sort
    # np.random.seed(draw(st.integers(min_value=0, max_value=int(1e6))))
    # gr.df = df.reindex(np.random.permutation(df.index.values))

    return gr


@st.composite
def dfs_no_min(draw):  # nosec
    df = draw(better_dfs_no_min)
    # strand = draw(use_strand)
    df.loc[:, "End"] += df.Start

    df.insert(3, "Name", "a")
    df.insert(4, "Score", 0)

    # stranded = draw(st.booleans())
    # if not strand:
    #     df = df.drop("Strand", axis=1)

    gr = PyRanges(df, int64=True)
    # gr = PyRanges(df)

    # do not sort like this, use pyranges sort
    # np.random.seed(draw(st.integers(min_value=0, max_value=int(1e6))))
    # gr.df = df.reindex(np.random.permutation(df.index.values))

    return gr


@st.composite
def dfs_min_with_id(draw):  # nosec
    df = draw(better_dfs_min)
    ids = df.Start
    # strand = draw(use_strand)
    df.loc[:, "End"] += df.Start

    df.insert(3, "Name", "a")
    df.insert(4, "Score", 0)
    df.insert(5, "ID", ids)

    # df.Start = df.Start.astype(np.int32)
    # df.End = df.End.astype(np.int32)
    # print(df.dtypes)
    # stranded = draw(st.booleans())
    # if not strand:
    #     df = df.drop("Strand", axis=1)

    gr = PyRanges(df)
    # print(gr)
    # raise
    # gr = PyRanges(df)

    # do not sort like this, use pyranges sort
    # np.random.seed(draw(st.integers(min_value=0, max_value=int(1e6))))
    # gr.df = df.reindex(np.random.permutation(df.index.values))

    return gr


@st.composite
def dfs_min_with_gene_id(draw):  # nosec
    df = draw(better_dfs_min)
    ids = df.Start
    # strand = draw(use_strand)
    df.loc[:, "End"] += df.Start

    df.insert(3, "Name", "a")
    df.insert(4, "Score", 0)
    df.insert(5, "gene_id", ids)
    df.insert(6, "exon_id", list(range(len(df))))

    # df.Start = df.Start.astype(np.int32)
    # df.End = df.End.astype(np.int32)
    # print(df.dtypes)
    # stranded = draw(st.booleans())
    # if not strand:
    #     df = df.drop("Strand", axis=1)

    gr = PyRanges(df)
    # print(gr)
    # raise
    # gr = PyRanges(df)

    # do not sort like this, use pyranges sort
    # np.random.seed(draw(st.integers(min_value=0, max_value=int(1e6))))
    # gr.df = df.reindex(np.random.permutation(df.index.values))

    return gr


@st.composite
def df_data(draw):
    df = draw(better_dfs_min)
    df.loc[:, "End"] += df.Start

    chromosome_type = draw(datatype)
    start_type = draw(datatype)
    end_type = draw(datatype)
    strand_type = draw(datatype)

    strand = strand_type(df.Strand)
    chromosome = chromosome_type(df.Chromosome)
    start = start_type(df.Start)
    end = end_type(df.End)

    return chromosome, start, end, strand


@st.composite
def dfs_min_single_chromosome(draw):
    df = draw(better_dfs_min_single_chromosome)
    df.loc[:, "End"] += df.Start
    df.insert(3, "Name", "a")
    df.insert(4, "Score", 0)

    return df


@st.composite
def genomicfeature(draw):
    dataset_name = draw(feature_data)
    print("dataset name " * 5, dataset_name)
    dataset = getattr(pr.data, dataset_name)()
    dataset = dataset[dataset.Feature.isin(["gene", "transcript", "exon"])]

    # subsetter = draw(arrays(np.bool, shape=len(dataset)))
    gene_ids = list(dataset.gene_id.drop_duplicates())
    genes = draw(st.lists(st.sampled_from(gene_ids), unique="True", min_size=1))
    dataset = dataset[dataset.gene_id.isin(genes)]

    return dataset


@st.composite
def selector(draw):
    df = draw(better_dfs_min)
    h = df.head(1)
    chromosome = h["Chromosome"].iloc[0]
    start = h["Start"].iloc[0]
    end = h["End"].iloc[0]
    strand = h["Strand"].iloc[0]

    chromosome_bool = draw(st.booleans())
    strand_bool = draw(st.booleans())
    start_bool = draw(st.booleans())
    end_bool = draw(st.booleans())

    _slice = {
        (True, True): slice(start, end),
        (True, False): slice(start, None),
        (False, True): slice(None, end),
        (False, False): slice(None, None),
    }[start_bool, end_bool]

    to_return = []
    if chromosome_bool:
        to_return.append(chromosome)
    if strand_bool:
        to_return.append(strand)
    if start_bool or end_bool:
        to_return.append(_slice)

    return to_return
