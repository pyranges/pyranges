import subprocess  # nosec
import tempfile
from io import StringIO

import numpy as np
import pandas as pd
import pytest
from hypothesis import reproduce_failure  # noqa: F401
from hypothesis import HealthCheck, given, settings
from natsort import natsorted  # type: ignore

import pyranges as pr
from tests.helpers import assert_df_equal
from tests.property_based.hypothesis_helper import deadline, df_data, dfs_min, dfs_min_with_id, max_examples, selector

# if environ.get("TRAVIS"):
#     max_examples = 100
#     deadline = None
# else:
#     max_examples = 1000
#     deadline = None

merge_command = "bedtools merge -o first,count -c 6,1 {} -i <(sort -k1,1 -k2,2n {})"


@pytest.mark.bedtools
@pytest.mark.parametrize("strand", [True, False])
@settings(
    max_examples=max_examples,
    deadline=deadline,
    suppress_health_check=HealthCheck.all(),
)
@given(gr=dfs_min())  # pylint: disable=no-value-for-parameter
def test_merge(gr, strand):
    bedtools_strand = {True: "-s", False: ""}[strand]

    print(gr)

    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = "{}/f1.bed".format(temp_dir)
        gr.df.to_csv(f1, sep="\t", header=False, index=False)

        cmd = merge_command.format(bedtools_strand, f1)
        print(cmd)

        # ignoring bandit security warning. All strings created by test suite
        result = subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()  # nosec  # nosec

        print("result" * 10)
        print(result)

        if not strand:
            print("if not " * 10)
            bedtools_df = pd.read_csv(
                StringIO(result),
                sep="\t",
                header=None,
                usecols=[0, 1, 2, 4],
                names="Chromosome Start End Count".split(),
                dtype={"Chromosome": "category"},
            )
        else:
            bedtools_df = pd.read_csv(
                StringIO(result),
                sep="\t",
                header=None,
                names="Chromosome Start End Strand Count".split(),
                dtype={"Chromosome": "category"},
            )

    print("bedtools_df\n", bedtools_df)
    result = gr.merge(strand=strand, count=True)
    print("result\n", result.df)

    if not bedtools_df.empty:
        # need to sort because bedtools sometimes gives the result in non-natsorted chromosome order!
        if result.stranded:
            assert_df_equal(
                result.df.sort_values("Chromosome Start Strand".split()),
                bedtools_df.sort_values("Chromosome Start Strand".split()),
            )
        else:
            assert_df_equal(
                result.df.sort_values("Chromosome Start".split()),
                bedtools_df.sort_values("Chromosome Start".split()),
            )
    else:
        assert bedtools_df.empty == result.df.empty


cluster_command = "bedtools cluster {} -i <(sort -k1,1 -k2,2n {})"


@pytest.mark.bedtools
@pytest.mark.parametrize("strand", [True, False])
@settings(
    max_examples=max_examples,
    deadline=deadline,
    suppress_health_check=HealthCheck.all(),
)
@given(gr=dfs_min())  # pylint: disable=no-value-for-parameter
def test_cluster(gr, strand):
    bedtools_strand = {True: "-s", False: ""}[strand]

    print(gr)

    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = "{}/f1.bed".format(temp_dir)
        gr.df.to_csv(f1, sep="\t", header=False, index=False)

        cmd = cluster_command.format(bedtools_strand, f1)
        print(cmd)

        # ignoring bandit security warning. All strings created by test suite
        result = subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()  # nosec  # nosec

        bedtools_df = pd.read_csv(
            StringIO(result),
            sep="\t",
            header=None,
            names="Chromosome Start End Name Score Strand Cluster".split(),
            dtype={"Chromosome": "category"},
        )

    print("bedtools_df\n", bedtools_df)

    print("gr\n", gr)
    result = gr.cluster(strand=strand)
    print("result\n", result[["Cluster"]])
    print("result\n", result.df)

    if not bedtools_df.empty:
        # need to sort because bedtools sometimes gives the result in non-natsorted chromosome order!
        if result.stranded:
            sort_values = "Chromosome Start Strand".split()
        else:
            sort_values = "Chromosome Start".split()

        result_df = result.df.sort_values(sort_values)
        print(bedtools_df)
        bedtools_df = bedtools_df.sort_values(sort_values)

        cluster_ids = {
            k: v
            for k, v in zip(
                result_df.Cluster.drop_duplicates(),
                bedtools_df.Cluster.drop_duplicates(),
            )
        }

        # bedtools gives different cluster ids than pyranges
        result_df.Cluster.replace(cluster_ids, inplace=True)

        bedtools_df.Cluster = bedtools_df.Cluster.astype("int32")
        assert_df_equal(result_df.drop("Cluster", axis=1), bedtools_df.drop("Cluster", axis=1))
    else:
        assert bedtools_df.empty == result.df.empty


@pytest.mark.parametrize("strand", [True, False])
@settings(
    max_examples=max_examples,
    deadline=deadline,
    suppress_health_check=HealthCheck.all(),
)
@given(gr=dfs_min_with_id())  # pylint: disable=no-value-for-parameter
def test_cluster_by(gr, strand):
    result = gr.cluster(by="ID", strand=strand).df
    print(result)
    df = gr.df

    if strand:
        groupby = ["Chromosome", "Strand", "ID"]
    else:
        groupby = ["Chromosome", "ID"]

    grs = []

    for _, gdf in natsorted(df.groupby(groupby)):
        grs.append(pr.PyRanges(gdf))

    clusters = [gr.cluster(strand=strand) for gr in grs]
    i = 1
    new_clusters = []
    for c in clusters:
        print("c")
        print(c)
        c.Cluster = i
        i += 1
        new_clusters.append(c)

    expected = pr.concat(new_clusters).df
    expected.loc[:, "Cluster"] = expected.Cluster.astype(np.int32)
    # expected = expected.drop_duplicates()

    print(expected)
    print(result)

    assert_df_equal(result.drop("Cluster", axis=1), expected.drop("Cluster", axis=1))


@pytest.mark.parametrize("strand", [True, False])
@settings(
    max_examples=max_examples,
    deadline=deadline,
    suppress_health_check=HealthCheck.all(),
)
@given(gr=dfs_min_with_id())  # pylint: disable=no-value-for-parameter
def test_merge_by(gr, strand):
    print(gr)
    result = gr.merge(by="ID").df.drop("ID", axis=1)

    df = gr.df

    grs = []
    for _, gdf in df.groupby("ID"):
        grs.append(pr.PyRanges(gdf))

    expected = pr.concat([gr.merge() for gr in grs]).df

    print(expected)
    print(result)

    assert_df_equal(result, expected)


makewindows_command = "bedtools makewindows -w 10 -b <(sort -k1,1 -k2,2n {})"


@pytest.mark.bedtools
@settings(
    max_examples=max_examples,
    print_blob=True,
    deadline=deadline,
    suppress_health_check=HealthCheck.all(),
)
@given(gr=dfs_min())  # pylint: disable=no-value-for-parameter
# @reproduce_failure('5.5.4', b'AXicY2RgYGAEIzgAsRkBAFsABg==')
def test_windows(gr):
    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = "{}/f1.bed".format(temp_dir)
        gr.df.to_csv(f1, sep="\t", header=False, index=False)

        cmd = makewindows_command.format(f1)
        print(cmd)

        # ignoring bandit security warning. All strings created by test suite
        result = subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()  # nosec  # nosec

        bedtools_df = pd.read_csv(
            StringIO(result),
            sep="\t",
            header=None,
            names="Chromosome Start End".split(),
            dtype={"Chromosome": "category"},
        )

    print("bedtools_df\n", bedtools_df)

    result = gr.window(10)["Chromosome Start End".split()].unstrand()
    print("result\n", result.df)

    if not bedtools_df.empty:
        assert_df_equal(result.df, bedtools_df)
    else:
        assert bedtools_df.empty == result.df.empty


@pytest.mark.parametrize("strand", [True, False])
@settings(
    max_examples=max_examples,
    deadline=deadline,
    suppress_health_check=HealthCheck.all(),
)
@given(gr=df_data())  # pylint: disable=no-value-for-parameter
def test_init(gr, strand):
    c, s, e, strands = gr

    if strand:
        pr.PyRanges(chromosomes=c, starts=s, ends=e, strands=strands)
    else:
        pr.PyRanges(chromosomes=c, starts=s, ends=e)


chipseq = pr.data.chipseq()


@settings(
    max_examples=max_examples,
    deadline=deadline,
    suppress_health_check=HealthCheck.all(),
)
@given(selector=selector())  # pylint: disable=no-value-for-parameter
def test_getitem(selector):
    if len(selector) == 3:
        a, b, c = selector
        chipseq[a, b, c]
    elif len(selector) == 2:
        a, b = selector
        chipseq[a, b]
    elif len(selector) == 1:
        a = selector[0]
        chipseq[a]
    elif len(selector) == 0:
        pass
    else:
        raise Exception("Should never happen")


@pytest.mark.bedtools
@settings(
    max_examples=max_examples,
    deadline=deadline,
    print_blob=True,
    suppress_health_check=HealthCheck.all(),
)
@given(gr=dfs_min())  # pylint: disable=no-value-for-parameter
def test_summary(gr):
    print(gr.to_example())
    # merely testing that it does not error
    # contents are just (pandas) dataframe.describe()
    gr.summary()
