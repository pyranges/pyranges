import pytest

from hypothesis import given, settings, HealthCheck
from hypothesis import reproduce_failure  # nosec

import tempfile
import subprocess  # nosec
from io import StringIO

import pandas as pd

from tests.helpers import assert_df_equal
from tests.hypothesis_helper import dfs_min

from os import environ

if environ.get("TRAVIS"):
    max_examples = 100
    deadline = None
else:
    max_examples = 100
    deadline = None

merge_command = "bedtools merge -o first -c 6 {} -i <(sort -k1,1 -k2,2n {})"


@pytest.mark.bedtools
@pytest.mark.parametrize("strand", [True, False])
@settings(
    max_examples=max_examples,
    deadline=deadline,
    suppress_health_check=HealthCheck.all())
@given(gr=dfs_min())  # nosec
def test_cluster(gr, strand):

    bedtools_strand = {True: "-s", False: ""}[strand]

    print(gr)

    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = "{}/f1.bed".format(temp_dir)
        gr.df.to_csv(f1, sep="\t", header=False, index=False)

        cmd = merge_command.format(bedtools_strand, f1)
        print(cmd)

        # ignoring bandit security warning. All strings created by test suite
        result = subprocess.check_output(  # nosec
            cmd, shell=True, executable="/bin/bash").decode()  # nosec

        if not strand:
            bedtools_df = pd.read_csv(
                StringIO(result),
                sep="\t",
                header=None,
                squeeze=True,
                usecols=[0, 1, 2],
                names="Chromosome Start End".split(),
                dtype={"Chromosome": "category"})
        else:
            bedtools_df = pd.read_csv(
                StringIO(result),
                sep="\t",
                header=None,
                squeeze=True,
                names="Chromosome Start End Strand".split(),
                dtype={"Chromosome": "category"})

    print("bedtools_df\n", bedtools_df)
    result = gr.merge(strand=strand)
    print("result\n", result.df)

    if not bedtools_df.empty:
        # need to sort because bedtools sometimes gives the result in non-natsorted chromosome order!
        if result.stranded:
            assert_df_equal(
                result.df.sort_values("Chromosome Start Strand".split()),
                bedtools_df.sort_values("Chromosome Start Strand".split()))
        else:
            assert_df_equal(
                result.df.sort_values("Chromosome Start".split()),
                bedtools_df.sort_values("Chromosome Start".split()))
    else:
        assert bedtools_df.empty == result.df.empty
