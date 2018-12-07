
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
        # need to sort because bedtools sometimes gives the result in non-natsorted chromosome order!
        if result.stranded:
            assert_df_equal(result.df.sort_values("Chromosome Start Strand".split()), bedtools_df.sort_values("Chromosome Start Strand".split()))
        else:
            assert_df_equal(result.df.sort_values("Chromosome Start".split()), bedtools_df.sort_values("Chromosome Start".split()))
    else:
        assert bedtools_df.empty == result.df.empty
