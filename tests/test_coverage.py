
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

# runlengths = data_frames(index=indexes(dtype=np.int64, min_size=1, unique=True),
#                          columns=[column("Runs", st.integers(min_value=1, max_value=int(1e7))),
#                                   # must have a min/max on floats because R S4vectors translates too big ones into inf.
#                                   # which is unequal to eg -1.79769e+308 so the tests fail
#                                   column("Values", st.integers(min_value=-int(1e7), max_value=int(1e7)))])
from tests.hypothesis_helper import dfs_min_single_chromosome

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


coverage_cmd = "Rscript --vanilla tests/compute_coverage.R {} {}"

@pytest.mark.r
@settings(max_examples=max_examples, deadline=deadline, timeout=unlimited, suppress_health_check=HealthCheck.all())
@given(df=dfs_min_single_chromosome())
def test_coverage(df):

    print("---" * 10)
    p = pr.PyRanges(df)
    print("pyranges\n", p)

    c = p.coverage(strand=False)["chr1"]

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
