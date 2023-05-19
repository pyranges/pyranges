import numpy as np
import pandas as pd
import pytest

from pyranges import PyRanges


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "bedtools: tests rely on",
    )

    config.addinivalue_line("markers", "explore: functionality not ready for prime-time")


@pytest.fixture(autouse=True, scope="session")
def numpy_seed():
    np.random.RandomState(42)


@pytest.fixture(autouse=True)
def add_np(doctest_namespace):
    doctest_namespace["np_seed"] = np.random.seed(0)


@pytest.fixture
def names():
    return "Chromosome  Start  End  Name Score Strand".split()


@pytest.fixture
def chip_10(names):
    df = pd.read_csv("tests/unit/chip_10.bed", header=None, names=names, sep="\t")

    gr = PyRanges(df)

    assert gr.stranded

    return gr


@pytest.fixture
def f1(names):
    df = pd.read_csv(
        "tests/unit/f1.bed",
        sep="\t",
        header=None,
        names="Chromosome  Start  End  Name Score Strand".split(),
    )

    return PyRanges(df)


@pytest.fixture
def f2(names):
    df = pd.read_csv("tests/unit/f2.bed", sep="\t", header=None, names=names)

    return PyRanges(df)


@pytest.fixture
def chromsizes():
    from io import StringIO

    df = pd.read_csv(
        StringIO(
            """
chr1                     249250621
chr2                     243199373
chr3                     198022430
chr4                     191154276
chr5                     180915260
chr6                     171115067
chr7                     159138663
chrX                     155270560
chr8                     146364022
chr9                     141213431
chr10                    135534747
chr11                    135006516
chr12                    133851895
chr13                    115169878
chr14                    107349540
chr15                    102531392
chr16                     90354753
chr17                     81195210
chr18                     78077248
chr20                     63025520
chrY                      59373566
chr19                     59128983
chr22                     51304566
chr21                     48129895"""
        ),
        sep=r"\s+",
        header=None,
        index_col=0,
    )

    df.insert(0, "Start", 0)
    df = df.reset_index()
    df.columns = ["Chromosome", "Start", "End"]

    return PyRanges(df)
