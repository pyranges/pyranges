"""Test that pyranges is able to download data from databases."""

import pytest
import pandas as pd

import pyranges.db as db


@pytest.fixture
def df_to_parse():

    return pd.read_csv(
        "tests/test_data/ucsc_df_to_parse.txt", header=0, index_col=0)


@pytest.mark.external()
def test_ucsc_genomes():

    df = db.ucsc.genomes()
    assert not df.isnull().sum().any()


@pytest.mark.external()
def test_ucsc_chromosome_sizes():

    df = db.ucsc.chromosome_sizes("hg19")
    assert not df.isnull().sum().any()


def test_ucsc_genes_parse(df_to_parse):

    gr = db.ucsc.parse_genes(df_to_parse)
    print(gr)


@pytest.mark.external()
def test_ucsc_genes(df_to_parse):

    gr = db.ucsc.genes_df(
        "hg19",
        'select chrom, txStart, txEnd, exonStarts, exonEnds, name, name2, strand from refGene limit 10;'
    )
    print(gr)


@pytest.mark.external()
@pytest.mark.ftp()
def test_gencode_genes():

    gr = db.gencode.genes("human", head=True)
    print(gr)
