from io import StringIO
import pandas as pd
import pyranges as pr
import pytest

from hypothesis import given, settings, HealthCheck  #, assume
from hypothesis import reproduce_failure  # pylint: disable=unused-import

from tests.hypothesis_helper import genomicfeature

from tests.hypothesis_helper import max_examples, deadline, slow_max_examples

by_to_id = {"gene": "gene_id", "transcript": "transcript_id"}

def compute_introns_single(df, by):
    id_column = by_to_id[by]
    g = df[df.Feature == by]
    x = pr.PyRanges(df[df.Feature == "exon"]).merge(by=id_column).df

    x = pd.concat([x, x.tail(1)], sort=False).reset_index(drop=True)
    cols = list(x.columns)
    original_cols = cols[::1]
    cols[2] = "Start"
    cols[1] = "End"
    x = x[cols]
    x.columns = original_cols
    x.loc[:, "Start"] = x.Start.shift()
    x.at[x.index[0], "Start"] = g.Start.iloc[0]
    x.at[x.index[-1], "End"] = g.End.iloc[0]
    x.loc[:, "Start"] = x.Start.astype(int)
    x = x[x.Start != x.End]

    return x


# def test_introns_single():

#     "Assert that our fast method of computing introns is the same as the slow, correct one in compute_introns_single"

#     gr = pr.data.gencode_gtf()[["gene_id", "Feature"]]
#     exons = gr[gr.Feature == "exon"].merge(by="gene_id")
#     exons.Feature = "exon"
#     exons = exons.df
#     df = pd.concat([gr[gr.Feature == "gene"].df, exons], sort=False)
#     print(df)

#     for gid, gdf in df.groupby("gene_id"):
#         print("-------" * 20)
#         print(gid)
#         print(gdf)
#         print("gdf", len(gdf))
#         expected = compute_introns_single(gdf)
#         print("expected", len(expected))
#         actual = pr.PyRanges(gdf).features.introns().df
#         print("actual", len(actual))
#         if actual.empty:
#             assert expected.empty
#             continue

#         assert list(expected.Start) == list(actual.Start)
#         assert list(expected.End) == list(actual.End)



def _introns_correct(full, genes, exons, introns, by):

    """Testing that introns:

    1: ends larger than starts
    2: the intersection of the computed introns and exons per gene are 0
    3: that the number of introns overlapping each gene is the same as number of introns per gene
    4 & 5: that the intron positions are the same as the ones computed with the slow, but correct algo"""
    id_column = by_to_id[by]
    if len(introns):
        assert (introns.Start < introns.End).all(), str(introns[(introns.Start >= introns.End)])

    expected_results = {}
    based_on = {}
    for gene_id, gdf in full.groupby(id_column): # #[full.gene_id.isin(["ENSG00000078808.16"])]
        expected_results[gene_id] = compute_introns_single(gdf, by)
        based_on[gene_id] = pr.PyRanges(gdf[gdf.Feature.isin([by, "exon"])]).df

    if not len(introns):
        for v in expected_results.values():
            assert v.empty
        return "okay"

    for gene_id, idf in introns.groupby(id_column):
        expected = expected_results[gene_id]
        # print("based_on", based_on[gene_id])
        exons = pr.PyRanges(based_on[gene_id]).subset(lambda df: df.Feature == "exon").merge(by=id_column)
        genes = pr.PyRanges(based_on[gene_id]).subset(lambda df: df.Feature == by)
        # print("exons", exons)
        # print("actual", idf)
        # print("expected", expected)
        _introns = pr.PyRanges(idf)
        assert len(exons.intersect(_introns)) == 0
        assert len(genes.intersect(_introns)) == len(_introns)
        assert list(idf.Start) == list(expected.Start), "not equal"
        assert list(idf.End) == list(expected.End), "not equal"

by = ["gene", "transcript"]



@pytest.mark.parametrize("by", by)
@settings(
    max_examples=slow_max_examples,
    deadline=deadline,
    suppress_health_check=HealthCheck.all())
@given(gr=genomicfeature())  # pylint: disable=no-value-for-parameter
def test_introns(gr, by):

    id_column = by_to_id[by]

    delim = "-" * 20
    print("\n" + delim +" new test " + delim)
    print(gr)

    genes = gr[gr.Feature.isin(["gene"])].df
    exons = gr[gr.Feature.isin(["exon"])].merge(by=id_column).df
    introns = gr.features.introns(by=by).df

    _introns_correct(gr.df, genes, exons, introns, by)

# @pytest.mark.parametrize("by", ["transcript"])
# def test_introns(by):

#     id_column = by_to_id[by]

#     gr = pr.data.gencode_gtf()[[id_column, "Feature"]].sort(id_column)
#     genes = gr[gr.Feature.isin(["gene"])].df
#     exons = gr[gr.Feature.isin(["exon"])].merge(by=id_column).df
#     introns = gr.features.introns(by=by).df

#     _introns_correct(gr.df, genes, exons, introns, by)



# @pytest.mark.parametrize("by", by)
# def test_introns_gencode(by):

#     gr = pr.data.gencode_gtf()
#     introns = gr.features.introns(by=by)


# @pytest.mark.parametrize("by", by)
# def test_introns_ensembl(by):

#     gr = pr.data.ensembl_gtf()
#     introns = gr.features.introns(by=by)

# @pytest.mark.parametrize("by", by)
# def test_introns_ucsc(by):

#     gr = pr.data.ucsc_bed()
#     introns = gr.features.introns(by=by)


# @pytest.fixture
# def simple_exons():

#     c = """Chromosome    Start      End Strand             gene_id Feature
# chr1  1471769  1472089      +  ENSG00000160072.19 exon
# chr1  1477274  1477350      +  ENSG00000160072.19 exon
# chr1  1478026  1478745      +  ENSG00000160072.19 exon
# chr1  1479049  1479108      +  ENSG00000160072.19 exon
# chr1  1480867  1480936      +  ENSG00000160072.19 exon
# chr1  1482138  1482303      +  ENSG00000160072.19 exon
# chr1  1482545  1482662      +  ENSG00000160072.19 exon
# chr1  1483485  1485171      +  ENSG00000160072.19 exon
# chr1  1485782  1486235      +  ENSG00000160072.19 exon
# chr1  1486544  1486668      +  ENSG00000160072.19 exon
# chr1  1487863  1487914      +  ENSG00000160072.19 exon
# chr1  1489204  1490671      +  ENSG00000160072.19 exon
# chr1  1495485  1497848      +  ENSG00000160072.19 exon"""

#     return pr.PyRanges(pd.read_csv(StringIO(c), sep=r"\s+"))


# @pytest.fixture
# def simple_gene():

#     c = """Chromosome    Start      End Strand             gene_id Feature
# chr1    1471769  1497848      +  ENSG00000160072.19  gene"""

#     return pr.PyRanges(pd.read_csv(StringIO(c), sep=r"\s+"))


# def test_simple_introns(simple_exons, simple_gene):

#     gr = pr.concat([simple_exons, simple_gene])
#     introns = gr.features.introns()
#     expected_starts = [1472089, 1477350, 1478745, 1479108, 1480936, 1482303, 1482662, 1485171, 1486235, 1486668, 1487914, 1490671]
#     expected_ends = [1477274, 1478026, 1479049, 1480867, 1482138, 1482545, 1483485, 1485782, 1486544, 1487863, 1489204, 1495485]
#     assert list(introns.Start) == expected_starts
#     assert list(introns.End) ==  expected_ends
#     assert len(introns.intersect(simple_exons)) == 0
