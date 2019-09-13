import pandas as pd
import pyranges as pr

from hypothesis import given, settings
from hypothesis import reproduce_failure  # pylint: disable=unused-import

from tests.hypothesis_helper import deadline, slow_max_examples

by_to_id = {"gene": "gene_id", "transcript": "transcript_id"}


def compute_introns_single(df, by):
    id_column = by_to_id[by]
    df = df.sort_values("Start End".split())
    g = df[df.Feature == by]
    x = pr.PyRanges(df[df.Feature == "exon"]).merge(by=id_column, strand=False)
    if not len(x) or not len(g):
        return
    x.Strand = "-"
    x = x.df

    print("g " * 100)
    print(g)
    print("x " * 100)
    print(x)

    if g.empty or x.empty:
        return pd.DataFrame()

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
    # if x.Start.iloc[0] == g.Start.iloc[0]:
    #     print("Start equals")
    # x = x.drop(x.head(1).index)
    # if x.End.iloc[-1] == g.End.iloc[0]:
    #     print("End equals")
    #     x = x.drop(x.head(1).index)
    # x = x[~(x.Start == g.Start.iloc[0])]
    # x = x[~(x.End == g.End.iloc[0])]
    x = x[x.Start != x.End]

    return x.sort_values("Start End".split())  #.reset_index(drop=True)


# actual   Chromosome     Start       End Strand
# 0       chr1  13163233  13164184      -
# 1       chr1  13164760  13165156      -
# 2       chr1  13165468  13167115      -
# 3       chr1  13167193  13194517      -
# 4       chr1  13194595  13196244      -
# 5       chr1  13196556  13196952      -
# 6       chr1  13197528  13198478      -
# expected   Chromosome     Start       End Strand
# 5       chr1  13163233  13164184      -
# 6       chr1  13164760  13165156      -
# 7       chr1  13165468  13167115      -
# 8       chr1  13167193  13199727      -
# 2       chr1  13196556  13196952      +
# 3       chr1  13197528  13198478      +
# 4       chr1  13199727  13161985      -


def _introns_correct(full, genes, exons, introns, by):
    """Testing that introns:

    1: ends larger than starts
    2: the intersection of the computed introns and exons per gene are 0
    3: that the number of introns overlapping each gene is the same as number of introns per gene
    4 & 5: that the intron positions are the same as the ones computed with the slow, but correct algo"""
    id_column = by_to_id[by]
    if len(introns):
        assert (introns.Start < introns.End).all(), str(
            introns[(introns.Start >= introns.End)])

    expected_results = {}
    based_on = {}
    for gene_id, gdf in full.groupby(
            id_column):  # #[full.gene_id.isin(["ENSG00000078808.16"])]
        # print("gdf " * 10)
        # print(gdf)
        if not len(gdf[gdf.Feature == "gene"]) or not len(
                gdf[gdf.Feature == "transcript"]):
            continue
        expected_results[gene_id] = compute_introns_single(gdf, by)

        based_on[gene_id] = pr.PyRanges(gdf[gdf.Feature.isin([by, "exon"])]).df

    if not len(introns):
        for v in expected_results.values():
            assert v.empty
        return  # test passed

    for gene_id, idf in introns.groupby(id_column):
        idf = idf.sort_values("Start End".split())
        if not gene_id in expected_results:
            continue
        expected = expected_results[gene_id]
        exons = pr.PyRanges(
            based_on[gene_id]).subset(lambda df: df.Feature == "exon").merge(
                by=id_column)
        genes = pr.PyRanges(
            based_on[gene_id]).subset(lambda df: df.Feature == by)
        print("exons", exons)
        print("based_on", based_on[gene_id])
        print("actual", idf["Chromosome Start End Strand".split()])
        print("expected", expected["Chromosome Start End Strand".split()])
        _introns = pr.PyRanges(idf)
        assert len(exons.intersect(_introns)) == 0
        assert len(genes.intersect(_introns)) == len(_introns)
        assert list(idf.Start) == list(expected.Start), "not equal"
        assert list(idf.End) == list(expected.End), "not equal"


by = ["gene", "transcript"]

# @pytest.mark.parametrize("by", by)
# @settings(
#     max_examples=slow_max_examples,
#     deadline=deadline,
#     suppress_health_check=HealthCheck.all())
# @given(gr=genomicfeature())  # pylint: disable=no-value-for-parameter
# # @reproduce_failure('4.15.0', b'AAIAAAABAP8A')
# # @reproduce_failure('4.15.0', b'AAIAAAABAOwBAO0A')
# def test_introns(gr, by):

#     id_column = by_to_id[by]

#     delim = "-" * 20
#     print("\n" + delim +" new test " + delim)

#     if not len(gr[gr.Feature == "gene"]) or not len(gr[gr.Feature == "exon"]):
#         pass

#     else:

#         genes = gr[gr.Feature.isin(["gene"])].df
#         exons = gr[gr.Feature.isin(["exon"])].merge(by=id_column).df
#         introns = gr.features.introns(by=by).df

#         print(genes)
#         print(exons)
#         print(introns)

#         _introns_correct(gr.df, genes, exons, introns, by)


def test_introns_single():

    "Assert that our fast method of computing introns is the same as the slow, correct one in compute_introns_single"

    gr = pr.data.gencode_gtf()[["gene_id", "Feature"]]
    exons = gr[gr.Feature == "exon"].merge(by="gene_id")
    exons.Feature = "exon"
    exons = exons.df
    df = pd.concat([gr[gr.Feature == "gene"].df, exons], sort=False)
    print(df)

    for gid, gdf in df.groupby("gene_id"):
        print("-------" * 20)
        print(gid)
        print(gdf)
        print("gdf", len(gdf))
        expected = compute_introns_single(gdf, by="gene")
        print("expected", len(expected))
        actual = pr.PyRanges(gdf).features.introns().df
        print("actual", len(actual))
        if actual.empty:
            assert expected.empty
            continue

        assert list(expected.Start) == list(actual.Start)
        assert list(expected.End) == list(actual.End)
