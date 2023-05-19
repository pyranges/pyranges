import pandas as pd

import pyranges as pr

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
    x = x[x.Start != x.End]

    return x.sort_values("Start End".split())  # .reset_index(drop=True)


def _introns_correct(full, genes, exons, introns, by):
    """Testing that introns:

    1: ends larger than starts
    2: the intersection of the computed introns and exons per gene are 0
    3: that the number of introns overlapping each gene is the same as number of introns per gene
    4 & 5: that the intron positions are the same as the ones computed with the slow, but correct algo
    """
    id_column = by_to_id[by]
    if len(introns):
        assert (introns.Start < introns.End).all(), str(introns[(introns.Start >= introns.End)])

    expected_results = {}
    based_on = {}
    for gene_id, gdf in full.groupby(id_column):  # #[full.gene_id.isin(["ENSG00000078808.16"])]
        # print("gdf " * 10)
        # print(gdf)
        if not len(gdf[gdf.Feature == "gene"]) or not len(gdf[gdf.Feature == "transcript"]):
            continue
        expected_results[gene_id] = compute_introns_single(gdf, by)

        based_on[gene_id] = pr.PyRanges(gdf[gdf.Feature.isin([by, "exon"])]).df

    if not len(introns):
        for v in expected_results.values():
            assert v.empty
        return  # test passed

    for gene_id, idf in introns.groupby(id_column):
        idf = idf.sort_values("Start End".split())
        if gene_id not in expected_results:
            continue
        expected = expected_results[gene_id]
        exons = pr.PyRanges(based_on[gene_id]).subset(lambda df: df.Feature == "exon").merge(by=id_column)
        genes = pr.PyRanges(based_on[gene_id]).subset(lambda df: df.Feature == by)
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
