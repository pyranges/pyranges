# pragma: no cover
import pandas as pd

import pyranges as pr
import numpy as np

from pyranges.multithreaded import pyrange_apply
from sorted_nearest.src.introns import find_introns


def _outside_bounds(df, kwargs):

    chromsizes = kwargs.get("chromsizes")
    size = chromsizes[df.Chromosome.iloc[0]]
    clip = kwargs.get("clip", False)

    ends_outside = df.End > size
    if not clip:  # i.e. remove
        df = df[~ends_outside]
    else:
        starts_inside = df.Start <= size
        df.loc[ends_outside & starts_inside, "End"] = size

    return df


def genome_bounds(gr, chromsizes, clip=False):

    if not isinstance(chromsizes, dict):
        chromsizes = {
            k: v
            for k, v in zip(chromsizes.Chromosome, chromsizes.End)
        }

    return gr.apply(_outside_bounds, chromsizes=chromsizes, clip=clip)


def random(n=1000, length=100, chromsizes=None, strand=True, int64=False):

    if chromsizes is None:
        chromsizes = pr.data.chromsizes()
        df = chromsizes.df
    elif isinstance(chromsizes, dict):
        df = pd.DataFrame({"Chromosome": list(chromsizes.keys()), "End": list(chromsizes.values())})
    else:
        df = chromsizes.df

    p = df.End / df.End.sum()

    n_per_chrom = pd.Series(np.random.choice(
        df.index, size=n, p=p)).value_counts(sort=False).to_frame()
    n_per_chrom.insert(1, "Chromosome", df.loc[n_per_chrom.index].Chromosome)
    n_per_chrom.columns = "Count Chromosome".split()

    random_dfs = []
    for _, (count, chrom) in n_per_chrom.iterrows():
        r = np.random.randint(
            0, df[df.Chromosome == chrom].End - length, size=count)
        _df = pd.DataFrame({
            "Chromosome": chrom,
            "Start": r,
            "End": r + length
        })
        random_dfs.append(_df)

    random_df = pd.concat(random_dfs)
    if strand:
        s = np.random.choice("+ -".split(), size=n)
        random_df.insert(3, "Strand", s)

    return pr.PyRanges(random_df, int64=int64)


def _last_tile(df, kwargs):
    # do not need copy, since it is only used internally by
    # tile_genome
    # df = df.copy()
    sizes = kwargs.get("sizes")
    size = sizes[df.Chromosome.iloc[0]].End.iloc[0]
    df.loc[df.tail(1).index, "End"] = size

    return df


def tile_genome(genome, tile_size, tile_last=False):

    if isinstance(genome, dict):
        chromosomes, ends = list(genome.keys()), list(genome.values())
        df = pd.DataFrame({"Chromosome": chromosomes, "Start": 0, "End": ends})
        genome = pr.PyRanges(df)

    gr = genome.tile(tile_size)

    if not tile_last:
        gr = gr.apply(_last_tile, sizes=genome)

    return gr


def _keep_transcript_with_most_exons(df):

    transcripts_with_most_exons = []

    for _, gdf in df.groupby("gene_id"):

        max_exon = gdf.exon_number.max()
        max_transcript = gdf.loc[gdf.exon_number == max_exon].Transcript.iloc[0]

        max_rows = gdf.loc[gdf.Transcript == max_transcript]

        transcripts_with_most_exons.append(max_rows)

    return pd.concat(transcripts_with_most_exons).reset_index(drop=True)



def filter_transcripts(df, keep="most_exons"):

    return _keep_transcript_with_most_exons(df)


def _tss(df, slack=0):

    intype = df.Start.dtype

    tss_pos = df.loc[df.Strand == "+"]

    tss_neg = df.loc[df.Strand == "-"].copy()

    # pd.options.mode.chained_assignment = None
    tss_neg.loc[:, "Start"] = tss_neg.End

    # pd.options.mode.chained_assignment = "warn"
    tss = pd.concat([tss_pos, tss_neg], sort=False)
    tss["End"] = tss.Start
    tss.End = tss.End + 1 + slack
    tss.Start = tss.Start - slack
    tss.loc[tss.Start < 0, "Start"] = 0

    tss.index = range(len(tss))

    tss[["Start", "End"]] = tss[["Start", "End"]].astype(intype)

    return tss


def _tes(df, slack=0):

    intype = df.Start.dtype
    # df = self.df

    tes_pos = df.loc[df.Strand == "+"]

    tes_neg = df.loc[df.Strand == "-"].copy()

    # pd.options.mode.chained_assignment = None
    tes_neg.loc[:, "Start"] = tes_neg.End

    # pd.options.mode.chained_assignment = "warn"
    tes = pd.concat([tes_pos, tes_neg], sort=False)
    tes["Start"] = tes.End
    tes.End = tes.End + 1 + slack
    tes.Start = tes.Start - slack
    tes.loc[tes.Start < 0, "Start"] = 0

    tes.index = range(len(tes))

    tes[["Start", "End"]] = tes[["Start", "End"]].astype(intype)

    return tes


by_to_id = {"gene": "gene_id", "transcript": "transcript_id"}


def _introns2(df, exons, kwargs):

    if df.empty or exons.empty:
        return None

    original_order = df.columns
    by = kwargs["by"]
    id_column = by_to_id[by]

    exons = exons[["Start", "End", id_column]]
    genes = df[["Start", "End", id_column]]
    exons.columns = ["Start", "End", "gene_id"]
    genes.columns = ["Start", "End", "gene_id"]

    intersection = pd.Series(
        np.intersect1d(exons["gene_id"], genes["gene_id"]))
    if len(intersection) == 0:
        return None

    exons = exons[exons["gene_id"].isin(intersection)].reset_index(
        drop=True).sort_values(["gene_id", "Start"])
    genes = genes[genes["gene_id"].isin(intersection)].reset_index(
        drop=True).sort_values(["gene_id", "Start"])
    df = df[df[id_column].isin(intersection)].reset_index(
        drop=True)  #.sort_values(id_column)

    assert len(genes) == len(
        genes.drop_duplicates("gene_id")
    ), "The {id_column}s need to be unique to compute the introns.".format(
        id_column=id_column)

    exon_ids = (exons["gene_id"].shift() != exons["gene_id"])
    gene_ids = pd.Series(range(1, len(genes) + 1))
    df.insert(0, "__temp__", gene_ids)

    if len(exons) > 1 and exons["gene_id"].iloc[0] == exons["gene_id"].iloc[1]:
        exon_ids.iloc[0] = False
        exon_ids = exon_ids.cumsum() + 1
    else:
        exon_ids = exon_ids.cumsum()

    assert (gene_ids == exon_ids.drop_duplicates().values).all()
    starts, ends, ids = find_introns(genes.Start.values, genes.End.values,
                                     gene_ids.values, exons.Start.values,
                                     exons.End.values, exon_ids.values)

    introns = pd.DataFrame(
        data={
            "Chromosome": df.Chromosome.iloc[0],
            "Start": starts,
            "End": ends,
            "gene_id": ids
        })

    vc = introns["gene_id"].value_counts(sort=False).to_frame().reset_index()
    vc.columns = ["gene_id", "counts"]

    genes_without_introns = pd.DataFrame(
        data={
            "gene_id": np.setdiff1d(gene_ids.values, vc.gene_id.values),
            "counts": 0
        })

    vc = pd.concat([vc, genes_without_introns]).sort_values("gene_id")

    original_ids = np.repeat(vc.gene_id, vc.counts).to_frame()
    original_ids = original_ids.merge(
        df[["__temp__", id_column]],
        right_on="__temp__",
        left_on="gene_id",
        suffixes=("_drop", ""))
    original_ids = original_ids.drop(
        ["__temp__"] +
        [c for c in original_ids.columns if c.endswith("_drop")],
        axis=1).sort_values("gene_id")
    introns.loc[:, "gene_id"] = original_ids[id_column].values
    introns = introns.merge(
        df, left_on="gene_id", right_on=id_column, suffixes=("", "_dropme"))
    introns = introns.drop(
        [c for c in introns.columns if c.endswith("_dropme")], axis=1)

    if introns.Feature.dtype.name == "category" and not "intron" in introns.Feature.cat.categories:
        introns.Feature.cat.add_categories(["intron"], inplace=True)
    introns.loc[:, "Feature"] = "intron"

    introns = introns[original_order]

    return introns


class GenomicFeaturesMethods():

    pr = None

    def __init__(self, pr):

        self.pr = pr

    def tss(self, slack=0):

        pr = self.pr

        if not pr.stranded:
            raise Exception(
                "Cannot compute TSSes or TESes without strand info. Perhaps use slack() instead?"
            )

        pr = pr[pr.Feature == "transcript"]
        pr = pr.apply(lambda df: _tss(df, slack))

        pr.Feature = "tss"

        return pr

    def tes(self, slack=0):

        pr = self.pr

        if not pr.stranded:
            raise Exception(
                "Cannot compute TSSes or TESes without strand info. Perhaps use slack() instead?"
            )

        pr = pr[pr.Feature == "transcript"]
        pr = pr.apply(lambda df: _tes(df, slack))

        pr.Feature = "tes"

        return pr


    def introns(self, by="gene"):

        kwargs = {"by": by}
        kwargs = pr.pyranges.fill_kwargs(kwargs)

        assert by in ["gene", "transcript"]

        id_column = by_to_id[by]
        gr = self.pr.sort(id_column)

        if not len(gr):
            return pr.PyRanges()

        exons = gr.subset(lambda df: df.Feature == "exon")
        exons = exons.merge(by=id_column)

        by_gr = gr.subset(lambda df: df.Feature == by)

        result = pyrange_apply(_introns2, by_gr, exons, **kwargs)

        return pr.PyRanges(result)
