import numpy as np
import pandas as pd
from sorted_nearest.src.introns import find_introns  # type: ignore

import pyranges as pr
from pyranges.multithreaded import pyrange_apply

__all__ = ["genome_bounds", "tile_genome", "GenomicFeaturesMethods"]


class GenomicFeaturesMethods:

    """Namespace for methods using feature information.

    Accessed through `gr.features`."""

    pr = None

    def __init__(self, pr):
        self.pr = pr

    def tss(self):
        """Return the transcription start sites.

        Returns the 5' for every interval with feature "transcript".

        See Also
        --------
        pyranges.genomicfeatures.GenomicFeaturesMethods.tes : return the transcription end sites

        Examples
        --------

        >>> gr = pr.data.ensembl_gtf()[["Source", "Feature"]]
        >>> gr
        +--------------+------------+--------------+-----------+-----------+--------------+
        | Chromosome   | Source     | Feature      | Start     | End       | Strand       |
        | (category)   | (object)   | (category)   | (int64)   | (int64)   | (category)   |
        |--------------+------------+--------------+-----------+-----------+--------------|
        | 1            | havana     | gene         | 11868     | 14409     | +            |
        | 1            | havana     | transcript   | 11868     | 14409     | +            |
        | 1            | havana     | exon         | 11868     | 12227     | +            |
        | 1            | havana     | exon         | 12612     | 12721     | +            |
        | ...          | ...        | ...          | ...       | ...       | ...          |
        | 1            | havana     | gene         | 1173055   | 1179555   | -            |
        | 1            | havana     | transcript   | 1173055   | 1179555   | -            |
        | 1            | havana     | exon         | 1179364   | 1179555   | -            |
        | 1            | havana     | exon         | 1173055   | 1176396   | -            |
        +--------------+------------+--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 2,446 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.features.tss()
        +--------------+------------+------------+-----------+-----------+--------------+
        | Chromosome   | Source     | Feature    | Start     | End       | Strand       |
        | (category)   | (object)   | (object)   | (int64)   | (int64)   | (category)   |
        |--------------+------------+------------+-----------+-----------+--------------|
        | 1            | havana     | tss        | 11868     | 11869     | +            |
        | 1            | havana     | tss        | 12009     | 12010     | +            |
        | 1            | havana     | tss        | 29553     | 29554     | +            |
        | 1            | havana     | tss        | 30266     | 30267     | +            |
        | ...          | ...        | ...        | ...       | ...       | ...          |
        | 1            | havana     | tss        | 1092812   | 1092813   | -            |
        | 1            | havana     | tss        | 1116086   | 1116087   | -            |
        | 1            | havana     | tss        | 1116088   | 1116089   | -            |
        | 1            | havana     | tss        | 1179554   | 1179555   | -            |
        +--------------+------------+------------+-----------+-----------+--------------+
        Stranded PyRanges object has 280 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        pr = self.pr

        if not pr.stranded:
            raise Exception(
                "Cannot compute TSSes or TESes without strand info. Perhaps use extend() or subsequence() or spliced_subsequence() instead?"
            )

        pr = pr[pr.Feature == "transcript"]
        pr = pr.apply(lambda df: _tss(df))

        pr.Feature = "tss"

        return pr

    def tes(self, slack=0):
        """Return the transcription end sites.

        Returns the 3' for every interval with feature "transcript".

        See Also
        --------
        pyranges.genomicfeatures.GenomicFeaturesMethods.tss : return the transcription start sites

        Examples
        --------

        >>> gr = pr.data.ensembl_gtf()[["Source", "Feature"]]
        >>> gr
        +--------------+------------+--------------+-----------+-----------+--------------+
        | Chromosome   | Source     | Feature      | Start     | End       | Strand       |
        | (category)   | (object)   | (category)   | (int64)   | (int64)   | (category)   |
        |--------------+------------+--------------+-----------+-----------+--------------|
        | 1            | havana     | gene         | 11868     | 14409     | +            |
        | 1            | havana     | transcript   | 11868     | 14409     | +            |
        | 1            | havana     | exon         | 11868     | 12227     | +            |
        | 1            | havana     | exon         | 12612     | 12721     | +            |
        | ...          | ...        | ...          | ...       | ...       | ...          |
        | 1            | havana     | gene         | 1173055   | 1179555   | -            |
        | 1            | havana     | transcript   | 1173055   | 1179555   | -            |
        | 1            | havana     | exon         | 1179364   | 1179555   | -            |
        | 1            | havana     | exon         | 1173055   | 1176396   | -            |
        +--------------+------------+--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 2,446 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.features.tes()
        +--------------+------------+------------+-----------+-----------+--------------+
        | Chromosome   | Source     | Feature    | Start     | End       | Strand       |
        | (category)   | (object)   | (object)   | (int64)   | (int64)   | (category)   |
        |--------------+------------+------------+-----------+-----------+--------------|
        | 1            | havana     | tes        | 14408     | 14409     | +            |
        | 1            | havana     | tes        | 13669     | 13670     | +            |
        | 1            | havana     | tes        | 31096     | 31097     | +            |
        | 1            | havana     | tes        | 31108     | 31109     | +            |
        | ...          | ...        | ...        | ...       | ...       | ...          |
        | 1            | havana     | tes        | 1090405   | 1090406   | -            |
        | 1            | havana     | tes        | 1091045   | 1091046   | -            |
        | 1            | havana     | tes        | 1091499   | 1091500   | -            |
        | 1            | havana     | tes        | 1173055   | 1173056   | -            |
        +--------------+------------+------------+-----------+-----------+--------------+
        Stranded PyRanges object has 280 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        pr = self.pr

        if not pr.stranded:
            raise Exception(
                "Cannot compute TSSes or TESes without strand info. Perhaps use extend() or subsequence() or spliced_subsequence() instead?"
            )

        pr = pr[pr.Feature == "transcript"]
        pr = pr.apply(lambda df: _tes(df))

        pr.Feature = "tes"

        return pr

    def introns(self, by="gene", nb_cpu=1):
        """Return the introns.

        Parameters
        ----------
        by : str, {"gene", "transcript"}, default "gene"
            Whether to find introns per gene or transcript.

        nb_cpu: int, default 1

            How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
            Will only lead to speedups on large datasets.

        See Also
        --------
        pyranges.genomicfeatures.GenomicFeaturesMethods.tss : return the transcription start sites

        Examples
        --------

        >>> gr = pr.data.ensembl_gtf()[["Feature", "gene_id", "transcript_id"]]
        >>> gr
        +--------------+--------------+-----------+-----------+--------------+-----------------+-----------------+
        | Chromosome   | Feature      | Start     | End       | Strand       | gene_id         | transcript_id   |
        | (category)   | (category)   | (int64)   | (int64)   | (category)   | (object)        | (object)        |
        |--------------+--------------+-----------+-----------+--------------+-----------------+-----------------|
        | 1            | gene         | 11868     | 14409     | +            | ENSG00000223972 | nan             |
        | 1            | transcript   | 11868     | 14409     | +            | ENSG00000223972 | ENST00000456328 |
        | 1            | exon         | 11868     | 12227     | +            | ENSG00000223972 | ENST00000456328 |
        | 1            | exon         | 12612     | 12721     | +            | ENSG00000223972 | ENST00000456328 |
        | ...          | ...          | ...       | ...       | ...          | ...             | ...             |
        | 1            | gene         | 1173055   | 1179555   | -            | ENSG00000205231 | nan             |
        | 1            | transcript   | 1173055   | 1179555   | -            | ENSG00000205231 | ENST00000379317 |
        | 1            | exon         | 1179364   | 1179555   | -            | ENSG00000205231 | ENST00000379317 |
        | 1            | exon         | 1173055   | 1176396   | -            | ENSG00000205231 | ENST00000379317 |
        +--------------+--------------+-----------+-----------+--------------+-----------------+-----------------+
        Stranded PyRanges object has 2,446 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.features.introns(by="gene")
        +--------------+------------+-----------+-----------+--------------+-----------------+-----------------+
        | Chromosome   | Feature    | Start     | End       | Strand       | gene_id         | transcript_id   |
        | (object)     | (object)   | (int64)   | (int64)   | (category)   | (object)        | (object)        |
        |--------------+------------+-----------+-----------+--------------+-----------------+-----------------|
        | 1            | intron     | 1173926   | 1174265   | +            | ENSG00000162571 | nan             |
        | 1            | intron     | 1174321   | 1174423   | +            | ENSG00000162571 | nan             |
        | 1            | intron     | 1174489   | 1174520   | +            | ENSG00000162571 | nan             |
        | 1            | intron     | 1175034   | 1179188   | +            | ENSG00000162571 | nan             |
        | ...          | ...        | ...       | ...       | ...          | ...             | ...             |
        | 1            | intron     | 874591    | 875046    | -            | ENSG00000283040 | nan             |
        | 1            | intron     | 875155    | 875525    | -            | ENSG00000283040 | nan             |
        | 1            | intron     | 875625    | 876526    | -            | ENSG00000283040 | nan             |
        | 1            | intron     | 876611    | 876754    | -            | ENSG00000283040 | nan             |
        +--------------+------------+-----------+-----------+--------------+-----------------+-----------------+
        Stranded PyRanges object has 311 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.features.introns(by="transcript")
        +--------------+------------+-----------+-----------+--------------+-----------------+-----------------+
        | Chromosome   | Feature    | Start     | End       | Strand       | gene_id         | transcript_id   |
        | (object)     | (object)   | (int64)   | (int64)   | (category)   | (object)        | (object)        |
        |--------------+------------+-----------+-----------+--------------+-----------------+-----------------|
        | 1            | intron     | 818202    | 818722    | +            | ENSG00000177757 | ENST00000326734 |
        | 1            | intron     | 960800    | 961292    | +            | ENSG00000187961 | ENST00000338591 |
        | 1            | intron     | 961552    | 961628    | +            | ENSG00000187961 | ENST00000338591 |
        | 1            | intron     | 961750    | 961825    | +            | ENSG00000187961 | ENST00000338591 |
        | ...          | ...        | ...       | ...       | ...          | ...             | ...             |
        | 1            | intron     | 732207    | 732980    | -            | ENSG00000230021 | ENST00000648019 |
        | 1            | intron     | 168165    | 169048    | -            | ENSG00000241860 | ENST00000655252 |
        | 1            | intron     | 165942    | 167958    | -            | ENSG00000241860 | ENST00000662089 |
        | 1            | intron     | 168165    | 169048    | -            | ENSG00000241860 | ENST00000662089 |
        +--------------+------------+-----------+-----------+--------------+-----------------+-----------------+
        Stranded PyRanges object has 1,043 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        kwargs = {"by": by, "nb_cpu": nb_cpu}
        kwargs = pr.pyranges_main.fill_kwargs(kwargs)

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


def _outside_bounds(df, **kwargs):
    df = df.copy()

    chromsizes = kwargs.get("chromsizes")

    if not isinstance(chromsizes, dict):
        size_df = chromsizes.df
        chromsizes = {k: v for k, v in zip(size_df.Chromosome, size_df.End)}

    size = int(chromsizes[df.Chromosome.iloc[0]])
    clip = kwargs.get("clip", False)
    only_right = kwargs.get("only_right", False)

    ends_outright = df.End > size
    if not only_right:
        starts_outleft = df.Start < 0

    if not clip:  # i.e. remove
        if only_right:
            df = df[~ends_outright]
        else:
            df = df[~ends_outright & ~starts_outleft]

    else:
        starts_outright = df.Start >= size

        if only_right:
            df.loc[ends_outright, "End"] = size

            # removing intervals completely out of bounds
            df = df[~starts_outright]

        else:
            ends_outleft = df.End <= 0

            df.loc[ends_outright, "End"] = size
            df.loc[starts_outleft, "Start"] = 0

            # removing intervals completely out of bounds:
            df = df[~starts_outright & ~ends_outleft]

    return df


def genome_bounds(gr, chromsizes, clip=False, only_right=False):
    """Remove or clip intervals outside of genome bounds.

    Parameters
    ----------

    gr : PyRanges

        Input intervals

    chromsizes : dict or PyRanges or pyfaidx.Fasta

        Dict or PyRanges describing the lengths of the chromosomes.
        pyfaidx.Fasta object is also accepted since it conveniently loads chromosome length

    clip : bool, default False

        Returns the portions of intervals within bounds,
        instead of dropping intervals entirely if they are even partially
        out of bounds

    only_right : bool, default False

        If True, remove or clip only intervals that are out-of-bounds on the right,
        and do not alter those out-of-bounds on the left (whose Start is < 0)


    Examples
    --------

    >>> d = {"Chromosome": [1, 1, 3], "Start": [1, 249250600, 5], "End": [2, 249250640, 7]}
    >>> gr = pr.from_dict(d)
    >>> gr
    +--------------+-----------+-----------+
    |   Chromosome |     Start |       End |
    |   (category) |   (int64) |   (int64) |
    |--------------+-----------+-----------|
    |            1 |         1 |         2 |
    |            1 | 249250600 | 249250640 |
    |            3 |         5 |         7 |
    +--------------+-----------+-----------+
    Unstranded PyRanges object has 3 rows and 3 columns from 2 chromosomes.
    For printing, the PyRanges was sorted on Chromosome.

    >>> chromsizes = {"1": 249250621, "3": 500}
    >>> chromsizes
    {'1': 249250621, '3': 500}

    >>> pr.gf.genome_bounds(gr, chromsizes)
    +--------------+-----------+-----------+
    |   Chromosome |     Start |       End |
    |   (category) |   (int64) |   (int64) |
    |--------------+-----------+-----------|
    |            1 |         1 |         2 |
    |            3 |         5 |         7 |
    +--------------+-----------+-----------+
    Unstranded PyRanges object has 2 rows and 3 columns from 2 chromosomes.
    For printing, the PyRanges was sorted on Chromosome.

    >>> pr.gf.genome_bounds(gr, chromsizes, clip=True)
    +--------------+-----------+-----------+
    |   Chromosome |     Start |       End |
    |   (category) |   (int64) |   (int64) |
    |--------------+-----------+-----------|
    |            1 |         1 |         2 |
    |            1 | 249250600 | 249250621 |
    |            3 |         5 |         7 |
    +--------------+-----------+-----------+
    Unstranded PyRanges object has 3 rows and 3 columns from 2 chromosomes.
    For printing, the PyRanges was sorted on Chromosome.

    >>> del chromsizes['3']
    >>> chromsizes
    {'1': 249250621}

    >>> pr.gf.genome_bounds(gr, chromsizes)
    Traceback (most recent call last):
    ...
    KeyError: '3'
    """

    if isinstance(chromsizes, pr.PyRanges):
        chromsizes = {k: v for k, v in zip(chromsizes.Chromosome, chromsizes.End)}

    elif isinstance(chromsizes, dict):
        pass

    else:
        try:
            import pyfaidx  # type: ignore

            if isinstance(chromsizes, pyfaidx.Fasta):
                chromsizes = {k: len(chromsizes[k]) for k in chromsizes.keys()}
        except ImportError:
            pass

    assert isinstance(
        chromsizes, dict
    ), "ERROR chromsizes must be a dictionary, or a PyRanges, or a pyfaidx.Fasta object"

    return gr.apply(_outside_bounds, chromsizes=chromsizes, clip=clip, only_right=only_right)


def _last_tile(df, **kwargs):
    # do not need copy, since it is only used internally by
    # tile_genome
    # df = df.copy()
    sizes = kwargs.get("sizes")
    size = sizes[df.Chromosome.iloc[0]].End.iloc[0]
    df.loc[df.tail(1).index, "End"] = size

    return df


def tile_genome(genome, tile_size, tile_last=False):
    """Create a tiled genome.

    Parameters
    ----------
    chromsizes : dict or PyRanges

        Dict or PyRanges describing the lengths of the chromosomes.

    tile_size : int
        Length of the tiles.

    tile_last : bool, default False

        Use genome length as end of last tile.

    See Also
    --------

    pyranges.PyRanges.tile : split intervals into adjacent non-overlapping tiles.

    Examples
    --------

    >>> chromsizes = pr.data.chromsizes()
    >>> chromsizes
    +--------------+-----------+-----------+
    | Chromosome   | Start     | End       |
    | (category)   | (int64)   | (int64)   |
    |--------------+-----------+-----------|
    | chr1         | 0         | 249250621 |
    | chr2         | 0         | 243199373 |
    | chr3         | 0         | 198022430 |
    | chr4         | 0         | 191154276 |
    | ...          | ...       | ...       |
    | chr22        | 0         | 51304566  |
    | chrM         | 0         | 16571     |
    | chrX         | 0         | 155270560 |
    | chrY         | 0         | 59373566  |
    +--------------+-----------+-----------+
    Unstranded PyRanges object has 25 rows and 3 columns from 25 chromosomes.
    For printing, the PyRanges was sorted on Chromosome.

    >>> pr.gf.tile_genome(chromsizes, int(1e6))
    +--------------+-----------+-----------+
    | Chromosome   | Start     | End       |
    | (category)   | (int64)   | (int64)   |
    |--------------+-----------+-----------|
    | chr1         | 0         | 1000000   |
    | chr1         | 1000000   | 2000000   |
    | chr1         | 2000000   | 3000000   |
    | chr1         | 3000000   | 4000000   |
    | ...          | ...       | ...       |
    | chrY         | 56000000  | 57000000  |
    | chrY         | 57000000  | 58000000  |
    | chrY         | 58000000  | 59000000  |
    | chrY         | 59000000  | 59373566  |
    +--------------+-----------+-----------+
    Unstranded PyRanges object has 3,114 rows and 3 columns from 25 chromosomes.
    For printing, the PyRanges was sorted on Chromosome.
    """

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
    tss_neg.loc[:, "Start"] = tss_neg.End - 1

    # pd.options.mode.chained_assignment = "warn"
    tss = pd.concat([tss_pos, tss_neg], sort=False)
    tss["End"] = tss.Start + 1
    tss.End = tss.End + slack
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
    tes_neg.loc[:, "End"] = tes_neg.Start + 1

    # pd.options.mode.chained_assignment = "warn"
    tes = pd.concat([tes_pos, tes_neg], sort=False)
    tes["Start"] = tes.End - 1
    tes.End = tes.End + slack
    tes.Start = tes.Start - slack
    tes.loc[tes.Start < 0, "Start"] = 0

    tes.index = range(len(tes))

    tes[["Start", "End"]] = tes[["Start", "End"]].astype(intype)

    return tes


by_to_id = {"gene": "gene_id", "transcript": "transcript_id"}


def _introns2(df, exons, **kwargs):
    """TODO: refactor"""

    if df.empty or exons.empty:
        return None

    original_order = df.columns
    by = kwargs["by"]
    id_column = by_to_id[by]

    exons = exons[["Start", "End", id_column]]
    genes = df[["Start", "End", id_column]]
    exons.columns = ["Start", "End", "by_id"]
    genes.columns = ["Start", "End", "by_id"]

    intersection = pd.Series(np.intersect1d(exons["by_id"], genes["by_id"]))
    if len(intersection) == 0:
        return None

    exons = exons[exons["by_id"].isin(intersection)].reset_index(drop=True).sort_values(["by_id", "Start"])
    genes = genes[genes["by_id"].isin(intersection)].reset_index(drop=True).sort_values(["by_id", "Start"])
    df = df[df[id_column].isin(intersection)].reset_index(drop=True)

    assert len(genes) == len(
        genes.drop_duplicates("by_id")
    ), "The {id_column}s need to be unique to compute the introns.".format(id_column=id_column)

    exon_ids = exons["by_id"].shift() != exons["by_id"]
    by_ids = pd.Series(range(1, len(genes) + 1))
    df.insert(0, "__temp__", by_ids)

    if len(exons) > 1 and exons["by_id"].iloc[0] == exons["by_id"].iloc[1]:
        exon_ids.iloc[0] = False
        exon_ids = exon_ids.cumsum() + 1
    else:
        exon_ids = exon_ids.cumsum()

    assert (by_ids == exon_ids.drop_duplicates().values).all()
    starts, ends, ids = find_introns(
        genes.Start.values,
        genes.End.values,
        by_ids.values,
        exons.Start.values,
        exons.End.values,
        exon_ids.values,
    )

    introns = pd.DataFrame(
        data={
            "Chromosome": df.Chromosome.iloc[0],
            "Start": starts,
            "End": ends,
            "by_id": ids,
        }
    )

    vc = introns["by_id"].value_counts(sort=False).to_frame().reset_index()
    vc.columns = ["by_id", "counts"]

    genes_without_introns = pd.DataFrame(data={"by_id": np.setdiff1d(by_ids.values, vc.by_id.values), "counts": 0})

    vc = pd.concat([vc, genes_without_introns]).sort_values("by_id")

    original_ids = np.repeat(vc.by_id, vc.counts).to_frame()
    original_ids = original_ids.merge(
        df[["__temp__", id_column]],
        right_on="__temp__",
        left_on="by_id",
        suffixes=("_drop", ""),
    )
    original_ids = original_ids.drop(
        ["__temp__"] + [c for c in original_ids.columns if c.endswith("_drop")], axis=1
    ).sort_values("by_id")
    introns.loc[:, "by_id"] = original_ids[id_column].values
    introns = introns.merge(df, left_on="by_id", right_on=id_column, suffixes=("", "_dropme"))
    introns = introns.drop([c for c in introns.columns if c.endswith("_dropme")], axis=1)

    if introns.Feature.dtype.name == "category" and "intron" not in introns.Feature.cat.categories:
        introns.Feature.cat.add_categories(["intron"])
    introns.loc[:, "Feature"] = "intron"

    introns = introns[original_order]

    return introns
