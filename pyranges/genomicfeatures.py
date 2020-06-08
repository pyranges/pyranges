import pandas as pd

import pyranges as pr
import numpy as np

from pyranges.multithreaded import pyrange_apply
from sorted_nearest.src.introns import find_introns

__all__ = ["genome_bounds", "tile_genome", "GenomicFeaturesMethods"]

class GenomicFeaturesMethods():

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

        >>> gr = pr.data.ensembl_gtf()
        >>> gr
        +--------------+------------+--------------+-----------+-----------+------------+--------------+------------+------------------------------------+-------+
        | Chromosome   | Source     | Feature      | Start     | End       | Score      | Strand       | Frame      | gene_biotype                       | +19   |
        | (category)   | (object)   | (category)   | (int32)   | (int32)   | (object)   | (category)   | (object)   | (object)                           | ...   |
        |--------------+------------+--------------+-----------+-----------+------------+--------------+------------+------------------------------------+-------|
        | 1            | havana     | gene         | 11868     | 14409     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
        | 1            | havana     | transcript   | 11868     | 14409     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
        | 1            | havana     | exon         | 11868     | 12227     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
        | 1            | havana     | exon         | 12612     | 12721     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
        | ...          | ...        | ...          | ...       | ...       | ...        | ...          | ...        | ...                                | ...   |
        | 1            | havana     | gene         | 1173055   | 1179555   | .          | -            | .          | lncRNA                             | ...   |
        | 1            | havana     | transcript   | 1173055   | 1179555   | .          | -            | .          | lncRNA                             | ...   |
        | 1            | havana     | exon         | 1179364   | 1179555   | .          | -            | .          | lncRNA                             | ...   |
        | 1            | havana     | exon         | 1173055   | 1176396   | .          | -            | .          | lncRNA                             | ...   |
        +--------------+------------+--------------+-----------+-----------+------------+--------------+------------+------------------------------------+-------+
        Stranded PyRanges object has 2,446 rows and 28 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        19 hidden columns: gene_id, gene_name, gene_source, gene_version, tag, transcript_biotype, transcript_id, transcript_name, transcript_source, transcript_support_level, ... (+ 9 more.)

        >>> gr.features.tss()
        +--------------+------------+------------+-----------+-----------+------------+--------------+------------+------------------------------------+-------+
        | Chromosome   | Source     | Feature    | Start     | End       | Score      | Strand       | Frame      | gene_biotype                       | +19   |
        | (category)   | (object)   | (object)   | (int32)   | (int32)   | (object)   | (category)   | (object)   | (object)                           | ...   |
        |--------------+------------+------------+-----------+-----------+------------+--------------+------------+------------------------------------+-------|
        | 1            | havana     | tss        | 11868     | 11869     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
        | 1            | havana     | tss        | 12009     | 12010     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
        | 1            | havana     | tss        | 29553     | 29554     | .          | +            | .          | lncRNA                             | ...   |
        | 1            | havana     | tss        | 30266     | 30267     | .          | +            | .          | lncRNA                             | ...   |
        | ...          | ...        | ...        | ...       | ...       | ...        | ...          | ...        | ...                                | ...   |
        | 1            | havana     | tss        | 1092813   | 1092814   | .          | -            | .          | protein_coding                     | ...   |
        | 1            | havana     | tss        | 1116087   | 1116088   | .          | -            | .          | protein_coding                     | ...   |
        | 1            | havana     | tss        | 1116089   | 1116090   | .          | -            | .          | protein_coding                     | ...   |
        | 1            | havana     | tss        | 1179555   | 1179556   | .          | -            | .          | lncRNA                             | ...   |
        +--------------+------------+------------+-----------+-----------+------------+--------------+------------+------------------------------------+-------+
        Stranded PyRanges object has 280 rows and 28 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        19 hidden columns: gene_id, gene_name, gene_source, gene_version, tag, transcript_biotype, transcript_id, transcript_name, transcript_source, transcript_support_level, ... (+ 9 more.)
        """

        pr = self.pr

        if not pr.stranded:
            raise Exception(
                "Cannot compute TSSes or TESes without strand info. Perhaps use slack() instead?"
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

        >>> gr = pr.data.ensembl_gtf()
        >>> gr
        +--------------+------------+--------------+-----------+-----------+------------+--------------+------------+------------------------------------+-------+
        | Chromosome   | Source     | Feature      | Start     | End       | Score      | Strand       | Frame      | gene_biotype                       | +19   |
        | (category)   | (object)   | (category)   | (int32)   | (int32)   | (object)   | (category)   | (object)   | (object)                           | ...   |
        |--------------+------------+--------------+-----------+-----------+------------+--------------+------------+------------------------------------+-------|
        | 1            | havana     | gene         | 11868     | 14409     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
        | 1            | havana     | transcript   | 11868     | 14409     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
        | 1            | havana     | exon         | 11868     | 12227     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
        | 1            | havana     | exon         | 12612     | 12721     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
        | ...          | ...        | ...          | ...       | ...       | ...        | ...          | ...        | ...                                | ...   |
        | 1            | havana     | gene         | 1173055   | 1179555   | .          | -            | .          | lncRNA                             | ...   |
        | 1            | havana     | transcript   | 1173055   | 1179555   | .          | -            | .          | lncRNA                             | ...   |
        | 1            | havana     | exon         | 1179364   | 1179555   | .          | -            | .          | lncRNA                             | ...   |
        | 1            | havana     | exon         | 1173055   | 1176396   | .          | -            | .          | lncRNA                             | ...   |
        +--------------+------------+--------------+-----------+-----------+------------+--------------+------------+------------------------------------+-------+
        Stranded PyRanges object has 2,446 rows and 28 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        19 hidden columns: gene_id, gene_name, gene_source, gene_version, tag, transcript_biotype, transcript_id, transcript_name, transcript_source, transcript_support_level, ... (+ 9 more.)

        >>> gr.features.tes()
        +--------------+------------+------------+-----------+-----------+------------+--------------+------------+------------------------------------+-------+
        | Chromosome   | Source     | Feature    | Start     | End       | Score      | Strand       | Frame      | gene_biotype                       | +19   |
        | (category)   | (object)   | (object)   | (int32)   | (int32)   | (object)   | (category)   | (object)   | (object)                           | ...   |
        |--------------+------------+------------+-----------+-----------+------------+--------------+------------+------------------------------------+-------|
        | 1            | havana     | tes        | 14409     | 14410     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
        | 1            | havana     | tes        | 13670     | 13671     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
        | 1            | havana     | tes        | 31097     | 31098     | .          | +            | .          | lncRNA                             | ...   |
        | 1            | havana     | tes        | 31109     | 31110     | .          | +            | .          | lncRNA                             | ...   |
        | ...          | ...        | ...        | ...       | ...       | ...        | ...          | ...        | ...                                | ...   |
        | 1            | havana     | tes        | 1092813   | 1092814   | .          | -            | .          | protein_coding                     | ...   |
        | 1            | havana     | tes        | 1116087   | 1116088   | .          | -            | .          | protein_coding                     | ...   |
        | 1            | havana     | tes        | 1116089   | 1116090   | .          | -            | .          | protein_coding                     | ...   |
        | 1            | havana     | tes        | 1179555   | 1179556   | .          | -            | .          | lncRNA                             | ...   |
        +--------------+------------+------------+-----------+-----------+------------+--------------+------------+------------------------------------+-------+
        Stranded PyRanges object has 280 rows and 28 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        19 hidden columns: gene_id, gene_name, gene_source, gene_version, tag, transcript_biotype, transcript_id, transcript_name, transcript_source, transcript_support_level, ... (+ 9 more.)
        """

        pr = self.pr

        if not pr.stranded:
            raise Exception(
                "Cannot compute TSSes or TESes without strand info. Perhaps use slack() instead?"
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

        >>> gr = pr.data.ensembl_gtf()
        >>> gr
        +--------------+------------+--------------+-----------+-----------+------------+--------------+------------+------------------------------------+-------+
        | Chromosome   | Source     | Feature      | Start     | End       | Score      | Strand       | Frame      | gene_biotype                       | +19   |
        | (category)   | (object)   | (category)   | (int32)   | (int32)   | (object)   | (category)   | (object)   | (object)                           | ...   |
        |--------------+------------+--------------+-----------+-----------+------------+--------------+------------+------------------------------------+-------|
        | 1            | havana     | gene         | 11868     | 14409     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
        | 1            | havana     | transcript   | 11868     | 14409     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
        | 1            | havana     | exon         | 11868     | 12227     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
        | 1            | havana     | exon         | 12612     | 12721     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
        | ...          | ...        | ...          | ...       | ...       | ...        | ...          | ...        | ...                                | ...   |
        | 1            | havana     | gene         | 1173055   | 1179555   | .          | -            | .          | lncRNA                             | ...   |
        | 1            | havana     | transcript   | 1173055   | 1179555   | .          | -            | .          | lncRNA                             | ...   |
        | 1            | havana     | exon         | 1179364   | 1179555   | .          | -            | .          | lncRNA                             | ...   |
        | 1            | havana     | exon         | 1173055   | 1176396   | .          | -            | .          | lncRNA                             | ...   |
        +--------------+------------+--------------+-----------+-----------+------------+--------------+------------+------------------------------------+-------+
        Stranded PyRanges object has 2,446 rows and 28 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        19 hidden columns: gene_id, gene_name, gene_source, gene_version, tag, transcript_biotype, transcript_id, transcript_name, transcript_source, transcript_support_level, ... (+ 9 more.)

        >>> gr.features.introns(by="gene")
        +--------------+----------------+------------+-----------+-----------+------------+--------------+------------+-------+
        | Chromosome   | Source         | Feature    | Start     | End       | Score      | Strand       | Frame      | +20   |
        | (object)     | (object)       | (object)   | (int32)   | (int32)   | (object)   | (category)   | (object)   | ...   |
        |--------------+----------------+------------+-----------+-----------+------------+--------------+------------+-------|
        | 1            | ensembl_havana | intron     | 1173926   | 1174265   | .          | +            | .          | ...   |
        | 1            | ensembl_havana | intron     | 1174321   | 1174423   | .          | +            | .          | ...   |
        | 1            | ensembl_havana | intron     | 1174489   | 1174520   | .          | +            | .          | ...   |
        | 1            | ensembl_havana | intron     | 1175034   | 1179188   | .          | +            | .          | ...   |
        | ...          | ...            | ...        | ...       | ...       | ...        | ...          | ...        | ...   |
        | 1            | havana         | intron     | 874591    | 875046    | .          | -            | .          | ...   |
        | 1            | havana         | intron     | 875155    | 875525    | .          | -            | .          | ...   |
        | 1            | havana         | intron     | 875625    | 876526    | .          | -            | .          | ...   |
        | 1            | havana         | intron     | 876611    | 876754    | .          | -            | .          | ...   |
        +--------------+----------------+------------+-----------+-----------+------------+--------------+------------+-------+
        Stranded PyRanges object has 311 rows and 28 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        20 hidden columns: gene_biotype, gene_id, gene_name, gene_source, gene_version, tag, transcript_biotype, transcript_id, transcript_name, ... (+ 11 more.)

        >>> gr.features.introns(by="transcript")
        +--------------+----------------+------------+-----------+-----------+------------+--------------+------------+----------------------------------+-------+
        | Chromosome   | Source         | Feature    | Start     | End       | Score      | Strand       | Frame      | gene_biotype                     | +19   |
        | (object)     | (object)       | (object)   | (int32)   | (int32)   | (object)   | (category)   | (object)   | (object)                         | ...   |
        |--------------+----------------+------------+-----------+-----------+------------+--------------+------------+----------------------------------+-------|
        | 1            | havana         | intron     | 818202    | 818722    | .          | +            | .          | lncRNA                           | ...   |
        | 1            | ensembl_havana | intron     | 960800    | 961292    | .          | +            | .          | protein_coding                   | ...   |
        | 1            | ensembl_havana | intron     | 961552    | 961628    | .          | +            | .          | protein_coding                   | ...   |
        | 1            | ensembl_havana | intron     | 961750    | 961825    | .          | +            | .          | protein_coding                   | ...   |
        | ...          | ...            | ...        | ...       | ...       | ...        | ...          | ...        | ...                              | ...   |
        | 1            | havana         | intron     | 732207    | 732980    | .          | -            | .          | transcribed_processed_pseudogene | ...   |
        | 1            | havana_tagene  | intron     | 168165    | 169048    | .          | -            | .          | lncRNA                           | ...   |
        | 1            | havana_tagene  | intron     | 165942    | 167958    | .          | -            | .          | lncRNA                           | ...   |
        | 1            | havana_tagene  | intron     | 168165    | 169048    | .          | -            | .          | lncRNA                           | ...   |
        +--------------+----------------+------------+-----------+-----------+------------+--------------+------------+----------------------------------+-------+
        Stranded PyRanges object has 1,043 rows and 28 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        19 hidden columns: gene_id, gene_name, gene_source, gene_version, tag, transcript_biotype, transcript_id, transcript_name, transcript_source, transcript_support_level, ... (+ 9 more.)
        """

        kwargs = {"by": by, "nb_cpu": nb_cpu}
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

def _outside_bounds(df, **kwargs):

    chromsizes = kwargs.get("chromsizes")

    if not isinstance(chromsizes, dict):
        size_df = chromsizes.df
        chromsizes = {k: v for k, v in zip(size_df.Chromosome, size_df.End)}

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

    """Remove or clip intervals outside of genome bounds.

    Parameters
    ----------
    chromsizes : dict or PyRanges

        Dict or PyRanges describing the lengths of the chromosomes.

    clip : bool, default False

        Part of interval within bounds.

    Examples
    --------

    >>> d = {"Chromosome": [1, 1, 3], "Start": [1, 249250600, 5], "End": [2, 249250640, 7]}
    >>> gr = pr.from_dict(d)
    >>> gr
    +--------------+-----------+-----------+
    |   Chromosome |     Start |       End |
    |   (category) |   (int32) |   (int32) |
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
    |   (category) |   (int32) |   (int32) |
    |--------------+-----------+-----------|
    |            1 |         1 |         2 |
    |            3 |         5 |         7 |
    +--------------+-----------+-----------+
    Unstranded PyRanges object has 2 rows and 3 columns from 2 chromosomes.
    For printing, the PyRanges was sorted on Chromosome.

    >>> pr.gf.genome_bounds(gr, chromsizes, clip=True)
    +--------------+-----------+-----------+
    |   Chromosome |     Start |       End |
    |   (category) |   (int32) |   (int32) |
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

    if not isinstance(chromsizes, dict):
        chromsizes = {
            k: v
            for k, v in zip(chromsizes.Chromosome, chromsizes.End)
        }

    return gr.apply(_outside_bounds, chromsizes=chromsizes, clip=clip)




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
    | (category)   | (int32)   | (int32)   |
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
    | (category)   | (int32)   | (int32)   |
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


def _introns2(df, exons, **kwargs):

    """TODO: refactor

    """

    if df.empty or exons.empty:
        return None

    original_order = df.columns
    by = kwargs["by"]
    id_column = by_to_id[by]

    exons = exons[["Start", "End", id_column]]
    genes = df[["Start", "End", id_column]]
    exons.columns = ["Start", "End", "by_id"]
    genes.columns = ["Start", "End", "by_id"]

    intersection = pd.Series(
        np.intersect1d(exons["by_id"], genes["by_id"]))
    if len(intersection) == 0:
        return None

    exons = exons[exons["by_id"].isin(intersection)].reset_index(
        drop=True).sort_values(["by_id", "Start"])
    genes = genes[genes["by_id"].isin(intersection)].reset_index(
        drop=True).sort_values(["by_id", "Start"])
    df = df[df[id_column].isin(intersection)].reset_index(
        drop=True)

    assert len(genes) == len(
        genes.drop_duplicates("by_id")
    ), "The {id_column}s need to be unique to compute the introns.".format(
        id_column=id_column)

    exon_ids = (exons["by_id"].shift() != exons["by_id"])
    by_ids = pd.Series(range(1, len(genes) + 1))
    df.insert(0, "__temp__", by_ids)

    if len(exons) > 1 and exons["by_id"].iloc[0] == exons["by_id"].iloc[1]:
        exon_ids.iloc[0] = False
        exon_ids = exon_ids.cumsum() + 1
    else:
        exon_ids = exon_ids.cumsum()

    assert (by_ids == exon_ids.drop_duplicates().values).all()
    starts, ends, ids = find_introns(genes.Start.values, genes.End.values,
                                     by_ids.values, exons.Start.values,
                                     exons.End.values, exon_ids.values)

    introns = pd.DataFrame(
        data={
            "Chromosome": df.Chromosome.iloc[0],
            "Start": starts,
            "End": ends,
            "by_id": ids
        })

    vc = introns["by_id"].value_counts(sort=False).to_frame().reset_index()
    vc.columns = ["by_id", "counts"]

    genes_without_introns = pd.DataFrame(
        data={
            "by_id": np.setdiff1d(by_ids.values, vc.by_id.values),
            "counts": 0
        })

    vc = pd.concat([vc, genes_without_introns]).sort_values("by_id")

    original_ids = np.repeat(vc.by_id, vc.counts).to_frame()
    original_ids = original_ids.merge(
        df[["__temp__", id_column]],
        right_on="__temp__",
        left_on="by_id",
        suffixes=("_drop", ""))
    original_ids = original_ids.drop(
        ["__temp__"] +
        [c for c in original_ids.columns if c.endswith("_drop")],
        axis=1).sort_values("by_id")
    introns.loc[:, "by_id"] = original_ids[id_column].values
    introns = introns.merge(
        df, left_on="by_id", right_on=id_column, suffixes=("", "_dropme"))
    introns = introns.drop(
        [c for c in introns.columns if c.endswith("_dropme")], axis=1)

    if introns.Feature.dtype.name == "category" and not "intron" in introns.Feature.cat.categories:
        introns.Feature.cat.add_categories(["intron"], inplace=True)
    introns.loc[:, "Feature"] = "intron"

    introns = introns[original_order]

    return introns
