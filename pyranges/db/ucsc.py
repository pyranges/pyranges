
import pandas as pd
import numpy as np

import pyranges as pr

import MySQLdb

from itertools import chain


def _exons(df):

    _df = df.drop(["Start", "End"], axis=1)

    cols = _df.columns.difference(['XS', 'XE'])
    exon_starts = _df['XS'].str.decode("utf-8").str.replace(",$", "").str.split(',')
    exon_ends = _df['XE'].str.decode("utf-8").str.replace(",$", "").str.split(',')

    exon_numbers = exon_starts.apply(lambda x: [i for i in range(1, len(x) + 1)])

    ls = exon_starts.str.len()
    exon_starts = [int(i) for i in chain.from_iterable(exon_starts.to_list())]
    exon_ends = [int(i) for i in chain.from_iterable(exon_ends.to_list())]
    exon_starts = pd.Series(exon_starts, dtype=np.int32)
    exon_ends = pd.Series(exon_ends, dtype=np.int32)
    exon_numbers = pd.Series(chain.from_iterable(exon_numbers.to_list()), dtype="category")
    _exons = (_df.loc[_df.index.repeat(ls), cols]
             .assign(Start=exon_starts, End=exon_ends, Feature="exon")
             .reset_index(drop=True)
             .assign(ExonNumber=exon_numbers))

    return _exons


def ucsc(genome, query):

    host = "genome-mysql.cse.ucsc.edu"
    user = "genome"
    db = genome

    conn = MySQLdb.Connection(host=host, user=user, db=db)

    df = pd.read_sql(query, conn)

    conn.close()

    return df


def genes(genome):

    query = 'select chrom, txStart, txEnd, exonStarts, exonEnds, name, name2, strand from refGene'

    df = ucsc(genome, query)

    df.columns = "Chromosome Start End XS XE TranscriptID GeneID Strand".split()

    df = df.astype({"Chromosome": "category", "TranscriptID": "category", "GeneID": "category", "Strand": "category", "Start": int, "End": int})

    exons = _exons(df)
    _df = (df.drop("XS XE".split(), axis=1)
           .assign(Feature="TranscriptID"))

    df = pd.concat([_df, exons], sort=False).sort_values("Chromosome Start End".split())["Chromosome Start End Strand Feature TranscriptID ExonNumber".split()]

    return pr.PyRanges(df)


def chromosome_lengths(genome):

    query = 'select chrom,size from chromInfo'

    df = ucsc(genome, query)

    return pd.Series(df.size, index=df.chrom)
