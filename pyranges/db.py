

import pandas as pd

import MySQLdb



def ucsc():

    conn = MySQLdb.Connection(host="genome-mysql.cse.ucsc.edu", user="genome", db="hg38")

    query = 'select chrom, txStart, txEnd, exonStarts, exonEnds, name, name2, strand from refGene'

    df = pd.read_sql(query, conn)

    df.columns = "Chromosome Start End XS XE TranscriptID GeneID Strand".split()

    cols = df.columns.difference(['XS', 'XE'])
    exon_starts = df['XS'].str.decode("utf-8").str.replace(",$", "").str.split(',')
    exon_ends = df['XE'].str.decode("utf-8").str.replace(",$", "").str.split(',')

    ls = exon_starts.str.len()
    from itertools import chain
    xs = list(chain.from_iterable(exon_starts))
    xe = list(chain.from_iterable(exon_ends))
    (df.loc[df.index.repeat(ls), cols]
         .assign(ExonStart=xs, ExonEnd=xe))

# pd.DataFrame({
#     'Name' : df['Name'].values.repeat(genres.str.len()),
#     'genres' : list(chain.from_iterable(genres.tolist()))
# })

    conn.close()

    return df


