
from __future__ import print_function

import sys

import pandas as pd
import numpy as np

import pkg_resources

from pyranges.pyranges import PyRanges
from pyranges import data

from pyrle import PyRles, Rle

from pyranges.version import __version__


def read_bed(f, output_df=False):

    columns = "Chromosome Start End Name Score Strand ThickStart ThickEnd ItemRGB BlockCount BlockSizes BlockStarts".split()

    df = pd.read_csv(f, dtype={"Chromosome": "category", "Strand": "category"}, header=None, sep="\t")

    df.columns = columns[:df.shape[1]]

    if df.Strand.str.contains("\.").any():
        print("Bed file contains '.' strands and is considered unstranded.", file=sys.stderr)

    if not output_df:
        return PyRanges(df)
    else:
        return df


def read_bam(f, output_df=False):

    import pysam

    samfile = pysam.AlignmentFile(f, "rb")

    chromosomes, starts, ends, names, sequences, strands = [], [], [], [], [], []
    for read in samfile.fetch():
        strand = '+'
        if(read.is_reverse):
            strand = '-'

        chromosomes.append(samfile.getrname(read.tid))
        starts.append(read.pos)
        ends.append(read.aend)
        sequences.append(read.seq)
        strands.append(strand)
        names.append(read.qname)

    df = pd.DataFrame({"Chromosome": chromosomes, "Start": starts, "End": ends,
                       "Sequence": sequences, "Strand": strands, "Name": names})["Chromosome Start End Sequence Name Strand".split()]

    if not output_df:
        return PyRanges(df)
    else:
        return df


def _fetch_gene_transcript_exon_id(attribute, annotation):

    no_quotes = attribute.str.replace('"', '').str.replace("'", "")

    df = no_quotes.str.extract("gene_id.?(.+?);(?:.*transcript_id.?(.+?);)?(?:.*exon_number.?(.+?);)?(?:.*exon_id.?(.+?);)?", expand=True)# .iloc[:, [1, 2, 3]]

    df.columns = "GeneID TranscriptID ExonNumber ExonID".split()

    # df.loc[:, "ExonNumber"] = df.ExonNumber.astype(int)

    if annotation == "ensembl":
        newdf = []
        for c in "GeneID TranscriptID ExonID".split():
            r = df[c].astype(str).str.extract('(\d+)').astype(float)
            newdf.append(r)

        newdf = pd.concat(newdf, axis=1)
        newdf.insert(2, "ExonNumber", df["ExonNumber"])
        df = newdf

    return df


def read_gtf(f, annotation=None, output_df=False):

    """seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
    # source - name of the program that generated this feature, or the data source (database or project name)
    feature - feature type name, e.g. Gene, Variation, Similarity
    start - Start position of the feature, with sequence numbering starting at 1.
    end - End position of the feature, with sequence numbering starting at 1.
    score - A floating point value.
    strand - defined as + (forward) or - (reverse).
    # frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature."""
    dtypes = {"Chromosome": "category", "Feature": "category", "Strand": "category"}

    df_iter = pd.read_csv(f, sep="\t", comment="#", usecols=[0, 2, 3, 4, 5, 6, 8], header=None, names="Chromosome Feature Start End Score Strand Attribute".split(), dtype=dtypes, chunksize=int(1e4))

    dfs = []
    for df in df_iter:
        # Since Start is 1-indexed
        df.Start -= 1

        if sum(df.Score == ".") == len(df):
            cols_to_concat = "Chromosome Start End Strand Feature".split()
        else:
            cols_to_concat = "Chromosome Start End Strand Feature Score".split()

        extract = _fetch_gene_transcript_exon_id(df.Attribute, annotation)
        extract.columns = "GeneID TranscriptID ExonNumber ExonID".split()

        extract.ExonNumber = extract.ExonNumber.astype(float)

        df = pd.concat([df[cols_to_concat],
                            extract], axis=1)

        dfs.append(df)

    df = pd.concat(dfs)

    if not output_df:
        return PyRanges(df)
    else:
        return df

def filehandle(f):

    if f.lower().endswith(".gz"):
        import gzip
        fh = gzip.open(f)
    else:
        fh = open(f)

    return fh

def guess_delim(f):

    fh = filehandle(f)

    # in case there is a header
    lines = [fh.readline() for i in range(0, 11)]

    fh.close()

    import re
    from collections import Counter

    counts = []
    for l in lines[1:]:
        no_floats = re.sub(r"\d*\.\d+", "", l)
        stripped_line = re.sub(r"[a-z0-9]+", "", no_floats)
        counts.append(Counter(stripped_line))

    chars_in_all_lines = set.intersection(*[set(c) for c in counts])

    # find the chars that have the same number of counts in all lines
    same_number_in_all_lines = []

    for c in chars_in_all_lines:
        for i in range(10):
            if i == 0:
                first = counts[i][c]
            else:
                if first != counts[i][c]:
                    break

            if i == 9:
                same_number_in_all_lines.append(c)

    # find the most common chars
    most_common = Counter({c: counts[0][c] for c in same_number_in_all_lines}).most_common()

    # several chars are equally common
    equally_common = [c for c in most_common if c[1] == most_common[0][1]]

    if len(equally_common) > 1:
        _equally_common = [c[0] for c in equally_common]
        if " " in _equally_common and "\t" in _equally_common: 
            guess = "\s+"
        else:
            import csv
            guess = csv.Sniffer("".join(lines)).sniff().delimiter
    else:
        guess = most_common[0][0]

    return guess



def guess_header(f, delim):

    df = pd.read_csv(f, sep=delim, nrows=10, header=None)
    df2 = pd.read_csv(f, sep=delim, nrows=10, header=0)

    if all(df.dtypes.values == df2.dtypes.values):
        return False

    return True


def guess_strand(df, number_unique):

    strand_cols = []
    for k, v in number_unique.items():

        # strand must be cateory or object
        if str(df[k].dtype) not in ["category", "object"]:
            continue

        # strand col has at most 3 values
        if v <= 3:
            # strand col can only contain "+, -, or ."
            if df[k].str.contains("\+|-|\.").all():
                strand_cols.append(k)

    guess = strand_cols[0]

    # bedpe for example
    if len(strand_cols) > 1:
        all_equal = True
        first = df[strand_cols[0]]
        for strand_col in strand_cols[1:]:
            if not (first == df[strand_col]).all():
                all_equal = False

        if not all_equal:
            print("More than one possible strand column found:", ", ".join(strand_cols),
                  "arbitrarily choosing:", guess)

    return guess


def guess_start_end(df, number_unique):

    # starts and ends should be int

    # should be at least 90% of original length, not too many equal starts/ends

    # all starts should be before all ends
    pass


def guess_chromosome(df, number_unique):

    # should be at most 10% of length of df

    # if integer must contain at most 256 values

    # if string/cateory ...

    pass


def parse_file(f):

    sep = guess_delim(f)
    header = guess_header(f, sep)

    df = pd.read_csv(f, sep=sep,
                     header={True: 0, False: None}[header])

    # all necessary columns found
    if header and df.columns.isin(["Chromosome", "Start", "End"]).sum() == 3:
        return df


    # helps find out which should be categorical and also which could be chromosome, strand or other
    number_unique = {c: len(df[c].unique()) for c in df}

    # guessing

    strand_guess = guess_strand(df, number_unique)
    start_guess, end_guess = guess_start_end(df, number_unique)
    

    return df


def guess_chromosome_col(df):

    # if df.columns.dtype == object:
    #     header_candidates = df.colums.str.lower().str.contains("chr|seq|scaffold")

    # no other possible types
    possible_cols = df.select_dtypes(["category", "int", "object"])

    



def guess_columns(f):

    df = parse_file(f)

    chromosome_idx = guess_chromosome_col(df)
