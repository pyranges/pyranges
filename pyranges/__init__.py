from __future__ import print_function

import pandas as pd

import pkg_resources

from pyranges.pyranges import PyRanges


def load_dataset(basename):

    full_path = pkg_resources.resource_filename("pyranges", "example_data/{}.bed".format(basename))

    df = pd.read_table(full_path, header=None,
                       names="Chromosome Start End Name Score Strand".split())

    return PyRanges(df)


def list_datasets():

    datasets = [f.replace(".bed", "") for f in pkg_resources.resource_listdir("pyranges", "example_data")]

    print(datasets)


def read_bed(f):

    columns = "Chromosome Start End Name Score Strand ThickStart ThickEnd ItemRGB BlockCount BlockSizes BlockStarts".split()

    df = pd.read_table(f)

    df.columns = columns[:df.shape[1]]

    if df.Strand.str.contains("\.").any():
        print("Bed file contains '.' strands and is considered unstranded.", file=sys.stderr)

    return PyRanges(df)


def read_bam(f):

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

    return PyRanges(df)


def _fetch_gene_transcript_exon_id(attribute):

    no_quotes = attribute.str.replace('"', '').str.replace("'", "")

    return no_quotes.str.extract("gene_id.?(.+?);(?:.*transcript_id.?(.+?);)?(?:.*exon_number.?(.+?);)?(?:.*exon_id.?(.+?);)?", expand=True)


def read_gtf(f):

    """seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
    source - name of the program that generated this feature, or the data source (database or project name)
    feature - feature type name, e.g. Gene, Variation, Similarity
    start - Start position of the feature, with sequence numbering starting at 1.
    end - End position of the feature, with sequence numbering starting at 1.
    score - A floating point value.
    strand - defined as + (forward) or - (reverse).
    frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature."""

    df = pd.read_table(f, sep="\t", comment="#", header=None, names="Chromosome Source Feature Start End Score Strand Frame Attribute".split())

    extract = _fetch_gene_transcript_exon_id(df.Attribute)
    extract.columns = "GeneID TranscriptID ExonNumber ExonID".split()

    extract.ExonNumber = extract.ExonNumber.astype(float)

    df = pd.concat([df["Chromosome Start End Strand Feature Source Score Frame".split()],
                        extract], axis=1)

    return PyRanges(df)
