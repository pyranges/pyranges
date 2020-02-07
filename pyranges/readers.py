from __future__ import print_function

import sys

import pandas as pd

from pyranges.pyranges import PyRanges


from pyranges.version import __version__


def read_bed(f, output_df=False, nrows=None):

    columns = "Chromosome Start End Name Score Strand ThickStart ThickEnd ItemRGB BlockCount BlockSizes BlockStarts".split(
    )

    if f.endswith(".gz"):
        import gzip
        first_start = gzip.open(f).readline().split()[1]
    else:
        first_start = open(f).readline().split()[1]

    header = None

    try:
        int(first_start)
    except ValueError:
        header = 0

    df = pd.read_csv(
        f,
        dtype={
            "Chromosome": "category",
            "Strand": "category"
        },
        nrows=nrows,
        header=header,
        sep="\t")

    df.columns = columns[:df.shape[1]]

    if not output_df:
        return PyRanges(df)
    else:
        return df


def read_bam(f, sparse=True, output_df=False, mapq=0, required_flag=0, filter_flag=1540):

    try:
        import bamread
    except ModuleNotFoundError as e:
        print("bamread must be installed to read bam. Use `conda install -c bioconda bamread` or `pip install bamread` to install it.")
        sys.exit(1)

    if sparse:
        df = bamread.read_bam(f, mapq, required_flag, filter_flag)
    else:
        try:
            df = bamread.read_bam_full(f, mapq, required_flag, filter_flag)
        except AttributeError:
            print("bamread version 0.0.6 or higher is required to read bam non-sparsely.")

    if output_df:
        return df
    else:
        return PyRanges(df)

    # return bamread.read_bam(f, mapq, required_flag, filter_flag)


def _fetch_gene_transcript_exon_id(attribute, annotation):

    no_quotes = attribute.str.replace('"', '').str.replace("'", "")

    df = no_quotes.str.extract(
        "gene_id.?(.+?);(?:.*transcript_id.?(.+?);)?(?:.*exon_number.?(.+?);)?(?:.*exon_id.?(.+?);)?",
        expand=True)  # .iloc[:, [1, 2, 3]]

    df.columns = "gene_id transcript_id exon_number exon_id".split()

    if annotation == "ensembl":
        newdf = []
        for c in "gene_id transcript_id exon_id".split():
            r = df[c].astype(str).str.extract('(\d+)').astype(float)
            newdf.append(r)

        newdf = pd.concat(newdf, axis=1)
        newdf.insert(2, "exon_number", df["exon_number"])
        df = newdf

    return df


def skiprows(f):

    try:
        import gzip
        fh = gzip.open(f)
        for i, l in enumerate(fh):
            if l.decode()[0] != "#":
                break
    except (OSError, TypeError):  # not a gzipped file, or StringIO
        fh = open(f)
        for i, l in enumerate(fh):
            if l[0] != "#":
                break

    fh.close()

    return i


def read_gtf(f, full=True, annotation=None, output_df=False, nrows=None, duplicate_attr=False):

    _skiprows = skiprows(f)

    if full:
        gr = read_gtf_full(f, annotation, output_df, nrows, _skiprows,
                           duplicate_attr=duplicate_attr)
    else:
        gr = read_gtf_restricted(f, annotation, output_df, nrows, _skiprows)

    return gr


def read_gtf_full(f, annotation=None, output_df=False, nrows=None, skiprows=0, duplicate_attr=False):
    """seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
    # source - name of the program that generated this feature, or the data source (database or project name)
    feature - feature type name, e.g. Gene, Variation, Similarity
    start - Start position of the feature, with sequence numbering starting at 1.
    end - End position of the feature, with sequence numbering starting at 1.
    score - A floating point value.
    strand - defined as + (forward) or - (reverse).
    # frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature."""
    dtypes = {
        "Chromosome": "category",
        "Feature": "category",
        "Strand": "category"
    }

    names = "Chromosome Source Feature Start End Score Strand Frame Attribute".split(
    )
    # names = "Chromosome Start End Score Strand Source Feature Frame Attribute".split()
    df_iter = pd.read_csv(
        f,
        sep="\t",
        header=None,
        names=names,
        dtype=dtypes,
        chunksize=int(1e5),
        skiprows=skiprows,
        nrows=nrows)

    _to_rows = to_rows_keep_duplicates if duplicate_attr else to_rows

    dfs = []
    for df in df_iter:
        extra = _to_rows(df.Attribute)
        df = df.drop("Attribute", axis=1)
        ndf = pd.concat([df, extra], axis=1, sort=False)
        dfs.append(ndf)

    df = pd.concat(dfs, sort=False)
    df.loc[:, "Start"] = df.Start - 1

    if not output_df:
        return PyRanges(df)
    else:
        return df


def to_rows(anno):
    rowdicts = []
    for l in anno:
        l = l.replace('"', '').replace(";", "").split()
        rowdicts.append({k: v for k, v in zip(*([iter(l)] * 2))})

    return pd.DataFrame.from_dict(rowdicts).set_index(anno.index)


def to_rows_keep_duplicates(anno):
    rowdicts = []
    for l in anno:
        rowdict = {}
        l = l.replace('"', '').replace(";", "").split()
        for k, v in zip(*([iter(l)] * 2)):
            if k not in rowdict:
                rowdict[k] = v
            elif k in rowdict and isinstance(rowdict[k], list):
                rowdict[k].append(v)
            else:
                rowdict[k] = [rowdict[k], v]

        rowdicts.append({
            k: ','.join(v) if isinstance(v, list) else v
            for k, v in rowdict.items()
        })

    return pd.DataFrame.from_dict(rowdicts).set_index(anno.index)


def read_gtf_restricted(f,
                        annotation=None,
                        output_df=False,
                        skiprows=0,
                        nrows=None):
    """seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
    # source - name of the program that generated this feature, or the data source (database or project name)
    feature - feature type name, e.g. Gene, Variation, Similarity
    start - Start position of the feature, with sequence numbering starting at 1.
    end - End position of the feature, with sequence numbering starting at 1.
    score - A floating point value.
    strand - defined as + (forward) or - (reverse).
    # frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature."""
    dtypes = {
        "Chromosome": "category",
        "Feature": "category",
        "Strand": "category"
    }

    df_iter = pd.read_csv(
        f,
        sep="\t",
        comment="#",
        usecols=[0, 2, 3, 4, 5, 6, 8],
        header=None,
        names="Chromosome Feature Start End Score Strand Attribute".split(),
        dtype=dtypes,
        chunksize=int(1e5),
        nrows=nrows)

    dfs = []
    for df in df_iter:
        # Since Start is 1-indexed
        df.Start -= 1

        if sum(df.Score == ".") == len(df):
            cols_to_concat = "Chromosome Start End Strand Feature".split()
        else:
            cols_to_concat = "Chromosome Start End Strand Feature Score".split(
            )

        extract = _fetch_gene_transcript_exon_id(df.Attribute, annotation)
        extract.columns = "gene_id transcript_id exon_number exon_id".split()

        extract.exon_number = extract.exon_number.astype(float)

        df = pd.concat([df[cols_to_concat], extract], axis=1, sort=False)

        dfs.append(df)

    df = pd.concat(dfs, sort=False)

    df.loc[:, "Start"] = df.Start - 1

    if not output_df:
        return PyRanges(df)
    else:
        return df


def to_rows_gff3(anno):
    rowdicts = []

    # anno = anno.str.replace(";$", "")

    for l in list(anno):
        l = (it.split("=") for it in l.split(";"))
        rowdicts.append({k: v for k, v in l})

    return pd.DataFrame.from_dict(rowdicts).set_index(anno.index)


def read_gff3(f, annotation=None, output_df=False, nrows=None, skiprows=0):
    """seqid - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seq ID must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
source - name of the program that generated this feature, or the data source (database or project name)
type - type of feature. Must be a term or accession from the SOFA sequence ontology
start - Start position of the feature, with sequence numbering starting at 1.
end - End position of the feature, with sequence numbering starting at 1.
score - A floating point value.
strand - defined as + (forward) or - (reverse).
phase - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
attributes - A semicolon-separated list of tag-value pairs, providing additional information about each feature. Some of these tags are predefined, e.g. ID, Name, Alias, Parent - see the GFF documentation for more details."""

    dtypes = {
        "Chromosome": "category",
        "Feature": "category",
        "Strand": "category"
    }

    names = "Chromosome Source Feature Start End Score Strand Frame Attribute".split(
    )

    df_iter = pd.read_csv(
        f,
        comment="#",
        sep="\t",
        header=None,
        names=names,
        dtype=dtypes,
        chunksize=int(1e5),
        skiprows=skiprows,
        nrows=nrows)

    dfs = []
    for df in df_iter:
        extra = to_rows_gff3(df.Attribute.astype(str))
        df = df.drop("Attribute", axis=1)
        ndf = pd.concat([df, extra], axis=1, sort=False)
        dfs.append(ndf)

    df = pd.concat(dfs, sort=False)

    df.loc[:, "Start"] = df.Start - 1

    if not output_df:
        return PyRanges(df)
    else:
        return df
