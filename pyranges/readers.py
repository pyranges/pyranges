from __future__ import print_function

import sys

import pandas as pd

from pyranges.pyranges import PyRanges

from io import StringIO

import pyranges as pr
from pyranges.version import __version__


def read_bed(f, as_df=False, nrows=None):

    """Return bed file as PyRanges.

    This is a reader for files that follow the bed format. They can have from
    3-12 columns which will be named like so:

    Chromosome Start End Name Score Strand ThickStart ThickEnd ItemRGB
    BlockCount BlockSizes BlockStarts

    Parameters
    ----------
    f : str

        Path to bed file

    as_df : bool, default False

        Whether to return as pandas DataFrame instead of PyRanges.

    nrows : int, default None

        Number of rows to return.

    Notes
    -----

    If you just want to create a PyRanges from a tab-delimited bed-like file,
    use `pr.PyRanges(pandas.read_table(f))` instead.

    Examples
    --------

    >>> path = pr.get_example_path("aorta.bed")
    >>> pr.read_bed(path, nrows=5)
    +--------------+-----------+-----------+------------+-----------+--------------+
    | Chromosome   |     Start |       End | Name       |     Score | Strand       |
    | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
    |--------------+-----------+-----------+------------+-----------+--------------|
    | chr1         |      9939 |     10138 | H3K27me3   |         7 | +            |
    | chr1         |      9953 |     10152 | H3K27me3   |         5 | +            |
    | chr1         |      9916 |     10115 | H3K27me3   |         5 | -            |
    | chr1         |      9951 |     10150 | H3K27me3   |         8 | -            |
    | chr1         |      9978 |     10177 | H3K27me3   |         7 | -            |
    +--------------+-----------+-----------+------------+-----------+--------------+
    Stranded PyRanges object has 5 rows and 6 columns from 1 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.

    >>> pr.read_bed(path, as_df=True, nrows=5)
      Chromosome  Start    End      Name  Score Strand
    0       chr1   9916  10115  H3K27me3      5      -
    1       chr1   9939  10138  H3K27me3      7      +
    2       chr1   9951  10150  H3K27me3      8      -
    3       chr1   9953  10152  H3K27me3      5      +
    4       chr1   9978  10177  H3K27me3      7      -
    """

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

    if not as_df:
        return PyRanges(df)
    else:
        return df


def read_bam(f, sparse=True, as_df=False, mapq=0, required_flag=0, filter_flag=1540):

    """Return bam file as PyRanges.

    Parameters
    ----------
    f : str

        Path to bam file

    sparse : bool, default True

        Whether to return only.

    as_df : bool, default False

        Whether to return as pandas DataFrame instead of PyRanges.

    mapq : int, default 0

        Minimum mapping quality score.

    required_flag : int, default 0

        Flags which must be present for the interval to be read.

    filter_flag : int, default 1540

        Ignore reads with these flags. Default 1540, which means that either
        the read is unmapped, the read failed vendor or platfrom quality
        checks, or the read is a PCR or optical duplicate.

    Notes
    -----

    This functionality requires the library `bamread`. It can be installed with
    `pip install bamread` or `conda install -c bioconda bamread`.

    Examples
    --------

    >>> path = pr.get_example_path("control.bam")
    >>> pr.read_bam(path)
    +--------------+-----------+-----------+--------------+------------+
    | Chromosome   | Start     | End       | Strand       | Flag       |
    | (category)   | (int32)   | (int32)   | (category)   | (uint16)   |
    |--------------+-----------+-----------+--------------+------------|
    | chr1         | 887771    | 887796    | +            | 16         |
    | chr1         | 994660    | 994685    | +            | 16         |
    | chr1         | 1770383   | 1770408   | +            | 16         |
    | chr1         | 1995141   | 1995166   | +            | 16         |
    | ...          | ...       | ...       | ...          | ...        |
    | chrY         | 57402214  | 57402239  | +            | 16         |
    | chrY         | 10643526  | 10643551  | -            | 0          |
    | chrY         | 11776321  | 11776346  | -            | 0          |
    | chrY         | 20557165  | 20557190  | -            | 0          |
    +--------------+-----------+-----------+--------------+------------+
    Stranded PyRanges object has 10,000 rows and 5 columns from 25 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.

    >>> pr.read_bam(path, sparse=False, as_df=True)
         Chromosome    Start      End Strand  Flag  QueryStart  QueryEnd Name Cigar Quality
    0          chr1   887771   887796      +    16           0        25   U0   25M    None
    1          chr1   994660   994685      +    16           0        25   U0   25M    None
    2          chr1  1041102  1041127      -     0           0        25   U0   25M    None
    3          chr1  1770383  1770408      +    16           0        25   U0   25M    None
    4          chr1  1995141  1995166      +    16           0        25   U0   25M    None
    ...         ...      ...      ...    ...   ...         ...       ...  ...   ...     ...
    9995       chrM     3654     3679      -     0           0        25   U0   25M    None
    9996       chrM     3900     3925      +    16           0        25   U0   25M    None
    9997       chrM    13006    13031      +    16           0        25   U0   25M    None
    9998       chrM    14257    14282      -     0           0        25   U0   25M    None
    9999       chrM    14257    14282      -     0           0        25   U0   25M    None
    <BLANKLINE>
    [10000 rows x 10 columns]
    """

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

    if as_df:
        return df
    else:
        return PyRanges(df)

    # return bamread.read_bam(f, mapq, required_flag, filter_flag)


def _fetch_gene_transcript_exon_id(attribute, annotation=None):

    no_quotes = attribute.str.replace('"', '').str.replace("'", "")

    df = no_quotes.str.extract(
        "gene_id.?(.+?);(?:.*transcript_id.?(.+?);)?(?:.*exon_number.?(.+?);)?(?:.*exon_id.?(.+?);)?",
        expand=True)  # .iloc[:, [1, 2, 3]]

    df.columns = "gene_id transcript_id exon_number exon_id".split()

    if annotation == "ensembl":
        newdf = []
        for c in "gene_id transcript_id exon_id".split():
            r = df[c].astype(str).str.extract(r'(\d+)').astype(float)
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


def read_gtf(f, full=True, as_df=False, nrows=None, duplicate_attr=False):

    """Read files in the Gene Transfer Format.

    Parameters
    ----------
    f : str

        Path to GTF file.

    as_df : bool, default False

        Whether to return as pandas DataFrame instead of PyRanges.

    nrows : int, default None

        Number of rows to read. Default None, i.e. all.

    duplicate_attr : bool, default False

        Whether to handle (potential) duplicate attributes or just keep last one.

    See Also
    --------

    pyranges.read_gff3 : read files in the General Feature Format

    Examples
    --------

    >>> path = pr.get_example_path("ensembl.gtf")
    >>> gr = pr.read_gtf(path)

    >>> # +--------------+------------+--------------+-----------+-----------+------------+--------------+------------+-----------------+----------------+-------+
    >>> # | Chromosome   | Source     | Feature      | Start     | End       | Score      | Strand       | Frame      | gene_id         | gene_version   | +18   |
    >>> # | (category)   | (object)   | (category)   | (int32)   | (int32)   | (object)   | (category)   | (object)   | (object)        | (object)       | ...   |
    >>> # |--------------+------------+--------------+-----------+-----------+------------+--------------+------------+-----------------+----------------+-------|
    >>> # | 1            | havana     | gene         | 11868     | 14409     | .          | +            | .          | ENSG00000223972 | 5              | ...   |
    >>> # | 1            | havana     | transcript   | 11868     | 14409     | .          | +            | .          | ENSG00000223972 | 5              | ...   |
    >>> # | 1            | havana     | exon         | 11868     | 12227     | .          | +            | .          | ENSG00000223972 | 5              | ...   |
    >>> # | 1            | havana     | exon         | 12612     | 12721     | .          | +            | .          | ENSG00000223972 | 5              | ...   |
    >>> # | ...          | ...        | ...          | ...       | ...       | ...        | ...          | ...        | ...             | ...            | ...   |
    >>> # | 1            | ensembl    | transcript   | 120724    | 133723    | .          | -            | .          | ENSG00000238009 | 6              | ...   |
    >>> # | 1            | ensembl    | exon         | 133373    | 133723    | .          | -            | .          | ENSG00000238009 | 6              | ...   |
    >>> # | 1            | ensembl    | exon         | 129054    | 129223    | .          | -            | .          | ENSG00000238009 | 6              | ...   |
    >>> # | 1            | ensembl    | exon         | 120873    | 120932    | .          | -            | .          | ENSG00000238009 | 6              | ...   |
    >>> # +--------------+------------+--------------+-----------+-----------+------------+--------------+------------+-----------------+----------------+-------+
    >>> # Stranded PyRanges object has 95 rows and 28 columns from 1 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    >>> # 18 hidden columns: gene_name, gene_source, gene_biotype, transcript_id, transcript_version, transcript_name, transcript_source, transcript_biotype, tag, transcript_support_level, ... (+ 8 more.)
    """

    _skiprows = skiprows(f)

    gr = read_gtf_full(f, as_df, nrows, _skiprows, duplicate_attr)

    return gr


def read_gtf_full(f, as_df=False, nrows=None, skiprows=0, duplicate_attr=False):

    dtypes = {
        "Chromosome": "category",
        "Feature": "category",
        "Strand": "category"
    }

    names = "Chromosome Source Feature Start End Score Strand Frame Attribute".split(
    )

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

    if not as_df:
        return PyRanges(df)
    else:
        return df


def to_rows(anno):
    rowdicts = []
    try:
        l = anno.head(1)
        for l in l:
            l.replace('"', '').replace(";", "").split()
    except AttributeError:
        raise Exception("Invalid attribute string: {l}. If the file is in GFF3 format, use pr.read_gff3 instead.".format(l=l))

    for l in anno:
        rowdicts.append({k: v
                         for k, v in [kv.replace('"', '').split(None, 1)
                                      # l[:-1] removes final ";" cheaply
                                      for kv in l[:-1].split("; ")]})

    # for l in anno:
    #     l = l.replace('"', '').replace(";", "").split()
    #     rowdicts.append({k: v for k, v in zip(*([iter(l)] * 2))})

    return pd.DataFrame.from_dict(rowdicts).set_index(anno.index)


def to_rows_keep_duplicates(anno):
    rowdicts = []
    for l in anno:
        rowdict = {}

        # l[:-1] removes final ";" cheaply
        for k, v in (kv.replace('"', '').split(None, 1) for kv in l[:-1].split("; ")):

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
                        as_df=False,
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

        extract = _fetch_gene_transcript_exon_id(df.Attribute)
        extract.columns = "gene_id transcript_id exon_number exon_id".split()

        extract.exon_number = extract.exon_number.astype(float)

        df = pd.concat([df[cols_to_concat], extract], axis=1, sort=False)

        dfs.append(df)

    df = pd.concat(dfs, sort=False)

    df.loc[:, "Start"] = df.Start - 1

    if not as_df:
        return PyRanges(df)
    else:
        return df


def to_rows_gff3(anno):
    rowdicts = []

    for l in list(anno):
        l = (it.split("=") for it in l.split(";"))
        rowdicts.append({k: v for k, v in l})

    return pd.DataFrame.from_dict(rowdicts).set_index(anno.index)


def read_gff3(f, annotation=None, as_df=False, nrows=None, skiprows=0):

    """Read files in the General Feature Format.

    Parameters
    ----------
    f : str

        Path to GFF file.

    as_df : bool, default False

        Whether to return as pandas DataFrame instead of PyRanges.

    nrows : int, default None

        Number of rows to read. Default None, i.e. all.

    See Also
    --------

    pyranges.read_gtf : read files in the Gene Transfer Format
    """


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

    if not as_df:
        return PyRanges(df)
    else:
        return df


def read_bigwig(f, as_df=False):

    try:
        import pyBigWig
    except ModuleNotFoundError:
        print("bwread must be installed to read bigwigs. Use `conda install -c bioconda bwread` or `pip install bwread` to install it.")
        import sys
        sys.exit(1)

    """Read bigwig files.

    Parameters
    ----------
    f : str

        Path to bw file.

    as_df : bool, default False

        Whether to return as pandas DataFrame instead of PyRanges.

    Examples
    --------

    >>> f = pr.get_example_path("bw.bw")
    >>> gr = pr.read_bigwig(f)
    >>> gr
    """

    from natsort import natsorted

    bw = pyBigWig.open(f)

    size = int(1e5)
    chromosomes = bw.chroms()

    dfs = {}

    for chromosome in natsorted(chromosomes):
        outstarts = []
        outends = []
        outvalues = []

        length = chromosomes[chromosome]

        starts = list(range(0, length, size))
        ends = list(range(size, length + size, size))
        ends[-1] = length
        for start, end in zip(starts, ends):
            intervals = bw.intervals(chromosome, start, end)
            if intervals is not None:
                for s, e, v in intervals:
                    outstarts.append(s)
                    outends.append(e)
                    outvalues.append(v)

        outstarts = pd.Series(outstarts)
        outends = pd.Series(outends)
        outvalues = pd.Series(outvalues)
        dfs[chromosome] = pd.DataFrame({"Chromosome": chromosome, "Start": outstarts, "End": outends, "Value": outvalues})

    return pr.PyRanges(dfs)
