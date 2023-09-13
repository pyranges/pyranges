import sys

import pandas as pd

import pyranges as pr  # noqa: F401


def get_sequence(gr, path=None, pyfaidx_fasta=None):
    """Get the sequence of the intervals from a fasta file

    Parameters
    ----------
    gr : PyRanges

        Coordinates.

    path : str

        Path to fasta file. It will be indexed using pyfaidx if an index is not found

    pyfaidx_fasta : pyfaidx.Fasta

        Alternative method to provide fasta target, as a pyfaidx.Fasta object


    Returns
    -------
    Series

        Sequences, one per interval.

    Note
    ----

    This function requires the library pyfaidx, it can be installed with
    ``conda install -c bioconda pyfaidx`` or ``pip install pyfaidx``.

    Sorting the PyRanges is likely to improve the speed.
    Intervals on the negative strand will be reverse complemented.

    Warning
    -------

    Note that the names in the fasta header and gr must be the same.

    See also
    --------
    get_transcript_sequence : obtain mRNA sequences, by joining exons belonging to the same transcript


    Examples
    --------

    >>> gr = pr.from_dict({"Chromosome": ["chr1", "chr1"],
    ...                    "Start": [5, 0], "End": [8, 5]})

    >>> gr
    +--------------+-----------+-----------+
    | Chromosome   |     Start |       End |
    | (category)   |   (int64) |   (int64) |
    |--------------+-----------+-----------|
    | chr1         |         5 |         8 |
    | chr1         |         0 |         5 |
    +--------------+-----------+-----------+
    Unstranded PyRanges object has 2 rows and 3 columns from 1 chromosomes.
    For printing, the PyRanges was sorted on Chromosome.

    >>> tmp_handle = open("temp.fasta", "w+")
    >>> _ = tmp_handle.write(">chr1\\n")
    >>> _ = tmp_handle.write("ATTACCAT\\n")
    >>> tmp_handle.close()

    >>> seq = pr.get_sequence(gr, "temp.fasta")

    >>> seq
    0      CAT
    1    ATTAC
    dtype: object

    >>> gr.seq = seq
    >>> gr
    +--------------+-----------+-----------+------------+
    | Chromosome   |     Start |       End | seq        |
    | (category)   |   (int64) |   (int64) | (object)   |
    |--------------+-----------+-----------+------------|
    | chr1         |         5 |         8 | CAT        |
    | chr1         |         0 |         5 | ATTAC      |
    +--------------+-----------+-----------+------------+
    Unstranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
    For printing, the PyRanges was sorted on Chromosome.

    """

    try:
        import pyfaidx  # type: ignore
    except ImportError:
        print(
            "pyfaidx must be installed to get fasta sequences. Use `conda install -c bioconda pyfaidx` or `pip install pyfaidx` to install it."
        )
        sys.exit(1)

    if pyfaidx_fasta is None:
        if path is None:
            raise Exception("ERROR get_sequence : you must provide a fasta path or pyfaidx_fasta object")
        pyfaidx_fasta = pyfaidx.Fasta(path, read_ahead=int(1e5))

    seqs = []
    for k, df in gr:
        if type(k) is tuple:  # input is Stranded
            _fasta = pyfaidx_fasta[k[0]]
            if k[1] == "-":
                for start, end in zip(df.Start, df.End):
                    seqs.append((-_fasta[start:end]).seq)  # reverse complement
            else:
                for start, end in zip(df.Start, df.End):
                    seqs.append(_fasta[start:end].seq)

        else:
            _fasta = pyfaidx_fasta[k]
            for start, end in zip(df.Start, df.End):
                seqs.append(_fasta[start:end].seq)

    return pd.concat([pd.Series(s) for s in seqs]).reset_index(drop=True)


def get_fasta(*args, **kwargs):
    """Deprecated: this function has been moved to Pyranges.get_sequence"""
    return get_sequence(*args, **kwargs)


def get_transcript_sequence(gr, group_by, path=None, pyfaidx_fasta=None):
    """Get the sequence of mRNAs, e.g. joining intervals corresponding to exons of the same transcript

    Parameters
    ----------
    gr : PyRanges

        Coordinates.

    group_by : str or list of str

        intervals are grouped by this/these ID column(s): these are exons belonging to same transcript

    path : str

        Path to fasta file. It will be indexed using pyfaidx if an index is not found

    pyfaidx_fasta : pyfaidx.Fasta

        Alternative method to provide fasta target, as a pyfaidx.Fasta object


    Returns
    -------
    DataFrame

        Pandas DataFrame with a column for Sequence, plus ID column(s) provided with "group_by"

    Note
    ----

    This function requires the library pyfaidx, it can be installed with
    ``conda install -c bioconda pyfaidx`` or ``pip install pyfaidx``.

    Sorting the PyRanges is likely to improve the speed.
    Intervals on the negative strand will be reverse complemented.

    Warning
    -------

    Note that the names in the fasta header and gr must be the same.

    See also
    --------
    get_sequence : obtain sequence of single intervals


    Examples
    --------

    >>> gr = pr.from_dict({"Chromosome": ['chr1', 'chr1', 'chr1'],
    ...                    "Start": [0, 9, 18], "End": [4, 13, 21],
    ...                    "Strand":['+', '-', '-'],
    ...                    "transcript": ['t1', 't2', 't2']})
    >>> gr
    +--------------+-----------+-----------+--------------+--------------+
    | Chromosome   |     Start |       End | Strand       | transcript   |
    | (category)   |   (int64) |   (int64) | (category)   | (object)     |
    |--------------+-----------+-----------+--------------+--------------|
    | chr1         |         0 |         4 | +            | t1           |
    | chr1         |         9 |        13 | -            | t2           |
    | chr1         |        18 |        21 | -            | t2           |
    +--------------+-----------+-----------+--------------+--------------+
    Stranded PyRanges object has 3 rows and 5 columns from 1 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.


    >>> tmp_handle = open("temp.fasta", "w+")
    >>> _ = tmp_handle.write(">chr1\\n")
    >>> _ = tmp_handle.write("AAACCCTTTGGGAAACCCTTTGGG\\n")
    >>> tmp_handle.close()

    >>> seq = pr.get_transcript_sequence(gr, path="temp.fasta", group_by='transcript')
    >>> seq
      transcript Sequence
    0         t1     AAAC
    1         t2  AAATCCC

    To write to a file in fasta format:
    # with open('outfile.fasta', 'w') as fw:
    #     nchars=60
    #     for row in seq.itertuples():
    #         s = '\\n'.join([ row.Sequence[i:i+nchars] for i in range(0, len(row.Sequence), nchars)])
    #         fw.write(f'>{row.transcript}\\n{s}\\n')
    """

    if gr.stranded:
        gr = gr.sort("5")
    else:
        gr = gr.sort()

    z = gr.df
    z["Sequence"] = get_sequence(gr, path=path, pyfaidx_fasta=pyfaidx_fasta)

    return z.groupby(group_by, as_index=False).agg({"Sequence": "".join})
