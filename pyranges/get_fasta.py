
import sys
import pyranges as pr
import pandas as pd

def get_fasta(gr, path=None, pyfaidx_fasta=None):
    """Get fasta sequence.

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

    Sorting the PyRanges is likely to improve the speed.
    Intervals on the negative strand will be reverse complemented.

    Warning
    -------

    Note that the names in the fasta header and gr must be the same.

    Examples
    --------

    >>> gr = pr.from_dict({"Chromosome": ["chr1", "chr1"],
    ...                    "Start": [5, 0], "End": [8, 5]})

    >>> gr
    +--------------+-----------+-----------+
    | Chromosome   |     Start |       End |
    | (category)   |   (int32) |   (int32) |
    |--------------+-----------+-----------|
    | chr1         |         5 |         8 |
    | chr1         |         0 |         5 |
    +--------------+-----------+-----------+
    Unstranded PyRanges object has 2 rows and 3 columns from 1 chromosomes.
    For printing, the PyRanges was sorted on Chromosome.

    >>> tmp_handle = open("temp.fasta", "w+")
    >>> _ = tmp_handle.write("> chr1\\n")
    >>> _ = tmp_handle.write("ATTACCAT")
    >>> tmp_handle.close()

    >>> seq = pr.get_fasta(gr, "temp.fasta")

    >>> seq
    0      CAT
    1    ATTAC
    dtype: object

    >>> gr.seq = seq
    >>> gr
    +--------------+-----------+-----------+------------+
    | Chromosome   |     Start |       End | seq        |
    | (category)   |   (int32) |   (int32) | (object)   |
    |--------------+-----------+-----------+------------|
    | chr1         |         5 |         8 | CAT        |
    | chr1         |         0 |         5 | ATTAC      |
    +--------------+-----------+-----------+------------+
    Unstranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
    For printing, the PyRanges was sorted on Chromosome.

    """

    try:
        import pyfaidx
    except ModuleNotFoundError as e:
        print("pyfaidx must be installed to get fasta sequences. Use `conda install -c bioconda pyfaidx` or `pip install pyfaidx` to install it.")
        sys.exit(1)

    if pyfaidx_fasta is None:
        pyfaidx_fasta = pyfaidx.Fasta(path, read_ahead=int(1e5))

    seqs = []
    for k, df in gr:
        if type(k) is tuple:  #input is Stranded
            _fasta = pyfaidx_fasta[k[0]]
            if k[1]=='-':
                for start, end in zip(df.Start, df.End):
                    seqs.append( (-_fasta[start:end]).seq ) # reverse complement
            else:
                for start, end in zip(df.Start, df.End):
                    seqs.append(_fasta[start:end].seq)

        else:        
            _fasta = pyfaidx_fasta[k]
            for start, end in zip(df.Start, df.End):
                seqs.append(_fasta[start:end].seq)

    return pd.concat([pd.Series(s) for s in seqs]).reset_index(drop=True)
