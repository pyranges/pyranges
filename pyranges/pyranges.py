"""Data structure for genomic intervals and their annotation."""

import pandas as pd
import numpy as np

from natsort import natsorted

import pyranges as pr

from pyranges.tostring2 import tostring

from pyranges.methods.intersection import _intersection, _overlap
from pyranges.multithreaded import pyrange_apply, pyrange_apply_single, pyrange_apply_chunks, _slack, _tes, _tss

__all__ = ["PyRanges"]

def fill_kwargs(kwargs):
    """Give the kwargs dict default options."""

    defaults = {
        "strandedness": None,
        "overlap": True,
        "how": None,
        "invert": None,
        "new_pos": None,
        "suffixes": ["_a", "_b"],
        "suffix": "_b",
        "sparse": {
            "self": False,
            "other": False
        }
    }

    defaults.update(kwargs)

    return defaults


class PyRanges():

    """Two-dimensional representation of genomic intervals and their annotations.

    A PyRanges object must have the columns Chromosome, Start and End. These
    describe the genomic position and function as implicit row labels. A Strand
    column is optional and adds strand information to the intervals. Any other
    columns are allowed and are considered metadata.

    Operations between PyRanges align intervals based on their position.

    If a PyRanges is built using the arguments chromosomes, starts, ends and
    optionally strands, all non-scalars must be of the same length.

    Parameters
    ----------
    df : pandas.DataFrame or dict of pandas.DataFrame, default None
        The data to be stored in the PyRanges.

    chromosomes : array-like or scalar value, default None
        The chromosome(s) in the PyRanges.

    starts : array-like, default None
        The start postions in the PyRanges.

    ends : array-like, default None
        The end postions in the PyRanges.

    strands : array-like or scalar value, default None
        The strands in the PyRanges.

    int64 : bool, default False
        Use np.int64 to represent starts and ends

    copy_df : bool, default True
        Copy input pandas.DataFrame

    See Also
    --------

    pyranges.read_bed: read bed-file into PyRanges
    pyranges.read_bam: read bam-file into PyRanges
    pyranges.read_gff: read gff-file into PyRanges
    pyranges.read_gtf: read gtf-file into PyRanges
    pyranges.from_dict: create PyRanges from dict of columns
    pyranges.from_string: create PyRanges from multiline string

    Notes
    -----

    A PyRanges object is represented internally as a dictionary efficiency. The keys are
    chromosomes or chromosome/strand tuples and the values are pandas DataFrames.

    Examples
    --------

    >>> pr.PyRanges()
    Empty PyRanges

    >>> pr.PyRanges(chromosomes="chr1", starts=(1, 5), ends=[3, 149],
    ...             strands=("+", "-"), int64=True)
    +--------------+-----------+-----------+--------------+
    | Chromosome   |     Start |       End | Strand       |
    | (category)   |   (int64) |   (int64) | (category)   |
    |--------------+-----------+-----------+--------------|
    | chr1         |         1 |         3 | +            |
    | chr1         |         5 |       149 | -            |
    +--------------+-----------+-----------+--------------+
    Stranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.

    >>> df = pd.DataFrame({"Chromosome": ["chr1", "chr2"], "Start": [100, 200],
    ...                    "End": [150, 201]})
    >>> df
      Chromosome  Start  End
    0       chr1    100  150
    1       chr2    200  201
    >>> pr.PyRanges(df)
    +--------------+-----------+-----------+
    | Chromosome   |     Start |       End |
    | (category)   |   (int32) |   (int32) |
    |--------------+-----------+-----------|
    | chr1         |       100 |       150 |
    | chr2         |       200 |       201 |
    +--------------+-----------+-----------+
    Unstranded PyRanges object has 2 rows and 3 columns from 2 chromosomes.
    For printing, the PyRanges was sorted on Chromosome.


    >>> gr = pr.from_dict({"Chromosome": [1, 1], "Strand": ["+", "-"], "Start": [1, 4], "End": [2, 27],
    ...                    "TP": [0, 1], "FP": [12, 11], "TN": [10, 9], "FN": [2, 3]})
    >>> gr
    +--------------+--------------+-----------+-----------+-----------+-----------+-----------+-----------+
    |   Chromosome | Strand       |     Start |       End |        TP |        FP |        TN |        FN |
    |   (category) | (category)   |   (int32) |   (int32) |   (int64) |   (int64) |   (int64) |   (int64) |
    |--------------+--------------+-----------+-----------+-----------+-----------+-----------+-----------|
    |            1 | +            |         1 |         2 |         0 |        12 |        10 |         2 |
    |            1 | -            |         4 |        27 |         1 |        11 |         9 |         3 |
    +--------------+--------------+-----------+-----------+-----------+-----------+-----------+-----------+
    Stranded PyRanges object has 2 rows and 8 columns from 1 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    dfs = None
    """Dict mapping chromosomes or chromosome/strand pairs to pandas DataFrames."""

    features = None
    """Namespace for genomic-features methods.

    See Also
    --------
    pyranges.genomicfeatures : namespace for feature-functionality
    pyranges.genomicfeatures.GenomicFeaturesMethods : namespace for feature-functionality
    """

    stats = None
    """Namespace for statistcal methods.

    See Also
    --------
    pyranges.statistics : namespace for statistics
    pyranges.stats.StatisticsMethods : namespace for statistics
    """

    def __init__(self,
                 df=None,
                 chromosomes=None,
                 starts=None,
                 ends=None,
                 strands=None,
                 int64=False,
                 copy_df=True):

        from pyranges.methods.init import _init

        if df is None and chromosomes is None:
            df = pd.DataFrame(columns="Chromosome Start End".split())

        _init(self, df, chromosomes, starts, ends, strands, int64, copy_df)


    def __array_ufunc__(self, *args, **kwargs):

        """Apply unary numpy-function.


        Apply function to all columns which are not index, i.e. Chromosome,
        Start, End nor Strand.

        Notes
        -----

        Function must produce a vector of equal length.

        Examples
        --------

        >>> gr = pr.from_dict({"Chromosome": [1, 2, 3], "Start": [1, 2, 3],
        ... "End": [2, 3, 4], "Score": [9, 16, 25], "Score2": [121, 144, 169],
        ... "Name": ["n1", "n2", "n3"]})
        >>> gr
        +--------------+-----------+-----------+-----------+-----------+------------+
        |   Chromosome |     Start |       End |     Score |    Score2 | Name       |
        |   (category) |   (int32) |   (int32) |   (int64) |   (int64) | (object)   |
        |--------------+-----------+-----------+-----------+-----------+------------|
        |            1 |         1 |         2 |         9 |       121 | n1         |
        |            2 |         2 |         3 |        16 |       144 | n2         |
        |            3 |         3 |         4 |        25 |       169 | n3         |
        +--------------+-----------+-----------+-----------+-----------+------------+
        Unstranded PyRanges object has 3 rows and 6 columns from 3 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> np.sqrt(gr)
        +--------------+-----------+-----------+-------------+-------------+------------+
        |   Chromosome |     Start |       End |       Score |      Score2 | Name       |
        |   (category) |   (int32) |   (int32) |   (float64) |   (float64) | (object)   |
        |--------------+-----------+-----------+-------------+-------------+------------|
        |            1 |         1 |         2 |           3 |          11 | n1         |
        |            2 |         2 |         3 |           4 |          12 | n2         |
        |            3 |         3 |         4 |           5 |          13 | n3         |
        +--------------+-----------+-----------+-------------+-------------+------------+
        Unstranded PyRanges object has 3 rows and 6 columns from 3 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        func, call, gr = args

        columns = list(gr.columns)
        non_index = [c for c in columns if c not in ["Chromosome", "Start", "End", "Strand"]]

        for chromosome, df in gr:
            subset = df.head(1)[non_index].select_dtypes(include=np.number).columns
            _v = getattr(func, call)(df[subset], **kwargs)
            # print(_v)
            # print(df[_c])
            df[subset] = _v


        return gr


        # self.apply()

    def __getattr__(self, name):

        """Return column.

        Parameters
        ----------
        name : str

            Column to return

        Returns
        -------
        pandas.Series

        Example
        -------

        >>> gr = pr.from_dict({"Chromosome": [1, 1, 1], "Start": [0, 100, 250], "End": [10, 125, 251]})
        >>> gr.Start
        0      0
        1    100
        2    250
        Name: Start, dtype: int32
        """

        from pyranges.methods.attr import _getattr

        return _getattr(self, name)

    def __setattr__(self, column_name, column):

        """Insert or update column.

        Parameters
        ----------
        column_name : str

            Name of column to update or insert.

        column : list, np.array or pd.Series

            Data to insert.

        Example
        -------

        >>> gr = pr.from_dict({"Chromosome": [1, 1, 1], "Start": [0, 100, 250], "End": [10, 125, 251]})
        >>> gr.Start = np.array([1, 1, 2])
        >>> gr
        +--------------+-----------+-----------+
        |   Chromosome |     Start |       End |
        |   (category) |   (int64) |   (int32) |
        |--------------+-----------+-----------|
        |            1 |         1 |        10 |
        |            1 |         1 |       125 |
        |            1 |         2 |       251 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 3 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        from pyranges.methods.attr import _setattr

        if column_name == "columns":
            dfs = {}
            for k, df in self:
                df.columns = column
                dfs[k] = df
            self.__dict__["dfs"] = dfs

        else:
            _setattr(self, column_name, column)

            if column_name in ["Start", "End"]:
                if self.dtypes["Start"] != self.dtypes["End"]:
                    print("Warning! Start and End columns now have different dtypes: {} and {}".format(
                        self.dtypes["Start"], self.dtypes["End"]))

    def __getitem__(self, val):

        """Fetch columns or subset on position.

        If a list is provided, the column(s) in the list is returned. This subsets on columns.

        If a numpy array is provided, it must be of type bool and the same length as the PyRanges.

        Otherwise, a subset of the rows is returned with the location info provided.

        Parameters
        ----------
        val : bool array/Series, tuple, list, str or slice

            Data to fetch.

        Examples
        --------

        >>> gr = pr.data.ensembl_gtf()
        >>> gr.columns
        Index(['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand',
            'Frame', 'gene_biotype', 'gene_id', 'gene_name', 'gene_source',
            'gene_version', 'tag', 'transcript_biotype', 'transcript_id',
            'transcript_name', 'transcript_source', 'transcript_support_level',
            'transcript_version', 'exon_id', 'exon_number', 'exon_version',
            '(assigned', 'previous', 'protein_id', 'protein_version', 'ccds_id'],
            dtype='object')

        >>> gr = gr[["Source", "Feature", "gene_id"]]
        >>> gr
        +--------------+------------+--------------+-----------+-----------+--------------+-----------------+
        | Chromosome   | Source     | Feature      | Start     | End       | Strand       | gene_id         |
        | (category)   | (object)   | (category)   | (int32)   | (int32)   | (category)   | (object)        |
        |--------------+------------+--------------+-----------+-----------+--------------+-----------------|
        | 1            | havana     | gene         | 11868     | 14409     | +            | ENSG00000223972 |
        | 1            | havana     | transcript   | 11868     | 14409     | +            | ENSG00000223972 |
        | 1            | havana     | exon         | 11868     | 12227     | +            | ENSG00000223972 |
        | 1            | havana     | exon         | 12612     | 12721     | +            | ENSG00000223972 |
        | ...          | ...        | ...          | ...       | ...       | ...          | ...             |
        | 1            | havana     | gene         | 1173055   | 1179555   | -            | ENSG00000205231 |
        | 1            | havana     | transcript   | 1173055   | 1179555   | -            | ENSG00000205231 |
        | 1            | havana     | exon         | 1179364   | 1179555   | -            | ENSG00000205231 |
        | 1            | havana     | exon         | 1173055   | 1176396   | -            | ENSG00000205231 |
        +--------------+------------+--------------+-----------+-----------+--------------+-----------------+
        Stranded PyRanges object has 2,446 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        Create boolean Series and use it to subset:

        >>> s = (gr.Feature == "gene") | (gr.gene_id == "ENSG00000223972")
        >>> gr[s]
        +--------------+----------------+--------------+-----------+-----------+--------------+-----------------+
        | Chromosome   | Source         | Feature      | Start     | End       | Strand       | gene_id         |
        | (category)   | (object)       | (category)   | (int32)   | (int32)   | (category)   | (object)        |
        |--------------+----------------+--------------+-----------+-----------+--------------+-----------------|
        | 1            | havana         | gene         | 11868     | 14409     | +            | ENSG00000223972 |
        | 1            | havana         | transcript   | 11868     | 14409     | +            | ENSG00000223972 |
        | 1            | havana         | exon         | 11868     | 12227     | +            | ENSG00000223972 |
        | 1            | havana         | exon         | 12612     | 12721     | +            | ENSG00000223972 |
        | ...          | ...            | ...          | ...       | ...       | ...          | ...             |
        | 1            | havana         | gene         | 1062207   | 1063288   | -            | ENSG00000273443 |
        | 1            | ensembl_havana | gene         | 1070966   | 1074306   | -            | ENSG00000237330 |
        | 1            | ensembl_havana | gene         | 1081817   | 1116361   | -            | ENSG00000131591 |
        | 1            | havana         | gene         | 1173055   | 1179555   | -            | ENSG00000205231 |
        +--------------+----------------+--------------+-----------+-----------+--------------+-----------------+
        Stranded PyRanges object has 95 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> cs = pr.data.chipseq()
        >>> cs[10000:100000]
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr2         |     33241 |     33266 | U0         |         0 | +            |
        | chr2         |     13611 |     13636 | U0         |         0 | -            |
        | chr2         |     32620 |     32645 | U0         |         0 | -            |
        | chr3         |     87179 |     87204 | U0         |         0 | +            |
        | chr4         |     45413 |     45438 | U0         |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 5 rows and 6 columns from 3 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> cs["chr1", "-"]
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   | Start     | End       | Name       | Score     | Strand       |
        | (category)   | (int32)   | (int32)   | (object)   | (int64)   | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         | 100079649 | 100079674 | U0         | 0         | -            |
        | chr1         | 223587418 | 223587443 | U0         | 0         | -            |
        | chr1         | 202450161 | 202450186 | U0         | 0         | -            |
        | chr1         | 156338310 | 156338335 | U0         | 0         | -            |
        | ...          | ...       | ...       | ...        | ...       | ...          |
        | chr1         | 203557775 | 203557800 | U0         | 0         | -            |
        | chr1         | 28114107  | 28114132  | U0         | 0         | -            |
        | chr1         | 21622765  | 21622790  | U0         | 0         | -            |
        | chr1         | 80668132  | 80668157  | U0         | 0         | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 437 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> cs["chr5", "-", 90000:]
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   | Start     | End       | Name       | Score     | Strand       |
        | (category)   | (int32)   | (int32)   | (object)   | (int64)   | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr5         | 399682    | 399707    | U0         | 0         | -            |
        | chr5         | 1847502   | 1847527   | U0         | 0         | -            |
        | chr5         | 5247533   | 5247558   | U0         | 0         | -            |
        | chr5         | 5300394   | 5300419   | U0         | 0         | -            |
        | ...          | ...       | ...       | ...        | ...       | ...          |
        | chr5         | 178786234 | 178786259 | U0         | 0         | -            |
        | chr5         | 179268931 | 179268956 | U0         | 0         | -            |
        | chr5         | 179289594 | 179289619 | U0         | 0         | -            |
        | chr5         | 180513795 | 180513820 | U0         | 0         | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 285 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> cs["chrM"]
        Empty PyRanges
        """

        from pyranges.methods.getitem import _getitem

        return _getitem(self, val)


    def __iter__(self):

        """Iterate over the keys and values.

        See Also
        --------
        pyranges.iter : iterate over multiple PyRanges

        Examples
        --------
        >>> gr = pr.from_dict({"Chromosome": [1, 1, 1], "Start": [0, 100, 250],
        ...                   "End": [10, 125, 251], "Strand": ["+", "+", "-"]})

        >>> for k, v in gr:
        ...     print(k)
        ...     print(v)
        ('1', '+')
        Chromosome  Start  End Strand
        0          1      0   10      +
        1          1    100  125      +
        ('1', '-')
        Chromosome  Start  End Strand
        2          1    250  251      -
        """

        return iter(self.items())


    def __len__(self):
        """Return the number of intervals in the PyRanges."""
        return sum([len(d) for d in self.values()])


    def __str__(self):

        """Return string representation."""

        return tostring(self)

    def __repr__(self):

        """Return REPL representation."""

        return str(self)

    def _repr_html_(self):

        """Return REPL HTML representation for Jupyter Noteboooks."""

        return self.df._repr_html_()

    def apply(self, f, strand=None, as_pyranges=True, nb_cpu=1, **kwargs):

        """Apply a function to the PyRanges.

        Parameters
        ----------
        f : function
            Function to apply on each DataFrame in a PyRanges

        strand : bool, default None, i.e. auto

            Whether to do operations on chromosome/strand pairs or chromosomes. If None, will use
            chromosome/strand pairs if the PyRanges is stranded.

        as_pyranges : bool, default True

            Whether to return as a PyRanges or dict. If `f` does not return a DataFrame valid for
            PyRanges, `as_pyranges` must be False.

        nb_cpu: int, default 1

            How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
            Will only lead to speedups on large datasets.

        **kwargs
            Additional keyword arguments to pass as keyword arguments to `f`

        Returns
        -------
        PyRanges or dict
            Result of applying f to each DataFrame in the PyRanges

        See also
        --------

        pyranges.PyRanges.apply_pair: apply a function to a pair of PyRanges
        pyranges.PyRanges.apply_chunks: apply a row-based function to a PyRanges in parallel

        Note
        ----

        This is the function used internally to carry out almost all unary PyRanges methods.

        Examples
        --------

        >>> gr = pr.from_dict({"Chromosome": [1, 1, 2, 2], "Strand": ["+", "+", "-", "+"],
        ...                    "Start": [1, 4, 2, 9], "End": [2, 27, 13, 10]})
        >>> gr
        +--------------+--------------+-----------+-----------+
        |   Chromosome | Strand       |     Start |       End |
        |   (category) | (category)   |   (int32) |   (int32) |
        |--------------+--------------+-----------+-----------|
        |            1 | +            |         1 |         2 |
        |            1 | +            |         4 |        27 |
        |            2 | +            |         9 |        10 |
        |            2 | -            |         2 |        13 |
        +--------------+--------------+-----------+-----------+
        Stranded PyRanges object has 4 rows and 4 columns from 2 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.apply(lambda df: len(df), as_pyranges=False)
        {('1', '+'): 2, ('2', '+'): 1, ('2', '-'): 1}

        >>> gr.apply(lambda df: len(df), as_pyranges=False, strand=False)
        {'1': 2, '2': 2}

        >>> def add_to_ends(df, **kwargs):
        ...     df.loc[:, "End"] = kwargs["slack"] + df.End
        ...     return df
        >>> gr.apply(add_to_ends, slack=500)
        +--------------+--------------+-----------+-----------+
        |   Chromosome | Strand       |     Start |       End |
        |   (category) | (category)   |   (int32) |   (int32) |
        |--------------+--------------+-----------+-----------|
        |            1 | +            |         1 |       502 |
        |            1 | +            |         4 |       527 |
        |            2 | +            |         9 |       510 |
        |            2 | -            |         2 |       513 |
        +--------------+--------------+-----------+-----------+
        Stranded PyRanges object has 4 rows and 4 columns from 2 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        if strand is None:
            strand = self.stranded

        kwargs.update({"strand": strand})
        kwargs.update(kwargs.get("kwargs", {}))
        kwargs = fill_kwargs(kwargs)

        result = pyrange_apply_single(f, self, **kwargs)

        if not as_pyranges:
            return result
        else:
            return PyRanges(result)


    def apply_chunks(self, f, as_pyranges=False, nb_cpu=1, **kwargs):

        """Apply a row-based function to arbitrary partitions of the PyRanges.

        apply_chunks speeds up the application of functions where the result is not affected by
        applying the function to ordered, non-overlapping splits of the data.

        Parameters
        ----------
        f : function
            Row-based or associative function to apply on the partitions.

        as_pyranges : bool, default False

            Whether to return as a PyRanges or dict.

        nb_cpu: int, default 1

            How many cpus to use. The data is split into nb_cpu partitions.

        **kwargs
            Additional keyword arguments to pass as keyword arguments to `f`

        Returns
        -------
        dict of lists
            Result of applying f to each partition of the DataFrames in the PyRanges.

        See also
        --------

        pyranges.PyRanges.apply_pair: apply a function to a pair of PyRanges
        pyranges.PyRanges.apply_chunks: apply a row-based function to a PyRanges in parallel

        Note
        ----

        apply_chunks will only lead to speedups on large datasets or slow-running functions. Using
        it with nb_cpu=1 is pointless; use apply instead.

        Examples
        --------

        >>> gr = pr.from_dict({"Chromosome": [1, 1, 1], "Start": [2, 3, 5], "End": [9, 4, 6]})
        >>> gr
        +--------------+-----------+-----------+
        |   Chromosome |     Start |       End |
        |   (category) |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        |            1 |         2 |         9 |
        |            1 |         3 |         4 |
        |            1 |         5 |         6 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 3 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.apply_chunks(
        ... lambda df, **kwargs: list(df.End + kwargs["add"]), nb_cpu=1, add=1000)
        {'1': [[1009, 1004, 1006]]}
        """

        kwargs.update(kwargs.get("kwargs", {}))
        kwargs = fill_kwargs(kwargs)

        result = pyrange_apply_chunks(f, self, as_pyranges, **kwargs)

        return result


    def apply_pair(self,
                   other,
                   f,
                   strandedness=None,
                   as_pyranges=True,
                   **kwargs):

        """Apply a function to a pair of PyRanges.

        The function is applied to each chromosome or chromosome/strand pair found in at least one
        of the PyRanges.

        Parameters
        ----------
        f : function
            Row-based or associative function to apply on the DataFrames.

        strandedness : {None, "same", "opposite", False}, default None, i.e. auto

            Whether to compare PyRanges on the same strand, the opposite or ignore strand
            information. The default, None, means use "same" if both PyRanges are strande,
            otherwise ignore the strand information.

        as_pyranges : bool, default False

            Whether to return as a PyRanges or dict. If `f` does not return a DataFrame valid for
            PyRanges, `as_pyranges` must be False.

        nb_cpu: int, default 1

            How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
            Will only lead to speedups on large datasets.

        **kwargs
            Additional keyword arguments to pass as keyword arguments to `f`

        Returns
        -------
        dict of lists
            Result of applying f to each partition of the DataFrames in the PyRanges.

        See also
        --------

        pyranges.PyRanges.apply_pair: apply a function to a pair of PyRanges
        pyranges.PyRanges.apply_chunks: apply a row-based function to a PyRanges in parallel
        pyranges.iter: iterate over two or more PyRanges

        Note
        ----

        This is the function used internally to carry out almost all comparison functions in
        PyRanges.

        Examples
        --------

        >>> gr = pr.data.chipseq()
        >>> gr2 = pr.data.chipseq_background()

        >>> gr.apply_pair(gr2, pr.methods.intersection._intersection) # same as gr.intersect(gr2)
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         | 226987603 | 226987617 | U0         |         0 | +            |
        | chr8         |  38747236 |  38747251 | U0         |         0 | -            |
        | chr15        |  26105515 |  26105518 | U0         |         0 | +            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 6 columns from 3 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> f1 = pr.data.f1()
        >>> f1
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         3 |         6 | interval1  |         0 | +            |
        | chr1         |         8 |         9 | interval3  |         0 | +            |
        | chr1         |         5 |         7 | interval2  |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> f2 = pr.data.f2()
        >>> f2
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         1 |         2 | a          |         0 | +            |
        | chr1         |         6 |         7 | b          |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> f1.apply_pair(f2, lambda df, df2: (len(df), len(df2)), as_pyranges=False)
        {('chr1', '+'): (2, 2), ('chr1', '-'): (1, 2)}

        """

        kwargs.update({"strandedness": strandedness})
        kwargs.update(kwargs.get("kwargs", {}))
        kwargs = fill_kwargs(kwargs)

        result = pyrange_apply(f, self, other, **kwargs)

        if not as_pyranges:
            return result
        else:
            return PyRanges(result)


    def as_df(self):

        """Return PyRanges as DataFrame.

        Returns
        -------
        DataFrame

            A DataFrame natural sorted on Chromosome and Strand. The ordering of rows within
            chromosomes and strands is preserved.

        See also
        --------

        PyRanges.df : Return PyRanges as DataFrame.

        Examples
        --------

        >>> gr = pr.from_dict({"Chromosome": [1, 1, 2, 2], "Start": [1, 2, 3, 9],
        ...                    "End": [3, 3, 10, 12], "Gene": ["A", "B", "C", "D"]})
        >>> gr
        +--------------+-----------+-----------+------------+
        |   Chromosome |     Start |       End | Gene       |
        |   (category) |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------|
        |            1 |         1 |         3 | A          |
        |            1 |         2 |         3 | B          |
        |            2 |         3 |        10 | C          |
        |            2 |         9 |        12 | D          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 4 rows and 4 columns from 2 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.as_df()
          Chromosome  Start  End Gene
        0          1      1    3    A
        1          1      2    3    B
        2          2      3   10    C
        3          2      9   12    D
        """

        if len(self) == 0:
            return pd.DataFrame()
        elif len(self) == 1:
            return self.values()[0]
        else:
            return pd.concat(self.values()).reset_index(drop=True)

    def assign(self, col, f, strand=None, nb_cpu=1, **kwargs):

        """Add or replace a column.

        Does not change the original PyRanges.

        Parameters
        ----------

        col : str

            Name of column.

        f : function
            Function to create new column.

        strand : bool, default None, i.e. auto

            Whether to do operations on chromosome/strand pairs or chromosomes. If None, will use
            chromosome/strand pairs if the PyRanges is stranded.

        nb_cpu: int, default 1

            How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
            Will only lead to speedups on large datasets.

        **kwargs
            Additional keyword arguments to pass as keyword arguments to `f`


        Returns
        -------
        PyRanges
            A copy of the PyRanges with the column inserted.

        Examples
        --------

        >>> gr = pr.from_dict({"Chromosome": [1, 1], "Start": [1, 2], "End": [3, 5],
        ... "Name": ["a", "b"]})
        >>> gr
        +--------------+-----------+-----------+------------+
        |   Chromosome |     Start |       End | Name       |
        |   (category) |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------|
        |            1 |         1 |         3 | a          |
        |            1 |         2 |         5 | b          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.assign("Blabla", lambda df: df.Chromosome.astype(str) + "_yadayada")
        +--------------+-----------+-----------+------------+------------+
        |   Chromosome |     Start |       End | Name       | Blabla     |
        |   (category) |   (int32) |   (int32) | (object)   | (object)   |
        |--------------+-----------+-----------+------------+------------|
        |            1 |         1 |         3 | a          | 1_yadayada |
        |            1 |         2 |         5 | b          | 1_yadayada |
        +--------------+-----------+-----------+------------+------------+
        Unstranded PyRanges object has 2 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        Note that assigning to an existing name replaces the column:

        >>> gr.assign("Name",
        ... lambda df, **kwargs: df.Start.astype(str) + kwargs["sep"] +
        ... df.Name.str.capitalize(), sep="_")
        +--------------+-----------+-----------+------------+
        |   Chromosome |     Start |       End | Name       |
        |   (category) |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------|
        |            1 |         1 |         3 | 1_A        |
        |            1 |         2 |         5 | 2_B        |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        self = self.copy()

        if strand is None:
            strand = self.stranded

        kwargs["strand"] = strand
        kwargs = fill_kwargs(kwargs)

        result = pyrange_apply_single(f, self, **kwargs)

        first_result = next(iter(result.values()))

        assert isinstance(
            first_result, pd.Series
        ), "result of assign function must be Series, but is {}".format(
            type(first_result))

        # do a deepcopy of object
        new_self = self.copy()
        new_self.__setattr__(col, result)

        return new_self


    @property
    def chromosomes(self):

        """Return chromosomes in natsorted order."""

        if self.stranded:
            return natsorted(set([k[0] for k in self.keys()]))
        else:
            return natsorted(set([k for k in self.keys()]))

    def cluster(self, strand=None, by=None, slack=0, count=False, nb_cpu=1):

        """Give overlapping intervals a common id.

        Parameters
        ----------
        strand : bool, default None, i.e. auto

            Whether to ignore strand information if PyRanges is stranded.

        by : str or list, default None

            Only intervals with an equal value in column(s) `by` are clustered.

        slack : int, default 0

            Consider intervals separated by less than `slack` to be in the same cluster. If `slack`
            is negative, intervals overlapping less than `slack` are not considered to be in the
            same cluster.

        nb_cpu: int, default 1

            How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
            Will only lead to speedups on large datasets.

        Returns
        -------
        PyRanges
            PyRanges with an ID-column "Cluster" added.

        See also
        --------

        PyRanges.merge: combine overlapping intervals into one

        Examples
        --------

        >>> gr = pr.from_dict({"Chromosome": [1, 1, 1, 1], "Start": [1, 2, 3, 9],
        ...                    "End": [3, 3, 10, 12], "Gene": [1, 2, 3, 3]})
        >>> gr
        +--------------+-----------+-----------+-----------+
        |   Chromosome |     Start |       End |      Gene |
        |   (category) |   (int32) |   (int32) |   (int64) |
        |--------------+-----------+-----------+-----------|
        |            1 |         1 |         3 |         1 |
        |            1 |         2 |         3 |         2 |
        |            1 |         3 |        10 |         3 |
        |            1 |         9 |        12 |         3 |
        +--------------+-----------+-----------+-----------+
        Unstranded PyRanges object has 4 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.cluster()
        +--------------+-----------+-----------+-----------+-----------+
        |   Chromosome |     Start |       End |      Gene |   Cluster |
        |   (category) |   (int32) |   (int32) |   (int64) |   (int32) |
        |--------------+-----------+-----------+-----------+-----------|
        |            1 |         1 |         3 |         1 |         1 |
        |            1 |         2 |         3 |         2 |         1 |
        |            1 |         3 |        10 |         3 |         1 |
        |            1 |         9 |        12 |         3 |         1 |
        +--------------+-----------+-----------+-----------+-----------+
        Unstranded PyRanges object has 4 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.cluster(by="Gene", count=True)
        +--------------+-----------+-----------+-----------+-----------+-----------+
        |   Chromosome |     Start |       End |      Gene |   Cluster |     Count |
        |   (category) |   (int32) |   (int32) |   (int64) |   (int32) |   (int64) |
        |--------------+-----------+-----------+-----------+-----------+-----------|
        |            1 |         1 |         3 |         1 |         1 |         1 |
        |            1 |         2 |         3 |         2 |         2 |         1 |
        |            1 |         3 |        10 |         3 |         3 |         2 |
        |            1 |         9 |        12 |         3 |         3 |         2 |
        +--------------+-----------+-----------+-----------+-----------+-----------+
        Unstranded PyRanges object has 4 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        Avoid clustering bookended intervals with slack=-1:

        >>> gr.cluster(slack=-1)
        +--------------+-----------+-----------+-----------+-----------+
        |   Chromosome |     Start |       End |      Gene |   Cluster |
        |   (category) |   (int32) |   (int32) |   (int64) |   (int32) |
        |--------------+-----------+-----------+-----------+-----------|
        |            1 |         1 |         3 |         1 |         1 |
        |            1 |         2 |         3 |         2 |         1 |
        |            1 |         3 |        10 |         3 |         2 |
        |            1 |         9 |        12 |         3 |         2 |
        +--------------+-----------+-----------+-----------+-----------+
        Unstranded PyRanges object has 4 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr2 = pr.data.ensembl_gtf()[["Feature", "Source"]]
        >>> gr2.cluster(by=["Feature", "Source"])
        +--------------+--------------+---------------+-----------+-----------+--------------+-----------+
        | Chromosome   | Feature      | Source        | Start     | End       | Strand       | Cluster   |
        | (category)   | (category)   | (object)      | (int32)   | (int32)   | (category)   | (int32)   |
        |--------------+--------------+---------------+-----------+-----------+--------------+-----------|
        | 1            | CDS          | ensembl       | 69090     | 70005     | +            | 1         |
        | 1            | CDS          | ensembl       | 925941    | 926013    | +            | 2         |
        | 1            | CDS          | ensembl       | 925941    | 926013    | +            | 2         |
        | 1            | CDS          | ensembl       | 925941    | 926013    | +            | 2         |
        | ...          | ...          | ...           | ...       | ...       | ...          | ...       |
        | 1            | transcript   | havana_tagene | 167128    | 169240    | -            | 1142      |
        | 1            | transcript   | mirbase       | 17368     | 17436     | -            | 1143      |
        | 1            | transcript   | mirbase       | 187890    | 187958    | -            | 1144      |
        | 1            | transcript   | mirbase       | 632324    | 632413    | -            | 1145      |
        +--------------+--------------+---------------+-----------+-----------+--------------+-----------+
        Stranded PyRanges object has 2,446 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        if strand is None:
            strand = self.stranded

        kwargs = {"strand": strand, "slack": slack, "count": count, "by": by}
        kwargs = fill_kwargs(kwargs)

        _stranded = self.stranded
        if not strand and _stranded:
            self.Strand2 = self.Strand
            self = self.unstrand()

        if not by:
            from pyranges.methods.cluster import _cluster
            df = pyrange_apply_single(_cluster, self, **kwargs)
        else:
            from pyranges.methods.cluster import _cluster_by
            kwargs["by"] = by
            df = pyrange_apply_single(_cluster_by, self, **kwargs)

        gr = PyRanges(df)

        # each chromosome got overlapping ids (0 to len). Need to make unique!
        new_dfs = {}
        first = True
        max_id = 0
        for k, v in gr.items():
            if first:
                max_id = v.Cluster.max()
                new_dfs[k] = v
                first = False
                continue

            v.loc[:, "Cluster"] += max_id
            max_id = v.Cluster.max()
            new_dfs[k] = v

        if not strand and _stranded:
            new_dfs = {
                k: d.rename(columns={"Strand2": "Strand"})
                for k, d in new_dfs.items()
            }

        self = PyRanges(new_dfs)

        return self

    def copy(self):

        """Make a deep copy of the PyRanges.

        Notes
        -----

        See the pandas docs for deep-copying caveats."""

        return self.apply(lambda df: df.copy(deep=True))

    @property
    def columns(self):
        """Return the column labels of the PyRanges.

        Returns
        -------
        pandas.Index

        See also
        --------

        PyRanges.chromosomes : return the chromosomes in the PyRanges

        Examples
        --------
        >>> f2 = pr.data.f2()
        >>> f2
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         1 |         2 | a          |         0 | +            |
        | chr1         |         6 |         7 | b          |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> f2.columns
        Index(['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand'], dtype='object')

        >>> f2.columns = f2.columns.str.replace("Sco|re", "NYAN")
        >>> f2
        +--------------+-----------+-----------+------------+------------+--------------+
        | Chromosome   |     Start |       End | Name       |   NYANNYAN | Strand       |
        | (category)   |   (int32) |   (int32) | (object)   |    (int64) | (category)   |
        |--------------+-----------+-----------+------------+------------+--------------|
        | chr1         |         1 |         2 | a          |          0 | +            |
        | chr1         |         6 |         7 | b          |          0 | -            |
        +--------------+-----------+-----------+------------+------------+--------------+
        Stranded PyRanges object has 2 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        if not len(self.values()):
            return []

        first = next(iter(self.values()))
        columns = first.columns

        return columns

    def count_overlaps(self, other, strandedness=None, keep_nonoverlapping=True, overlap_col="NumberOverlaps"):

        """Count number of overlaps per interval.

        Count how many intervals in self overlap with those in other.

        Parameters
        ----------
        strandedness : {"same", "opposite", None, False}, default None, i.e. auto

            Whether to perform the operation on the same, opposite or no strand. Use False to
            ignore the strand. None means use "same" if both PyRanges are stranded, otherwise
            ignore.

        keep_nonoverlapping : bool, default True

            Keep intervals without overlaps.

        overlap_col : str, default "NumberOverlaps"

            Name of column with overlap counts.

        nb_cpu : int, default 1

            How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
            Will only lead to speedups on large datasets.

        Returns
        -------
        PyRanges
            PyRanges with a column of overlaps added.

        See also
        --------

        PyRanges.coverage: find coverage of PyRanges
        pyranges.count_overlaps: count overlaps from multiple PyRanges

        Examples
        --------
        >>> f1 = pr.data.f1().drop()
        >>> f1
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int32) |   (int32) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         3 |         6 | +            |
        | chr1         |         8 |         9 | +            |
        | chr1         |         5 |         7 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        >>> f2 = pr.data.f2().drop()
        >>> f2
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int32) |   (int32) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         1 |         2 | +            |
        | chr1         |         6 |         7 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> f1.count_overlaps(f2, overlap_col="Count")
        +--------------+-----------+-----------+--------------+-----------+
        | Chromosome   |     Start |       End | Strand       |     Count |
        | (category)   |   (int32) |   (int32) | (category)   |   (int64) |
        |--------------+-----------+-----------+--------------+-----------|
        | chr1         |         3 |         6 | +            |         0 |
        | chr1         |         8 |         9 | +            |         0 |
        | chr1         |         5 |         7 | -            |         1 |
        +--------------+-----------+-----------+--------------+-----------+
        Stranded PyRanges object has 3 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        kwargs = {"strandedness": strandedness, "keep_nonoverlapping": keep_nonoverlapping,
                  "overlap_col": overlap_col}
        kwargs = fill_kwargs(kwargs)

        from pyranges.methods.coverage import _number_overlapping
        counts = pyrange_apply(_number_overlapping, self, other, **kwargs)

        return pr.PyRanges(counts)


    def coverage(self, other, strandedness=None, keep_nonoverlapping=True, overlap_col="NumberOverlaps", fraction_col="FractionOverlaps", nb_cpu=1):

        """Count number of overlaps and their fraction per interval.

        Count how many intervals in self overlap with those in other.

        Parameters
        ----------
        strandedness : {"same", "opposite", None, False}, default None, i.e. auto

            Whether to perform the operation on the same, opposite or no strand. Use False to
            ignore the strand. None means use "same" if both PyRanges are stranded, otherwise
            ignore.

        keep_nonoverlapping : bool, default True

            Keep intervals without overlaps.

        overlap_col : str, default "NumberOverlaps"

            Name of column with overlap counts.

        fraction_col : str, default "FractionOverlaps"

            Name of column with fraction of counts.

        nb_cpu: int, default 1

            How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
            Will only lead to speedups on large datasets.

        Returns
        -------
        PyRanges
            PyRanges with a column of overlaps added.

        See also
        --------

        pyranges.count_overlaps: count overlaps from multiple PyRanges

        Examples
        --------
        >>> f1 = pr.from_dict({"Chromosome": [1, 1, 1], "Start": [3, 8, 5],
        ...                    "End": [6,  9, 7]})
        >>> f1
        +--------------+-----------+-----------+
        |   Chromosome |     Start |       End |
        |   (category) |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        |            1 |         3 |         6 |
        |            1 |         8 |         9 |
        |            1 |         5 |         7 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 3 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        >>> f2 = pr.from_dict({"Chromosome": [1, 1], "Start": [1, 6],
        ...                    "End": [2, 7]})
        >>> f2
        +--------------+-----------+-----------+
        |   Chromosome |     Start |       End |
        |   (category) |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        |            1 |         1 |         2 |
        |            1 |         6 |         7 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 2 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> f1.coverage(f2, overlap_col="C", fraction_col="F")
        +--------------+-----------+-----------+-----------+-------------+
        |   Chromosome |     Start |       End |         C |           F |
        |   (category) |   (int32) |   (int32) |   (int64) |   (float64) |
        |--------------+-----------+-----------+-----------+-------------|
        |            1 |         3 |         6 |         0 |         0   |
        |            1 |         8 |         9 |         0 |         0   |
        |            1 |         5 |         7 |         1 |         0.5 |
        +--------------+-----------+-----------+-----------+-------------+
        Unstranded PyRanges object has 3 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        kwargs = {"strandedness": strandedness, "keep_nonoverlapping": keep_nonoverlapping,
                  "overlap_col": overlap_col, "fraction_col": fraction_col, "nb_cpu": nb_cpu}
        kwargs = fill_kwargs(kwargs)

        counts = self.count_overlaps(other, keep_nonoverlapping=True, overlap_col=overlap_col, strandedness=strandedness)

        strand = True if kwargs["strandedness"] else False
        other = other.merge(count=True, strand=strand)

        from pyranges.methods.coverage import _coverage

        counts = pr.PyRanges(pyrange_apply(_coverage, counts, other, **kwargs))

        return counts


    @property
    def df(self):

        """Return PyRanges as DataFrame.

        See also
        --------

        PyRanges.as_df : return PyRanges as DataFrame."""

        return self.as_df()

    def drop(self, drop=None, like=None):
        """Drop column(s).

        If no arguments are given, all the columns except Chromosome, Start, End and Strand are
        dropped.

        Parameters
        ----------

        drop : str or list, default None

            Columns to drop.

        like : str, default None

            Regex-string matching columns to drop. Matches with Chromosome, Start, End or Strand
            are ignored.

        See also
        --------

        PyRanges.unstrand : drop strand information

        Examples
        --------

        >>> gr = pr.from_dict({"Chromosome": [1, 1], "Start": [1, 4], "End": [5, 6],
        ...                    "Strand": ["+", "-"], "Count": [1, 2],
        ...                    "Type": ["exon", "exon"]})
        >>> gr
        +--------------+-----------+-----------+--------------+-----------+------------+
        |   Chromosome |     Start |       End | Strand       |     Count | Type       |
        |   (category) |   (int32) |   (int32) | (category)   |   (int64) | (object)   |
        |--------------+-----------+-----------+--------------+-----------+------------|
        |            1 |         1 |         5 | +            |         1 | exon       |
        |            1 |         4 |         6 | -            |         2 | exon       |
        +--------------+-----------+-----------+--------------+-----------+------------+
        Stranded PyRanges object has 2 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.drop()
        +--------------+-----------+-----------+--------------+
        |   Chromosome |     Start |       End | Strand       |
        |   (category) |   (int32) |   (int32) | (category)   |
        |--------------+-----------+-----------+--------------|
        |            1 |         1 |         5 | +            |
        |            1 |         4 |         6 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        Matches with position-columns are ignored:

        >>> gr.drop(like="Chromosome|Strand")
        +--------------+-----------+-----------+--------------+-----------+------------+
        |   Chromosome |     Start |       End | Strand       |     Count | Type       |
        |   (category) |   (int32) |   (int32) | (category)   |   (int64) | (object)   |
        |--------------+-----------+-----------+--------------+-----------+------------|
        |            1 |         1 |         5 | +            |         1 | exon       |
        |            1 |         4 |         6 | -            |         2 | exon       |
        +--------------+-----------+-----------+--------------+-----------+------------+
        Stranded PyRanges object has 2 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.drop(like="e$")
        +--------------+-----------+-----------+--------------+-----------+
        |   Chromosome |     Start |       End | Strand       |     Count |
        |   (category) |   (int32) |   (int32) | (category)   |   (int64) |
        |--------------+-----------+-----------+--------------+-----------|
        |            1 |         1 |         5 | +            |         1 |
        |            1 |         4 |         6 | -            |         2 |
        +--------------+-----------+-----------+--------------+-----------+
        Stranded PyRanges object has 2 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        from pyranges.methods.drop import _drop
        return _drop(self, drop, like)

    def drop_duplicate_positions(self, strand=None, keep="first"):

        """Return PyRanges with duplicate postion rows removed.

        Parameters
        ----------

        strand : bool, default None, i.e. auto

            Whether to take strand-information into account when considering duplicates.

        keep : {"first", "last", False}

            Whether to keep first, last or drop all duplicates.

        Examples
        --------

        >>> gr = pr.from_string('''Chromosome Start End Strand Name
        ... 1 1 2 + A
        ... 1 1 2 - B
        ... 1 1 2 + Z''')
        >>> gr
        +--------------+-----------+-----------+--------------+------------+
        |   Chromosome |     Start |       End | Strand       | Name       |
        |   (category) |   (int32) |   (int32) | (category)   | (object)   |
        |--------------+-----------+-----------+--------------+------------|
        |            1 |         1 |         2 | +            | A          |
        |            1 |         1 |         2 | +            | Z          |
        |            1 |         1 |         2 | -            | B          |
        +--------------+-----------+-----------+--------------+------------+
        Stranded PyRanges object has 3 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.drop_duplicate_positions()
        +--------------+-----------+-----------+--------------+------------+
        |   Chromosome |     Start |       End | Strand       | Name       |
        |   (category) |   (int32) |   (int32) | (category)   | (object)   |
        |--------------+-----------+-----------+--------------+------------|
        |            1 |         1 |         2 | +            | A          |
        |            1 |         1 |         2 | -            | B          |
        +--------------+-----------+-----------+--------------+------------+
        Stranded PyRanges object has 2 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.drop_duplicate_positions(keep="last")
        +--------------+-----------+-----------+--------------+------------+
        |   Chromosome |     Start |       End | Strand       | Name       |
        |   (category) |   (int32) |   (int32) | (category)   | (object)   |
        |--------------+-----------+-----------+--------------+------------|
        |            1 |         1 |         2 | +            | Z          |
        |            1 |         1 |         2 | -            | B          |
        +--------------+-----------+-----------+--------------+------------+
        Stranded PyRanges object has 2 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        Note that the reverse strand is considered to be behind the forward strand:

        >>> gr.drop_duplicate_positions(keep="last", strand=False)
        +--------------+-----------+-----------+--------------+------------+
        |   Chromosome |     Start |       End | Strand       | Name       |
        |   (category) |   (int32) |   (int32) | (category)   | (object)   |
        |--------------+-----------+-----------+--------------+------------|
        |            1 |         1 |         2 | -            | B          |
        +--------------+-----------+-----------+--------------+------------+
        Stranded PyRanges object has 1 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.drop_duplicate_positions(keep=False, strand=False)
        Empty PyRanges
        """

        from pyranges.methods.drop_duplicates import _drop_duplicate_positions
        if strand is None:
            strand = self.stranded

        kwargs = {}
        kwargs["sparse"] = {"self": False}
        kwargs["keep"] = keep
        kwargs = fill_kwargs(kwargs)
        kwargs["strand"] = strand and self.stranded
        return PyRanges(
            pyrange_apply_single(_drop_duplicate_positions, self, **kwargs))

    @property
    def dtypes(self):
        """Return the dtypes of the PyRanges.

        Examples
        --------

        >>> gr = pr.data.chipseq()
        >>> gr
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   | Start     | End       | Name       | Score     | Strand       |
        | (category)   | (int32)   | (int32)   | (object)   | (int64)   | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         | 212609534 | 212609559 | U0         | 0         | +            |
        | chr1         | 169887529 | 169887554 | U0         | 0         | +            |
        | chr1         | 216711011 | 216711036 | U0         | 0         | +            |
        | chr1         | 144227079 | 144227104 | U0         | 0         | +            |
        | ...          | ...       | ...       | ...        | ...       | ...          |
        | chrY         | 15224235  | 15224260  | U0         | 0         | -            |
        | chrY         | 13517892  | 13517917  | U0         | 0         | -            |
        | chrY         | 8010951   | 8010976   | U0         | 0         | -            |
        | chrY         | 7405376   | 7405401   | U0         | 0         | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 10,000 rows and 6 columns from 24 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.dtypes
        Chromosome    category
        Start            int32
        End              int32
        Name            object
        Score            int64
        Strand        category
        dtype: object
        """

        df = next(iter(self.dfs.values()))

        return df.dtypes

    @property
    def empty(self):

        """Indicate whether PyRanges is empty."""

        return len(self) == 0



    # # TODO: use subtract code here instead, easier
    # def no_overlap(self, other, **kwargs):

    #     kwargs = fill_kwargs(kwargs)
    #     kwargs["invert"] = True
    #     kwargs["sparse"] = {"self": False, "other": True}

    #     # if kwargs["strandedness"] in ["same", "opposite"]:
    #     #     kwargs["strandedness"] = {
    #     #         "same": "opposite",
    #     #         "opposite": "same"
    #     #     }[kwargs["strandedness"]]
    #     dfs = pyrange_apply(_overlap, self, other, **kwargs)

    #     return PyRanges(dfs)


    # @profile



    def five_end(self):

        """Return the five prime end of intervals.

        The five prime end is the start of a forward strand or the end of a reverse strand.

        Returns
        -------
        PyRanges

            PyRanges with the five prime ends

        Notes
        -----

        Requires the PyRanges to be stranded.

        See Also
        --------

        PyRanges.three_end : return the 3' end

        Examples
        --------

        >>> gr = pr.from_dict({'Chromosome': ['chr1', 'chr1'], 'Start': [3, 5], 'End': [9, 7],
        ...                    'Strand': ["+", "-"]})
        >>> gr
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int32) |   (int32) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         3 |         9 | +            |
        | chr1         |         5 |         7 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.five_end()
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int32) |   (int32) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         3 |         4 | +            |
        | chr1         |         7 |         8 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        assert self.stranded, "Need stranded pyrange to find 5'."
        kwargs = fill_kwargs({"strand": self.stranded})
        return PyRanges(
            pyrange_apply_single(_tss, self, **kwargs))

    def head(self, n=8):

        """Return the n first rows.

        Parameters
        ----------

        n : int, default 8

            Return n rows.

        Returns
        -------
        PyRanges

            PyRanges with the n first rows.

        See Also
        --------

        PyRanges.tail : return the last rows
        PyRanges.sample : return random rows

        Examples
        --------

        >>> gr = pr.data.chipseq()
        >>> gr
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   | Start     | End       | Name       | Score     | Strand       |
        | (category)   | (int32)   | (int32)   | (object)   | (int64)   | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         | 212609534 | 212609559 | U0         | 0         | +            |
        | chr1         | 169887529 | 169887554 | U0         | 0         | +            |
        | chr1         | 216711011 | 216711036 | U0         | 0         | +            |
        | chr1         | 144227079 | 144227104 | U0         | 0         | +            |
        | ...          | ...       | ...       | ...        | ...       | ...          |
        | chrY         | 15224235  | 15224260  | U0         | 0         | -            |
        | chrY         | 13517892  | 13517917  | U0         | 0         | -            |
        | chrY         | 8010951   | 8010976   | U0         | 0         | -            |
        | chrY         | 7405376   | 7405401   | U0         | 0         | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 10,000 rows and 6 columns from 24 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.head(3)
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         | 212609534 | 212609559 | U0         |         0 | +            |
        | chr1         | 169887529 | 169887554 | U0         |         0 | +            |
        | chr1         | 216711011 | 216711036 | U0         |         0 | +            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        subsetter = np.zeros(len(self), dtype=np.bool)
        subsetter[:n] = True
        return self[subsetter]


    def insert(self, other, loc=None):

        """Add one or more columns to the PyRanges.

        Parameters
        ----------
        other : Series, DataFrame or dict
            Data to insert into the PyRanges. `other` must have the same number of rows as the PyRanges.

        loc : int, default None, i.e. after last column of PyRanges.
            Insertion index.

        Returns
        -------
        PyRanges
            A copy of the PyRanges with the column(s) inserted starting at `loc`.

        Note
        ----

        If a Series, or a dict of Series is used, the Series must have a name.

        Examples
        --------

        >>> gr = pr.from_dict({"Chromosome": ["L", "E", "E", "T"], "Start": [1, 1, 2, 3], "End": [5, 8, 13, 21]})
        >>> gr
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        | E            |         1 |         8 |
        | E            |         2 |        13 |
        | L            |         1 |         5 |
        | T            |         3 |        21 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 4 rows and 3 columns from 3 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> s = pd.Series(data = [1, 3, 3, 7], name="Column")
        >>> gr.insert(s)
        +--------------+-----------+-----------+-----------+
        | Chromosome   |     Start |       End |    Column |
        | (category)   |   (int32) |   (int32) |   (int64) |
        |--------------+-----------+-----------+-----------|
        | E            |         1 |         8 |         1 |
        | E            |         2 |        13 |         3 |
        | L            |         1 |         5 |         3 |
        | T            |         3 |        21 |         7 |
        +--------------+-----------+-----------+-----------+
        Unstranded PyRanges object has 4 rows and 4 columns from 3 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> df = pd.DataFrame({"NY": s, "AN": s})
        >>> df
           NY  AN
        0   1   1
        1   3   3
        2   3   3
        3   7   7

        Note that the original PyRanges was not affected by previously inserting Column:

        >>> gr.insert(df, 1)
        +--------------+-----------+-----------+-----------+-----------+
        | Chromosome   |        NY |        AN |     Start |       End |
        | (category)   |   (int64) |   (int64) |   (int32) |   (int32) |
        |--------------+-----------+-----------+-----------+-----------|
        | E            |         1 |         1 |         1 |         8 |
        | E            |         3 |         3 |         2 |        13 |
        | L            |         3 |         3 |         1 |         5 |
        | T            |         7 |         7 |         3 |        21 |
        +--------------+-----------+-----------+-----------+-----------+
        Unstranded PyRanges object has 4 rows and 5 columns from 3 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> arbitrary_result = gr.apply(
        ... lambda df: pd.Series(df.Start + df.End, name="Hi!"), as_pyranges=False)
        >>> arbitrary_result
        {'E': 1     9
        2    15
        Name: Hi!, dtype: int32, 'L': 0    6
        Name: Hi!, dtype: int32, 'T': 3    24
        Name: Hi!, dtype: int32}

        >>> gr.insert(arbitrary_result)
        +--------------+-----------+-----------+-----------+
        | Chromosome   |     Start |       End |       Hi! |
        | (category)   |   (int32) |   (int32) |   (int32) |
        |--------------+-----------+-----------+-----------|
        | E            |         1 |         8 |         9 |
        | E            |         2 |        13 |        15 |
        | L            |         1 |         5 |         6 |
        | T            |         3 |        21 |        24 |
        +--------------+-----------+-----------+-----------+
        Unstranded PyRanges object has 4 rows and 4 columns from 3 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        if loc is None:
            loc = len(self.columns)

        self = self.copy()

        from pyranges.methods.attr import _setattr

        if isinstance(other, (pd.Series, pd.DataFrame)):
            assert len(other) == len(self), "Pandas Series or DataFrame must be same length as PyRanges!"

            if isinstance(other, pd.Series):
                if not other.name:
                    raise Exception("Series must have a name!")

                _setattr(self, other.name, other, loc)

            if isinstance(other, pd.DataFrame):
                for c in other:
                    _setattr(self, c, other[c], loc)
                    loc += 1

        elif isinstance(other, dict) and other:

            first = next(iter(other.values()))
            is_dataframe = isinstance(first, pd.DataFrame)
            if is_dataframe:
                columns = first.columns

                ds = []
                for c in columns:
                    ds.append({k: v[c] for k, v in other.items()})

                for c, d in zip(columns, ds):
                    _setattr(self, str(c), d, loc)
                    loc += 1
            else:
                if not first.name:
                    raise Exception("Series must have a name!")

                d = {k: v for k, v in other.items()}
                _setattr(self, first.name, d, loc)

        return self


    def intersect(self, other, strandedness=None, how=None, invert=False, nb_cpu=1):

        """Return overlapping subintervals.

        Returns the segments of the intervals in self which overlap with those in other.

        Parameters
        ----------
        other : PyRanges

            PyRanges to intersect.

        strandedness : {None, "same", "opposite", False}, default None, i.e. auto

            Whether to compare PyRanges on the same strand, the opposite or ignore strand
            information. The default, None, means use "same" if both PyRanges are strande,
            otherwise ignore the strand information.

        how : {None, "first", "last", "containment"}, default None, i.e. all

            What intervals to report. By default reports all overlapping intervals. "containment"
            reports intervals where the overlapping is contained within it.

        invert : bool, default False

            Whether to return the intervals without overlaps.

        nb_cpu: int, default 1

            How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
            Will only lead to speedups on large datasets.

        Returns
        -------
        PyRanges

            A PyRanges with overlapping subintervals.

        See also
        --------

        PyRanges.set_intersect : set-intersect PyRanges
        PyRanges.overlap : report overlapping intervals

        Examples
        --------

        >>> gr = pr.from_dict({"Chromosome": ["chr1"] * 3, "Start": [1, 4, 10],
        ...                    "End": [3, 9, 11], "ID": ["a", "b", "c"]})
        >>> gr
        +--------------+-----------+-----------+------------+
        |   Chromosome |     Start |       End | ID         |
        |   (category) |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------|
        |         chr1 |         1 |         3 | a          |
        |         chr1 |         4 |         9 | b          |
        |         chr1 |        10 |        11 | c          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr2 = pr.from_dict({"Chromosome": ["chr1"] * 3, "Start": [2, 2, 9], "End": [3, 9, 10]})
        >>> gr2
        +--------------+-----------+-----------+
        |   Chromosome |     Start |       End |
        |   (category) |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        |         chr1 |         2 |         3 |
        |         chr1 |         2 |         9 |
        |         chr1 |         9 |        10 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 3 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.intersect(gr2)
        +--------------+-----------+-----------+------------+
        |   Chromosome |     Start |       End | ID         |
        |   (category) |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------|
        |         chr1 |         2 |         3 | a          |
        |         chr1 |         2 |         3 | a          |
        |         chr1 |         4 |         9 | b          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.intersect(gr2, how="first")
        +--------------+-----------+-----------+------------+
        |   Chromosome |     Start |       End | ID         |
        |   (category) |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------|
        |         chr1 |         2 |         3 | a          |
        |         chr1 |         4 |         9 | b          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.intersect(gr2, how="containment")
        +--------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | ID         |
        | (category)   |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------|
        | chr1         |         4 |         9 | b          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 1 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        kwargs = {"how": how, "strandedness": strandedness, "nb_cpu": nb_cpu}
        kwargs = fill_kwargs(kwargs)
        kwargs["sparse"] = {"self": False, "other": True}

        if len(self) == 0:
            return self

        if invert:
            self.__ix__ = np.arange(len(self))

        dfs = pyrange_apply(_intersection, self, other, **kwargs)
        result = pr.PyRanges(dfs)

        if invert:
            found_idxs = getattr(result, "__ix__", [])
            result = self[~self.__ix__.isin(found_idxs)]
            result = result.drop("__ix__")

        return result

    def items(self):

        """Return the pairs of keys and DataFrames.

        Returns
        -------
        dict

            The dict mapping keys to DataFrames in the PyRanges.

        See Also
        --------

        PyRanges.chromosomes : return the chromosomes
        PyRanges.keys : return the keys
        PyRanges.values : return the DataFrames in the PyRanges

        Examples
        --------

        >>> gr = pr.data.f1()
        >>> gr.items()
        [(('chr1', '+'),   Chromosome  Start  End       Name  Score Strand
        0       chr1      3    6  interval1      0      +
        2       chr1      8    9  interval3      0      +), (('chr1', '-'),   Chromosome  Start  End       Name  Score Strand
        1       chr1      5    7  interval2      0      -)]
        """

        return natsorted([(k, df) for (k, df) in self.dfs.items()])

    def join(self, other, strandedness=None, how=None, report_overlap=False, slack=0, suffix="_b", nb_cpu=1, apply_strand_suffix=None):

        """Join PyRanges on genomic location.

        Parameters
        ----------
        other : PyRanges

            PyRanges to join.

        strandedness : {None, "same", "opposite", False}, default None, i.e. auto

            Whether to compare PyRanges on the same strand, the opposite or ignore strand
            information. The default, None, means use "same" if both PyRanges are strande,
            otherwise ignore the strand information.

        how : {None, "left", "right"}, default None, i.e. "inner"

            How to handle intervals without overlap. None means only keep overlapping intervals.
            "left" keeps all intervals in self, "right" keeps all intervals in other.

        report_overlap : bool, default False

            Report amount of overlap in base pairs.

        slack : int, default 0

            Lengthen intervals in self before joining.

        suffix : str or tuple, default "_b"

            Suffix to give overlapping columns in other.

        apply_strand_suffix : bool, default None

            If first pyranges is unstranded, but the second is not, the first will be given a strand column.
            apply_strand_suffix makes the added strand column a regular data column instead by adding a suffix.

        nb_cpu: int, default 1

            How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
            Will only lead to speedups on large datasets.

        Returns
        -------
        PyRanges

            A PyRanges appended with columns of another.

        Notes
        -----

        The chromosome from other will never be reported as it is always the same as in self.

        As pandas did not have NaN for non-float datatypes until recently, "left" and "right" join
        give non-overlapping rows the value -1 to avoid promoting columns to object. This will
        change to NaN in a future version as general NaN becomes stable in pandas.

        See also
        --------

        PyRanges.new_position : give joined PyRanges new coordinates

        Examples
        --------

        >>> f1 = pr.from_dict({'Chromosome': ['chr1', 'chr1', 'chr1'], 'Start': [3, 8, 5],
        ...                    'End': [6, 9, 7], 'Name': ['interval1', 'interval3', 'interval2']})
        >>> f1
        +--------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | Name       |
        | (category)   |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------|
        | chr1         |         3 |         6 | interval1  |
        | chr1         |         8 |         9 | interval3  |
        | chr1         |         5 |         7 | interval2  |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> f2 = pr.from_dict({'Chromosome': ['chr1', 'chr1'], 'Start': [1, 6],
        ...                    'End': [2, 7], 'Name': ['a', 'b']})
        >>> f2
        +--------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | Name       |
        | (category)   |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------|
        | chr1         |         1 |         2 | a          |
        | chr1         |         6 |         7 | b          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> f1.join(f2)
        +--------------+-----------+-----------+------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | Name       |   Start_b |     End_b | Name_b     |
        | (category)   |   (int32) |   (int32) | (object)   |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------+-----------+-----------+------------|
        | chr1         |         5 |         7 | interval2  |         6 |         7 | b          |
        +--------------+-----------+-----------+------------+-----------+-----------+------------+
        Unstranded PyRanges object has 1 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> f1.join(f2, how="right")
        +--------------+-----------+-----------+------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | Name       |   Start_b |     End_b | Name_b     |
        | (category)   |   (int32) |   (int32) | (object)   |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------+-----------+-----------+------------|
        | chr1         |         5 |         7 | interval2  |         6 |         7 | b          |
        | chr1         |        -1 |        -1 | -1         |         1 |         2 | a          |
        +--------------+-----------+-----------+------------+-----------+-----------+------------+
        Unstranded PyRanges object has 2 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        With slack 1, bookended features are joined (see row 1):

        >>> f1.join(f2, slack=1)
        +--------------+-----------+-----------+------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | Name       |   Start_b |     End_b | Name_b     |
        | (category)   |   (int32) |   (int32) | (object)   |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------+-----------+-----------+------------|
        | chr1         |         3 |         6 | interval1  |         6 |         7 | b          |
        | chr1         |         5 |         7 | interval2  |         6 |         7 | b          |
        +--------------+-----------+-----------+------------+-----------+-----------+------------+
        Unstranded PyRanges object has 2 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        from pyranges.methods.join import _write_both

        kwargs = {"strandedness": strandedness, "how": how, "report_overlap":report_overlap, "suffix": suffix, "nb_cpu": nb_cpu, "apply_strand_suffix": apply_strand_suffix}
        if slack:
            self = self.copy()
            self.Start__slack = self.Start
            self.End__slack = self.End

            self = self.slack(slack)

        if "suffix" in kwargs and isinstance(kwargs["suffix"], str):
            suffixes = "", kwargs["suffix"]
            kwargs["suffixes"] = suffixes

        kwargs = fill_kwargs(kwargs)

        how = kwargs.get("how")

        if how in ["left", "outer"]:
            kwargs["example_header_other"] = other.head(1).df
        if how in ["right", "outer"]:
            kwargs["example_header_self"] = self.head(1).df

        dfs = pyrange_apply(_write_both, self, other, **kwargs)
        gr = PyRanges(dfs)

        if slack and len(gr) > 0:
            gr.Start = gr.Start__slack
            gr.End = gr.End__slack
            gr = gr.drop(like="(Start|End).*__slack")

        if not self.stranded and other.stranded:
            if apply_strand_suffix is None:
                import sys
                print("join: Strand data from other will be added as strand data to self.\nIf this is undesired use the flag apply_strand_suffix=False.\nTo turn off the warning set apply_strand_suffix to True or False.", file=sys.stderr)
            elif apply_strand_suffix:
                gr.columns = gr.columns.str.replace("Strand", "Strand" + kwargs["suffix"])

        return gr

    def keys(self):

        """Return the keys.

        Returns
        -------

        Returns the keys (chromosomes or chromosome/strand pairs) as strings or tuples of strings
        in natsorted order.

        See Also
        --------

        PyRanges.chromosomes : return the chromosomes

        Examples
        --------

        >>> gr = pr.data.chipseq()
        >>> gr.keys()
        [('chr1', '+'), ('chr1', '-'), ('chr2', '+'), ('chr2', '-'), ('chr3', '+'), ('chr3', '-'), ('chr4', '+'), ('chr4', '-'), ('chr5', '+'), ('chr5', '-'), ('chr6', '+'), ('chr6', '-'), ('chr7', '+'), ('chr7', '-'), ('chr8', '+'), ('chr8', '-'), ('chr9', '+'), ('chr9', '-'), ('chr10', '+'), ('chr10', '-'), ('chr11', '+'), ('chr11', '-'), ('chr12', '+'), ('chr12', '-'), ('chr13', '+'), ('chr13', '-'), ('chr14', '+'), ('chr14', '-'), ('chr15', '+'), ('chr15', '-'), ('chr16', '+'), ('chr16', '-'), ('chr17', '+'), ('chr17', '-'), ('chr18', '+'), ('chr18', '-'), ('chr19', '+'), ('chr19', '-'), ('chr20', '+'), ('chr20', '-'), ('chr21', '+'), ('chr21', '-'), ('chr22', '+'), ('chr22', '-'), ('chrX', '+'), ('chrX', '-'), ('chrY', '+'), ('chrY', '-')]
        >>> gr.unstrand().keys()
        ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

        """

        return natsorted(self.dfs.keys())

    def k_nearest(self, other, k=1, ties=None, strandedness=None, overlap=True, how=None, suffix="_b", nb_cpu=1, apply_strand_suffix=None):

        """Find k nearest intervals.

        Parameters
        ----------
        other : PyRanges

            PyRanges to find nearest interval in.

        k : int or list/array/Series of int

            Number of closest to return. If iterable, must be same length as PyRanges.

        ties : {None, "first", "last", "different"}, default None

            How to resolve ties, i.e. closest intervals with equal distance. None means that the k nearest intervals are kept.
            "first" means that the first tie is kept, "last" meanst that the last is kept.
            "different" means that all nearest intervals with the k unique nearest distances are kept.

        strandedness : {None, "same", "opposite", False}, default None, i.e. auto

            Whether to compare PyRanges on the same strand, the opposite or ignore strand
            information. The default, None, means use "same" if both PyRanges are stranded,
            otherwise ignore the strand information.

        overlap : bool, default True

            Whether to include overlaps.

        how : {None, "upstream", "downstream"}, default None, i.e. both directions

            Whether to only look for nearest in one direction. Always with respect to the PyRanges
            it is called on.

        suffix : str, default "_b"

            Suffix to give columns with shared name in other.

        apply_strand_suffix : bool, default None

            If first pyranges is unstranded, but the second is not, the first will be given a strand column.
            apply_strand_suffix makes the added strand column a regular data column instead by adding a suffix.


        nb_cpu: int, default 1

            How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
            Will only lead to speedups on large datasets.

        Returns
        -------
        PyRanges

            A PyRanges with columns of nearest interval horizontally appended.

        Notes
        -----

        nearest also exists, and is more performant.

        See also
        --------

        PyRanges.new_position : give joined PyRanges new coordinates
        PyRanges.nearest : find nearest intervals

        Examples
        --------

        >>> f1 = pr.from_dict({'Chromosome': ['chr1', 'chr1', 'chr1'], 'Start': [3, 8, 5],
        ...                    'End': [6, 9, 7], 'Strand': ['+', '+', '-']})
        >>> f1
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int32) |   (int32) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         3 |         6 | +            |
        | chr1         |         8 |         9 | +            |
        | chr1         |         5 |         7 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> f2 = pr.from_dict({'Chromosome': ['chr1', 'chr1'], 'Start': [1, 6],
        ...                    'End': [2, 7], 'Strand': ['+', '-']})
        >>> f2
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int32) |   (int32) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         1 |         2 | +            |
        | chr1         |         6 |         7 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> f1.k_nearest(f2, k=2)
        +--------------+-----------+-----------+--------------+-----------+-----------+--------------+------------+
        | Chromosome   |     Start |       End | Strand       |   Start_b |     End_b | Strand_b     |   Distance |
        | (category)   |   (int32) |   (int32) | (category)   |   (int32) |   (int32) | (category)   |    (int32) |
        |--------------+-----------+-----------+--------------+-----------+-----------+--------------+------------|
        | chr1         |         3 |         6 | +            |         6 |         7 | -            |          1 |
        | chr1         |         3 |         6 | +            |         1 |         2 | +            |         -2 |
        | chr1         |         8 |         9 | +            |         6 |         7 | -            |         -2 |
        | chr1         |         8 |         9 | +            |         1 |         2 | +            |         -7 |
        | chr1         |         5 |         7 | -            |         6 |         7 | -            |          0 |
        | chr1         |         5 |         7 | -            |         1 |         2 | +            |          4 |
        +--------------+-----------+-----------+--------------+-----------+-----------+--------------+------------+
        Stranded PyRanges object has 6 rows and 8 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> f1.k_nearest(f2, how="upstream", k=2)
        +--------------+-----------+-----------+--------------+-----------+-----------+--------------+------------+
        | Chromosome   |     Start |       End | Strand       |   Start_b |     End_b | Strand_b     |   Distance |
        | (category)   |   (int32) |   (int32) | (category)   |   (int32) |   (int32) | (category)   |    (int32) |
        |--------------+-----------+-----------+--------------+-----------+-----------+--------------+------------|
        | chr1         |         3 |         6 | +            |         1 |         2 | +            |         -2 |
        | chr1         |         8 |         9 | +            |         6 |         7 | -            |         -2 |
        | chr1         |         8 |         9 | +            |         1 |         2 | +            |         -7 |
        | chr1         |         5 |         7 | -            |         6 |         7 | -            |          0 |
        +--------------+-----------+-----------+--------------+-----------+-----------+--------------+------------+
        Stranded PyRanges object has 4 rows and 8 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> f1.k_nearest(f2, k=[1, 2, 1])
        +--------------+-----------+-----------+--------------+-----------+-----------+--------------+------------+
        | Chromosome   |     Start |       End | Strand       |   Start_b |     End_b | Strand_b     |   Distance |
        | (category)   |   (int32) |   (int32) | (category)   |   (int32) |   (int32) | (category)   |    (int32) |
        |--------------+-----------+-----------+--------------+-----------+-----------+--------------+------------|
        | chr1         |         3 |         6 | +            |         6 |         7 | -            |          1 |
        | chr1         |         8 |         9 | +            |         6 |         7 | -            |         -2 |
        | chr1         |         8 |         9 | +            |         1 |         2 | +            |         -7 |
        | chr1         |         5 |         7 | -            |         6 |         7 | -            |          0 |
        +--------------+-----------+-----------+--------------+-----------+-----------+--------------+------------+
        Stranded PyRanges object has 4 rows and 8 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> d1 = {"Chromosome": [1], "Start": [5], "End": [6]}
        >>> d2 = {"Chromosome": 1, "Start": [1] * 2 + [5] * 2 + [9] * 2,
        ...       "End": [3] * 2 + [7] * 2 + [11] * 2, "ID": range(6)}
        >>> gr, gr2 = pr.from_dict(d1), pr.from_dict(d2)

        >>> gr
        +--------------+-----------+-----------+
        |   Chromosome |     Start |       End |
        |   (category) |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        |            1 |         5 |         6 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 1 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr2
        +--------------+-----------+-----------+-----------+
        |   Chromosome |     Start |       End |        ID |
        |   (category) |   (int32) |   (int32) |   (int64) |
        |--------------+-----------+-----------+-----------|
        |            1 |         1 |         3 |         0 |
        |            1 |         1 |         3 |         1 |
        |            1 |         5 |         7 |         2 |
        |            1 |         5 |         7 |         3 |
        |            1 |         9 |        11 |         4 |
        |            1 |         9 |        11 |         5 |
        +--------------+-----------+-----------+-----------+
        Unstranded PyRanges object has 6 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.k_nearest(gr2, k=2)
        +--------------+-----------+-----------+-----------+-----------+-----------+------------+
        |   Chromosome |     Start |       End |   Start_b |     End_b |        ID |   Distance |
        |   (category) |   (int32) |   (int32) |   (int32) |   (int32) |   (int64) |    (int64) |
        |--------------+-----------+-----------+-----------+-----------+-----------+------------|
        |            1 |         5 |         6 |         5 |         7 |         2 |          0 |
        |            1 |         5 |         6 |         5 |         7 |         3 |          0 |
        +--------------+-----------+-----------+-----------+-----------+-----------+------------+
        Unstranded PyRanges object has 2 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.k_nearest(gr2, k=2, ties="different")
        +--------------+-----------+-----------+-----------+-----------+-----------+------------+
        |   Chromosome |     Start |       End |   Start_b |     End_b |        ID |   Distance |
        |   (category) |   (int32) |   (int32) |   (int32) |   (int32) |   (int64) |    (int64) |
        |--------------+-----------+-----------+-----------+-----------+-----------+------------|
        |            1 |         5 |         6 |         5 |         7 |         2 |          0 |
        |            1 |         5 |         6 |         5 |         7 |         3 |          0 |
        |            1 |         5 |         6 |         1 |         3 |         1 |         -3 |
        |            1 |         5 |         6 |         1 |         3 |         0 |         -3 |
        +--------------+-----------+-----------+-----------+-----------+-----------+------------+
        Unstranded PyRanges object has 4 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.k_nearest(gr2, k=3, ties="first")
        +--------------+-----------+-----------+-----------+-----------+-----------+------------+
        |   Chromosome |     Start |       End |   Start_b |     End_b |        ID |   Distance |
        |   (category) |   (int32) |   (int32) |   (int32) |   (int32) |   (int64) |    (int64) |
        |--------------+-----------+-----------+-----------+-----------+-----------+------------|
        |            1 |         5 |         6 |         5 |         7 |         2 |          0 |
        |            1 |         5 |         6 |         1 |         3 |         1 |         -3 |
        |            1 |         5 |         6 |         9 |        11 |         4 |          4 |
        +--------------+-----------+-----------+-----------+-----------+-----------+------------+
        Unstranded PyRanges object has 3 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.k_nearest(gr2, k=1, overlap=False)
        +--------------+-----------+-----------+-----------+-----------+-----------+------------+
        |   Chromosome |     Start |       End |   Start_b |     End_b |        ID |   Distance |
        |   (category) |   (int32) |   (int32) |   (int32) |   (int32) |   (int64) |    (int32) |
        |--------------+-----------+-----------+-----------+-----------+-----------+------------|
        |            1 |         5 |         6 |         1 |         3 |         1 |         -3 |
        |            1 |         5 |         6 |         1 |         3 |         0 |         -3 |
        +--------------+-----------+-----------+-----------+-----------+-----------+------------+
        Unstranded PyRanges object has 2 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        from pyranges.methods.k_nearest import _nearest
        from sorted_nearest import get_all_ties, get_different_ties

        kwargs = {"strandedness": strandedness, "how": how, "overlap": overlap, "nb_cpu": nb_cpu,
                  "k": k, "ties": ties}
        kwargs = fill_kwargs(kwargs)
        kwargs["stranded"] = self.stranded and other.stranded

        overlap = kwargs.get("overlap", True)
        ties = kwargs.get("ties", False)

        self = self.copy()

        try: # if k is a Series
            k = k.values
        except:
            pass

        # how many to nearest to find; might be different for each
        self.__k__ = k
        # give each their own unique ID
        self.__IX__ = np.arange(len(self))

        dfs = pyrange_apply(_nearest, self, other, **kwargs)
        nearest = PyRanges(dfs)

        if not overlap:
            result = nearest
        else:
            from collections import defaultdict
            overlap_how = defaultdict(lambda: None, {"first": "first", "last": "last"})[kwargs.get("ties")]
            overlaps = self.join(other, strandedness=strandedness, how=overlap_how, nb_cpu=nb_cpu, apply_strand_suffix=apply_strand_suffix)
            overlaps.Distance = 0
            result = pr.concat([overlaps, nearest])

        if not len(result):
            return pr.PyRanges()
        new_result = {}
        if ties in ["first", "last"]:
            for c, df in result:
                df = df.sort_values(["__IX__", "Distance"])
                grpby = df.groupby("__k__", sort=False)
                dfs = []
                for k, kdf in grpby:
                    grpby2 = kdf.groupby("__IX__", sort=False)
                    _df = grpby2.head(k)
                    dfs.append(_df)

                if dfs:
                    new_result[c] = pd.concat(dfs)

        elif ties == "different" or not ties:
            for c, df in result:

                if df.empty:
                    continue
                dfs = []

                df = df.sort_values(["__IX__", "Distance"])
                grpby = df.groupby("__k__", sort=False)

                for k, kdf in grpby:
                    if ties:
                        lx = get_different_ties(kdf.index.values, kdf.__IX__.values, kdf.Distance.astype(np.int64).values, k)
                        _df = kdf.reindex(lx)
                    else:
                        lx = get_all_ties(kdf.index.values, kdf.__IX__.values, kdf.Distance.astype(np.int64).values, k)
                        _df = kdf.reindex(lx)
                        _df = _df.groupby("__IX__").head(k)
                    dfs.append(_df)

                if dfs:
                    new_result[c] = pd.concat(dfs)

        result = pr.PyRanges(new_result)

        if not result.__IX__.is_monotonic:
            result = result.sort("__IX__")

        result = result.drop(like="__IX__|__k__")

        self = self.drop(like="__k__|__IX__")

        def prev_to_neg(df, **kwargs):

            strand = df.Strand.iloc[0] if "Strand" in df else "+"

            suffix = kwargs["suffix"]

            bools = df["End" + suffix] < df.Start
            if not strand == "+":
                bools = ~bools

            df.loc[bools, "Distance"] = -df.loc[bools, "Distance"]
            return df

        result = result.apply(prev_to_neg, suffix=kwargs["suffix"])

        if not self.stranded and other.stranded:
            if apply_strand_suffix is None:
                import sys
                print("join: Strand data from other will be added as strand data to self.\nIf this is undesired use the flag apply_strand_suffix=False.\nTo turn off the warning set apply_strand_suffix to True or False.", file=sys.stderr)
            elif apply_strand_suffix:
                result.columns = result.columns.str.replace("Strand", "Strand" + kwargs["suffix"])

        return result



    @property
    def length(self):

        """Return the total length of the intervals.

        See Also
        --------

        PyRanges.lengths : return the intervals lengths

        Examples
        --------

        >>> gr = pr.data.f1()
        >>> gr
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         3 |         6 | interval1  |         0 | +            |
        | chr1         |         8 |         9 | interval3  |         0 | +            |
        | chr1         |         5 |         7 | interval2  |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.length
        6

        To find the length of the genome covered by the intervals, use merge first:

        >>> gr.merge(strand=False).length
        5
        """

        return int(self.lengths(as_dict=False).sum())


    def lengths(self, as_dict=False):

        """Return the length of each interval.

        Parameters
        ----------

        as_dict : bool, default False

            Whether to return lengths as Series or dict of Series per key.

        Returns
        -------
        Series or dict of Series with the lengths of each interval.

        See Also
        --------

        PyRanges.lengths : return the intervals lengths

        Examples
        --------

        >>> gr = pr.data.f1()
        >>> gr
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         3 |         6 | interval1  |         0 | +            |
        | chr1         |         8 |         9 | interval3  |         0 | +            |
        | chr1         |         5 |         7 | interval2  |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.lengths()
        0    3
        1    1
        2    2
        dtype: int32

        To find the length of the genome covered by the intervals, use merge first:

        >>> gr.Length = gr.lengths()
        >>> gr
        +--------------+-----------+-----------+------------+-----------+--------------+-----------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |    Length |
        | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |   (int32) |
        |--------------+-----------+-----------+------------+-----------+--------------+-----------|
        | chr1         |         3 |         6 | interval1  |         0 | +            |         3 |
        | chr1         |         8 |         9 | interval3  |         0 | +            |         1 |
        | chr1         |         5 |         7 | interval2  |         0 | -            |         2 |
        +--------------+-----------+-----------+------------+-----------+--------------+-----------+
        Stranded PyRanges object has 3 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        if as_dict:
            if not len(self):
                return {}
            lengths = {}
            for k, df in self.items():
                lengths[k] = df.End - df.Start

            return lengths
        else:
            _lengths = []
            if not len(self):
                return np.array(_lengths, dtype=int)
            for _, df in self:
                lengths = df.End - df.Start
                _lengths.append(lengths)

            return pd.concat(_lengths).reset_index(drop=True)

    def max_disjoint(self, strand=None, slack=0, **kwargs):

        """Find the maximal disjoint set of intervals.

        Parameters
        ----------
        strand : bool, default None, i.e. auto

            Find the max disjoint set separately for each strand.

        slack : int, default 0

            Consider intervals within a distance of slack to be overlapping.

        Returns
        -------
        PyRanges

            PyRanges with maximal disjoint set of intervals.

        Examples
        --------
        >>> gr = pr.data.f1()
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         3 |         6 | interval1  |         0 | +            |
        | chr1         |         8 |         9 | interval3  |         0 | +            |
        | chr1         |         5 |         7 | interval2  |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.max_disjoint(strand=False)
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         3 |         6 | interval1  |         0 | +            |
        | chr1         |         8 |         9 | interval3  |         0 | +            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        if strand is None:
            strand = self.stranded

        kwargs = {"strand": strand, "slack": slack}
        kwargs = fill_kwargs(kwargs)

        from pyranges.methods.max_disjoint import _max_disjoint
        df = pyrange_apply_single(_max_disjoint, self, **kwargs)

        return pr.PyRanges(df)

    def merge(self, strand=None, count=False, count_col="Count", by=None, slack=0):

        """Merge overlapping intervals into one.

        Parameters
        ----------
        strand : bool, default None, i.e. auto

            Only merge intervals on same strand.

        count : bool, default False

            Count intervals in each superinterval.

        count_col : str, default "Count"

            Name of column with counts.

        by : str or list of str, default None

            Only merge intervals with equal values in these columns.

        slack : int, default 0

            Allow this many nucleotides between each interval to merge.

        Returns
        -------
        PyRanges

            PyRanges with superintervals.

        Notes
        -----

        To avoid losing metadata, use cluster instead. If you want to perform a reduction function
        on the metadata, use pandas groupby.

        See Also
        --------

        PyRanges.cluster : annotate overlapping intervals with common ID

        Examples
        --------

        >>> gr = pr.data.ensembl_gtf()[["Feature", "gene_name"]]
        >>> gr
        +--------------+--------------+-----------+-----------+--------------+-------------+
        | Chromosome   | Feature      | Start     | End       | Strand       | gene_name   |
        | (category)   | (category)   | (int32)   | (int32)   | (category)   | (object)    |
        |--------------+--------------+-----------+-----------+--------------+-------------|
        | 1            | gene         | 11868     | 14409     | +            | DDX11L1     |
        | 1            | transcript   | 11868     | 14409     | +            | DDX11L1     |
        | 1            | exon         | 11868     | 12227     | +            | DDX11L1     |
        | 1            | exon         | 12612     | 12721     | +            | DDX11L1     |
        | ...          | ...          | ...       | ...       | ...          | ...         |
        | 1            | gene         | 1173055   | 1179555   | -            | TTLL10-AS1  |
        | 1            | transcript   | 1173055   | 1179555   | -            | TTLL10-AS1  |
        | 1            | exon         | 1179364   | 1179555   | -            | TTLL10-AS1  |
        | 1            | exon         | 1173055   | 1176396   | -            | TTLL10-AS1  |
        +--------------+--------------+-----------+-----------+--------------+-------------+
        Stranded PyRanges object has 2,446 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.merge(count=True, count_col="Count")
        +--------------+-----------+-----------+--------------+-----------+
        | Chromosome   | Start     | End       | Strand       | Count     |
        | (category)   | (int32)   | (int32)   | (category)   | (int32)   |
        |--------------+-----------+-----------+--------------+-----------|
        | 1            | 11868     | 14409     | +            | 12        |
        | 1            | 29553     | 31109     | +            | 11        |
        | 1            | 52472     | 53312     | +            | 3         |
        | 1            | 57597     | 64116     | +            | 7         |
        | ...          | ...       | ...       | ...          | ...       |
        | 1            | 1062207   | 1063288   | -            | 4         |
        | 1            | 1070966   | 1074306   | -            | 10        |
        | 1            | 1081817   | 1116361   | -            | 319       |
        | 1            | 1173055   | 1179555   | -            | 4         |
        +--------------+-----------+-----------+--------------+-----------+
        Stranded PyRanges object has 62 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.merge(by="Feature", count=True)
        +--------------+-----------+-----------+--------------+--------------+-----------+
        | Chromosome   | Start     | End       | Strand       | Feature      | Count     |
        | (category)   | (int32)   | (int32)   | (category)   | (category)   | (int32)   |
        |--------------+-----------+-----------+--------------+--------------+-----------|
        | 1            | 65564     | 65573     | +            | CDS          | 1         |
        | 1            | 69036     | 70005     | +            | CDS          | 2         |
        | 1            | 924431    | 924948    | +            | CDS          | 1         |
        | 1            | 925921    | 926013    | +            | CDS          | 11        |
        | ...          | ...       | ...       | ...          | ...          | ...       |
        | 1            | 1062207   | 1063288   | -            | transcript   | 1         |
        | 1            | 1070966   | 1074306   | -            | transcript   | 1         |
        | 1            | 1081817   | 1116361   | -            | transcript   | 19        |
        | 1            | 1173055   | 1179555   | -            | transcript   | 1         |
        +--------------+-----------+-----------+--------------+--------------+-----------+
        Stranded PyRanges object has 748 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.merge(by=["Feature", "gene_name"], count=True)
        +--------------+-----------+-----------+--------------+--------------+-------------+-----------+
        | Chromosome   | Start     | End       | Strand       | Feature      | gene_name   | Count     |
        | (category)   | (int32)   | (int32)   | (category)   | (category)   | (object)    | (int32)   |
        |--------------+-----------+-----------+--------------+--------------+-------------+-----------|
        | 1            | 1020172   | 1020373   | +            | CDS          | AGRN        | 1         |
        | 1            | 1022200   | 1022462   | +            | CDS          | AGRN        | 2         |
        | 1            | 1034555   | 1034703   | +            | CDS          | AGRN        | 2         |
        | 1            | 1035276   | 1035324   | +            | CDS          | AGRN        | 4         |
        | ...          | ...       | ...       | ...          | ...          | ...         | ...       |
        | 1            | 347981    | 348366    | -            | transcript   | RPL23AP24   | 1         |
        | 1            | 1173055   | 1179555   | -            | transcript   | TTLL10-AS1  | 1         |
        | 1            | 14403     | 29570     | -            | transcript   | WASH7P      | 1         |
        | 1            | 185216    | 195411    | -            | transcript   | WASH9P      | 1         |
        +--------------+-----------+-----------+--------------+--------------+-------------+-----------+
        Stranded PyRanges object has 807 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        if strand is None:
            strand = self.stranded

        kwargs = {"strand": strand, "count": count, "by": by, "count_col": count_col, "slack": slack}

        if not kwargs["by"]:
            kwargs["sparse"] = {"self": True}
            from pyranges.methods.merge import _merge
            df = pyrange_apply_single(_merge, self, **kwargs)
        else:
            kwargs["sparse"] = {"self": False}
            from pyranges.methods.merge import _merge_by
            df = pyrange_apply_single(_merge_by, self, **kwargs)

        return PyRanges(df)

    def mp(self, n=8, formatting=None):

        """Merge location and print.

        See Also
        --------

        PyRanges.print : print PyRanges."""

        print(tostring(self, n=n, merge_position=True, formatting=formatting))

    def mpc(self, n=8, formatting=None):

        """Merge location, print and return self.

        See Also
        --------

        PyRanges.print : print PyRanges."""

        print(tostring(self, n=n, merge_position=True, formatting=formatting))

        return self

    def msp(self, n=30, formatting=None):

        """Sort on location, merge location info and print.

        See Also
        --------

        PyRanges.print : print PyRanges."""

        print(
            tostring(
                self,
                n=n,
                merge_position=True,
                sort=True,
                formatting=formatting))


    def mspc(self, n=30, formatting=None):

        """Sort on location, merge location, print and return self.

        See Also
        --------

        PyRanges.print : print PyRanges."""

        print(
            tostring(
                self,
                n=n,
                merge_position=True,
                sort=True,
                formatting=formatting))

        return self


    def nearest(self, other, strandedness=None, overlap=True, how=None, suffix="_b", nb_cpu=1, apply_strand_suffix=None):

        """Find closest interval.

        Parameters
        ----------
        other : PyRanges

            PyRanges to find nearest interval in.

        strandedness : {None, "same", "opposite", False}, default None, i.e. auto

            Whether to compare PyRanges on the same strand, the opposite or ignore strand
            information. The default, None, means use "same" if both PyRanges are strande,
            otherwise ignore the strand information.

        overlap : bool, default True

            Whether to include overlaps.

        how : {None, "upstream", "downstream"}, default None, i.e. both directions

            Whether to only look for nearest in one direction. Always with respect to the PyRanges
            it is called on.

        suffix : str, default "_b"

            Suffix to give columns with shared name in other.

        apply_strand_suffix : bool, default None

            If first pyranges is unstranded, but the second is not, the first will be given the strand column of the second.
            apply_strand_suffix makes the added strand column a regular data column instead by adding a suffix.

        nb_cpu: int, default 1

            How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
            Will only lead to speedups on large datasets.

        Returns
        -------
        PyRanges

            A PyRanges with columns representing nearest interval horizontally appended.

        Notes
        -----

        A k_nearest also exists, but is less performant.

        See also
        --------

        PyRanges.new_position : give joined PyRanges new coordinates
        PyRanges.k_nearest : find k nearest intervals

        Examples
        --------

        >>> f1 = pr.from_dict({'Chromosome': ['chr1', 'chr1', 'chr1'], 'Start': [3, 8, 5],
        ...                    'End': [6, 9, 7], 'Strand': ['+', '+', '-']})
        >>> f1
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int32) |   (int32) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         3 |         6 | +            |
        | chr1         |         8 |         9 | +            |
        | chr1         |         5 |         7 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> f2 = pr.from_dict({'Chromosome': ['chr1', 'chr1'], 'Start': [1, 6],
        ...                    'End': [2, 7], 'Strand': ['+', '-']})
        >>> f2
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int32) |   (int32) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         1 |         2 | +            |
        | chr1         |         6 |         7 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> f1.nearest(f2)
        +--------------+-----------+-----------+--------------+-----------+-----------+--------------+------------+
        | Chromosome   |     Start |       End | Strand       |   Start_b |     End_b | Strand_b     |   Distance |
        | (category)   |   (int32) |   (int32) | (category)   |   (int32) |   (int32) | (category)   |    (int64) |
        |--------------+-----------+-----------+--------------+-----------+-----------+--------------+------------|
        | chr1         |         3 |         6 | +            |         6 |         7 | -            |          1 |
        | chr1         |         8 |         9 | +            |         6 |         7 | -            |          2 |
        | chr1         |         5 |         7 | -            |         6 |         7 | -            |          0 |
        +--------------+-----------+-----------+--------------+-----------+-----------+--------------+------------+
        Stranded PyRanges object has 3 rows and 8 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> f1.nearest(f2, how="upstream")
        +--------------+-----------+-----------+--------------+-----------+-----------+--------------+------------+
        | Chromosome   |     Start |       End | Strand       |   Start_b |     End_b | Strand_b     |   Distance |
        | (category)   |   (int32) |   (int32) | (category)   |   (int32) |   (int32) | (category)   |    (int64) |
        |--------------+-----------+-----------+--------------+-----------+-----------+--------------+------------|
        | chr1         |         3 |         6 | +            |         1 |         2 | +            |          2 |
        | chr1         |         8 |         9 | +            |         6 |         7 | -            |          2 |
        | chr1         |         5 |         7 | -            |         6 |         7 | -            |          0 |
        +--------------+-----------+-----------+--------------+-----------+-----------+--------------+------------+
        Stranded PyRanges object has 3 rows and 8 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        from pyranges.methods.nearest import _nearest

        kwargs = {"strandedness": strandedness, "how": how, "overlap": overlap, "nb_cpu": nb_cpu, "suffix": suffix, "apply_strand_suffix": apply_strand_suffix}
        kwargs = fill_kwargs(kwargs)
        if kwargs.get("how") in "upstream downstream".split():
            assert other.stranded, "If doing upstream or downstream nearest, other pyranges must be stranded"

        dfs = pyrange_apply(_nearest, self, other, **kwargs)
        gr = PyRanges(dfs)

        if not self.stranded and other.stranded:
            if apply_strand_suffix is None:
                import sys
                print("join: Strand data from other will be added as strand data to self.\nIf this is undesired use the flag apply_strand_suffix=False.\nTo turn off the warning set apply_strand_suffix to True or False.", file=sys.stderr)
            elif apply_strand_suffix:
                gr.columns = gr.columns.str.replace("Strand", "Strand" + kwargs["suffix"])

        return gr


    def new_position(self, new_pos, columns=None):

        """Give new position.

        The operation join produces a PyRanges with two pairs of start coordinates and two pairs of
        end coordinates. This operation uses these to give the PyRanges a new position.

        Parameters
        ----------
        new_pos : {"union", "intersection", "swap"}

           Change of coordinates.

        suffixes : tuple of str, default ("", "_b")

           Suffixes of coordinate-columns to switch.

        columns : tuple of str, default None, i.e. auto

           The name of the coordinate columns. By default uses the two first columns containing
           "Start" and the two first columns containing "End".

        See Also
        --------

        PyRanges.join : combine two PyRanges horizontally with SQL-style joins.

        Returns
        -------
        PyRanges

            PyRanges with new coordinates.

        Examples
        --------

        >>> gr = pr.from_dict({'Chromosome': ['chr1', 'chr1', 'chr1'],
        ...                    'Start': [3, 8, 5], 'End': [6, 9, 7]})
        >>> gr
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        | chr1         |         3 |         6 |
        | chr1         |         8 |         9 |
        | chr1         |         5 |         7 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 3 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr2 = pr.from_dict({'Chromosome': ['chr1', 'chr1'], 'Start': [1, 6],
        ...                     'End': [4, 7]})
        >>> gr2
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        | chr1         |         1 |         4 |
        | chr1         |         6 |         7 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 2 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> j = gr.join(gr2)
        >>> j
        +--------------+-----------+-----------+-----------+-----------+
        | Chromosome   |     Start |       End |   Start_b |     End_b |
        | (category)   |   (int32) |   (int32) |   (int32) |   (int32) |
        |--------------+-----------+-----------+-----------+-----------|
        | chr1         |         3 |         6 |         1 |         4 |
        | chr1         |         5 |         7 |         6 |         7 |
        +--------------+-----------+-----------+-----------+-----------+
        Unstranded PyRanges object has 2 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> j.new_position("swap")
        +--------------+-----------+-----------+-----------+-----------+
        | Chromosome   |     Start |       End |   Start_b |     End_b |
        | (category)   |   (int32) |   (int32) |   (int32) |   (int32) |
        |--------------+-----------+-----------+-----------+-----------|
        | chr1         |         1 |         4 |         3 |         6 |
        | chr1         |         6 |         7 |         5 |         7 |
        +--------------+-----------+-----------+-----------+-----------+
        Unstranded PyRanges object has 2 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> j.new_position("union").mp()
        +--------------------+-----------+-----------+
        | - Position -       |   Start_b |     End_b |
        | (Multiple types)   |   (int32) |   (int32) |
        |--------------------+-----------+-----------|
        | chr1 1-6           |         1 |         4 |
        | chr1 5-7           |         6 |         7 |
        +--------------------+-----------+-----------+
        Unstranded PyRanges object has 2 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> j.new_position("intersection").mp()
        +--------------------+-----------+-----------+
        | - Position -       |   Start_b |     End_b |
        | (Multiple types)   |   (int32) |   (int32) |
        |--------------------+-----------+-----------|
        | chr1 1-4           |         1 |         4 |
        | chr1 6-7           |         6 |         7 |
        +--------------------+-----------+-----------+
        Unstranded PyRanges object has 2 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> j2 = pr.from_dict({"Chromosome": [1], "Start": [3],
        ...                   "End": [4], "A": [1], "B": [3], "C": [2], "D": [5]})
        >>> j2
        +--------------+-----------+-----------+-----------+-----------+-----------+-----------+
        |   Chromosome |     Start |       End |         A |         B |         C |         D |
        |   (category) |   (int32) |   (int32) |   (int64) |   (int64) |   (int64) |   (int64) |
        |--------------+-----------+-----------+-----------+-----------+-----------+-----------|
        |            1 |         3 |         4 |         1 |         3 |         2 |         5 |
        +--------------+-----------+-----------+-----------+-----------+-----------+-----------+
        Unstranded PyRanges object has 1 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> j2.new_position("intersection", ("A", "B", "C", "D"))
        +--------------+-----------+-----------+-----------+-----------+-----------+-----------+
        |   Chromosome |     Start |       End |         A |         B |         C |         D |
        |   (category) |   (int32) |   (int32) |   (int64) |   (int64) |   (int64) |   (int64) |
        |--------------+-----------+-----------+-----------+-----------+-----------+-----------|
        |            1 |         2 |         3 |         1 |         3 |         2 |         5 |
        +--------------+-----------+-----------+-----------+-----------+-----------+-----------+
        Unstranded PyRanges object has 1 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        """

        from pyranges.methods.new_position import _new_position

        if self.empty:
            return self

        kwargs = {"strand": None}
        kwargs["sparse"] = {"self": False}
        kwargs["new_pos"] = new_pos

        if columns is None:
            start1, start2 = self.columns[self.columns.str.contains("Start")][:2]
            end1, end2 = self.columns[self.columns.str.contains("End")][:2]
            columns = (start1, end1, start2, end2)

        kwargs["columns"] = columns

        kwargs = fill_kwargs(kwargs)

        dfs = pyrange_apply_single(_new_position, self, **kwargs)

        return pr.PyRanges(dfs)


    def overlap(self, other, strandedness=None, how="first", invert=False, nb_cpu=1):

        """Return overlapping intervals.

        Returns the intervals in self which overlap with those in other.

        Parameters
        ----------
        other : PyRanges

            PyRanges to find overlaps with.

        strandedness : {None, "same", "opposite", False}, default None, i.e. auto

            Whether to compare PyRanges on the same strand, the opposite or ignore strand
            information. The default, None, means use "same" if both PyRanges are strande,
            otherwise ignore the strand information.

        how : {"first", "containment", False, None}, default "first"

            What intervals to report. By default reports every interval in self with overlap once.
            "containment" reports all intervals where the overlapping is contained within it.

        invert : bool, default False

            Whether to return the intervals without overlaps.

        nb_cpu: int, default 1

            How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
            Will only lead to speedups on large datasets.

        Returns
        -------
        PyRanges

            A PyRanges with overlapping intervals.

        See also
        --------

        PyRanges.intersect : report overlapping subintervals
        PyRanges.set_intersect : set-intersect PyRanges

        Examples
        --------

        >>> gr = pr.from_dict({"Chromosome": ["chr1"] * 3, "Start": [1, 4, 10],
        ...                    "End": [3, 9, 11], "ID": ["a", "b", "c"]})
        >>> gr
        +--------------+-----------+-----------+------------+
        |   Chromosome |     Start |       End | ID         |
        |   (category) |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------|
        |         chr1 |         1 |         3 | a          |
        |         chr1 |         4 |         9 | b          |
        |         chr1 |        10 |        11 | c          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr2 = pr.from_dict({"Chromosome": ["chr1"] * 3, "Start": [2, 2, 9], "End": [3, 9, 10]})
        >>> gr2
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        | chr1         |         2 |         3 |
        | chr1         |         2 |         9 |
        | chr1         |         9 |        10 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 3 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.overlap(gr2)
        +--------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | ID         |
        | (category)   |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------|
        | chr1         |         1 |         3 | a          |
        | chr1         |         4 |         9 | b          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.overlap(gr2, how=None)
        +--------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | ID         |
        | (category)   |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------|
        | chr1         |         1 |         3 | a          |
        | chr1         |         1 |         3 | a          |
        | chr1         |         4 |         9 | b          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.overlap(gr2, how="containment")
        +--------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | ID         |
        | (category)   |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------|
        | chr1         |         4 |         9 | b          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 1 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.overlap(gr2, invert=True)
        +--------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | ID         |
        | (category)   |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------|
        | chr1         |        10 |        11 | c          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 1 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        kwargs = {"strandedness": strandedness, "nb_cpu": nb_cpu}
        kwargs["sparse"] = {"self": False, "other": True}
        kwargs["how"] = how
        kwargs["invert"] = invert
        kwargs = fill_kwargs(kwargs)

        if len(self) == 0:
            return self

        if invert:
            self = self.copy()
            self.__ix__ = np.arange(len(self))

        dfs = pyrange_apply(_overlap, self, other, **kwargs)
        result = pr.PyRanges(dfs)

        if invert:
            found_idxs = getattr(result, "__ix__", [])
            result = self[~self.__ix__.isin(found_idxs)]
            result = result.drop("__ix__")

        return result

    def pc(self, n=8, formatting=None):

        """Print and return self.

        See Also
        --------

        PyRanges.print : print PyRanges."""

        print(tostring(self, n=n, formatting=formatting))

        return self

    def print(self, n=8, merge_position=False, sort=False, formatting=None, chain=False):

        """Print the PyRanges.

        Parameters
        ----------

        n : int, default 8

            The number of rows to print.

        merge_postion : bool, default False

            Print location in same column to save screen space.

        sort : bool or str, default False

            Sort the PyRanges before printing. Will print chromosomsomes or strands interleaved on
            sort columns.

        formatting : dict, default None

            Formatting options per column.

        chain : False

            Return the PyRanges. Useful to print intermediate results in call chains.

        See Also
        --------

        PyRanges.pc : print chain
        PyRanges.sp : sort print
        PyRanges.mp : merge print
        PyRanges.spc : sort print chain
        PyRanges.mpc : merge print chain
        PyRanges.msp : merge sort print
        PyRanges.mspc : merge sort print chain
        PyRanges.rp : raw print dictionary of DataFrames

        Examples
        --------

        >>> d = {'Chromosome': ['chr1', 'chr1', 'chr1'], 'Start': [3, 8, 5000],
        ...      'End': [6, 9, 7000], 'Name': ['i1', 'i3', 'i2'],
        ...      'Score': [1.1, 2.3987, 5.9999995], 'Strand': ['+', '+', '-']}
        >>> gr = pr.from_dict(d)
        >>> gr
        +--------------+-----------+-----------+------------+-------------+--------------+
        | Chromosome   |     Start |       End | Name       |       Score | Strand       |
        | (category)   |   (int32) |   (int32) | (object)   |   (float64) | (category)   |
        |--------------+-----------+-----------+------------+-------------+--------------|
        | chr1         |         3 |         6 | i1         |      1.1    | +            |
        | chr1         |         8 |         9 | i3         |      2.3987 | +            |
        | chr1         |      5000 |      7000 | i2         |      6      | -            |
        +--------------+-----------+-----------+------------+-------------+--------------+
        Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.print(formatting={"Start": "{:,}", "Score": "{:.2f}"})
        +--------------+-----------+-----------+------------+-------------+--------------+
        | Chromosome   |     Start |       End | Name       |       Score | Strand       |
        | (category)   |   (int32) |   (int32) | (object)   |   (float64) | (category)   |
        |--------------+-----------+-----------+------------+-------------+--------------|
        | chr1         |         3 |         6 | i1         |         1.1 | +            |
        | chr1         |         8 |         9 | i3         |         2.4 | +            |
        | chr1         |     5,000 |      7000 | i2         |         6   | -            |
        +--------------+-----------+-----------+------------+-------------+--------------+
        Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.print(merge_position=True) # gr.mp()
        +--------------------+------------+-------------+
        | - Position -       | Name       |       Score |
        | (Multiple types)   | (object)   |   (float64) |
        |--------------------+------------+-------------|
        | chr1 3-6 +         | i1         |      1.1    |
        | chr1 8-9 +         | i3         |      2.3987 |
        | chr1 5000-7000 -   | i2         |      6      |
        +--------------------+------------+-------------+
        Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> chipseq = pr.data.chipseq()
        >>> chipseq
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   | Start     | End       | Name       | Score     | Strand       |
        | (category)   | (int32)   | (int32)   | (object)   | (int64)   | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         | 212609534 | 212609559 | U0         | 0         | +            |
        | chr1         | 169887529 | 169887554 | U0         | 0         | +            |
        | chr1         | 216711011 | 216711036 | U0         | 0         | +            |
        | chr1         | 144227079 | 144227104 | U0         | 0         | +            |
        | ...          | ...       | ...       | ...        | ...       | ...          |
        | chrY         | 15224235  | 15224260  | U0         | 0         | -            |
        | chrY         | 13517892  | 13517917  | U0         | 0         | -            |
        | chrY         | 8010951   | 8010976   | U0         | 0         | -            |
        | chrY         | 7405376   | 7405401   | U0         | 0         | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 10,000 rows and 6 columns from 24 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        To interleave strands in output, use print with `sort=True`:

        >>> chipseq.print(sort=True, n=20) # chipseq.sp()
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   | Start     | End       | Name       | Score     | Strand       |
        | (category)   | (int32)   | (int32)   | (object)   | (int64)   | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         | 1325303   | 1325328   | U0         | 0         | -            |
        | chr1         | 1541598   | 1541623   | U0         | 0         | +            |
        | chr1         | 1599121   | 1599146   | U0         | 0         | +            |
        | chr1         | 1820285   | 1820310   | U0         | 0         | -            |
        | chr1         | 2448322   | 2448347   | U0         | 0         | -            |
        | chr1         | 3046141   | 3046166   | U0         | 0         | -            |
        | chr1         | 3437168   | 3437193   | U0         | 0         | -            |
        | chr1         | 3504032   | 3504057   | U0         | 0         | +            |
        | chr1         | 3637087   | 3637112   | U0         | 0         | -            |
        | chr1         | 3681903   | 3681928   | U0         | 0         | -            |
        | ...          | ...       | ...       | ...        | ...       | ...          |
        | chrY         | 15224235  | 15224260  | U0         | 0         | -            |
        | chrY         | 15548022  | 15548047  | U0         | 0         | +            |
        | chrY         | 16045242  | 16045267  | U0         | 0         | -            |
        | chrY         | 16495497  | 16495522  | U0         | 0         | -            |
        | chrY         | 21559181  | 21559206  | U0         | 0         | +            |
        | chrY         | 21707662  | 21707687  | U0         | 0         | -            |
        | chrY         | 21751211  | 21751236  | U0         | 0         | -            |
        | chrY         | 21910706  | 21910731  | U0         | 0         | -            |
        | chrY         | 22054002  | 22054027  | U0         | 0         | -            |
        | chrY         | 22210637  | 22210662  | U0         | 0         | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 10,000 rows and 6 columns from 24 chromosomes.
        For printing, the PyRanges was sorted on Chromosome, Start, End and Strand.

        >>> pr.data.chromsizes().print()
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
        """

        s = tostring(
            self,
            n=n,
            merge_position=merge_position,
            sort=sort,
            formatting=formatting)

        print(s)

        if chain:
            return self

    def rp(self):

        """Print dict of DataFrames.

        See Also
        --------

        PyRanges.print : print PyRanges."""

        print(self.dfs)



    def rpc(self):

        """Print dict of DataFrames and return self.

        See Also
        --------

        PyRanges.print : print PyRanges."""

        print(self.dfs)

        return self

    def sample(self, n=8, replace=False):
        """Subsample arbitrary rows of PyRanges.

        If n is larger than length of PyRanges, replace must be True.

        Parameters
        ----------
        n : int, default 8

            Number of rows to return

        replace : bool, False

            Reuse rows.

        Examples
        --------

        >>> gr = pr.data.chipseq()
        >>> np.random.seed(0)
        >>> gr.sample(n=3)
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr2         |  76564764 |  76564789 | U0         |         0 | +            |
        | chr3         | 185739979 | 185740004 | U0         |         0 | -            |
        | chr20        |  40373657 |  40373682 | U0         |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 6 columns from 3 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.sample(10001)
        Traceback (most recent call last):
        ...
        ValueError: Cannot take a larger sample than population when 'replace=False'
        """
        sample = np.random.choice(len(self), size=n, replace=False)
        subsetter = np.zeros(len(self), dtype=np.bool)
        subsetter[sample] = True
        return self[subsetter]


    def set_intersect(self, other, strandedness=None, how=None, new_pos=False, nb_cpu=1):

        """Return set-theoretical intersection.

        Like intersect, but both PyRanges are merged first.

        Parameters
        ----------
        other : PyRanges

            PyRanges to set-intersect.

        strandedness : {None, "same", "opposite", False}, default None, i.e. auto

            Whether to compare PyRanges on the same strand, the opposite or ignore strand
            information. The default, None, means use "same" if both PyRanges are strande,
            otherwise ignore the strand information.

        how : {None, "first", "last", "containment"}, default None, i.e. all

            What intervals to report. By default reports all overlapping intervals. "containment"
            reports intervals where the overlapping is contained within it.

        nb_cpu: int, default 1

            How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
            Will only lead to speedups on large datasets.

        Returns
        -------
        PyRanges

            A PyRanges with overlapping subintervals.

        See also
        --------

        PyRanges.intersect : find overlapping subintervals
        PyRanges.overlap : report overlapping intervals

        Examples
        --------

        >>> gr = pr.from_dict({"Chromosome": ["chr1"] * 3, "Start": [1, 4, 10],
        ...                    "End": [3, 9, 11], "ID": ["a", "b", "c"]})
        >>> gr
        +--------------+-----------+-----------+------------+
        |   Chromosome |     Start |       End | ID         |
        |   (category) |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------|
        |         chr1 |         1 |         3 | a          |
        |         chr1 |         4 |         9 | b          |
        |         chr1 |        10 |        11 | c          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr2 = pr.from_dict({"Chromosome": ["chr1"] * 3, "Start": [2, 2, 9], "End": [3, 9, 10]})
        >>> gr2
        +--------------+-----------+-----------+
        |   Chromosome |     Start |       End |
        |   (category) |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        |         chr1 |         2 |         3 |
        |         chr1 |         2 |         9 |
        |         chr1 |         9 |        10 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 3 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.set_intersect(gr2)
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        | chr1         |         2 |         3 |
        | chr1         |         4 |         9 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 2 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        In this simple unstranded case, this is the same as the below:

        >>> gr.merge().intersect(gr2.merge())
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        | chr1         |         2 |         3 |
        | chr1         |         4 |         9 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 2 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.set_intersect(gr2, how="containment")
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        | chr1         |         4 |         9 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 1 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        kwargs = {"strandedness": strandedness, "how": how, "nb_cpu": nb_cpu, "new_pos": new_pos}
        kwargs = fill_kwargs(kwargs)
        strand = True if strandedness else False
        self_clusters = self.merge(strand=strand)
        other_clusters = other.merge(strand=strand)
        dfs = pyrange_apply(_intersection, self_clusters, other_clusters,
                            **kwargs)

        return PyRanges(dfs)

    def set_union(self, other, strandedness=None, nb_cpu=1):

        """Return set-theoretical union.

        Parameters
        ----------
        other : PyRanges

            PyRanges to do union with.

        strandedness : {None, "same", "opposite", False}, default None, i.e. auto

            Whether to compare PyRanges on the same strand, the opposite or ignore strand
            information. The default, None, means use "same" if both PyRanges are strande,
            otherwise ignore the strand information.

        nb_cpu: int, default 1

            How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
            Will only lead to speedups on large datasets.

        Returns
        -------
        PyRanges

            A PyRanges with the union of intervals.

        See also
        --------

        PyRanges.set_intersect : set-theoretical intersection
        PyRanges.overlap : report overlapping intervals

        Examples
        --------

        >>> gr = pr.from_dict({"Chromosome": ["chr1"] * 3, "Start": [1, 4, 10],
        ...                    "End": [3, 9, 11], "ID": ["a", "b", "c"]})
        >>> gr
        +--------------+-----------+-----------+------------+
        |   Chromosome |     Start |       End | ID         |
        |   (category) |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------|
        |         chr1 |         1 |         3 | a          |
        |         chr1 |         4 |         9 | b          |
        |         chr1 |        10 |        11 | c          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr2 = pr.from_dict({"Chromosome": ["chr1"] * 3, "Start": [2, 2, 9], "End": [3, 9, 10]})
        >>> gr2
        +--------------+-----------+-----------+
        |   Chromosome |     Start |       End |
        |   (category) |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        |         chr1 |         2 |         3 |
        |         chr1 |         2 |         9 |
        |         chr1 |         9 |        10 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 3 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.set_union(gr2)
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        | chr1         |         1 |        11 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 1 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        if self.empty and other.empty:
            return pr.PyRanges()

        strand = True if strandedness else False

        if not strand:
            self = self.unstrand()
            other = other.unstrand()

        if strandedness == "opposite" and len(other):
            other = other.copy()
            other.Strand = other.Strand.replace({"+": "-", "-": "+"})


        gr = pr.concat([self, other], strand)

        gr = gr.merge(strand=strand)

        return gr


    def sort(self, by=None, nb_cpu=1):

        """Sort by position or columns.

        Parameters
        ----------
        by : str or list of str, default None

            Columns to sort by. Default is Start and End.

        nb_cpu: int, default 1

            How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
            Will only lead to speedups on large datasets.

        Note
        ----

        Since a PyRanges contains multiple DataFrames, the sorting only happens within dataframes.

        Returns
        -------
        PyRanges

            Sorted PyRanges

        See Also
        --------

        pyranges.multioverlap : find overlaps with multiple PyRanges

        Examples
        --------

        >>> gr = pr.data.f1()
        >>> gr
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         3 |         6 | interval1  |         0 | +            |
        | chr1         |         8 |         9 | interval3  |         0 | +            |
        | chr1         |         5 |         7 | interval2  |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.split(between=True)
        +--------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | Strand     |
        | (object)     |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------|
        | chr1         |         3 |         6 | +          |
        | chr1         |         6 |         8 | +          |
        | chr1         |         8 |         9 | +          |
        | chr1         |         5 |         7 | -          |
        +--------------+-----------+-----------+------------+
        Stranded PyRanges object has 4 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.split(strand=False)
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (object)     |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        | chr1         |         3 |         5 |
        | chr1         |         5 |         6 |
        | chr1         |         6 |         7 |
        | chr1         |         8 |         9 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 4 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.split(strand=False, between=True)
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (object)     |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        | chr1         |         3 |         5 |
        | chr1         |         5 |         6 |
        | chr1         |         6 |         7 |
        | chr1         |         7 |         8 |
        | chr1         |         8 |         9 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 5 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        """

        from pyranges.methods.sort import _sort
        kwargs = {"strand": self.stranded}
        kwargs["sparse"] = {"self": False}
        if by:
            kwargs["by"] = by
        kwargs = fill_kwargs(kwargs)
        return PyRanges(
            pyrange_apply_single(_sort, self, **kwargs))


    def sp(self, n=30, formatting=None):

        """Sort on location and print.

        See Also
        --------

        PyRanges.print : print PyRanges."""

        print(tostring(self, n=n, sort=True, formatting=formatting))

    def spc(self, n=30, formatting=None):

        """Sort on location, print and return self.

        See Also
        --------

        PyRanges.print : print PyRanges."""

        print(tostring(self, n=n, sort=True, formatting=formatting))

        return self


    def slack(self, slack):

        """Extend the intervals from the ends.

        Parameters
        ----------

        slack : int or dict of ints with "3" and/or "5" as keys.

            The number of nucleotides to extend the ends with.

        Examples
        --------

        >>> d = {'Chromosome': ['chr1', 'chr1', 'chr1'], 'Start': [3, 8, 5], 'End': [6, 9, 7],
        ...      'Strand': ['+', '+', '-']}
        >>> gr = pr.from_dict(d)
        >>> gr
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int32) |   (int32) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         3 |         6 | +            |
        | chr1         |         8 |         9 | +            |
        | chr1         |         5 |         7 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.slack(4)
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int32) |   (int32) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         0 |        10 | +            |
        | chr1         |         4 |        13 | +            |
        | chr1         |         1 |        11 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.slack({"3": 1})
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int32) |   (int32) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         3 |         7 | +            |
        | chr1         |         8 |        10 | +            |
        | chr1         |         4 |         7 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.slack({"3": 1, "5": 2})
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int32) |   (int32) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         1 |         7 | +            |
        | chr1         |         6 |        10 | +            |
        | chr1         |         4 |         9 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.slack(-1)
        Traceback (most recent call last):
        ...
        AssertionError: Some intervals are negative or zero length after applying slack!
        """

        if isinstance(slack, dict):
            assert self.stranded, "PyRanges must be stranded to add 5/3-end specific slack."

        kwargs = fill_kwargs({"slack": slack, "strand": self.stranded})

        prg = PyRanges(
            pyrange_apply_single(_slack, self, **kwargs))

        return prg


    def split(self, strand=None, between=False, nb_cpu=1):

        """Split into non-overlapping intervals.

        Parameters
        ----------
        strand : bool, default None, i.e. auto

            Whether to ignore strand information if PyRanges is stranded.

        between : bool, default False

            Include lengths between intervals.

        nb_cpu: int, default 1

            How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
            Will only lead to speedups on large datasets.

        Returns
        -------
        PyRanges

            PyRanges with intervals split at overlap points.

        See Also
        --------

        pyranges.multioverlap : find overlaps with multiple PyRanges

        Examples
        --------

        >>> d = {'Chromosome': ['chr1', 'chr1', 'chr1', 'chr1'], 'Start': [3, 5, 5, 11],
        ...       'End': [6, 9, 7, 12], 'Strand': ['+', '+', '-', '-']}
        >>> gr = pr.from_dict(d)
        >>> gr
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int32) |   (int32) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         3 |         6 | +            |
        | chr1         |         5 |         9 | +            |
        | chr1         |         5 |         7 | -            |
        | chr1         |        11 |        12 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 4 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.split()
        +--------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | Strand     |
        | (object)     |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------|
        | chr1         |         3 |         5 | +          |
        | chr1         |         5 |         6 | +          |
        | chr1         |         6 |         9 | +          |
        | chr1         |         5 |         7 | -          |
        | chr1         |        11 |        12 | -          |
        +--------------+-----------+-----------+------------+
        Stranded PyRanges object has 5 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.split(between=True)
        +--------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | Strand     |
        | (object)     |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------|
        | chr1         |         3 |         5 | +          |
        | chr1         |         5 |         6 | +          |
        | chr1         |         6 |         9 | +          |
        | chr1         |         5 |         7 | -          |
        | chr1         |         7 |        11 | -          |
        | chr1         |        11 |        12 | -          |
        +--------------+-----------+-----------+------------+
        Stranded PyRanges object has 6 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.split(strand=False)
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (object)     |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        | chr1         |         3 |         5 |
        | chr1         |         5 |         6 |
        | chr1         |         6 |         7 |
        | chr1         |         7 |         9 |
        | chr1         |        11 |        12 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 5 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.split(strand=False, between=True)
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (object)     |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        | chr1         |         3 |         5 |
        | chr1         |         5 |         6 |
        | chr1         |         6 |         7 |
        | chr1         |         7 |         9 |
        | chr1         |         9 |        11 |
        | chr1         |        11 |        12 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 6 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        if strand is None:
            strand = self.stranded

        kwargs = fill_kwargs({"strand": strand})

        from pyranges.methods.split import _split
        df = pyrange_apply_single(_split, self, **kwargs)

        split = pr.PyRanges(df)
        if not between:
            strandedness = "same" if strand else False
            split = split.overlap(self, strandedness=strandedness)

        return split

    @property
    def stranded(self):
        """Whether PyRanges has (valid) strand info.

        Note
        ----

        A PyRanges can have invalid values in the Strand-column. It is not considered stranded.

        See Also
        --------

        PyRanges.strands : return the strands

        Examples
        --------

        >>> d =  {'Chromosome': ['chr1', 'chr1'], 'Start': [1, 6],
        ...       'End': [5, 8], 'Strand': ['+', '.']}
        >>> gr = pr.from_dict(d)
        >>> gr
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int32) |   (int32) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         1 |         5 | +            |
        | chr1         |         6 |         8 | .            |
        +--------------+-----------+-----------+--------------+
        Unstranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        Considered unstranded due to these Strand values: '.'

        >>> gr.stranded
        False

        >>> "Strand" in gr.columns
        True
        """
        keys = self.keys()

        if not len(keys):
            # so that stranded ops work with empty dataframes
            return True

        key = keys[0]

        return isinstance(key, tuple)

    @property
    def strands(self):

        """Return strands.

        Notes
        -----

        If the strand-column contains an invalid value, [] is returned.

        See Also
        --------

        PyRanges.stranded : whether has valid strand info

        Examples
        --------
        >>> d =  {'Chromosome': ['chr1', 'chr1'], 'Start': [1, 6],
        ...       'End': [5, 8], 'Strand': ['+', '.']}
        >>> gr = pr.from_dict(d)
        >>> gr
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int32) |   (int32) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         1 |         5 | +            |
        | chr1         |         6 |         8 | .            |
        +--------------+-----------+-----------+--------------+
        Unstranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        Considered unstranded due to these Strand values: '.'

        >>> gr.strands
        []

        >>> gr.Strand.drop_duplicates().to_list()
        ['+', '.']

        >>> gr.Strand = ["+", "-"]
        >>> gr.strands
        ['+', '-']
        """

        if not self.stranded:
            return []

        return natsorted(set([k[1] for k in self.keys()]))


    def subset(self, f, strand=None, **kwargs):

        """Return a subset of the rows.

        Parameters
        ----------
        f : function
            Function which returns boolean Series equal to length of df.

        strand : bool, default None, i.e. auto

            Whether to do operations on chromosome/strand pairs or chromosomes. If None, will use
            chromosome/strand pairs if the PyRanges is stranded.

        nb_cpu : int, default 1

            How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
            Will only lead to speedups on large datasets.

        **kwargs
            Additional keyword arguments to pass as keyword arguments to `f`

        Notes
        -----

        PyRanges can also be subsetted directly with a boolean Series. This function is slightly
        faster, but more cumbersome.

        Returns
        -------
        PyRanges

            PyRanges subset on rows.

        Examples
        --------

        >>> gr = pr.data.f1()
        >>> gr
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         3 |         6 | interval1  |         0 | +            |
        | chr1         |         8 |         9 | interval3  |         0 | +            |
        | chr1         |         5 |         7 | interval2  |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.subset(lambda df: df.Start > 4)
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         8 |         9 | interval3  |         0 | +            |
        | chr1         |         5 |         7 | interval2  |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        Also possible:

        >>> gr[gr.Start > 4]
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         8 |         9 | interval3  |         0 | +            |
        | chr1         |         5 |         7 | interval2  |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        kwargs = fill_kwargs(kwargs)

        if strand is None:
            strand = self.stranded

        if self.stranded and not strand:
            self = self.unstrand()


        kwargs.update({"strand": strand})

        result = pyrange_apply_single(f, self, **kwargs)

        if not result:
            return pr.PyRanges()

        first_result = next(iter(result.values()))

        assert first_result.dtype == bool, "result of subset function must be bool, but is {}".format(
            first_result.dtype)

        return self[result]


    def subtract(self, other, strandedness=None, nb_cpu=1):

        """Subtract intervals.

        Parameters
        ----------
        strandedness : {None, "same", "opposite", False}, default None, i.e. auto

            Whether to compare PyRanges on the same strand, the opposite or ignore strand
            information. The default, None, means use "same" if both PyRanges are strande,
            otherwise ignore the strand information.

        nb_cpu: int, default 1

            How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
            Will only lead to speedups on large datasets.

        See Also
        --------
        pyranges.PyRanges.overlap : use with invert=True to return all intervals without overlap

        Examples
        --------

        >>> gr = pr.from_dict({"Chromosome": ["chr1"] * 3, "Start": [1, 4, 10],
        ...                    "End": [3, 9, 11], "ID": ["a", "b", "c"]})
        >>> gr2 = pr.from_dict({"Chromosome": ["chr1"] * 3, "Start": [2, 2, 9], "End": [3, 9, 10]})
        >>> gr
        +--------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | ID         |
        | (category)   |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------|
        | chr1         |         1 |         3 | a          |
        | chr1         |         4 |         9 | b          |
        | chr1         |        10 |        11 | c          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr2
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        | chr1         |         2 |         3 |
        | chr1         |         2 |         9 |
        | chr1         |         9 |        10 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 3 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.subtract(gr2)
        +--------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | ID         |
        | (category)   |   (int32) |   (int32) | (object)   |
        |--------------+-----------+-----------+------------|
        | chr1         |         1 |         2 | a          |
        | chr1         |        10 |        11 | c          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        from pyranges.methods.subtraction import _subtraction

        kwargs = {"strandedness": strandedness}
        kwargs["sparse"] = {"self": False, "other": True}
        kwargs = fill_kwargs(kwargs)

        strand = True if strandedness else False
        other_clusters = other.merge(strand=strand)

        self = self.count_overlaps(other_clusters, strandedness=strandedness, overlap_col="__num__")

        result = pyrange_apply(_subtraction, self, other_clusters, **kwargs)

        self = self.drop("__num__")

        return PyRanges(result).drop("__num__")


    def summary(self, to_stdout=True, return_df=False):

        """Return info.

        Count refers to the number of intervals, the rest to the lengths.

        The column "pyrange" describes the data as is. "coverage_forward" and "coverage_reverse"
        describe the data after strand-specific merging of overlapping intervals.
        "coverage_unstranded" describes the data after merging, without considering the strands.

        The row "count" is the number of intervals and "sum" is their total length. The rest describe the lengths of the
        intervals.

        Parameters
        ----------
        to_stdout : bool, default True

            Print summary.

        return_df : bool, default False

            Return df with summary.

        Returns
        -------
            None or DataFrame with summary.


        Examples
        --------

        >>> gr = pr.data.ensembl_gtf()[["Feature", "gene_id"]]
        >>> gr
        +--------------+--------------+-----------+-----------+--------------+-----------------+
        | Chromosome   | Feature      | Start     | End       | Strand       | gene_id         |
        | (category)   | (category)   | (int32)   | (int32)   | (category)   | (object)        |
        |--------------+--------------+-----------+-----------+--------------+-----------------|
        | 1            | gene         | 11868     | 14409     | +            | ENSG00000223972 |
        | 1            | transcript   | 11868     | 14409     | +            | ENSG00000223972 |
        | 1            | exon         | 11868     | 12227     | +            | ENSG00000223972 |
        | 1            | exon         | 12612     | 12721     | +            | ENSG00000223972 |
        | ...          | ...          | ...       | ...       | ...          | ...             |
        | 1            | gene         | 1173055   | 1179555   | -            | ENSG00000205231 |
        | 1            | transcript   | 1173055   | 1179555   | -            | ENSG00000205231 |
        | 1            | exon         | 1179364   | 1179555   | -            | ENSG00000205231 |
        | 1            | exon         | 1173055   | 1176396   | -            | ENSG00000205231 |
        +--------------+--------------+-----------+-----------+--------------+-----------------+
        Stranded PyRanges object has 2,446 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.summary()
        +-------+------------------+--------------------+--------------------+-----------------------+
        |       |          pyrange |   coverage_forward |   coverage_reverse |   coverage_unstranded |
        |-------+------------------+--------------------+--------------------+-----------------------|
        | count |   2446           |               39   |               23   |                  32   |
        | mean  |   2291.92        |             7058.1 |            30078.6 |               27704.2 |
        | std   |  11906.9         |            10322.3 |            59467.7 |               67026.9 |
        | min   |      1           |               83   |              154   |                  83   |
        | 25%   |     90           |             1051   |             1204   |                1155   |
        | 50%   |    138           |             2541   |             6500   |                6343   |
        | 75%   |    382.25        |             7168   |            23778   |               20650.8 |
        | max   | 241726           |            43065   |           241726   |              291164   |
        | sum   |      5.60603e+06 |           275266   |           691807   |              886534   |
        +-------+------------------+--------------------+--------------------+-----------------------+

        >>> gr.summary(return_df=True, to_stdout=False)
                    pyrange  coverage_forward  coverage_reverse  coverage_unstranded
        count  2.446000e+03         39.000000         23.000000            32.000000
        mean   2.291918e+03       7058.102564      30078.565217         27704.187500
        std    1.190685e+04      10322.309347      59467.695265         67026.868647
        min    1.000000e+00         83.000000        154.000000            83.000000
        25%    9.000000e+01       1051.000000       1204.000000          1155.000000
        50%    1.380000e+02       2541.000000       6500.000000          6343.000000
        75%    3.822500e+02       7168.000000      23778.000000         20650.750000
        max    2.417260e+05      43065.000000     241726.000000        291164.000000
        sum    5.606031e+06     275266.000000     691807.000000        886534.000000
        """

        from pyranges.methods.summary import _summary

        return _summary(self, to_stdout, return_df)


    def tail(self, n=8):

        """Return the n last rows.

        Parameters
        ----------

        n : int, default 8

            Return n rows.

        Returns
        -------
        PyRanges

            PyRanges with the n last rows.

        See Also
        --------

        PyRanges.head : return the first rows
        PyRanges.sample : return random rows

        Examples
        --------

        >>> gr = pr.data.chipseq()
        >>> gr
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   | Start     | End       | Name       | Score     | Strand       |
        | (category)   | (int32)   | (int32)   | (object)   | (int64)   | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         | 212609534 | 212609559 | U0         | 0         | +            |
        | chr1         | 169887529 | 169887554 | U0         | 0         | +            |
        | chr1         | 216711011 | 216711036 | U0         | 0         | +            |
        | chr1         | 144227079 | 144227104 | U0         | 0         | +            |
        | ...          | ...       | ...       | ...        | ...       | ...          |
        | chrY         | 15224235  | 15224260  | U0         | 0         | -            |
        | chrY         | 13517892  | 13517917  | U0         | 0         | -            |
        | chrY         | 8010951   | 8010976   | U0         | 0         | -            |
        | chrY         | 7405376   | 7405401   | U0         | 0         | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 10,000 rows and 6 columns from 24 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.tail(3)
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chrY         |  13517892 |  13517917 | U0         |         0 | -            |
        | chrY         |   8010951 |   8010976 | U0         |         0 | -            |
        | chrY         |   7405376 |   7405401 | U0         |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        subsetter = np.zeros(len(self), dtype=np.bool)
        subsetter[(len(self) - n):] = True
        return self[subsetter]


    def tile(self, tile_size, overlap=False, strand=None, nb_cpu=1):

        """Return overlapping genomic tiles.

        The genome is divided into bookended tiles of length `tile_size` and one is returned per
        overlapping interval.

        Parameters
        ----------
        tile_size : int
            Length of the tiles.

        overlap : bool, default False

            Add column of nucleotide overlap to each tile.

        strand : bool, default None, i.e. auto

            Whether to do operations on chromosome/strand pairs or chromosomes. If None, will use
            chromosome/strand pairs if the PyRanges is stranded.

        nb_cpu: int, default 1

            How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
            Will only lead to speedups on large datasets.

        **kwargs
            Additional keyword arguments to pass as keyword arguments to `f`

        Returns
        -------
        PyRanges

            Tiled PyRanges.

        See also
        --------

        pyranges.PyRanges.window : divide intervals into windows

        Examples
        --------

        >>> gr = pr.data.ensembl_gtf()[["Feature", "gene_name"]]
        >>> gr
        +--------------+--------------+-----------+-----------+--------------+-------------+
        | Chromosome   | Feature      | Start     | End       | Strand       | gene_name   |
        | (category)   | (category)   | (int32)   | (int32)   | (category)   | (object)    |
        |--------------+--------------+-----------+-----------+--------------+-------------|
        | 1            | gene         | 11868     | 14409     | +            | DDX11L1     |
        | 1            | transcript   | 11868     | 14409     | +            | DDX11L1     |
        | 1            | exon         | 11868     | 12227     | +            | DDX11L1     |
        | 1            | exon         | 12612     | 12721     | +            | DDX11L1     |
        | ...          | ...          | ...       | ...       | ...          | ...         |
        | 1            | gene         | 1173055   | 1179555   | -            | TTLL10-AS1  |
        | 1            | transcript   | 1173055   | 1179555   | -            | TTLL10-AS1  |
        | 1            | exon         | 1179364   | 1179555   | -            | TTLL10-AS1  |
        | 1            | exon         | 1173055   | 1176396   | -            | TTLL10-AS1  |
        +--------------+--------------+-----------+-----------+--------------+-------------+
        Stranded PyRanges object has 2,446 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.tile(200)
        +--------------+--------------+-----------+-----------+--------------+-------------+
        | Chromosome   | Feature      | Start     | End       | Strand       | gene_name   |
        | (category)   | (category)   | (int32)   | (int32)   | (category)   | (object)    |
        |--------------+--------------+-----------+-----------+--------------+-------------|
        | 1            | gene         | 11800     | 12000     | +            | DDX11L1     |
        | 1            | gene         | 12000     | 12200     | +            | DDX11L1     |
        | 1            | gene         | 12200     | 12400     | +            | DDX11L1     |
        | 1            | gene         | 12400     | 12600     | +            | DDX11L1     |
        | ...          | ...          | ...       | ...       | ...          | ...         |
        | 1            | exon         | 1175600   | 1175800   | -            | TTLL10-AS1  |
        | 1            | exon         | 1175800   | 1176000   | -            | TTLL10-AS1  |
        | 1            | exon         | 1176000   | 1176200   | -            | TTLL10-AS1  |
        | 1            | exon         | 1176200   | 1176400   | -            | TTLL10-AS1  |
        +--------------+--------------+-----------+-----------+--------------+-------------+
        Stranded PyRanges object has 30,538 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.tile(100, overlap=True)
        +--------------+--------------+-----------+-----------+--------------+-------------+---------------+
        | Chromosome   | Feature      | Start     | End       | Strand       | gene_name   | TileOverlap   |
        | (category)   | (category)   | (int32)   | (int32)   | (category)   | (object)    | (int32)       |
        |--------------+--------------+-----------+-----------+--------------+-------------+---------------|
        | 1            | gene         | 11800     | 11900     | +            | DDX11L1     | 32            |
        | 1            | gene         | 11900     | 12000     | +            | DDX11L1     | 100           |
        | 1            | gene         | 12000     | 12100     | +            | DDX11L1     | 100           |
        | 1            | gene         | 12100     | 12200     | +            | DDX11L1     | 100           |
        | ...          | ...          | ...       | ...       | ...          | ...         | ...           |
        | 1            | exon         | 1176000   | 1176100   | -            | TTLL10-AS1  | 100           |
        | 1            | exon         | 1176100   | 1176200   | -            | TTLL10-AS1  | 100           |
        | 1            | exon         | 1176200   | 1176300   | -            | TTLL10-AS1  | 100           |
        | 1            | exon         | 1176300   | 1176400   | -            | TTLL10-AS1  | 96            |
        +--------------+--------------+-----------+-----------+--------------+-------------+---------------+
        Stranded PyRanges object has 58,516 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        from pyranges.methods.windows import _tiles

        if strand is None:
            strand = self.stranded

        kwargs = {"strand": strand, "overlap": overlap}
        kwargs["sparse"] = {"self": False}
        kwargs["tile_size"] = tile_size

        df = pyrange_apply_single(_tiles, self, **kwargs)

        return PyRanges(df)

    def to_example(self, n=10):

        """Return as dict.

        Used for easily creating examples for copy and pasting.

        Parameters
        ----------
        n : int, default 10
            Number of rows. Half is taken from the start, the other half from the end.

        See Also
        --------

        PyRanges.from_dict : create PyRanges from dict

        Examples
        --------

        >>> gr = pr.data.chipseq()
        >>> gr
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   | Start     | End       | Name       | Score     | Strand       |
        | (category)   | (int32)   | (int32)   | (object)   | (int64)   | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         | 212609534 | 212609559 | U0         | 0         | +            |
        | chr1         | 169887529 | 169887554 | U0         | 0         | +            |
        | chr1         | 216711011 | 216711036 | U0         | 0         | +            |
        | chr1         | 144227079 | 144227104 | U0         | 0         | +            |
        | ...          | ...       | ...       | ...        | ...       | ...          |
        | chrY         | 15224235  | 15224260  | U0         | 0         | -            |
        | chrY         | 13517892  | 13517917  | U0         | 0         | -            |
        | chrY         | 8010951   | 8010976   | U0         | 0         | -            |
        | chrY         | 7405376   | 7405401   | U0         | 0         | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 10,000 rows and 6 columns from 24 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> d = gr.to_example(n=4)
        >>> d
        {'Chromosome': ['chr1', 'chr1', 'chrY', 'chrY'], 'Start': [212609534, 169887529, 8010951, 7405376], 'End': [212609559, 169887554, 8010976, 7405401], 'Name': ['U0', 'U0', 'U0', 'U0'], 'Score': [0, 0, 0, 0], 'Strand': ['+', '+', '-', '-']}
        >>> pr.from_dict(d)
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int32) |   (int32) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         | 212609534 | 212609559 | U0         |         0 | +            |
        | chr1         | 169887529 | 169887554 | U0         |         0 | +            |
        | chrY         |   8010951 |   8010976 | U0         |         0 | -            |
        | chrY         |   7405376 |   7405401 | U0         |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 4 rows and 6 columns from 2 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        nrows_half = int(min(n, len(self))/2)

        if n < len(self):
            first = self.head(nrows_half)
            last = self.tail(nrows_half)
            example = pr.concat([first, last])
        else:
            example = self

        d = {c: list(getattr(example, c)) for c in example.columns}

        return d

    def three_end(self):

        """Return the 3'-end.

        The 3'-end is the start of intervals on the reverse strand and the end of intervals on the
        forward strand.

        Returns
        -------
        PyRanges
            PyRanges with the 3'.

        See Also
        --------
        PyRanges.five_end : return the five prime end

        Examples
        --------

        >>> d =  {'Chromosome': ['chr1', 'chr1'], 'Start': [1, 6],
        ...       'End': [5, 8], 'Strand': ['+', '-']}
        >>> gr = pr.from_dict(d)
        >>> gr
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int32) |   (int32) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         1 |         5 | +            |
        | chr1         |         6 |         8 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.three_end()
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int32) |   (int32) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         5 |         6 | +            |
        | chr1         |         6 |         7 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        """

        assert self.stranded, "Need stranded pyrange to find 3'."
        kwargs = fill_kwargs({"strand": True})
        return PyRanges(
            pyrange_apply_single(_tes, self, **kwargs))


#     def to_bam(self, path=None, header=None, chromosome_sizes=None, chain=False):

#         r"""Write to bam.

#         Parameters
#         ----------
#         path : str, default None
#             Where to write. If None, returns string representation.

#         keep : bool, default True

#             Whether to keep all columns, not just Chromosome, Start, End,
#             Name, Score, Strand when writing.

#         compression : str, compression type to use, by default infer based on extension.
#             See pandas.DataFree.to_csv for more info.

#         header : dict

#             Header to use in the bamfile. See the pysam docs for how it should look.
#             Or use the header attribute from another pyasam.AlignmentFile.

#         chromosome_sizes : PyRanges or dict

#             If dict: map of chromosome names to chromosome length.

#         chain : bool, default False
#             Whether to return the PyRanges after writing.

#         Note
#         ----

#         The following pyranges columns are used when writing:

#         Chromosome, Start, End, Strand, MapQ, Flag, QueryStart, QueryEnd, Name, Cigar, Quality

#         Examples
#         --------

#         >>> header = {"SQ": [{"SN": 1, "LN": 249250621}]}

#         >>> c = '''Name	Flag Chromosome Start End MapQ Cigar QuerySequence Quality
# read1	115	1	142618765 142618790	255	25M	CGACCCACTCCGCCATTTTCATCCG	IIGIIIHIGIIFIIIIIIIGIGIII	NM:i:0ZP:i:65536	ZL:i:25
# read2	115	1	142618765 142618790  255	25M	CGACCCACTCCGCCATTTTCATCCG	IIGIIIHIGIIFIIIIIIIGIGIII	NM:i:0ZP:i:214748	ZL:i:25
# read3	115	1	142618765 142618790  255	25M	CGACCCACTCCGCCATTTTCATCCG	IIGIIIHIGIIFIIIIIIIGIGIII	NM:i:0ZP:i:2147484	ZL:i:25
# read4	115	1	142618765 142618790  255	25M	CGACCCACTCCGCCATTTTCATCCG	IIGIIIHIGIIFIIIIIIIGIGIII	NM:i:0ZP:i:2147483647	ZL:i:25
# read5	115	1	142618765 142618790  255	25M	CGACCCACTCCGCCATTTTCATCCG	IIGIIIHIGIIFIIIIIIIGIGIII	NM:i:0ZP:i:-65536	ZL:i:25
# read6	115	1	142618765 142618790  255	25M	CGACCCACTCCGCCATTTTCATCCG	IIGIIIHIGIIFIIIIIIIGIGIII	NM:i:0ZP:i:-214748	ZL:i:25
# read7	115	1	142618765 142618790  255	25M	CGACCCACTCCGCCATTTTCATCCG	IIGIIIHIGIIFIIIIIIIGIGIII	NM:i:0ZP:i:-2147484	ZL:i:25
# read8	115	1	142618765 142618790  255	25M	CGACCCACTCCGCCATTTTCATCCG	IIGIIIHIGIIFIIIIIIIGIGIII	NM:i:0ZP:i:-2147483647	ZL:i:25'''

#         >>>
#         """




    def to_bed(self, path=None, keep=True, compression="infer", chain=False):
        r"""Write to bed.

        Parameters
        ----------
        path : str, default None
            Where to write. If None, returns string representation.

        keep : bool, default True

            Whether to keep all columns, not just Chromosome, Start, End,
            Name, Score, Strand when writing.

        compression : str, compression type to use, by default infer based on extension.
            See pandas.DataFree.to_csv for more info.

        chain : bool, default False
            Whether to return the PyRanges after writing.

        Examples
        --------

        >>> d =  {'Chromosome': ['chr1', 'chr1'], 'Start': [1, 6],
        ...       'End': [5, 8], 'Strand': ['+', '-'], "Gene": [1, 2]}
        >>> gr = pr.from_dict(d)
        >>> gr
        +--------------+-----------+-----------+--------------+-----------+
        | Chromosome   |     Start |       End | Strand       |      Gene |
        | (category)   |   (int32) |   (int32) | (category)   |   (int64) |
        |--------------+-----------+-----------+--------------+-----------|
        | chr1         |         1 |         5 | +            |         1 |
        | chr1         |         6 |         8 | -            |         2 |
        +--------------+-----------+-----------+--------------+-----------+
        Stranded PyRanges object has 2 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> print(gr.to_bed())
        chr1	1	5	.	.	+	1
        chr1	6	8	.	.	-	2
        <BLANKLINE>

        Does not include noncanonical bed-column `Gene`:

        >>> print(gr.to_bed(keep=False))
        chr1	1	5	.	.	+
        chr1	6	8	.	.	-
        <BLANKLINE>

        >>> gr.to_bed("test.bed", chain=True)
        +--------------+-----------+-----------+--------------+-----------+
        | Chromosome   |     Start |       End | Strand       |      Gene |
        | (category)   |   (int32) |   (int32) | (category)   |   (int64) |
        |--------------+-----------+-----------+--------------+-----------|
        | chr1         |         1 |         5 | +            |         1 |
        | chr1         |         6 |         8 | -            |         2 |
        +--------------+-----------+-----------+--------------+-----------+
        Stranded PyRanges object has 2 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> open("test.bed").readlines()
        ['chr1\t1\t5\t.\t.\t+\t1\n', 'chr1\t6\t8\t.\t.\t-\t2\n']
        """
        from pyranges.out import _to_bed

        result = _to_bed(self, path, keep=keep, compression=compression)

        if path and chain:
            return self
        else:
            return result

    def to_bigwig(self, path=None, chromosome_sizes=None, rpm=True, divide=None, value_col=None, dryrun=False, chain=False):

        """Write regular or value coverage to bigwig.

        Note
        ----

        To create one bigwig per strand, subset the PyRanges first.

        Parameters
        ----------
        path : str

            Where to write bigwig.

        chromosome_sizes : PyRanges or dict

            If dict: map of chromosome names to chromosome length.

        rpm : True

            Whether to normalize data by dividing by total number of intervals and multiplying by
            1e6.

        divide : bool, default False

            (Only useful with value_col) Divide value coverage by regular coverage and take log2.

        value_col : str, default None

            Name of column to compute coverage of.

        dryrun : bool, default False

            Return data that would be written without writing bigwigs.

        chain : bool, default False
            Whether to return the PyRanges after writing.

        Note
        ----

        Requires pybigwig to be installed.

        If you require more control over the normalization process, use pyranges.to_bigwig()

        See Also
        --------
        pyranges.to_bigwig : write pandas DataFrame to bigwig.

        Examples
        --------

        >>> d =  {'Chromosome': ['chr1', 'chr1', 'chr1'], 'Start': [1, 4, 6],
        ...       'End': [7, 8, 10], 'Strand': ['+', '-', '-'],
        ...       'Value': [10, 20, 30]}
        >>> gr = pr.from_dict(d)
        >>> gr
        +--------------+-----------+-----------+--------------+-----------+
        | Chromosome   |     Start |       End | Strand       |     Value |
        | (category)   |   (int32) |   (int32) | (category)   |   (int64) |
        |--------------+-----------+-----------+--------------+-----------|
        | chr1         |         1 |         7 | +            |        10 |
        | chr1         |         4 |         8 | -            |        20 |
        | chr1         |         6 |        10 | -            |        30 |
        +--------------+-----------+-----------+--------------+-----------+
        Stranded PyRanges object has 3 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.to_bigwig(dryrun=True, rpm=False)
        +--------------+-----------+-----------+-------------+
        | Chromosome   |     Start |       End |       Score |
        | (category)   |   (int32) |   (int32) |   (float64) |
        |--------------+-----------+-----------+-------------|
        | chr1         |         1 |         4 |           1 |
        | chr1         |         4 |         6 |           2 |
        | chr1         |         6 |         7 |           3 |
        | chr1         |         7 |         8 |           2 |
        | chr1         |         8 |        10 |           1 |
        +--------------+-----------+-----------+-------------+
        Unstranded PyRanges object has 5 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.to_bigwig(dryrun=True, rpm=False, value_col="Value")
        +--------------+-----------+-----------+-------------+
        | Chromosome   |     Start |       End |       Score |
        | (category)   |   (int32) |   (int32) |   (float64) |
        |--------------+-----------+-----------+-------------|
        | chr1         |         1 |         4 |          10 |
        | chr1         |         4 |         6 |          30 |
        | chr1         |         6 |         7 |          60 |
        | chr1         |         7 |         8 |          50 |
        | chr1         |         8 |        10 |          30 |
        +--------------+-----------+-----------+-------------+
        Unstranded PyRanges object has 5 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.to_bigwig(dryrun=True, rpm=False, value_col="Value", divide=True)
        +--------------+-----------+-----------+-------------+
        | Chromosome   |     Start |       End |       Score |
        | (category)   |   (int32) |   (int32) |   (float64) |
        |--------------+-----------+-----------+-------------|
        | chr1         |         0 |         1 |   nan       |
        | chr1         |         1 |         4 |     3.32193 |
        | chr1         |         4 |         6 |     3.90689 |
        | chr1         |         6 |         7 |     4.32193 |
        | chr1         |         7 |         8 |     4.64386 |
        | chr1         |         8 |        10 |     4.90689 |
        +--------------+-----------+-----------+-------------+
        Unstranded PyRanges object has 6 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        from pyranges.out import _to_bigwig

        if chromosome_sizes is None:
            chromosome_sizes = pr.data.chromsizes()

        result = _to_bigwig(self, path, chromosome_sizes, rpm, divide, value_col, dryrun)

        if dryrun:
            return result

        if chain:
            return self
        else:
            pass

    def to_csv(self, path=None, sep=",", header=True, compression="infer", chain=False):

        r"""Write to comma- or other value-separated file.

        Parameters
        ----------
        path : str, default None, i.e. return string representation.

            Where to write file.

        sep : str, default ","

            String of length 1. Field delimiter for the output file.

        header : bool, default True

            Write out the column names.

        compression : {‘infer’, ‘gzip’, ‘bz2’, ‘zip’, ‘xz’, None}, default "infer"

            Which compression to use. Uses file extension to infer by default.

        chain: bool, default False

            Whether to return the PyRanges after writing.

        Examples
        --------

        >>> d = {"Chromosome": [1] * 3, "Start": [1, 3, 5], "End": [4, 6, 9], "Feature": ["gene", "exon", "exon"]}
        >>> gr = pr.from_dict(d)
        >>> print(gr.to_csv(sep="\t"))
        Chromosome	Start	End	Feature
        1	1	4	gene
        1	3	6	exon
        1	5	9	exon
        <BLANKLINE>
        """

        from pyranges.out import _to_csv
        result = _to_csv(
            self, path, sep=sep, header=header, compression=compression)
        if path and chain:
            return self
        else:
            return result

    def to_gff3(self, path=None, compression="infer", chain=False):

        """Write to General Feature Format.

        Parameters
        ----------
        path : str, default None, i.e. return string representation.

            Where to write file.

        compression : {‘infer’, ‘gzip’, ‘bz2’, ‘zip’, ‘xz’, None}, default "infer"

            Which compression to use. Uses file extension to infer by default.

        chain: bool, default False

            Whether to return the PyRanges after writing.

        Notes
        -----

        GTF uses a different naming-convention for columns than PyRanges.
        This is the mapping between column names:

        ``{"seqname": "Chromosome", "source": "Source", "type": "Feature", "start": "Start", "end": "End", "score": "Score", "strand": "Strand", "phase": "Frame", "attributes": "Attribute"}``

        All other columns are appended as a field in the attribute string.

        Nonexisting columns will be added with a '.' to represent the missing values.

        See Also
        --------
        pyranges.read_gff3 : read GFF3 files
        pyranges.to_gtf : write to GTF format

        Examples
        --------

        >>> d = {"Chromosome": [1] * 3, "Start": [1, 3, 5], "End": [4, 6, 9], "Feature": ["gene", "exon", "exon"]}
        >>> gr = pr.from_dict(d)
        >>> print(gr.to_gff3())
        1	.	gene	2	4	.	.	.
        1	.	exon	4	6	.	.	.
        1	.	exon	6	9	.	.	.
        <BLANKLINE>

        >>> gr.Gene = [1, 2, 3]
        >>> gr.function = ["a b", "c", "def"]
        >>> print(gr.to_gff3())
        1	.	gene	2	4	.	.	.	Gene=1;function=a b
        1	.	exon	4	6	.	.	.	Gene=2;function=c
        1	.	exon	6	9	.	.	.	Gene=3;function=def
        <BLANKLINE>
        """

        from pyranges.out import _to_gff3

        result = _to_gff3(self, path, compression=compression)

        if path and chain:
            return self
        else:
            return result

    def to_gtf(self, path=None, compression="infer", chain=False):

        """Write to Gene Transfer Format.

        Parameters
        ----------
        path : str, default None, i.e. return string representation.

            Where to write file.

        compression : {‘infer’, ‘gzip’, ‘bz2’, ‘zip’, ‘xz’, None}, default "infer"

            Which compression to use. Uses file extension to infer by default.

        chain: bool, default False

            Whether to return the PyRanges after writing.

        Notes
        -----

        GTF uses a different naming-convention for columns than PyRanges.
        This is the mapping between column names:

        ``{"seqname": "Chromosome", "source": "Source", "feature": "Feature", "start": "Start", "end": "End", "score": "Score", "strand": "Strand", "frame": "Frame", "attribute": "Attribute"}``

        All other columns are appended as a field in the attribute string.

        Nonexisting columns will be added with a '.' to represent the missing values.

        See Also
        --------
        pyranges.read_gtf : read GTF files
        pyranges.to_gff3 : write to GFF3 format

        Examples
        --------

        >>> d = {"Chromosome": [1] * 3, "Start": [1, 3, 5], "End": [4, 6, 9], "Feature": ["gene", "exon", "exon"]}
        >>> gr = pr.from_dict(d)
        >>> print(gr.to_gtf())
        1	.	gene	2	4	.	.	.
        1	.	exon	4	6	.	.	.
        1	.	exon	6	9	.	.	.
        <BLANKLINE>

        >>> gr.name = ["Tim", "Eric", "Endre"]
        >>> gr.prices = ["Cheap", "Premium", "Fine European"]
        >>> print(gr.to_gtf())
        1	.	gene	2	4	.	.	.	name "Tim"; prices "Cheap";
        1	.	exon	4	6	.	.	.	name "Eric"; prices "Premium";
        1	.	exon	6	9	.	.	.	name "Endre"; prices "Fine European";
        <BLANKLINE>
        """

        from pyranges.out import _to_gtf

        result = _to_gtf(self, path, compression=compression)

        if path and chain:
            return self
        else:
            return result


    def to_rle(self, value_col=None, strand=None, rpm=False, nb_cpu=1):

        """Return as RleDict.

        Create collection of Rles representing the coverage or other numerical value.

        Parameters
        ----------
        value_col : str, default None
            Numerical column to create RleDict from.

        strand : bool, default None, i.e. auto
            Whether to treat strands serparately.

        rpm : bool, default False
            Normalize by multiplying with `1e6/(number_intervals)`.

        nb_cpu : int, default 1
            How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
            Will only lead to speedups on large datasets.

        Returns
        -------
        pyrle.RleDict

            Rle with coverage or other info from the PyRanges.

        Examples
        --------

        >>> d = {'Chromosome': ['chr1', 'chr1', 'chr1'], 'Start': [3, 8, 5],
        ...      'End': [6, 9, 7], 'Score': [0.1, 5, 3.14], 'Strand': ['+', '+', '-']}
        >>> gr = pr.from_dict(d)
        >>> gr.to_rle()
        chr1 +
        --
        +--------+-----+-----+-----+-----+
        | Runs   | 3   | 3   | 2   | 1   |
        |--------+-----+-----+-----+-----|
        | Values | 0.0 | 1.0 | 0.0 | 1.0 |
        +--------+-----+-----+-----+-----+
        Rle of length 9 containing 4 elements (avg. length 2.25)
        <BLANKLINE>
        chr1 -
        --
        +--------+-----+-----+
        | Runs   | 5   | 2   |
        |--------+-----+-----|
        | Values | 0.0 | 1.0 |
        +--------+-----+-----+
        Rle of length 7 containing 2 elements (avg. length 3.5)
        RleDict object with 2 chromosomes/strand pairs.

        >>> gr.to_rle(value_col="Score")
        chr1 +
        --
        +--------+-----+-----+-----+-----+
        | Runs   | 3   | 3   | 2   | 1   |
        |--------+-----+-----+-----+-----|
        | Values | 0.0 | 0.1 | 0.0 | 5.0 |
        +--------+-----+-----+-----+-----+
        Rle of length 9 containing 4 elements (avg. length 2.25)
        <BLANKLINE>
        chr1 -
        --
        +--------+-----+------+
        | Runs   | 5   | 2    |
        |--------+-----+------|
        | Values | 0.0 | 3.14 |
        +--------+-----+------+
        Rle of length 7 containing 2 elements (avg. length 3.5)
        RleDict object with 2 chromosomes/strand pairs.

        >>> gr.to_rle(value_col="Score", strand=False)
        chr1
        +--------+-----+-----+------+------+-----+-----+
        | Runs   | 3   | 2   | 1    | 1    | 1   | 1   |
        |--------+-----+-----+------+------+-----+-----|
        | Values | 0.0 | 0.1 | 3.24 | 3.14 | 0.0 | 5.0 |
        +--------+-----+-----+------+------+-----+-----+
        Rle of length 9 containing 6 elements (avg. length 1.5)
        Unstranded RleDict object with 1 chromosome.

        >>> gr.to_rle(rpm=True)
        chr1 +
        --
        +--------+-----+-------------------+-----+-------------------+
        | Runs   | 3   | 3                 | 2   | 1                 |
        |--------+-----+-------------------+-----+-------------------|
        | Values | 0.0 | 333333.3333333333 | 0.0 | 333333.3333333333 |
        +--------+-----+-------------------+-----+-------------------+
        Rle of length 9 containing 4 elements (avg. length 2.25)
        <BLANKLINE>
        chr1 -
        --
        +--------+-----+-------------------+
        | Runs   | 5   | 2                 |
        |--------+-----+-------------------|
        | Values | 0.0 | 333333.3333333333 |
        +--------+-----+-------------------+
        Rle of length 7 containing 2 elements (avg. length 3.5)
        RleDict object with 2 chromosomes/strand pairs.
        """

        if strand is None:
            strand = self.stranded

        from pyranges.methods.to_rle import _to_rle

        return _to_rle(self, value_col, strand=strand, rpm=rpm, nb_cpu=nb_cpu)


    def unstrand(self):

        """Remove strand.

        Note
        ----

        Removes Strand column even if PyRanges is not stranded.

        See Also
        --------

        PyRanges.stranded : whether PyRanges contains valid strand info.

        Examples
        --------

        >>> d =  {'Chromosome': ['chr1', 'chr1'], 'Start': [1, 6],
        ...       'End': [5, 8], 'Strand': ['+', '-']}
        >>> gr = pr.from_dict(d)
        >>> gr
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int32) |   (int32) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         1 |         5 | +            |
        | chr1         |         6 |         8 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.unstrand()
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        | chr1         |         1 |         5 |
        | chr1         |         6 |         8 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 2 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        if not self.stranded and "Strand" in self.columns:
            return self.drop("Strand")
        elif not self.stranded:
            return self

        gr = pr.concat([self["+"], self["-"]])

        gr = gr.apply(lambda df: df.drop("Strand", axis=1).reset_index(drop=
                                                                       True))

        return pr.PyRanges(gr.dfs)


    def values(self):
        """Return the underlying DataFrames."""

        return [df for k, df in self.items() if not df.empty]

    def window(self, window_size, strand=None):

        """Return overlapping genomic windows.

        Windows of length `window_size` are returned.

        Parameters
        ----------
        window_size : int
            Length of the windows.

        strand : bool, default None, i.e. auto

            Whether to do operations on chromosome/strand pairs or chromosomes. If None, will use
            chromosome/strand pairs if the PyRanges is stranded.

        nb_cpu: int, default 1

            How many cpus to use. Can at most use 1 per chromosome or chromosome/strand tuple.
            Will only lead to speedups on large datasets.

        **kwargs
            Additional keyword arguments to pass as keyword arguments to `f`

        Returns
        -------
        PyRanges

            Tiled PyRanges.

        See also
        --------

        pyranges.PyRanges.tile : divide intervals into adjacent tiles.

        Examples
        --------

        >>> gr = pr.from_dict({"Chromosome": [1], "Start": [895], "End": [1259]})
        >>> gr
        +--------------+-----------+-----------+
        |   Chromosome |     Start |       End |
        |   (category) |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        |            1 |       895 |      1259 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 1 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.window(200)
        +--------------+-----------+-----------+
        |   Chromosome |     Start |       End |
        |   (category) |   (int32) |   (int32) |
        |--------------+-----------+-----------|
        |            1 |       895 |      1095 |
        |            1 |      1095 |      1259 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 2 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr = pr.data.ensembl_gtf()[["Feature", "gene_name"]]
        >>> gr
        +--------------+--------------+-----------+-----------+--------------+-------------+
        | Chromosome   | Feature      | Start     | End       | Strand       | gene_name   |
        | (category)   | (category)   | (int32)   | (int32)   | (category)   | (object)    |
        |--------------+--------------+-----------+-----------+--------------+-------------|
        | 1            | gene         | 11868     | 14409     | +            | DDX11L1     |
        | 1            | transcript   | 11868     | 14409     | +            | DDX11L1     |
        | 1            | exon         | 11868     | 12227     | +            | DDX11L1     |
        | 1            | exon         | 12612     | 12721     | +            | DDX11L1     |
        | ...          | ...          | ...       | ...       | ...          | ...         |
        | 1            | gene         | 1173055   | 1179555   | -            | TTLL10-AS1  |
        | 1            | transcript   | 1173055   | 1179555   | -            | TTLL10-AS1  |
        | 1            | exon         | 1179364   | 1179555   | -            | TTLL10-AS1  |
        | 1            | exon         | 1173055   | 1176396   | -            | TTLL10-AS1  |
        +--------------+--------------+-----------+-----------+--------------+-------------+
        Stranded PyRanges object has 2,446 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.window(1000)
        +--------------+--------------+-----------+-----------+--------------+-------------+
        | Chromosome   | Feature      | Start     | End       | Strand       | gene_name   |
        | (category)   | (category)   | (int32)   | (int32)   | (category)   | (object)    |
        |--------------+--------------+-----------+-----------+--------------+-------------|
        | 1            | gene         | 11868     | 12868     | +            | DDX11L1     |
        | 1            | gene         | 12868     | 13868     | +            | DDX11L1     |
        | 1            | gene         | 13868     | 14409     | +            | DDX11L1     |
        | 1            | transcript   | 11868     | 12868     | +            | DDX11L1     |
        | ...          | ...          | ...       | ...       | ...          | ...         |
        | 1            | exon         | 1173055   | 1174055   | -            | TTLL10-AS1  |
        | 1            | exon         | 1174055   | 1175055   | -            | TTLL10-AS1  |
        | 1            | exon         | 1175055   | 1176055   | -            | TTLL10-AS1  |
        | 1            | exon         | 1176055   | 1176396   | -            | TTLL10-AS1  |
        +--------------+--------------+-----------+-----------+--------------+-------------+
        Stranded PyRanges object has 7,516 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        from pyranges.methods.windows import _windows

        if strand is None:
            strand = self.stranded

        kwargs = {"strand": strand}
        kwargs["sparse"] = {"self": False}
        kwargs["window_size"] = window_size

        df = pyrange_apply_single(_windows, self, **kwargs)

        return PyRanges(df)


    def __getstate__(self):
        return self.dfs

    def __setstate__(self, d):
        self.__dict__["dfs"] = d
