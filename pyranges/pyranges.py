import pandas as pd
from natsort import natsorted

import pyranges as pr

from pyranges.tostring import tostring

from pyranges.methods.intersection import _intersection, _overlap
from pyranges.multithreaded import pyrange_apply, pyrange_apply_single, _slack, _tes, _tss


def fill_kwargs(kwargs):

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

    dfs = None
    gf = None

    def __init__(self,
                 df=None,
                 chromosomes=None,
                 starts=None,
                 ends=None,
                 strands=None,
                 int64=None):

        from pyranges.methods.init import _init

        if df is None and chromosomes is None:
            df = pd.DataFrame(columns="Chromosome Start End".split())

        _init(self, df, chromosomes, starts, ends, strands, int64)

    def __len__(self):
        return sum([len(d) for d in self.values()])

    def __call__(self, f, col=None, strand=None, as_pyranges=True):
        """Apply a function to each chromosome or chromosome/strand pair.

        Example:
            # add 5 to the Score column
            gr.apply(lambda df: df.Score + 5)
            # subset on rows whose sequence does not contain GC.
            gr.apply(lambda df: ~df.Sequence.str.contains("GC"))

        Args:
            f (function): a function which takes a pandas dataframe and modifies it.
            col (string): if the function returns a series, add it to the df with name col (or update
                the col, if it exists). Default: None.
            strand (bool or None): whether or not to do the operation on chromosome/strand pairs. None
                means auto, i.e. use strand if it exists in the PyRanges. Default: None.
            subset (bool): if True, and the return value is a boolean series, use it to subset the
                data. Default: True
            as_pyranges (bool): whether to return as a PyRanges or as a dict. Default: True

        Returns:
            A PyRanges or a dict of objects.
        """

        from pyranges.methods.call import _call

        return _call(
            self, f, col=col, strand=None, subset=True, as_pyranges=True)

    def __getattr__(self, name):

        from pyranges.methods.attr import _getattr

        return _getattr(self, name)

    def __setattr__(self, column_name, column):

        from pyranges.methods.attr import _setattr

        _setattr(self, column_name, column)

    # def __eq__(self, other):

    #     return self.df.equals(other.df)

    def __getitem__(self, val):

        from pyranges.methods.getitem import _getitem

        return _getitem(self, val)

    def __str__(self):

        return tostring(self)

    def __repr__(self):

        return str(self)

    def __iter__(self):

        return iter(self.dfs.items())

    def eval(self, eval_cmd, strand=True, as_pyranges=True, **kwargs):

        f = lambda df: eval(eval_cmd)

        kwargs = fill_kwargs(kwargs)

        result = pyrange_apply_single(f, self, strand, kwargs)

        if not as_pyranges:
            return result
        else:
            return PyRanges(result)

    def overlap(self, other, **kwargs):

        kwargs["sparse"] = {"self": False, "other": True}
        kwargs = fill_kwargs(kwargs)

        dfs = pyrange_apply(_overlap, self, other, **kwargs)

        return PyRanges(dfs)

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

    def nearest(self, other, **kwargs):

        from pyranges.methods.nearest import _nearest

        kwargs = fill_kwargs(kwargs)
        if kwargs.get("how") in "upstream downstream".split():
            assert other.stranded, "If doing upstream or downstream nearest, other pyranges must be stranded"

        dfs = pyrange_apply(_nearest, self, other, **kwargs)

        return PyRanges(dfs)

    def intersect(self, other, **kwargs):

        kwargs = fill_kwargs(kwargs)
        kwargs["sparse"] = {"self": False, "other": True}

        dfs = pyrange_apply(_intersection, self, other, **kwargs)

        return PyRanges(dfs)

    def set_intersect(self, other, **kwargs):

        kwargs = fill_kwargs(kwargs)
        strandedness = kwargs["strandedness"]
        strand = True if strandedness else False
        self_clusters = self.merge(strand=strand, **kwargs)
        other_clusters = other.merge(strand=strand, **kwargs)
        dfs = pyrange_apply(_intersection, self_clusters, other_clusters,
                            **kwargs)

        return PyRanges(dfs)

    def set_union(self, other, **kwargs):

        kwargs = fill_kwargs(kwargs)
        strandedness = kwargs["strandedness"]
        strand = True if strandedness else False

        if not strand:
            self = self.unstrand()
            other = other.unstrand()

        # from pydbg import dbg
        # dbg(self)
        # dbg(other)
        gr = pr.concat([self, other], strand)

        # from pydbg import dbg
        # dbg(gr)

        gr = gr.merge(strand=strand, **kwargs)
        # dbg(gr)

        return gr

    def subtract(self, other, **kwargs):

        from pyranges.methods.subtraction import _subtraction

        kwargs["sparse"] = {"self": False, "other": True}
        kwargs = fill_kwargs(kwargs)
        strandedness = kwargs["strandedness"]

        strand = True if strandedness else False
        other_clusters = other.merge(strand=strand, **kwargs)
        result = pyrange_apply(_subtraction, self, other_clusters, **kwargs)

        return PyRanges(result)

    def join(self, other, **kwargs):

        from pyranges.methods.join import _write_both

        kwargs = fill_kwargs(kwargs)
        dfs = pyrange_apply(_write_both, self, other, **kwargs)

        return PyRanges(dfs)

    def cluster(self, strand=False, **kwargs):

        kwargs = fill_kwargs(kwargs)

        # TODO: implement in Cython. Will lead to 2X speedup.

        mr = self.merge(strand=strand)
        # from pydbg import dbg
        # dbg(mr)

        clusters = list(range(1, len(mr) + 1))
        # dbg(clusters)
        mr.Cluster = clusters
        # dbg(mr)
        j = self.join(mr, how="first", strandedness="same" if strand else None)
        # dbg(j.df)
        j = j.drop(drop="(Start|End|Strand)_b")
        # dbg(j.df)

        return j

    def windows(self, window_size, strand=None, **kwargs):

        from pyranges.methods.windows import _windows

        sparse = True
        if kwargs.get("keep_metadata"):
            sparse = False

        kwargs["sparse"] = {"self": sparse}
        kwargs["window_size"] = window_size

        df = pyrange_apply_single(_windows, self, strand, kwargs)

        return PyRanges(df)

    def merge(self, strand=None, **kwargs):

        from pyranges.methods.merge import _merge

        kwargs["sparse"] = {"self": True}
        df = pyrange_apply_single(_merge, self, strand, kwargs)

        return PyRanges(df)

    def subset(self, function, strand=False, **kwargs):

        kwargs = fill_kwargs(kwargs)

        if strand is None:
            strand = self.stranded

        if self.stranded and not strand:
            self = self.unstrand()

        result = pyrange_apply_single(function, self, strand, kwargs)

        if not result:
            return pr.PyRanges()

        first_result = next(iter(result.values()))

        assert first_result.dtype == bool, "result of subset function must be bool, but is {}".format(
            first_result.dtype)

        return self[result]

    def assign(self, col, function, strand=False, **kwargs):

        kwargs = fill_kwargs(kwargs)

        result = pyrange_apply_single(function, self, strand, kwargs)

        first_result = next(iter(result.values()))

        assert type(
            first_result
        ) == pd.Series, "result of assign function must be Series, but is {}".format(
            type(first_result))

        # do a deepcopy of object
        new_self = pr.PyRanges({k: v.copy() for k, v in self.items()})
        new_self.__setattr__(col, result)

        return new_self

    # def groupby():

    # df2 = df.assign(Group=df.groupby("CpG").ngroup()).sort_values("Group")

    def coverage(self, value_col=None, strand=False, rpm=False):

        from pyranges.methods.coverage import _coverage

        return _coverage(self, value_col, strand=strand, rpm=rpm)

    def apply(self, f, strand=False, as_pyranges=True, kwargs=None):

        if not kwargs:
            kwargs = {}
        kwargs = fill_kwargs(kwargs)

        result = pyrange_apply_single(f, self, strand, kwargs)

        if not as_pyranges:
            return result
        else:
            return PyRanges(result)

    def apply_pair(self, other, f, kwargs=None, strand=False,
                   as_pyranges=True):

        if kwargs is None:
            kwargs = {}
        kwargs = fill_kwargs(kwargs)

        result = pyrange_apply(f, self, other, **kwargs)

        if not as_pyranges:
            return result
        else:
            return PyRanges(result)

    def slack(self, slack):

        kwargs = fill_kwargs({"slack": slack})
        prg = PyRanges(
            pyrange_apply_single(_slack, self, self.stranded, kwargs))

        return prg

    def tssify(self, slack=0):

        kwargs = fill_kwargs({"slack": slack})
        return PyRanges(
            pyrange_apply_single(_tss, self, self.stranded, kwargs))

    def sort(self, columns=("Start", "End"), **kwargs):
        from pyranges.methods.sort import _sort
        kwargs["sparse"] = {"self": False}
        kwargs = fill_kwargs(kwargs)
        return PyRanges(
            pyrange_apply_single(_sort, self, self.stranded, kwargs))

    def tesify(self, slack=0):

        kwargs = fill_kwargs({"slack": slack})
        return PyRanges(
            pyrange_apply_single(_tes, self, self.stranded, kwargs))

    def keys(self):
        return natsorted(self.dfs.keys())

    @property
    def columns(self):
        """Return the list of column names in the dataframes."""

        columns = [list(df.columns) for df in self.values()]
        assert all([c == columns[0] for c in columns[1:]])
        if columns:
            return columns[0]
        else:
            return []

    def drop(self, drop=None, keep=None, drop_strand=False):
        """Drop column(s) from the PyRanges object.

        If no arguments are given, all the columns except Chromosome, Start, End and Strand are
        dropped. To drop Strand, the drop_strand argument needs to be given.

        Args:
            drop (None, iterable or str): An iterable of columns to drop or a string containing a
            substring/regex of the columns to drop.
            keep (None, iterable or str): An iterable of columns to drop or a string containing a
            substring/regex of the columns not to drop.
            drop_strand (bool): Whether or not to drop the Strand column
        """
        from pyranges.methods.drop import _drop
        return _drop(self, drop, keep, drop_strand)

    @property
    def stranded(self):
        keys = self.keys()
        if not len(keys):
            # so that stranded ops work with empty dataframes
            return True

        key = keys[0]

        return isinstance(key, tuple)

    @property
    def strands(self):

        if not self.stranded:
            raise Exception("PyRanges not stranded!")

        return natsorted(set([k[1] for k in self.keys()]))

    @property
    def chromosomes(self):

        if self.stranded:
            return natsorted(set([k[0] for k in self.keys()]))
        else:
            return natsorted(set([k for k in self.keys()]))

    def items(self):
        # from pydbg import dbg
        # dbg(self.dfs.items())

        return natsorted([(k, df) for (k, df) in self.dfs.items()])

    def values(self):

        return [df for k, df in self.items() if not df.empty]

    def unstrand(self):

        # from pydbg import dbg
        # dbg(self.stranded)

        if not self.stranded:
            return self

        gr = pr.concat([self["+"], self["-"]])
        # dbg(gr)

        gr = gr.drop("Strand", drop_strand=True)
        # dbg(gr)
        return gr

    @property
    def df(self):

        return self.as_df()

    @property
    def empty(self):
        return len(self) == 0

    def as_df(self):

        if len(self) == 0:
            return pd.DataFrame()
        elif len(self) == 1:
            return self.values()[0]
        else:
            return pd.concat(self.values())

    def lengths(self):

        lengths = {}
        for k, df in self.items():
            lengths[k] = df.End - df.Start

        return lengths

    # def midpoints(self):

    #     midpoints = {}
    #     for k, df in self.items():
    #         midpoints[k] = (df.End + df.Start) / 2

    #     return midpoints

    def summary(self):

        from pyranges.methods.summary import _summary

        return _summary(self)

    def to_csv(self, path=None, sep=",", header=False):
        from pyranges.out import _to_csv

        return _to_csv(self, path, sep=sep, header=header)

    def to_bed(self, path=None, keep=True):
        from pyranges.out import _to_bed

        return _to_bed(self, path, keep=keep)

    def to_gtf(self, path=None):
        from pyranges.out import _to_gtf

        return _to_gtf(self, path)

    def to_bigwig(self, path, chromosome_sizes, rpm=True):
        from pyranges.out import _to_bigwig

        _to_bigwig(self, path, chromosome_sizes, rpm)
