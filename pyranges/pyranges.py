import pandas as pd
from natsort import natsorted

import pyranges as pr

import pyranges.raymock as ray

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
                 copy_df=False,
                 extended=False):

        from pyranges.methods.init import _init

        _init(self, df, chromosomes, starts, ends, strands, copy_df, extended)

    def __len__(self):
        return sum([len(d) for d in self.values()])

    def __call__(self, f, strand=False):
        return self.apply(f)

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

        gr = pr.concat([self, other])
        gr = gr.merge(strand=strand, **kwargs)

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

    def merge(self, strand=None, **kwargs):

        from pyranges.methods.merge import _merge

        kwargs["sparse"] = {"self": True}
        df = pyrange_apply_single(_merge, self, strand, kwargs)

        return PyRanges(df)

    def coverage(self, value_col=None, strand=False, rpm=False):

        from pyranges.methods.coverage import _coverage

        return _coverage(self, value_col, strand=strand, rpm=rpm)

    def apply(self, f, strand=False, as_pyranges=True, kwargs=None):

        if not kwargs:
            kwargs = {}
        kwargs = fill_kwargs(kwargs)

        f = ray.remote(f)

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

        f = ray.remote(f)

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
    def stranded(self):
        if not len(self.keys()):
            return True

        return len(list(self.keys())[0]) == 2

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

        return natsorted([(k, df) for (k, df) in self.dfs.items()])

    def values(self):

        return [df for k, df in self.items() if not df.empty]

    @property
    def df(self):

        return self.as_df()

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
