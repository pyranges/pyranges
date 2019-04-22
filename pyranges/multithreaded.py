import pandas as pd

import pyranges as pr

from natsort import natsorted

import sys


def merge_dfs(df1, df2):

    if not df1.empty and not df2.empty:
        return pd.concat([df1, df2], sort=False)

    elif df1.empty and df2.empty:
        # can this happen?
        return None
    elif df1.empty:
        return df2
    else:
        return df1


def process_results(results, keys):

    results_dict = {k: r for k, r in zip(keys, results) if r is not None}

    try:
        first_item = next(iter(results_dict.values()))
    except StopIteration:  # empty collection
        return results_dict

    if not isinstance(first_item, pd.DataFrame):
        return results_dict

    to_delete = []
    # to ensure no duplicate indexes and no empty dataframes
    for k in results_dict:
        if results_dict[k] is None or results_dict[k].empty:
            to_delete.append(k)
        else:
            # pandas might make a df that is not always C-contiguous
            # copying fixes this
            # TODO: better to only fix columns that are not C-contiguous?
            results_dict[k] = results_dict[k].copy(deep=True)
            results_dict[k].index = range(len(results_dict[k]))

    for k in to_delete:
        del results_dict[k]

    return results_dict


def call_f(f, df, odf, kwargs):

    import inspect
    nparams = len(inspect.signature(f).parameters)

    if nparams == 3:
        return f.remote(df, odf, kwargs)
    else:
        return f.remote(df, odf)


def make_sparse(df):

    if "Strand" in df:
        cols = "Chromosome Start End Strand".split()
    else:
        cols = "Chromosome Start End".split()

    return df[cols]


def make_binary_sparse(kwargs, df, odf):

    sparse = kwargs.get("sparse")

    if not sparse:
        return df, odf

    if sparse.get("self"):

        df = make_sparse(df)

    if sparse.get("other"):

        odf = make_sparse(odf)

    return df, odf


def make_unary_sparse(kwargs, df):

    sparse = kwargs.get("sparse").get("self")

    if sparse:
        df = make_sparse(df)

    return df


def ray_initialized():
    def test_function():
        pass

    try:
        test_function = ray.remote(test_function)
    except Exception as e:
        if type(e) == NameError:
            return False

        raise e

    try:
        test_function.remote()
    except Exception as e:
        if "RayConnectionError" in str(type(e)):
            return True
        else:
            raise e


def get_multithreaded_funcs(function):

    if ray_initialized():
        _merge_dfs = ray.remote(merge_dfs)
        get = ray.get
        function = ray.remote(function)
    else:
        _merge_dfs = lambda: "dummy value"
        _merge_dfs.remote = merge_dfs
        get = lambda x: x
        function.remote = function

    return function, get, _merge_dfs


def pyrange_apply(function, self, other, **kwargs):

    function, get, _merge_dfs = get_multithreaded_funcs(function)

    strandedness = kwargs["strandedness"]

    other_strand = {"+": "-", "-": "+"}
    same_strand = {"+": "+", "-": "-"}

    if strandedness == "opposite":
        strand_dict = other_strand
    else:
        strand_dict = same_strand

    assert strandedness in ["same", "opposite", False, None]

    if strandedness:
        assert self.stranded and other.stranded, \
            "Can only do stranded operations when both PyRanges contain strand info"

    results = []

    items = natsorted(self.dfs.items())
    keys = natsorted(self.dfs.keys())

    if strandedness:

        for (c, s), df in items:

            os = strand_dict[s]

            if not (c, os) in other.keys() or len(other[c, os].values()) == 0:
                odf = pd.DataFrame(columns="Chromosome Start End".split())
            else:
                odf = other[c, os].values()[0]

            df, odf = make_binary_sparse(kwargs, df, odf)

            result = call_f(function, df, odf, kwargs)
            results.append(result)

    else:

        if self.stranded and not other.stranded:

            for (c, s), df in items:

                if not c in other.chromosomes:
                    odf = pd.DataFrame(columns="Chromosome Start End".split())
                else:
                    odf = other.dfs[c]

                df, odf = make_binary_sparse(kwargs, df, odf)
                result = call_f(function, df, odf, kwargs)
                results.append(result)

        elif not self.stranded and other.stranded:

            for c, df in items:

                if not c in other.chromosomes:
                    odf = pd.DataFrame(columns="Chromosome Start End".split())
                else:
                    odf1 = other[c, "+"].df
                    odf2 = other[c, "-"].df

                    odf = _merge_dfs.remote(odf1, odf2)

                df, odf = make_binary_sparse(kwargs, df, odf)

                result = call_f(function, df, odf, kwargs)
                results.append(result)

        elif self.stranded and other.stranded:

            for (c, s), df in self.items():

                if not c in other.chromosomes:
                    odfs = pr.PyRanges(
                        pd.DataFrame(columns="Chromosome Start End".split()))
                else:
                    odfs = other[c].values()

                # from pydbg import dbg
                # dbg(odfs)

                if len(odfs) == 2:
                    odf = _merge_dfs.remote(*odfs)
                elif len(odfs) == 1:
                    odf = odfs[0]
                else:
                    odf = pd.DataFrame(columns="Chromosome Start End".split())

                df, odf = make_binary_sparse(kwargs, df, odf)

                # dbg(df)
                # dbg(odf)

                result = call_f(function, df, odf, kwargs)
                results.append(result)

        else:

            for c, df in items:
                if not c in other.chromosomes:
                    odf = pd.DataFrame(columns="Chromosome Start End".split())
                else:
                    odf = other.dfs[c]

                df, odf = make_binary_sparse(kwargs, df, odf)

                result = call_f(function, df, odf, kwargs)
                results.append(result)

    results = get(results)

    results = process_results(results, keys)

    return results


def call_f_single(f, df, kwargs):

    import inspect
    nparams = len(inspect.signature(f).parameters)

    if nparams == 2:
        return f.remote(df, kwargs)
    else:
        return f.remote(df)


def pyrange_apply_single(function, self, strand, kwargs):

    function, get, _merge_dfs = get_multithreaded_funcs(function)

    if strand:
        assert self.stranded, \
            "Can only do stranded operation when PyRange contains strand info"

    results = []

    if strand:

        for (c, s), df in self.items():

            kwargs["chromosome"] = c
            _strand = s
            kwargs["strand"] = _strand

            df = make_unary_sparse(kwargs, df)
            result = call_f_single(function, df, kwargs)
            results.append(result)

        keys = self.keys()

    elif not self.stranded:

        keys = []
        for c, df in self.items():

            kwargs["chromosome"] = c

            df = make_unary_sparse(kwargs, df)
            result = call_f_single(function, df, kwargs)
            results.append(result)
            keys.append(c)

    else:

        keys = []
        for c in self.chromosomes:

            kwargs["chromosome"] = c

            dfs = self[c]

            if len(dfs.keys()) == 2:
                df1, df2 = dfs.values()
                # merge strands
                df1 = _merge_dfs.remote(df1, df2)
            else:
                df1 = dfs.values()[0]

            df1 = make_unary_sparse(kwargs, df1)
            result = call_f_single(function, df1, kwargs)
            results.append(result)
            keys.append(c)

    results = get(results)

    results = process_results(results, keys)

    return results


def _lengths(df):

    lengths = df.End - df.Start

    return lengths


def _tss(df, kwargs):

    slack = kwargs["slack"]

    tss_pos = df.loc[df.Strand == "+"]

    tss_neg = df.loc[df.Strand == "-"].copy()

    # pd.options.mode.chained_assignment = None
    tss_neg.loc[:, "Start"] = tss_neg.End

    # pd.options.mode.chained_assignment = "warn"
    tss = pd.concat([tss_pos, tss_neg], sort=False)
    tss["End"] = tss.Start
    tss.End = tss.End + 1 + slack
    tss.Start = tss.Start - slack
    tss.loc[tss.Start < 0, "Start"] = 0

    return tss.reindex(df.index)


def _tes(df, kwargs):

    slack = kwargs["slack"]

    tes_pos = df.loc[df.Strand == "+"]

    tes_neg = df.loc[df.Strand == "-"].copy()

    # pd.options.mode.chained_assignment = None
    tes_neg.loc[:, "Start"] = tes_neg.End

    # pd.options.mode.chained_assignment = "warn"
    tes = pd.concat([tes_pos, tes_neg], sort=False)
    tes["Start"] = tes.End
    tes.End = tes.End + 1 + slack
    tes.Start = tes.Start - slack
    tes.loc[tes.Start < 0, "Start"] = 0

    return tes.reindex(df.index)


def _slack(df, kwargs):

    slack = kwargs["slack"]
    # df = self.df.copy()
    df.Start = df.Start - slack
    df.loc[df.Start < 0, "Start"] = 0
    df.End = df.End + slack

    return df
