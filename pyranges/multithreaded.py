import numpy as np

import pandas as pd

import pyranges as pr

from natsort import natsorted

import os


def get_n_args(f):

    import inspect
    nparams = len(inspect.signature(f).parameters)
    return nparams


def call_f(f, nparams, df, odf, kwargs):

    if nparams == 3:
        return f.remote(df, odf, **kwargs)
    else:
        return f.remote(df, odf)


def call_f_single(f, nparams, df, **kwargs):

    if nparams == 2:
        return f.remote(df, **kwargs)
    else:
        return f.remote(df)


class suppress_stdout_stderr(object):
    '''
    A context manager for doing a "deep suppression" of stdout and stderr in
    Python, i.e. will suppress all print, even if the print originates in a
    compiled C/Fortran sub-function.
       This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).

    '''

    def __init__(self):
        # Open a pair of null files
        self.null_fds = [os.open(os.devnull, os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = (os.dup(1), os.dup(2))

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0], 1)
        os.dup2(self.null_fds[1], 2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0], 1)
        os.dup2(self.save_fds[1], 2)
        # Close the null files
        os.close(self.null_fds[0])
        os.close(self.null_fds[1])


def merge_dfs(df1, df2):

    if not df1.empty and not df2.empty:
        return pd.concat([df1, df2], sort=False).reset_index(drop=True)

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
        if isinstance(e, NameError):
            return False

        raise e

    try:
        test_function.remote()
    except Exception as e:
        if "RayConnectionError" in str(type(e)):
            return True
        else:
            raise e


def get_multithreaded_funcs(function, nb_cpu):

    if nb_cpu > 1:
        import ray
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

    nparams = get_n_args(function)
    nb_cpu = kwargs.get("nb_cpu", 1)

    if nb_cpu > 1:
        import ray
        with suppress_stdout_stderr():
            ray.init(num_cpus=nb_cpu, ignore_reinit_error=True)

    function, get, _merge_dfs = get_multithreaded_funcs(
        function, nb_cpu=nb_cpu)

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

    dummy = pd.DataFrame(columns="Chromosome Start End".split())

    other_chromosomes = other.chromosomes
    other_dfs = other.dfs

    if strandedness:

        for (c, s), df in items:

            os = strand_dict[s]

            if not (c, os) in other.keys() or len(other[c, os].values()) == 0:
                odf = dummy
            else:
                odf = other[c, os].values()[0]

            df, odf = make_binary_sparse(kwargs, df, odf)
            result = call_f(function, nparams, df, odf, kwargs)

            results.append(result)

    else:

        if self.stranded and not other.stranded:

            for (c, s), df in items:

                if not c in other_chromosomes:
                    odf = dummy
                else:
                    odf = other_dfs[c]

                df, odf = make_binary_sparse(kwargs, df, odf)
                result = call_f(function, nparams, df, odf, kwargs)
                results.append(result)

        elif not self.stranded and other.stranded:

            for c, df in items:

                if not c in other_chromosomes:
                    odf = dummy
                else:
                    odf1 = other_dfs.get((c, "+"), dummy)
                    odf2 = other_dfs.get((c, "-"), dummy)

                    odf = _merge_dfs.remote(odf1, odf2)

                df, odf = make_binary_sparse(kwargs, df, odf)

                result = call_f(function, nparams, df, odf, kwargs)
                results.append(result)

        elif self.stranded and other.stranded:

            for (c, s), df in self.items():

                if not c in other_chromosomes:
                    odfs = pr.PyRanges(dummy)
                else:
                    odfp = other_dfs.get((c, "+"), dummy)
                    odfm = other_dfs.get((c, "-"), dummy)

                    odfs = [odfp, odfm]

                if len(odfs) == 2:
                    odf = _merge_dfs.remote(*odfs)
                elif len(odfs) == 1:
                    odf = odfs[0]
                else:
                    odf = dummy

                df, odf = make_binary_sparse(kwargs, df, odf)

                result = call_f(function, nparams, df, odf, kwargs)
                results.append(result)

        else:

            for c, df in items:
                if not c in other_chromosomes:
                    odf = dummy
                else:
                    odf = other_dfs[c]

                df, odf = make_binary_sparse(kwargs, df, odf)

                result = call_f(function, nparams, df, odf, kwargs)
                results.append(result)

    results = get(results)

    results = process_results(results, keys)

    if nb_cpu > 1:
        ray.shutdown()

    return results


def pyrange_apply_single(function, self, **kwargs):

    nparams = get_n_args(function)
    nb_cpu = kwargs.get("nb_cpu", 1)
    strand = kwargs["strand"]

    if nb_cpu > 1:
        import ray
        with suppress_stdout_stderr():
            ray.init(num_cpus=nb_cpu, ignore_reinit_error=True)

    function, get, _merge_dfs = get_multithreaded_funcs(
        function, nb_cpu=nb_cpu)

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
            result = call_f_single(function, nparams, df, **kwargs)
            results.append(result)

        keys = self.keys()

    elif not self.stranded:

        keys = []
        for c, df in self.items():

            kwargs["chromosome"] = c

            df = make_unary_sparse(kwargs, df)
            result = call_f_single(function, nparams, df, **kwargs)
            results.append(result)
            keys.append(c)

    else:

        keys = []
        for c in self.chromosomes:

            kwargs["chromosome"] = c

            dfs = self[c]

            if len(dfs.keys()) == 2:
                df, df2 = dfs.values()
                # merge strands
                df = _merge_dfs.remote(df, df2)
            else:
                df = dfs.values()[0]

            df = make_unary_sparse(kwargs, df)
            result = call_f_single(function, nparams, df, **kwargs)
            results.append(result)
            keys.append(c)

    results = get(results)

    if nb_cpu > 1:
        ray.shutdown()

    results = process_results(results, keys)

    return results


def _lengths(df):

    lengths = df.End - df.Start

    return lengths


def _tss(df, **kwargs):

    df = df.copy(deep=True)
    dtype = df.dtypes["Start"]
    slack = kwargs.get("slack", 0)

    starts = np.where(df.Strand == "+", df.Start, df.End - 1)
    ends = starts + slack + 1
    starts = starts - slack
    starts = np.where(starts < 0, 0, starts)

    df.loc[:, "Start"] = starts.astype(dtype)
    df.loc[:, "End"] = ends.astype(dtype)

    return df

def _tes(df, **kwargs):

    df = df.copy(deep=True)
    dtype = df.dtypes["Start"]
    slack = kwargs.get("slack", 0)

    starts = np.where(df.Strand == "+", df.End - 1, df.Start)
    ends = starts + 1 + slack
    starts = starts - slack
    starts = np.where(starts < 0, 0, starts)

    df.loc[:, "Start"] = starts.astype(dtype)
    df.loc[:, "End"] = ends.astype(dtype)

    return df


def _slack(df, **kwargs):

    df = df.copy()
    dtype = df.Start.dtype
    slack = kwargs["slack"]

    assert isinstance(
        slack,
        (int, dict)), "Slack parameter must be integer or dict, is {}".format(
            type(slack))

    if isinstance(slack, int):
        df.loc[:, "Start"] = df.Start - slack
        df.loc[df.Start < 0, "Start"] = 0
        df.End = df.End + slack
    else:
        strand = df.Strand.iloc[0]
        slack_dict = slack
        five_end_slack = slack.get("5")
        three_end_slack = slack.get("3")

        if five_end_slack:
            df.loc[df.Strand == "+", "Start"] -= five_end_slack
            df.loc[df.Strand == "-", "End"] += five_end_slack

        if three_end_slack:
            df.loc[df.Strand == "-", "Start"] -= three_end_slack
            df.loc[df.Strand == "+", "End"] += three_end_slack

    df = df.astype({"Start": dtype, "End": dtype})

    assert (df.Start < df.End).all(), "Some intervals are negative or zero length after applying slack!"

    return df


def pyrange_apply_chunks(function, self, as_pyranges, **kwargs):

    nparams = get_n_args(function)
    nb_cpu = kwargs.get("nb_cpu", 1)
    if nb_cpu > 1:
        import ray
        with suppress_stdout_stderr():
            ray.init(num_cpus=nb_cpu, ignore_reinit_error=True)

    function, get, _merge_dfs = get_multithreaded_funcs(
        function, nb_cpu=nb_cpu)

    keys = []
    lengths = []
    results = []
    for k, v in self.items():
        dfs = np.array_split(v, nb_cpu)
        lengths.append(len(dfs))
        results.extend(
            [call_f_single(function, nparams, df, **kwargs) for df in dfs])
        keys.append(k)

    results = get(results)

    _results = []
    start = 0
    for _, length in zip(keys, lengths):
        end = start + length
        _r = results[start:end]

        if as_pyranges:
            _results.append(pd.concat(_r))
        else:
            _results.append(_r)

        start = end

    results = _results
    if nb_cpu > 1:
        ray.shutdown()

    results = process_results(results, keys)

    return results
