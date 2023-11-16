from typing import TYPE_CHECKING, Any, Callable, Dict, List, Tuple, Union

import numpy as np
import pandas as pd
from natsort import natsorted  # type: ignore
from pandas.core.frame import DataFrame

if TYPE_CHECKING:
    from pyranges.pyranges_main import PyRanges


def merge_dfs(df1: DataFrame, df2: DataFrame) -> DataFrame:
    if not df1.empty and not df2.empty:
        return pd.concat([df1, df2], sort=False).reset_index(drop=True)

    elif df1.empty and df2.empty:
        # can this happen?
        return pd.DataFrame()
    elif df1.empty:
        return df2
    else:
        return df1


def process_results(results: List[Any], keys: Union[List[str], List[Tuple[str, str]]]) -> dict:
    results_dict = {k: r for k, r in zip(keys, results) if r is not None}

    try:
        next(iter(results_dict.values()))
    except StopIteration:  # empty collection
        return results_dict

    # An arbitrary operation might make the keys in the dict and df out of sync.
    # This fixes that by having the PyRanges initializer find the correct keys again..
    try:
        if all(isinstance(v, pd.DataFrame) for v in results_dict.values()):
            df = pd.concat(results_dict.values())
            import pyranges as pr

            _results_dict = pr.PyRanges(df).dfs
        else:
            return results_dict
    except (ValueError, TypeError):
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

    return _results_dict


def make_sparse(df: DataFrame) -> DataFrame:
    if "Strand" in df:
        cols = "Chromosome Start End Strand".split()
    else:
        cols = "Chromosome Start End".split()

    return df[cols]


def make_binary_sparse(kwargs: Dict[str, Any], df: DataFrame, odf: DataFrame) -> Tuple[DataFrame, DataFrame]:
    sparse = kwargs.get("sparse")

    if not sparse:
        return df, odf

    if sparse.get("self"):
        df = make_sparse(df)

    if sparse.get("other"):
        odf = make_sparse(odf)

    return df, odf


def make_unary_sparse(kwargs: Dict[str, Any], df: DataFrame) -> DataFrame:
    sparse = kwargs.get("sparse", {}).get("self")

    return make_sparse(df) if sparse else df


def pyrange_apply(
    function: Callable, self: "PyRanges", other: "PyRanges", **kwargs
) -> Union[Dict[Tuple[str, str], Any], Dict[str, Any]]:
    strandedness = kwargs["strandedness"]

    other_strand = {"+": "-", "-": "+"}
    same_strand = {"+": "+", "-": "-"}

    if strandedness == "opposite":
        strand_dict = other_strand
    else:
        strand_dict = same_strand

    assert strandedness in ["same", "opposite", False, None]

    if strandedness:
        assert (
            self.stranded and other.stranded
        ), "Can only do stranded operations when both PyRanges contain strand info"

    results = []

    keys = natsorted(self.dfs.keys())

    dummy = pd.DataFrame(columns="Chromosome Start End".split())

    other_chromosomes = other.chromosomes
    other_dfs = other.dfs

    if strandedness:
        for c, s in self.chromosomes_and_strands:
            df = self._dfs_with_strand[c, s]
            os = strand_dict[s]

            if not (c, os) in other.keys() or len(other[c, os].values()) == 0:
                odf = dummy
            else:
                odf = other[c, os].values()[0]

            df, odf = make_binary_sparse(kwargs, df, odf)

            try:
                result = function(df, odf, **kwargs)
            except TypeError:
                result = function(df, odf)

            results.append(result)

    else:
        if self.stranded and not other.stranded:
            for (c, s), df in self._dfs_with_strand.items():
                if c not in other_chromosomes:
                    odf = dummy
                else:
                    odf = other._dfs_without_strands[c]

                df, odf = make_binary_sparse(kwargs, df, odf)

                try:
                    result = function(df, odf, **kwargs)
                except TypeError:
                    result = function(df, odf)
                results.append(result)

        elif not self.stranded and other.stranded:
            for c, df in self._dfs_without_strand.items():
                if c not in other_chromosomes:
                    odf = dummy
                else:
                    odf1 = other._dfs_with_strand.get((c, "+"), dummy)
                    odf2 = other._dfs_with_strand.get((c, "-"), dummy)

                    odf = merge_dfs(odf1, odf2)

                df, odf = make_binary_sparse(kwargs, df, odf)

                try:
                    result = function(df, odf, **kwargs)
                except TypeError:
                    result = function(df, odf)
                results.append(result)

        elif self.stranded and other.stranded:
            for (c, s), df in self._dfs_with_strand.items():
                if c not in other_chromosomes:
                    odfs = [dummy]
                else:
                    odfp = other_dfs.get((c, "+"), dummy)  # type: ignore
                    odfm = other_dfs.get((c, "-"), dummy)  # type: ignore

                    odfs = [odfp, odfm]

                if len(odfs) == 2:
                    odf = merge_dfs(*odfs)
                elif len(odfs) == 1:
                    odf = odfs[0]
                else:
                    odf = dummy

                df, odf = make_binary_sparse(kwargs, df, odf)

                try:
                    result = function(df, odf, **kwargs)
                except TypeError:
                    result = function(df, odf)
                results.append(result)

        else:
            for c, df in self._dfs_without_strand.items():
                if c not in other_chromosomes:
                    odf = dummy
                else:
                    odf = other._dfs_without_strand[c]

                df, odf = make_binary_sparse(kwargs, df, odf)

                try:
                    result = function(df, odf, **kwargs)
                except TypeError:
                    result = function(df, odf)
                results.append(result)

    return process_results(results, keys)


def pyrange_apply_single(function: Callable, self: "PyRanges", **kwargs) -> Any:
    strand = kwargs["strand"]

    if strand:
        assert self.stranded, "Can only do stranded operation when PyRange contains strand info"

    results = []

    keys: List = []
    if strand:
        for (c, s), df in self._dfs_with_strand.items():  # type: ignore
            kwargs["chromosome"] = c
            _strand = s
            kwargs["strand"] = _strand

            try:
                result = function(df, **kwargs)
            except TypeError:
                result = function(df)
            results.append(result)

        keys = self.keys()

    elif not self.stranded:
        for c, df in self._dfs_without_strand.items():
            kwargs["chromosome"] = c
            assert isinstance(c, str)

            try:
                result = function(df, **kwargs)
            except TypeError:
                result = function(df)

            results.append(result)
            keys.append(c)

    else:
        for c in self.chromosomes:
            assert isinstance(c, str)
            kwargs["chromosome"] = c

            dfs = self[c]

            if len(dfs.keys()) == 2:
                df, df2 = dfs.values()
                # merge strands
                df = merge_dfs(df, df2)
            else:
                df = dfs.values()[0]

            keys.append(c)

            try:
                result = function(df, **kwargs)
            except TypeError:
                result = function(df)

            results.append(result)

    return process_results(results, keys)


def _lengths(df):
    lengths = df.End - df.Start

    return lengths


def _tss(df: DataFrame, **kwargs) -> DataFrame:
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


def _tes(df: DataFrame, **kwargs) -> DataFrame:
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


def _extend(df: DataFrame, **kwargs) -> DataFrame:
    df = df.copy()
    dtype = df.Start.dtype
    slack = kwargs["ext"]

    assert isinstance(slack, (int, dict)), "Extend parameter must be integer or dict, is {}".format(type(slack))

    if isinstance(slack, int):
        df.loc[:, "Start"] = df.Start - slack
        df.loc[df.Start < 0, "Start"] = 0
        df.End = df.End + slack
    else:
        strand = df.Strand.iloc[0]
        five_end_slack = slack.get("5")
        three_end_slack = slack.get("3")

        if five_end_slack and strand == "+":
            df.loc[:, "Start"] -= five_end_slack
        elif five_end_slack and strand == "-":
            df.loc[:, "End"] += five_end_slack

        if three_end_slack and strand == "-":
            df.loc[:, "Start"] -= three_end_slack
        elif three_end_slack and strand == "+":
            df.loc[:, "End"] += three_end_slack

    df = df.astype({"Start": dtype, "End": dtype})

    assert (df.Start < df.End).all(), "Some intervals are negative or zero length after applying extend!"

    return df


def _extend_grp(df: pd.DataFrame, **kwargs):
    df = df.copy()
    dtype = df.Start.dtype
    slack = kwargs["ext"]
    by = kwargs["group_by"]
    g = df.groupby(by)

    assert isinstance(slack, (int, dict)), "Extend parameter must be integer or dict, is {}".format(type(slack))

    minstarts_pos = g.Start.idxmin()
    maxends_pos = g.End.idxmax()

    if isinstance(slack, int):
        df.loc[minstarts_pos, "Start"] = df.Start - slack
        df.loc[df.Start < 0, "Start"] = 0
        df.loc[maxends_pos, "End"] = df.End + slack

    else:
        strand = df.Strand.iloc[0]
        five_end_slack = slack.get("5")
        three_end_slack = slack.get("3")

        if five_end_slack and strand == "+":
            df.loc[minstarts_pos, "Start"] -= five_end_slack
        elif five_end_slack and strand == "-":
            df.loc[maxends_pos, "End"] += five_end_slack

        if three_end_slack and strand == "-":
            df.loc[minstarts_pos, "Start"] -= three_end_slack
        elif three_end_slack and strand == "+":
            df.loc[maxends_pos, "End"] += three_end_slack

    df = df.astype({"Start": dtype, "End": dtype})

    assert (df.Start < df.End).all(), "Some intervals are negative or zero length after applying extend!"

    return df
