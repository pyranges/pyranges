from typing import Any, Dict, List, Tuple, Union

import pandas as pd
from ncls import NCLS  # type: ignore
from numpy import int64
from pandas.core.frame import DataFrame

from pyranges.pyranges_main import PyRanges


def create_ncls(df: DataFrame) -> NCLS:
    return NCLS(df.Start.values, df.End.values, df.index.values)


def find_overlaps(df: DataFrame, start: int, end: Union[int64, int]) -> List[Union[int, Any]]:
    n = create_ncls(df)

    idxes = []
    for r in n.find_overlap(start, end):
        idxes.append(r[2])

    return idxes


def get_slice(self: PyRanges, val: slice) -> Union[Dict[str, DataFrame], Dict[Tuple[str, str], DataFrame]]:
    # 100:999

    if self.stranded:
        sd = {}
        for sk, sdf in self._dfs_with_strand.items():
            start = val.start or 0
            stop = val.stop or max(sdf.End.max(), start)
            idxes = find_overlaps(sdf, start, stop)
            sd[sk] = sdf.reindex(idxes)
        return sd
    else:
        d = {}
        for k, df in self._dfs_without_strand.items():
            start = val.start or 0
            stop = val.stop or max(df.End.max(), start)
            idxes = find_overlaps(df, start, stop)
            d[k] = df.reindex(idxes)
        return d


def get_string(self: PyRanges, val: str) -> Union[Dict[Tuple[str, str], DataFrame], Dict[str, DataFrame]]:
    if val in self.chromosomes:
        if self.stranded:
            return {k: df for k, df in self._dfs_with_strand.items() if k[0] == val}
        else:
            return {val: df for k, df in self._dfs_without_strand.items() if k == val}
    elif val in "+ -".split():
        return {k: v for k, v in self._dfs_with_strand.items() if k[1] == val}
    else:
        d: Dict[str, DataFrame] = {}
        return d


def get_2_tuple(
    self: PyRanges, first: str, second: Union[str, slice]
) -> Union[Dict[str, DataFrame], Dict[Tuple[str, str], DataFrame]]:
    if isinstance(first, str) and first in "+-" and isinstance(second, slice):
        return get_strand_and_slice(self, strand=first, loc=second)
    if isinstance(first, (int, str)) and isinstance(second, str):
        return get_chromosome_and_strand(self, chromosome=first, strand=second)
    if isinstance(first, (int, str)) and isinstance(second, slice):
        return get_chromosome_and_slice(self, chromosome=first, loc=second)
    else:
        raise TypeError(f"Incorrect types: {type(first)}, {type(second)}")


def get_chromosome_and_slice(
    self: PyRanges, chromosome: str, loc: slice
) -> Union[Dict[str, DataFrame], Dict[Tuple[str, str], DataFrame]]:
    if chromosome in self.chromosomes:
        start = loc.start or 0
        if self.stranded:
            dfs = [df for (c, _), df in self._dfs_with_strand.items() if c == chromosome]
        else:
            dfs = [df for c, df in self._dfs_without_strand.items() if c == chromosome]
        max_end = max([df.End.max() for df in dfs])

        # in case 1:None
        stop = loc.stop or max(max_end, start)

        out_dfs = [df[find_overlaps(df, start, stop)] for df in dfs]

    return PyRanges(pd.concat(out_dfs)).dfs


def get_strand_and_slice(self: PyRanges, strand: str, loc: slice) -> Dict[Tuple[str, str], DataFrame]:
    # "+", 5:10
    start = loc.start or 0

    dfs = [df for (c, s), df in self._dfs_with_strand.items() if s == strand]
    max_end = max([df.End.max() for df in dfs])

    stop = loc.stop or max(max_end, start)

    out_dfs = [df[find_overlaps(df, start, stop)] for df in dfs]

    return {k: v for k, v in PyRanges(pd.concat(out_dfs))._dfs_with_strand.items()}


# "chr1", "+"
def get_chromosome_and_strand(
    self: PyRanges, chromosome: Union[int, str], strand: str
) -> Dict[Tuple[str, str], DataFrame]:
    return {k: df for k, df in self._dfs_with_strand.items() if k == (chromosome, strand)}


def get_chromosome_strand_loc(
    self: PyRanges, chromosome: str, strand: str, loc: slice
) -> Dict[Tuple[str, str], DataFrame]:
    # "chr1", "+", 5:10
    start = loc.start or 0

    if strand not in "+ -".split():
        raise Exception("Strand '{}' invalid.".format(strand))

    r = self[chromosome, strand].values()
    if len(r):
        df = r[0]
    else:
        return {}

    max_end = df.End.max()

    stop = loc.stop or max(max_end, start)

    idxes = find_overlaps(df, start, stop)
    return {(chromosome, strand): df.reindex(idxes)}


def get_booldict(self, df):
    _overlapping = set(self.dfs.keys()).intersection(set(df.keys()))

    new_dfs = {}
    for k in _overlapping:
        new_dfs[k] = self.dfs[k][df[k]]

    return new_dfs
