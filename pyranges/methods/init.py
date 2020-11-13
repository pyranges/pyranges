import sys
import numpy as np
import pandas as pd
from natsort import natsorted

from pyranges.statistics import StatisticsMethods
from pyranges.genomicfeatures import GenomicFeaturesMethods
from pyranges import PyRanges
from pyranges.helpers import single_value_key, get_key_from_df


def set_dtypes(df, int64):

    # if extended is None:
    #     extended = False if df.Start.dtype == np.int32 else True

    if not int64:
        dtypes = {
            "Start": np.int32,
            "End": np.int32,
            "Chromosome": "category",
            "Strand": "category",
        }
    else:
        dtypes = {
            "Start": np.int64,
            "End": np.int64,
            "Chromosome": "category",
            "Strand": "category",
        }

    if "Strand" not in df:
        del dtypes["Strand"]

    # need to ascertain that object columns do not consist of multiple types
    # https://github.com/biocore-ntnu/epic2/issues/32
    for column in "Chromosome Strand".split():
        if column not in df:
            continue
        df[column] = df[column].astype(str)

    for col, dtype in dtypes.items():
        if df[col].dtype.name != dtype:
            df[col] = df[col].astype(dtype)
    return df


def create_df_dict(df, stranded):

    chrs = df.Chromosome.cat.remove_unused_categories()

    df["Chromosome"] = chrs

    if stranded:
        grpby_key = "Chromosome Strand".split()
        df["Strand"] = df.Strand.cat.remove_unused_categories()
    else:
        grpby_key = "Chromosome"

    return {k: v for k, v in df.groupby(grpby_key)}


def create_pyranges_df(chromosomes, starts, ends, strands=None):

    if isinstance(chromosomes, str) or isinstance(chromosomes, int):
        chromosomes = pd.Series([chromosomes] * len(starts), dtype="category")

    if strands is not None:

        if isinstance(strands, str):
            strands = pd.Series([strands] * len(starts), dtype="category")

        columns = [chromosomes, starts, ends, strands]
        lengths = list(str(len(s)) for s in columns)
        assert (
            len(set(lengths)) == 1
        ), "chromosomes, starts, ends and strands must be of equal length. But are {}".format(
            ", ".join(lengths)
        )
        colnames = "Chromosome Start End Strand".split()
    else:
        columns = [chromosomes, starts, ends]
        lengths = list(str(len(s)) for s in columns)
        assert (
            len(set(lengths)) == 1
        ), "chromosomes, starts and ends must be of equal length. But are {}".format(
            ", ".join(lengths)
        )
        colnames = "Chromosome Start End".split()

    idx = range(len(starts))
    series_to_concat = []
    for s in columns:
        if isinstance(s, pd.Series):
            s = pd.Series(s.values, index=idx)
        else:
            s = pd.Series(s, index=idx)

        series_to_concat.append(s)

    df = pd.concat(series_to_concat, axis=1)
    df.columns = colnames

    return df


def check_strandedness(df):
    """Check whether strand contains '.'"""

    if "Strand" not in df:
        return False

    contains_more_than_plus_minus_in_strand_col = False

    if str(df.Strand.dtype) == "category" and (
        set(df.Strand.cat.categories) - set("+-")
    ):
        contains_more_than_plus_minus_in_strand_col = True
    elif not ((df.Strand == "+") | (df.Strand == "-")).all():
        contains_more_than_plus_minus_in_strand_col = True

        # if contains_more_than_plus_minus_in_strand_col:
        #     logging.warning("Strand contained more symbols than '+' or '-'. Not supported (yet) in PyRanges.")

    return not contains_more_than_plus_minus_in_strand_col


def _init(
    self,
    df=None,
    chromosomes=None,
    starts=None,
    ends=None,
    strands=None,
    int64=False,
    copy_df=True,
):
    # TODO: add categorize argument with dict of args to categorize?

    if isinstance(df, PyRanges):
        raise Exception("Object is already a PyRange.")

    if isinstance(df, pd.DataFrame):
        assert all(
            c in df for c in "Chromosome Start End".split()
        ), "The dataframe does not have all the columns Chromosome, Start and End."
        if copy_df:
            df = df.copy()

    if df is False or df is None:
        df = create_pyranges_df(chromosomes, starts, ends, strands)

    if isinstance(df, pd.DataFrame):

        df = df.reset_index(drop=True)

        stranded = check_strandedness(df)

        df = set_dtypes(df, int64)

        self.__dict__["dfs"] = create_df_dict(df, stranded)

    # df is actually dict of dfs
    else:

        empty_removed = {k: v.copy() for k, v in df.items() if not v.empty}

        _single_value_key = True
        _key_same = True
        _strand_valid = True
        _has_strand = True
        for key, df in empty_removed.items():

            _key = get_key_from_df(df)
            _single_value_key = single_value_key(df) and _single_value_key
            _key_same = (_key == key) and _key_same

            if isinstance(_key, tuple):
                _strand_valid = _strand_valid and (_key[1] in ["+", "-"])
            else:
                _has_strand = False

        if not all([_single_value_key, _key_same, _strand_valid]):
            df = pd.concat(empty_removed.values()).reset_index(drop=True)

            if _has_strand and _strand_valid:
                empty_removed = df.groupby(["Chromosome", "Strand"])
            else:
                empty_removed = df.groupby("Chromosome")

            empty_removed = {k: v for (k, v) in empty_removed}

        self.__dict__["dfs"] = empty_removed

    self.__dict__["features"] = GenomicFeaturesMethods(self)
    self.__dict__["stats"] = StatisticsMethods(self)
