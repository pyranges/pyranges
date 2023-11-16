import numpy as np
import pandas as pd

from pyranges import PyRanges
from pyranges.genomicfeatures import GenomicFeaturesMethods
from pyranges.statistics import StatisticsMethods


def set_dtypes(df):
    from pandas.api.types import CategoricalDtype

    dtypes = {
        "Start": np.int64,
        "End": np.int64,
        "Chromosome": "category",
    }

    if "Strand" in df:
        dtypes["Strand"] = CategoricalDtype(
            categories=df["Strand"].drop_duplicates().to_list()
        )

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
    else:
        grpby_key = "Chromosome"

    return {k: v for k, v in df.groupby(grpby_key)}


def check_strandedness(df):
    """Check whether strand contains '.'"""

    if "Strand" not in df:
        return False

    contains_more_than_plus_minus_in_strand_col = False

    if str(df.Strand.dtype) == "category" and (set(df.Strand.cat.categories) - set("+-")):
        contains_more_than_plus_minus_in_strand_col = True
    elif not ((df.Strand == "+") | (df.Strand == "-")).all():
        contains_more_than_plus_minus_in_strand_col = True

        # if contains_more_than_plus_minus_in_strand_col:
        #     logging.warning("Strand contained more symbols than '+' or '-'. Not supported (yet) in PyRanges.")

    return not contains_more_than_plus_minus_in_strand_col


def _init(self, df: pd.DataFrame) -> None:
    if isinstance(df, PyRanges):
        raise Exception("Object is already a PyRange.")

    df = df.copy()

    df = df.reset_index(drop=True)

    stranded = check_strandedness(df)

    df = set_dtypes(df)

    self.__dict__["dfs"] = create_df_dict(df, stranded)

    self.__dict__["features"] = GenomicFeaturesMethods(self)
    self.__dict__["stats"] = StatisticsMethods(self)
