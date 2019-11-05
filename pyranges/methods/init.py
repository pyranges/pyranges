import sys
import numpy as np
import pandas as pd
from natsort import natsorted

from pyranges.statistics import StatisticsMethods
from pyranges.genomicfeatures.genomicfeatures import GenomicFeaturesMethods
from pyranges import PyRanges


def set_dtypes(df, int64):

    # if extended is None:
    #     extended = False if df.Start.dtype == np.int32 else True

    if not int64:
        dtypes = {
            "Start": np.int32,
            "End": np.int32,
            "Chromosome": "category",
            "Strand": "category"
        }
    else:
        dtypes = {
            "Start": np.int64,
            "End": np.int64,
            "Chromosome": "category",
            "Strand": "category"
        }

    if not "Strand" in df:
        del dtypes["Strand"]

    # need to ascertain that object columns do not consist of multiple types
    # https://github.com/biocore-ntnu/epic2/issues/32
    for column in "Chromosome Strand".split():
        if not column in df:
            continue

        if df[column].dtype == object and len(
                df[column].apply(type).drop_duplicates()) > 1:
            df[column] = df[column].astype(str)
        elif df[column].dtype != object:
            df[column] = df[column].astype(str)

    for col, dtype in dtypes.items():

        if df[col].dtype.name != dtype:

            df[col] = df[col].astype(dtype)

    return df


def create_df_dict(df, stranded):

    chrs = df.Chromosome.cat.remove_unused_categories()

    # with Debug(locals()):
    #     f = df.head(1).head(1)

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
        assert len(
            set(lengths)
        ) == 1, "chromosomes, starts, ends and strands must be of equal length. But are {}".format(
            ", ".join(lengths))
        colnames = "Chromosome Start End Strand".split()
    else:
        columns = [chromosomes, starts, ends]
        lengths = list(str(len(s)) for s in columns)
        assert len(
            set(lengths)
        ) == 1, "chromosomes, starts and ends must be of equal length. But are {}".format(
            ", ".join(lengths))
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
            set(df.Strand.cat.categories) - set("+-")):
        contains_more_than_plus_minus_in_strand_col = True
    elif not ((df.Strand == "+") | (df.Strand == "-")).all():
        contains_more_than_plus_minus_in_strand_col = True

        # if contains_more_than_plus_minus_in_strand_col:
        #     logging.warning("Strand contained more symbols than '+' or '-'. Not supported (yet) in PyRanges.")

    return not contains_more_than_plus_minus_in_strand_col


def _init(self,
          df=None,
          chromosomes=None,
          starts=None,
          ends=None,
          strands=None,
          int64=False,
          copy_df=True):
    # TODO: add categorize argument with dict of args to categorize?

    if isinstance(df, PyRanges):
        raise Exception("Object is already a PyRange.")

    if isinstance(df, pd.DataFrame):
        assert all(
            c in df for c in "Chromosome Start End".split()
        ), "The dataframe does not have all the columns Chromosome, Start and End."
        df = df.copy()

    if df is False or df is None:
        df = create_pyranges_df(chromosomes, starts, ends, strands)

    if isinstance(df, pd.DataFrame):

        stranded = check_strandedness(df)

        df = set_dtypes(df, int64)

        # below is not a good idea! then gr["chr1"] might change the dtypes of a gr!
        # elif isinstance(df, dict):
        #     df = {k: set_dtypes(v, extended) for k, v in df.items()}

        self.__dict__["dfs"] = create_df_dict(df, stranded)

    # df is actually dict of dfs
    else:

        empty_removed = {k: v.copy() for k, v in df.items() if not v.empty}
        if empty_removed:
            first_key, first_df = list(empty_removed.items())[0]
            # from pydbg import dbg;
            # dbg(first_df)
            stranded = "Strand" in first_df

            all_strands_valid = True
            if stranded:
                all_strands_valid = all([
                    len(set(df.Strand.drop_duplicates()) -
                        set(["+", "-"])) == 0 for df in empty_removed.values()
                ])

            assert all(c in first_df for c in "Chromosome Start End".split(
            )), "Columns Chromosome, Start and End must be in the dataframe!"

            # if not has strand key, but is stranded, need to add strand key
            has_strand_key = isinstance(first_key, tuple)
            if not has_strand_key and stranded and all_strands_valid:
                new_dfs = {}
                for k, v in empty_removed.items():
                    for s, sdf in v.groupby("Strand"):
                        new_dfs[k, s] = sdf
                empty_removed = new_dfs

            # need to merge strand keys if not strands valid anymore
            elif has_strand_key and (not all_strands_valid or not stranded):
                new_dfs = {}
                cs = set([k[0] for k in empty_removed.keys()])
                for c in natsorted(cs):
                    dfs = [empty_removed.get((c, s)) for s in "+-"]
                    new_dfs[c] = pd.concat(
                        [df for df in dfs if not df is None]).reset_index(drop=True)
                empty_removed = new_dfs

        self.__dict__["dfs"] = empty_removed

    self.__dict__["features"] = GenomicFeaturesMethods(self)
    self.__dict__["stats"] = StatisticsMethods(self)
