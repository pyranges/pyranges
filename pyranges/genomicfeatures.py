# pragma: no cover
import pandas as pd

import pyranges as pr


def pyrange_or_df(func):
    def extension(self, other, *args, **kwargs):
        df = func(self, other, *args, **kwargs)

        if kwargs.get("df"):
            return df

        return pr.PyRanges(df)

    return extension


def pyrange_or_df_single(func):

    # from pyranges import PyRanges
    def extension(self, *args, **kwargs):
        df = func(self, *args, **kwargs)

        if kwargs.get("df"):
            return df

        return pr.PyRanges(df)

    return extension


def _keep_transcript_with_most_exons(df):

    transcripts_with_most_exons = []

    for _, gdf in df.groupby("GeneID"):

        max_exon = gdf.ExonNumber.max()
        max_transcript = gdf.loc[gdf.ExonNumber == max_exon].Transcript.iloc[0]

        max_rows = gdf.loc[gdf.Transcript == max_transcript]

        transcripts_with_most_exons.append(max_rows)

    return pd.concat(transcripts_with_most_exons).reset_index(drop=True)


def tss_or_tes(df, which, slack=0):

    assert which in "tes tss".split()

    if "Feature" not in df:
        raise Exception("No Feature information in object.")

    _df = df[df.Feature == "transcript"]

    if which == "tss":
        _df = _tss(_df, slack)
    elif which == "tes":
        _df = _tes(_df, slack)

    return _df


def filter_transcripts(df, keep="most_exons"):

    return _keep_transcript_with_most_exons(df)


def _tss(df, slack=0, drop_duplicates=True):

    # try:
    #     df = self.df
    # except:
    #     df = self

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

    if drop_duplicates:
        tss = tss.drop_duplicates("Chromosome Start End".split())

    tss.index = range(len(tss))

    return tss


def _tes(df, slack=0, drop_duplicates=True):

    # df = self.df

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

    if drop_duplicates:
        tes = tes.drop_duplicates("Chromosome Start End".split())

    tes.index = range(len(tes))

    return tes


class GenomicFeaturesMethods():

    pr = None

    def __init__(self, pr):

        self.pr = pr

    @pyrange_or_df_single
    def tss(self, transcripts="all", drop_duplicate_tss=True, slack=0):

        pr = self.pr

        df = pr.df

        if transcripts == "all":
            pass
        elif transcripts == "most_exons":
            df = _keep_transcript_with_most_exons(df)

        if not pr.stranded:
            raise Exception(
                "Cannot compute TSSes or TESes without strand info. Perhaps use slack() instead?"
            )
        df = tss_or_tes(df, "tss", slack)

        if drop_duplicate_tss:
            df = df.drop_duplicates("Chromosome Start End".split())

        df = df.drop(["ExonNumber", "ExonID"], 1)

        return df

    @pyrange_or_df_single
    def tes(self, transcripts="all", drop_duplicate_tss=True, slack=0):

        pr = self.pr

        df = pr.df
        df = df.drop("ExonNumber", 1)

        if transcripts == "all":
            pass
        elif transcripts == "most_exons":
            df = _keep_transcript_with_most_exons(df)

        if not pr.stranded:
            raise Exception(
                "Cannot compute TSSes or TESes without strand info. Perhaps use slack() instead?"
            )
        df = tss_or_tes(df, "tes", slack)

        if drop_duplicate_tss:
            df = df.drop_duplicates("Chromosome Start End".split())

        df = df.drop(["ExonNumber", "ExonID"], 1)

        return df
