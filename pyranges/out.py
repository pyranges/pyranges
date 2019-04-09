import pandas as pd
import csv

from natsort import natsorted


def _fill_missing(df, all_columns):

    columns = list(df.columns)

    if not df.get(all_columns) is None:
        outdf = df.get(all_columns)
    else:
        missing = set(all_columns) - set(columns)
        missing_idx = {all_columns.index(m): m for m in missing}
        not_missing = set(columns).intersection(set(all_columns))
        not_missing_ordered = sorted(not_missing, key=all_columns.index)
        outdf = df.get(not_missing_ordered)

        for idx, missing in sorted(missing_idx.items()):
            outdf.insert(idx, missing, ".")

    return outdf


def _bed(df, all_cols=False):

    all_columns = "Chromosome Start End Name Score Strand".split()

    outdf = _fill_missing(df, all_columns)

    noncanonical = set(df.columns) - set(all_columns)

    if all_cols:
        return pd.concat([outdf, df.get(noncanonical)], axis=1)
    else:
        return outdf


def _gtf(df):

    all_columns = "Chromosome   Source   Feature    Start     End       Score    Strand   Frame".split(
    )
    columns = list(df.columns)

    outdf = _fill_missing(df, all_columns)

    # gotten all needed columns, need to join the rest
    rest = set(df.columns) - set(all_columns)
    rest = sorted(rest, key=columns.index)
    rest_df = df.get(rest).copy()
    for c in rest_df:
        col = rest_df[c]
        isnull = col.isnull()
        col = col.astype(str).str.replace("nan", "")
        # dbg(col.head())
        new_val = c + ' "' + col + '";'
        # dbg(new_val)
        rest_df.loc[:, c] = rest_df[c].astype(str)
        rest_df.loc[~isnull, c] = new_val
        rest_df.loc[isnull, c] = ""

    attribute = rest_df.apply(lambda r: " ".join([v for v in r if v]), axis=1)
    outdf.insert(outdf.shape[1], "Attribute", attribute)

    return outdf


class OutMethods():

    pr = None

    def __init__(self, pr):

        self.pr = pr

    def gtf(self, path=None):

        gr = self.pr

        outdfs = [_gtf(v) for k, v in sorted(gr.dfs.items())]

        if path:
            mode = "w+"
            for outdf in outdfs:
                outdf.to_csv(
                    path,
                    index=False,
                    header=False,
                    mode=mode,
                    sep="\t",
                    quoting=csv.QUOTE_NONE)
                mode = "a"
        else:
            return "".join([
                outdf.to_csv(
                    index=False,
                    header=False,
                    sep="\t",
                    quoting=csv.QUOTE_NONE) for outdf in outdfs
            ])

    def csv(self, path=None, sep=","):

        gr = self.pr

        if path:
            mode = "w+"
            for _, outdf in natsorted(gr.dfs.items()):
                outdf.to_csv(
                    path,
                    index=False,
                    header=False,
                    mode=mode,
                    sep=sep,
                    quoting=csv.QUOTE_NONE)
                mode = "a"
        else:
            return "".join([
                outdf.to_csv(
                    index=False, header=False, sep=sep, quoting=csv.QUOTE_NONE)
                for _, outdf in sorted(gr.dfs.items())
            ])

    def bed(self, path=None, sep="\t", all_cols=False):

        gr = self.pr

        outdfs = natsorted(gr.dfs.items())
        outdfs = [_bed(df, all_cols) for _, df in outdfs]

        if path:
            mode = "w+"
            for outdf in outdfs:
                outdf.to_csv(
                    path,
                    index=False,
                    header=False,
                    mode=mode,
                    sep="\t",
                    quoting=csv.QUOTE_NONE)
                mode = "a"

        else:
            return "".join([
                outdf.to_csv(
                    index=False,
                    header=False,
                    sep="\t",
                    quoting=csv.QUOTE_NONE)
                for _, outdf in sorted(gr.dfs.items())
            ])

    def bigwig(self, path, chromosome_sizes, rpm=True):

        gr = self.pr.coverage(rpm=rpm, strand=False).to_ranges()

        unique_chromosomes = gr.chromosomes

        def subset(df):
            return df[['Chromosome', 'Start', 'End', 'Score']]

        gr = gr(subset, strand=False)
        df = gr.sort(strand=False).df

        import pyBigWig

        header = [(c, int(chromosome_sizes[c])) for c in unique_chromosomes]

        bw = pyBigWig.open(path, "w")
        bw.addHeader(header)

        chromosomes = df.Chromosome.tolist()
        starts = df.Start.tolist()
        ends = df.End.tolist()
        values = df.Score.tolist()

        bw.addEntries(chromosomes, starts, ends=ends, values=values)
