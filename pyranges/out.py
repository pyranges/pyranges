import csv

import numpy as np
import pandas as pd
from natsort import natsorted  # type: ignore

_gtf_columns = {
    "seqname": "Chromosome",
    "source": "Source",
    "feature": "Feature",
    "start": "Start",
    "end": "End",
    "score": "Score",
    "strand": "Strand",
    "frame": "Frame",
    # "attribute": "Attribute"  # filled with all others columns
}

_gff3_columns = _gtf_columns.copy()
_gff3_columns["phase"] = _gff3_columns.pop("frame")

_ordered_gtf_columns = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute",
]
_ordered_gff3_columns = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attribute",
]


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


def _bed(df, keep):
    all_columns = "Chromosome Start End Name Score Strand".split()

    outdf = _fill_missing(df, all_columns)

    noncanonical = set(df.columns) - set(all_columns)
    noncanonical = [c for c in df.columns if c in noncanonical]

    if keep:
        return pd.concat([outdf, df.get(noncanonical)], axis=1)
    else:
        return outdf


def _gtf(df, mapping):
    pr_col2gff_col = {v: k for k, v in mapping.items()}

    df = df.rename(columns=pr_col2gff_col)  # copying here
    df.loc[:, "start"] = df.start + 1
    all_columns = _ordered_gtf_columns[:-1]
    columns = list(df.columns)

    outdf = _fill_missing(df, all_columns)

    if "attribute" in df.columns:
        attribute = mapping["attribute"] + ' "' + df.attribute + '";'
    else:
        # gotten all needed columns, need to join the rest
        rest = set(df.columns) - set(all_columns)
        rest = sorted(rest, key=columns.index)
        rest_df = df.get(rest).copy()
        for c in rest_df:
            col = rest_df[c]
            isnull = col.isnull()
            col = col.astype(str).str.replace("nan", "")
            new_val = c + ' "' + col + '";'
            rest_df.loc[:, c] = rest_df[c].astype(str)
            rest_df.loc[~isnull, c] = new_val
            rest_df.loc[isnull, c] = ""

        attribute = rest_df.apply(lambda r: " ".join([v for v in r if v]), axis=1)
    outdf.insert(outdf.shape[1], "attribute", attribute)

    return outdf


def _to_gtf(self, path=None, compression="infer", map_cols=None):
    mapping = _gtf_columns.copy()
    if map_cols:
        mapping.update(map_cols)

    gr = self

    outdfs = [_gtf(v, mapping) for k, v in sorted(gr.dfs.items())]

    if path:
        mode = "w+"
        for outdf in outdfs:
            outdf.to_csv(
                path,
                index=False,
                header=False,
                compression=compression,
                mode=mode,
                sep="\t",
                quoting=csv.QUOTE_NONE,
            )
            mode = "a"
    else:
        return "".join([outdf.to_csv(index=False, header=False, sep="\t", quoting=csv.QUOTE_NONE) for outdf in outdfs])


def _to_csv(self, path=None, sep=",", header=True, compression="infer"):
    gr = self

    if path:
        mode = "w+"
        for _, outdf in natsorted(gr.dfs.items()):
            outdf.to_csv(
                path,
                index=False,
                compression=compression,
                header=header,
                mode=mode,
                sep=sep,
                quoting=csv.QUOTE_NONE,
            )
            mode = "a"
            header = False
    else:
        return "".join(
            [
                outdf.to_csv(index=False, header=header, sep=sep, quoting=csv.QUOTE_NONE)
                for _, outdf in sorted(gr.dfs.items())
            ]
        )


def _to_bed(self, path=None, sep="\t", keep=True, compression="infer"):
    gr = self

    outdfs = natsorted(gr.dfs.items())
    outdfs = [_bed(df, keep) for _, df in outdfs]

    if path:
        mode = "w+"
        for outdf in outdfs:
            outdf.to_csv(
                path,
                index=False,
                header=False,
                compression=compression,
                mode=mode,
                sep="\t",
                quoting=csv.QUOTE_NONE,
            )
            mode = "a"

    else:
        res = "".join([outdf.to_csv(index=False, header=False, sep="\t", quoting=csv.QUOTE_NONE) for outdf in outdfs])
        return res


def _to_bigwig(self, path, chromosome_sizes, rpm=True, divide=False, value_col=None, dryrun=False):
    try:
        import pyBigWig  # type: ignore
    except ModuleNotFoundError:
        print(
            "pybigwig must be installed to create bigwigs. Use `conda install -c bioconda pybigwig` or `pip install pybigwig` to install it."
        )
        import sys

        sys.exit(1)

    if not divide:
        gr = self.to_rle(rpm=rpm, strand=False, value_col=value_col).to_ranges()
    else:
        gr = self.to_rle(rpm=rpm, strand=False, value_col=value_col)
        divide_by = self.to_rle(rpm=rpm, strand=False)
        c = gr / divide_by
        new_pyrles = {}
        for k, v in c.items():
            v.values = np.log2(v.values)
            v.defragment()
            new_pyrles[k] = v

        gr = c.defragment().to_ranges()

    unique_chromosomes = gr.chromosomes

    subset = ["Chromosome", "Start", "End", "Score"]

    gr = gr[subset].unstrand()

    gr = gr.sort()

    if dryrun:
        return gr

    if not isinstance(chromosome_sizes, dict):
        size_df = chromosome_sizes.df
        chromosome_sizes = {k: v for k, v in zip(size_df.Chromosome, size_df.End)}

    header = [(c, int(chromosome_sizes[c])) for c in unique_chromosomes]

    bw = pyBigWig.open(path, "w")
    bw.addHeader(header)

    for chromosome, df in gr:
        chromosomes = df.Chromosome.tolist()
        starts = df.Start.tolist()
        ends = df.End.tolist()
        values = df.Score.tolist()

        bw.addEntries(chromosomes, starts, ends=ends, values=values)


def _to_gff3(self, path=None, compression="infer", map_cols=None):
    mapping = _gff3_columns.copy()
    if map_cols:
        mapping.update(map_cols)

    gr = self

    outdfs = [_gff3(v, mapping) for k, v in sorted(gr.dfs.items())]

    if path:
        mode = "w+"
        for outdf in outdfs:
            outdf.to_csv(
                path,
                index=False,
                header=False,
                compression=compression,
                mode=mode,
                sep="\t",
                quoting=csv.QUOTE_NONE,
            )
            mode = "a"
    else:
        return "".join(
            [outdf.to_csv(index=False, header=False, sep="\t", quoting=csv.QUOTE_MINIMAL) for outdf in outdfs]
        )


def _gff3(df, mapping):
    pr_col2gff_col = {v: k for k, v in mapping.items()}

    df = df.rename(columns=pr_col2gff_col)  # copying here
    df.loc[:, "start"] = df.start + 1
    all_columns = _ordered_gff3_columns[:-1]
    columns = list(df.columns)

    outdf = _fill_missing(df, all_columns)

    if "attribute" in mapping:
        attribute_name = mapping["attribute"]
        attribute_value = df.attribute.iloc[0]
        attribute = f"{attribute_name}={attribute_value}"
    else:
        # gotten all needed columns, need to join the rest
        rest = set(df.columns) - set(all_columns)
        rest = sorted(rest, key=columns.index)
        rest_df = df.get(rest).copy()
        total_cols = rest_df.shape[1]
        for i, c in enumerate(rest_df, 1):
            col = rest_df[c]
            isnull = col.isnull()
            col = col.astype(str).str.replace("nan", "")
            if i != total_cols:
                new_val = c + "=" + col + ";"
            else:
                new_val = c + "=" + col
            rest_df.loc[:, c] = rest_df[c].astype(str)
            rest_df.loc[~isnull, c] = new_val
            rest_df.loc[isnull, c] = ""

        attribute = rest_df.apply(lambda r: "".join([v for v in r if v]), axis=1).str.replace(";$", "", regex=True)
    outdf.insert(outdf.shape[1], "attribute", attribute)

    return outdf


# def _to_bam(df, filename, header=None, chromsizes=None):

#     def _header_from_chromsizes(chromsizes):

#         try:
#             chromsizes = {k: v for k, v in zip(chromsizes.Chromosome, chromsizes.End)}
#         except:
#             pass

#         chromosomes = []
#         for chromosome, length in chromsizes.items():
#             chromosomes.append({"LN": length, "SN": chromosome})

#         return {"SQ": chromosomes}

#     if chromsizes:
#         header = _header_from_chromsizes(chromsizes)


#     import sys

#     try:
#         import bamread
#     except ModuleNotFoundError as e:
#         print("bamread must be installed to write bam. Use `conda install -c bioconda bamread` or `pip install bamread` to install it.")
#         sys.exit(1)
