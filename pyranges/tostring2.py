import functools
import os
import shutil
from typing import Optional

import natsort  # type: ignore
import pandas as pd

sort_cols = "Start End".split()

GITHUB_ACTIONS = os.environ.get("GITHUB_ACTIONS", False)


def _get_stranded_f(self, half_entries, f, sort=False):
    counter = 0
    dfs = []

    chromosomes = self.chromosomes

    if f == "tail":
        chromosomes = reversed(chromosomes)

    default = pd.DataFrame(columns=self.columns)
    for chromosome in chromosomes:
        plus = self.dfs.get((chromosome, "+"), default)
        minus = self.dfs.get((chromosome, "-"), default)

        if sort:
            plus = plus.sort_values(sort_cols)
            minus = minus.sort_values(sort_cols)

        plus = getattr(plus, f)(half_entries)
        minus = getattr(minus, f)(half_entries)

        df = pd.concat([plus, minus])
        if sort:
            df = df.sort_values(sort_cols)

        counter += len(df)

        dfs.append(df)

        if counter >= half_entries:
            break

    df = pd.concat(dfs)
    # got twice as many entries as needed before sort. Halve here:
    df = getattr(df, f)(half_entries)

    # dfs = {df.Chromosome.iloc[0]: df for df in}
    df = df.reset_index(drop=True)
    df = df.reindex(index=natsort.order_by_index(df.index, natsort.index_natsorted(zip(df.Chromosome))))

    return df


def _get_unstranded_f(self, half_entries, f, sort=False):
    chromosomes = self.chromosomes

    if f == "tail":
        chromosomes = reversed(chromosomes)

    default = pd.DataFrame(columns=self.columns)

    counter = 0
    dfs = []
    for chromosome in chromosomes:
        cdf = self.dfs.get((chromosome), default)
        cdf = getattr(cdf, f)(half_entries)

        if sort:
            cdf = cdf.sort_values(sort_cols)

        dfs.append(cdf)
        counter += len(cdf)

        if counter >= half_entries:
            break

    df = pd.concat(dfs)

    if f == "tail" and len(df.Chromosome.drop_duplicates()) > 1:
        df = df.iloc[::-1]

    return df


def _get_df(self, n, sort):
    half_entries = int(n / 2)

    if len(self) <= n:
        df = self.df.astype(object)
        if sort:
            df = df.sort_values(sort_cols)
    else:
        if self.stranded:
            top = _get_stranded_f(self, half_entries, "head", sort)
            bottom = _get_stranded_f(self, half_entries, "tail", sort)
        else:
            top = _get_unstranded_f(self, half_entries, "head", sort)
            bottom = _get_unstranded_f(self, half_entries, "tail", sort)

        middle = top.head(1)
        # dot_dot_line.loc[:, :] = "..."
        df = pd.concat([top, middle, bottom]).astype(object)

    # df = df.reset_index(drop=True)
    # df = df.reindex(index=natsort.order_by_index(df.index, natsort.index_natsorted(zip(df.Chromosome))))

    return df


def show_pos_merge_position(df):
    # all_dots = df.Start == "..."

    cols_to_drop = "Chromosome Start End".split()
    if "Strand" in df:
        pos = (
            df.Chromosome.astype(str)
            + " "
            + df.Start.astype(str)
            + "-"
            + df.End.astype(str)
            + " "
            + df.Strand.astype(str)
        )
        cols_to_drop.append("Strand")
    else:
        pos = df.Chromosome.astype(str) + " " + df.Start.astype(str) + "-" + df.End.astype(str)

    df = df.drop(cols_to_drop, axis=1)
    df.insert(0, "- Position -", pos)

    # df.loc[all_dots, :] = "..."

    return df


def get_columns_dtypes(self):
    _df = next(iter(self.dfs.values()))
    dtypes = [
        str(d)
        # ["\n(" + str(d) + ")"
        for d in _df.dtypes
    ]  # ["\n(" + d + ")" for d in _df.dtypes]

    columns = list(_df.columns)

    return {c: d for c, d in zip(columns, dtypes)}


def build_header(columns_dtypes):
    header = []
    for c, d in columns_dtypes.items():
        cd = "".join([str(c), "\n(", d, ")"])
        header.append(cd)

    return header


def add_hidden_col_dotdot(df, n_hidden_cols):
    ddd = pd.Series("...", index=df.index)
    ddd.name = "+{}\n...".format(n_hidden_cols)
    df = pd.concat([df, ddd], axis=1)

    return df


def _grow_string_representation(df, columns_dtypes, terminal_width: Optional[int] = None):
    from tabulate import tabulate

    _terminal_width = shutil.get_terminal_size().columns if terminal_width is None else terminal_width
    magic_number = 10  # length of '| ...   |' to append if there are hidden columns

    if len(columns_dtypes) < 15:
        header = build_header(columns_dtypes)
        str_repr = tabulate(df, headers=header, tablefmt="psql", showindex=False)

        table_width = len(str_repr.split("\n", 1)[0])

        if table_width <= _terminal_width:
            return str_repr, []

    header = build_header({k: columns_dtypes[k] for k in columns_dtypes})
    original_header = list(columns_dtypes)
    df.columns = header

    # know that any pyrange will have at least three columns
    build_df = df.get(list(df.columns[:3]))

    total_columns = len(df.columns)

    i = 3
    for i, c in enumerate(df.columns[3:], 3):
        new_build_df = pd.concat([build_df, df[c]], axis=1)

        new_str_repr = tabulate(
            new_build_df, headers=list(new_build_df.columns), tablefmt="psql", showindex=False  # type: ignore
        )

        table_width = len(new_str_repr.split("\n", 1)[0])

        if table_width >= _terminal_width - magic_number:
            break
        else:
            str_repr = new_str_repr
            build_df = new_build_df

    if i < total_columns:
        new_build_df = add_hidden_col_dotdot(build_df, len(original_header[i:]))
        str_repr = tabulate(
            new_build_df, headers=list(new_build_df.columns), tablefmt="psql", showindex=False  # type: ignore
        )

    return str_repr, original_header[i:]


# ugly hack to make doctests pass in github actions
grow_string_representation = functools.partial(
    _grow_string_representation,
    terminal_width=250 if GITHUB_ACTIONS else None,
)


def untraditional_strand_info(self, str_repr_width):
    _ustr = ""
    if "Strand" in self.columns and not self.stranded:
        strands = []
        for _, df in self:
            strands.extend(list(df.Strand.drop_duplicates()))

        untraditional_strands = ["'" + str(s) + "'" for s in set(strands) - set("+-")]
        n_untraditional_strands = len(untraditional_strands)

        if n_untraditional_strands:
            ustr = "Considered unstranded due to these Strand values: {}"
            for i in range(n_untraditional_strands + 1):
                _ustr = ustr.format(", ".join(untraditional_strands[: i + 1]))
                if len(_ustr) > str_repr_width - 20:
                    break

            if i < n_untraditional_strands - 1:
                untraditional_strands = untraditional_strands[: i + 1]
                untraditional_strands.append("...")
                _ustr = ustr.format(", ".join(untraditional_strands)) + " (+ {} more.)".format(
                    n_untraditional_strands - i - 1
                )

    return _ustr


def hidden_columns_info(hidden_columns, str_repr_width):
    n_hidden_cols = len(hidden_columns)
    _hstr = ""
    if n_hidden_cols:
        hstr = str(n_hidden_cols) + " hidden columns: {}"
        for i in range(n_hidden_cols + 1):
            _hstr = hstr.format(", ".join(hidden_columns[:i]))
            if len(_hstr) > str_repr_width:
                break

        if i < n_hidden_cols:
            hidden_columns = hidden_columns[:i]
            hidden_columns.append("...")

            _hstr = hstr.format(", ".join(hidden_columns)) + " (+ {} more.)".format(n_hidden_cols - i)

    return _hstr


def add_text_to_str_repr(self, str_repr, hidden_columns, sort):
    n_intervals = len(self)
    n_chromosomes = len(self.chromosomes)

    stranded = "Stranded" if self.stranded else "Unstranded"
    str1 = "{} PyRanges object has {:,} rows and {:,} columns from {} chromosomes.".format(
        stranded, n_intervals, len(self.columns), n_chromosomes
    )

    str_repr_width = len(str_repr.split("\n", 1)[0])

    hstr = hidden_columns_info(hidden_columns, str_repr_width)

    ustr = untraditional_strand_info(self, str_repr_width)

    order = {
        (True, True): "Chromosome, Start, End and Strand.",
        (True, False): "Chromosome, Start, End and Strand.",
        (False, False): "Chromosome.",
        (False, True): "Chromosome and Strand.",
    }[sort, self.stranded]

    order = "For printing, the PyRanges was sorted on " + order

    str_repr = "\n".join([s for s in [str_repr, str1, order, ustr, hstr] if s])

    return str_repr


def tostring(self, n=8, merge_position=False, formatting=None, sort=False):
    if len(self) == 0:
        return "Empty PyRanges"

    df = _get_df(self, n, sort)

    columns_dtypes = get_columns_dtypes(self)

    if formatting:
        for k, v in formatting.items():
            df[k] = df[k].map(v.format)

    if merge_position:
        df = show_pos_merge_position(df)
        for k in "Chromosome Start End Strand".split():
            if k in columns_dtypes:
                del columns_dtypes[k]

        _columns_dtypes = {}
        _columns_dtypes["- Position -"] = "Multiple types"
        for k, v in columns_dtypes.items():
            _columns_dtypes[k] = v
        columns_dtypes = _columns_dtypes

    df = df.astype(object).reset_index(drop=True)
    if len(self) > n:
        middle = int(n / 2)
        df.loc[middle, :] = "..."

    str_repr, hidden_columns = grow_string_representation(df, columns_dtypes)

    str_repr = add_text_to_str_repr(self, str_repr, hidden_columns, sort)

    return str_repr
