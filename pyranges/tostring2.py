import pandas as pd

sort_cols = "Start End".split()


def _get_stranded_f(self, half_entries, f):

    counter = 0
    dfs = []

    chromosomes = self.chromosomes

    if f == "tail":
        chromosomes = reversed(chromosomes)

    for chromosome in chromosomes:
        plus = self.dfs.get((chromosome, "+"))
        minus = self.dfs.get((chromosome, "-"))

        plus = getattr(plus, f)(half_entries)
        minus = getattr(minus, f)(half_entries)

        df = pd.concat([plus, minus]).sort_values(sort_cols)
        counter += len(df)

        dfs.append(df)

        if counter >= half_entries:
            break

    df = pd.concat(dfs)
    df = getattr(df, f)(half_entries)

    return df


def _get_unstranded_f(self, half_entries, f):

    counter = 0
    dfs = []
    for chromosome, cdf in self:

        cdf = getattr(cdf, f)(half_entries)

        dfs.append(cdf)
        counter += len(cdf)

        if counter >= half_entries:
            break

    df = pd.concat(dfs)
    df = getattr(df, f)(half_entries)

    return df


def _get_df(self, n):

    half_entries = int(n / 2)

    if len(self) <= n:
        df = self.df.astype(object)
    else:
        if self.stranded:
            top = _get_stranded_f(self, half_entries, "head")
            bottom = _get_stranded_f(self, half_entries, "tail")
        else:
            top = _get_unstranded_f(self, half_entries, "head")
            bottom = _get_unstranded_f(self, half_entries, "tail")

        dot_dot_line = top.head(1).astype(object)
        dot_dot_line.loc[:, :] = "..."
        df = pd.concat([top, dot_dot_line, bottom]).astype(object)

    return df


def show_pos_merge_position(df):

    all_dots = df.Start == "..."

    cols_to_drop = "Chromosome Start End".split()
    if "Strand" in df:
        pos = df.Chromosome.astype(str) + " " + df.Start.astype(
            str) + "-" + df.End.astype(str) + " " + df.Strand.astype(str)
        cols_to_drop.append("Strand")
    else:
        pos = df.Chromosome.astype(str) + " " + df.Start.astype(
            str) + "-" + df.End.astype(str)

    df = df.drop(cols_to_drop, axis=1)
    df.insert(0, "- Position -", pos)

    df.loc[all_dots, :] = "..."

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
        cd = "".join([c, "\n(", d, ")"])
        header.append(cd)

    return header


def grow_string_representation(df, columns_dtypes):

    from tabulate import tabulate
    import shutil

    terminal_width = shutil.get_terminal_size().columns
    magic_number = 9  # length of '| ...   |' to append if there are hidden columns

    if len(columns_dtypes) < 15:

        header = build_header(columns_dtypes)
        str_repr = tabulate(
            df, headers=header, tablefmt='psql', showindex=False)

        table_width = len(str_repr.split("\n", 1)[0])

        if table_width <= terminal_width:
            return str_repr, []

    columns, dtypes = [], []

    never_add = "Chromosome Start End Strand".split()
    build_df = df.get(never_add)
    header = build_header({k: columns_dtypes[k] for k in columns_dtypes})
    hidden_columns = []
    total_columns = len(df.columns)
    for i, c in enumerate(df.columns):

        if c in never_add:
            continue

        columns.append(c)
        dtype = columns_dtypes[c]
        dtypes.append()

        _header = c + "\n(" + str(dtype) + ")"
        header.append(_header)

        new_build_df = pd.concat([build_df, df[c]], axis=1)

        new_str_repr = tabulate(
            new_build_df, headers=header, tablefmt='psql', showindex=False)

        table_width = len(new_str_repr.split("\n", 1)[0])

        if table_width >= terminal_width - magic_number:
            hidden_columns = df.columns[i - 1:]
            columns = columns[:-1]
            dtypes = dtypes[:-1]
            break
        else:
            str_repr = new_str_repr
            build_df = new_build_df

    return str_repr, build_df, hidden_columns


def add_text_to_str_repr(self, build_df, hidden_cols):

    _df = next(iter(self.dfs.values()))
    n_intervals = len(self)
    n_chromosomes = len(self.chromosomes)

    n_hidden_cols = len(hidden_cols)
    stranded = "Stranded" if self.stranded else "Unstranded"
    str1 = "{} PyRanges object has {:,} rows and {:,} columns from {} chromosomes.".format(
        stranded, n_intervals, len(self.columns), n_chromosomes)

    if n_hidden_cols:
        # add ... as last col
        ddd = pd.Series("...", index=build_df.index)
        ddd.name = "...\n..."
        df = pd.concat([df, ddd], axis=1)

        # add hidden col info
        columns = list(df.columns)
        first_hidden_idx = columns.index(first_hidden)
        str2 = "Hidden columns: {}".format(", ".join(
            columns[first_hidden_idx:first_hidden_idx + 10]))
        if (n_hidden_cols - 10) > 0:
            str2 += "... (+ {} more.)".format(n_hidden_cols - 10)
            str_repr = "\n".join([str_repr, str1, str2])
        else:
            str_repr = "\n".join([str_repr, str1, str2])
    else:
        str_repr = "\n".join([str_repr, str1])

    if "Strand" in self.columns and not self.stranded:
        strands = []
        for _, df in self:
            strands.extend(list(df.Strand.drop_duplicates()))

        untraditional_strands = set(strands) - set("+-")
        more_than_10 = ", ..." if len(untraditional_strands) > 9 else ""
        from itertools import islice
        untraditional_strands = islice(untraditional_strands, 10)
        str_repr += "\n(Considered unstranded due to the values {}{} being present in the Strand column.)".format(
            ", ".join((str(x) for x in untraditional_strands)), more_than_10)


def tostring(self, n=8, merge_position=False):

    if len(self) == 0:
        return "Empty PyRanges"

    df = _get_df(self, n)

    columns_dtypes = get_columns_dtypes(self)

    if merge_position:
        df = show_pos_merge_position(df)
        columns_dtypes["- Position -"] = "Multiple types"
        [
            columns_dtypes.pop(k, None)
            for k in "Chromosome Start End Strand".split()
        ]

    # print(columns_dtypes)

    str_repr, columns = grow_string_representation(df, columns_dtypes)

    return str_repr


if __name__ == "__main__":

    from pyranges.tostring2 import _get_stranded_f, _get_unstranded_f
    from pyranges.tostring2 import tostring
    import pyranges as pr
    gr = pr.data.chipseq()
    df = gr.df
    _get_stranded_f(gr, 2, "head")
    _get_unstranded_f(gr.unstrand(), 6, "tail")

    tostring(gr)
