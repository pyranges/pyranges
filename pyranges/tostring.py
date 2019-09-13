import shutil

import pandas as pd

from tabulate import tabulate


def get_terminal_size():

    return shutil.get_terminal_size().columns


def show_pos_merge_position(self, df, first_df):

    if "Strand" in df:
        pos = df.Chromosome.astype(str) + " " + df.Start.astype(
            str) + "-" + df.End.astype(str) + " " + df.Strand.astype(str)
    else:
        pos = df.Chromosome.astype(str) + " " + df.Start.astype(
            str) + "-" + df.End.astype(str)

    return pos


def increase_string_width(self, df, first_df, merge_position):

    # TODO: increase until no longer fits the screen
    # instead of decreasing. Takes long time with 60k+ cols

    stranded = self.stranded
    n_intervals = len(self)
    n_chromosomes = len(self.chromosomes)

    terminal_width = get_terminal_size()

    columns = []
    all_dtypes = first_df.dtypes
    if not merge_position:
        if "Strand" in df:
            # if the df is stranded, always include enough cols to show strand
            for c in df.columns:
                columns.append(c)
                if c == "Strand":
                    break
            build_df = df.get(columns)
            never_add = columns[:]
        else:
            columns = "Chromosome Start End".split()
            build_df = df.get(columns)
            never_add = "Chromosome Start End Strand".split()

        dtypes = []
        for c in columns:
            dtype = all_dtypes[c]
            dtypes.append(dtype)
    else:
        build_df = show_pos_merge_position(self, df, first_df)
        columns = ["-Position-"]
        dtypes = ["Multiple types"]
        never_add = "Chromosome Start End Strand".split()

    h = [c + "\n(" + str(t) + ")" for c, t in zip(columns, list(dtypes))]

    str_repr = tabulate(
        df.get(columns), headers=h, tablefmt='psql', showindex=False)
    table_width = len(str_repr.split("\n", 1)[0])

    for _, c in enumerate(df.columns):
        if c in never_add:
            continue

        columns.append(c)
        t = all_dtypes[c]
        dtypes.append(t)
        _h = c + "\n(" + str(t) + ")"
        h.append(_h)

        new_build_df = pd.concat([build_df, df[c]], axis=1)

        new_str_repr = tabulate(
            new_build_df, headers=h, tablefmt='psql', showindex=False)

        table_width = len(new_str_repr.split("\n", 1)[0])
        # print("1", table_width, get_terminal_size())
        if table_width >= terminal_width:
            first_hidden = columns[-1]
            columns = columns[:-1]
            dtypes = dtypes[:-1]
            break
        else:
            str_repr = new_str_repr
            build_df = new_build_df
            # print("2", table_width, get_terminal_size())

    hidden_cols = set(df.columns) - (set(columns).union(never_add))
    n_hidden_cols = len(hidden_cols)
    stranded = "Stranded" if self.stranded else "Unstranded"
    str1 = "{} PyRanges object has {:,} rows and {:,} columns from {} chromosomes.".format(
        stranded, n_intervals, len(self.columns), n_chromosomes)

    if n_hidden_cols:

        # try to add ... as last col
        ddd = pd.Series("...", index=build_df.index)
        ddd.name = "...\n..."
        h[-1] = ddd.name

        new_build_df = pd.concat([build_df, ddd], axis=1)
        new_str_repr = tabulate(
            new_build_df, headers=h, tablefmt='psql', showindex=False)
        table_width = len(new_str_repr.split("\n", 1)[0])

        if table_width <= terminal_width:
            str_repr = new_str_repr
        else:
            # need to remove last columns, because adding ... made it too wide
            columns = list(build_df.columns)
            first_hidden = columns[-1]
            build_df = build_df.drop(columns[-1], axis=1)
            ddd = pd.Series("...", index=build_df.index)
            ddd.name = "..."
            h[-2] = ddd.name
            new_build_df = pd.concat([build_df, ddd], axis=1)
            str_repr = tabulate(
                new_build_df, headers=h, tablefmt='psql', showindex=False)

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

    return str_repr


def tostring(self, n=8, merge_position=False):

    entries = n
    half_entries = int(entries / 2)

    # TODO: test this f
    if len(self) == 0:
        return "Empty PyRanges"

    # keys = natsorted(list(self.dfs.keys()))

    if len(self.keys()) == 1:

        first_key = self.keys()[0]
        # last_key = list(self.dfs.keys())[-1]
        first_df = self.dfs[first_key]

        # last_df = ray.get(self.dfs[last_key]).tail(half_entries)
        h = first_df.head(half_entries).astype(object)
        m = first_df.head(1).astype(object)
        t = first_df.tail(half_entries).astype(object)
        m.loc[:, :] = "..."

        if len(self) > entries:
            s = pd.concat([h, m, t])
        elif len(self) == entries:
            s = pd.concat([h, t])
        else:  # < entries
            s = first_df.copy()
    elif len(self) <= entries:
        s = self.df.copy()
        s = s.astype(object)
        first_df = s.copy()
    else:
        keys = self.keys()
        first_key = keys[0]
        last_key = keys[-1]
        # first_key = self.keys[0]
        # last_key = self.keys[-1]
        first_df = self.dfs[first_key].head(half_entries)
        last_df = self.dfs[last_key].tail(half_entries)

        # last_df = self.dfs[list(self.dfs.keys())[-1]].tail(half_entries)

        h = first_df.head(half_entries).astype(object)
        m = first_df.head(1).astype(object)
        t = last_df.head(half_entries).astype(object)
        m.loc[:, :] = "..."
        # m.index = ["..."]
        # print((len(h) + len(t)) < entries, len(self) >= entries)
        if (len(h) + len(t)) < entries:

            keys_covered = set()
            # iterate from front until have three
            heads = []
            hl = 0
            for k in keys:
                keys_covered.add(k)
                h = self.dfs[k].head(half_entries)
                first_df = h
                hl += len(h)
                heads.append(h)
                if hl >= half_entries:
                    break

            tails = []
            tl = 0

            for k in keys[::-1]:
                if k in keys_covered:
                    continue

                t = self.dfs[k].tail(half_entries)
                tl += len(t)
                tails.append(t)
                if tl >= half_entries:
                    break
            # iterate from back until have three

            h = pd.concat(heads).head(half_entries).astype(object)
            if tails:
                t = pd.concat(tails).tail(half_entries).astype(object)
                if len(h) + len(t) >= entries:
                    m = h.head(1).astype(object)
                    m.loc[:, :] = "..."
                    s = pd.concat([h, m, t])
                else:
                    s = pd.concat([h, t])
            else:
                s = h

        elif len(h) + len(t) == entries:
            m.loc[:, :] = "..."
            s = pd.concat([h, m, t])
        else:
            s = pd.concat([h, t])

    str_repr = increase_string_width(self, s, first_df, merge_position)

    return str_repr


def sort_tostring(self, n=30, merge_position=False):

    sort_cols = "Start End".split()
    if self.stranded:
        sort_cols.append("Strand")

    first_df = self.values()[0]

    if len(self) <= n:
        df = self.df.sort_values(sort_cols)
    else:
        dfs = []
        half_n = int(n / 2)
        current_len = 0

        for chromosome in self.chromosomes:

            df = self[chromosome].df.sort_values(sort_cols)
            current_len += len(df)
            dfs.append(df.head(half_n))

            if current_len > half_n:
                break

        head = pd.concat(dfs)

        current_len = 0
        dfs = []
        for chromosome in reversed(self.chromosomes):

            df = self[chromosome].df.sort_values(sort_cols)
            current_len += len(df)
            dfs.append(df.tail(half_n))

            if current_len > half_n:
                break

        tail = pd.concat(dfs)

        middle = tail.head(1).copy().astype(object)
        middle.loc[:, :] = "..."

        df = pd.concat([head, middle, tail])

    str_repr = increase_string_width(self, df, first_df, merge_position)

    return str_repr
