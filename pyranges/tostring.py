import shutil

import pyranges as pr
import pandas as pd

from tabulate import tabulate


def get_terminal_size():

    return shutil.get_terminal_size().columns

def show_pos_merge_position(self, df, first_df):

    if self.stranded:
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
        if stranded:
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

    str_repr = tabulate(df.get(columns), headers=h, tablefmt='psql', showindex=False)
    table_width = len(str_repr.split("\n", 1)[0])

    for i, c in enumerate(df.columns):
        if c in never_add:
            continue

        columns.append(c)
        t = all_dtypes[c]
        dtypes.append(t)
        _h = c + "\n(" + str(t) + ")"
        h.append(_h)

        new_build_df =  pd.concat([build_df, df[c]], axis=1)

        new_str_repr = tabulate(new_build_df, headers=h, tablefmt='psql', showindex=False)

        table_width = len(str_repr.split("\n", 1)[0])
        if table_width > terminal_width:
            break

        str_repr = new_str_repr
        build_df = new_build_df

    hidden_cols = set(df.columns) - (set(columns).union(never_add))
    n_hidden_cols = len(hidden_cols)
    str1 = "PyRanges object has {} sequences from {} chromosomes.".format(
        n_intervals, n_chromosomes)

    if n_hidden_cols:
        str2 = "Hidden columns: {}".format(", ".join(hidden_cols))
        if (n_hidden_cols - 10) > 0:
            str3 = "(+ {} more.)".format(n_hidden_cols - 10)
            str_repr = "\n".join([str_repr, str1, str2, str3])
        else:
            str_repr = "\n".join([str_repr, str1, str2])
    else:
        str_repr = "\n".join([str_repr, str1])

    return str_repr


def tostring(self, n=8, merge_position=False):

    entries = n
    half_entries = int(entries / 2)

    # TODO: test this f
    # print("in str")
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
        else:
            s = h
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
        df = self.df.sort(sort_cols)
    else:
        dfs = []
        half_n = int(n/2)
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
