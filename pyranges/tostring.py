import shutil

import pandas as pd

from tabulate import tabulate


def get_terminal_size():

    return shutil.get_terminal_size().columns


def reduce_string_width(str_repr, s, h, n_intervals, n_chromosomes):

    terminal_width = get_terminal_size()

    hidden_cols = []
    header = str_repr.split("\n", 2)[1]

    while len(header) > terminal_width:
        header, hidden = header.rsplit("|", 1)
        hidden = hidden.strip()
        if hidden:
            hidden_cols.append(hidden.strip())

    columns = header.replace("|", "").split()

    str_repr = tabulate(
        s.get(columns), headers=h, tablefmt='psql', showindex=False)

    str_repr += "\nPyRanges object has {} sequences from {} chromosomes.".format(
        n_intervals, n_chromosomes)
    if hidden_cols:
        str_repr += "\nHidden columns: {}".format(", ".join(hidden_cols))

    return str_repr


def tostring(self):
    # TODO: test this f
    # print("in str")
    if len(self) == 0:
        return "Empty PyRanges"

    # keys = natsorted(list(self.dfs.keys()))
    if len(self.keys()) == 1:

        first_key = self.keys()[0]
        # last_key = list(self.dfs.keys())[-1]
        first_df = self.dfs[first_key]

        # last_df = ray.get(self.dfs[last_key]).tail(3)
        h = first_df.head(3).astype(object)
        m = first_df.head(1).astype(object)
        t = first_df.tail(3).astype(object)
        m.loc[:, :] = "..."

        if len(self) > 6:
            s = pd.concat([h, m, t])
        elif len(self) == 6:
            s = pd.concat([h, t])
        else:
            s = h
    else:
        keys = self.keys()
        first_key = keys[0]
        last_key = keys[-1]
        # first_key = self.keys[0]
        # last_key = self.keys[-1]
        first_df = self.dfs[first_key].head(3)
        last_df = self.dfs[last_key].tail(3)
        # last_df = self.dfs[list(self.dfs.keys())[-1]].tail(3)

        h = first_df.head(3).astype(object)
        m = first_df.head(1).astype(object)
        t = last_df.head(3).astype(object)
        m.loc[:, :] = "..."
        # m.index = ["..."]
        # print((len(h) + len(t)) < 6, len(self) >= 6)
        if (len(h) + len(t)) < 6:

            keys_covered = set()
            # iterate from front until have three
            heads = []
            hl = 0
            for k in keys:
                keys_covered.add(k)
                h = self.dfs[k].head(3)
                first_df = h
                hl += len(h)
                heads.append(h)
                if hl >= 3:
                    break

            tails = []
            tl = 0

            for k in keys[::-1]:
                if k in keys_covered:
                    continue

                t = self.dfs[k].tail(3)
                tl += len(t)
                tails.append(t)
                if tl >= 3:
                    break
            # iterate from back until have three

            h = pd.concat(heads).head(3).astype(object)
            if tails:
                t = pd.concat(tails).tail(3).astype(object)
                if len(h) + len(t) > 6:
                    m = h.head(1).astype(object)
                    m.loc[:, :] = "..."
                    s = pd.concat([h, m, t])
                else:
                    s = pd.concat([h, t])
            else:
                s = h

        elif len(h) + len(t) == 6:
            m.loc[:, :] = "..."
            s = pd.concat([h, m, t])
        else:
            s = pd.concat([h, t])

    if False:  # make setting
        if self.stranded:
            pos = s.Chromosome.astype(str) + " " + s.Start.astype(
                str) + "-" + s.End.astype(str) + " " + s.Strand.astype(str)
            s = s.drop("Chromosome Start End Strand".split(), axis=1)
            first_df = first_df.drop(
                "Chromosome Start End Strand".split(), axis=1)
        else:
            pos = s.Chromosome.astype(str) + " " + s.Start.astype(
                str) + "-" + s.End.astype(str)
            s = s.drop("Chromosome Start End".split(), axis=1)
            first_df = first_df.drop("Chromosome Start End".split(), axis=1)

        s.insert(0, "Position", pos)
        h = [
            c + "\n(" + str(t) + ")"
            for c, t in zip(s.columns, ["multiple types"] +
                            list(first_df.dtypes))
        ]
    else:
        dtypes = []
        for col, dtype in zip(s.columns, first_df.dtypes):
            # if str(dtype) == "category":

            #     dtype = first_df[col].cat.codes.dtype

            # dtype = str(dtype).replace("float", "f_").replace("int", "i_")
            dtypes.append(dtype)

        h = [c + "\n(" + str(t) + ")" for c, t in zip(s.columns, list(dtypes))]

    str_repr = tabulate(s, headers=h, tablefmt='psql', showindex=False) + \
                                    "\nPyRanges object has {} sequences from {} chromosomes.".format(len(self), len(self.chromosomes))

    str_repr = reduce_string_width(str_repr, s, h, len(self),
                                   len(self.chromosomes))

    return str_repr
