import pyranges as pr
import numpy as np
import pandas as pd


def _setattr(self, column_name, column, pos=False):

    if not len(self):
        return

    isiterable = isinstance(column, (list, pd.Series, np.ndarray))
    isdict = isinstance(column, dict)

    if isiterable:
        if not len(self) == len(column):
            raise Exception("DataFrame and column must be same length.")

    already_exists = column_name in self.columns

    if isinstance(column, pd.Series):
        column = column.values

    if already_exists:
        pos = list(self.values()[0].columns).index(column_name)
    elif not pos:
        pos = self.values()[0].shape[1]

    start_length, end_length = 0, 0

    dfs = {}
    for k, df in self.items():

        end_length += len(df)

        if already_exists:
            df = df.drop(column_name, axis=1)

        if isiterable:
            df.insert(pos, column_name, column[start_length:end_length])
        elif isdict:
            if isinstance(column[k], pd.Series):
                _column = column[k].values
            else:
                _column = column[k]

            df.insert(pos, column_name, _column)
        else:
            df.insert(pos, column_name, column)

        start_length = end_length

        dfs[k] = df

    if column_name not in ["Chromosome", "Strand"]:
        self.__dict__["dfs"] = dfs
    else:
        int64 = True if self.dtypes["Start"] == np.int64 else False
        self.__dict__["dfs"] = pr.PyRanges(pr.PyRanges(dfs).df, int64=int64).dfs # will merge the dfs, then split on keys again to ensure they are correct


def _getattr(self, name):

    if name in self.columns:
        return pd.concat([df[name] for df in self.values()])
    else:
        raise AttributeError("PyRanges object has no attribute", name)
