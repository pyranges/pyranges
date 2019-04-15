import numpy as np
import pandas as pd


def _setattr(self, column_name, column):

    if column_name in "Chromosome Strand".split():
        raise Exception("The columns Chromosome and Strand can not be reset.")

    isiterable = isinstance(column, list) or isinstance(
        column, pd.Series) or isinstance(column, np.ndarray)
    isdict = isinstance(column, dict)

    if isiterable:
        if not len(self) == len(column):
            raise Exception("DataFrame and column must be same length.")

    already_exists = column_name in self.values()[0]

    if already_exists:
        pos = list(self.values()[0].columns).index(column_name)
    else:
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
            df.insert(pos, column_name, column[k])
        else:
            df.insert(pos, column_name, column)

        start_length = end_length

        dfs[k] = df

    self.__dict__["dfs"] = dfs


def _getattr(self, name):
    if name in self.values()[0]:
        return pd.concat([df[name] for df in self.values()])
    else:
        raise Exception("PyRanges object has no attribute", name)
