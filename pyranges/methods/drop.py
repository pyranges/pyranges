import re
from collections.abc import Iterable


def _drop(self, drop=None):
    columns = self.columns
    if "Strand" in columns:
        self = self.unstrand()
        columns = [c for c in columns if c != "Strand"]

    if self.stranded:
        always_keep = "Chromosome Start End Strand".split()
    else:
        always_keep = "Chromosome Start End".split()

    _to_drop = []

    if not drop:
        _to_drop = set(columns) - set(always_keep)
    elif isinstance(drop, str):
        _to_drop = [drop]
    elif isinstance(drop, Iterable) or isinstance(drop, list):
        _to_drop = drop
    else:
        raise Exception("Not valid subsetters!")

    _to_drop = set(_to_drop) - set(always_keep)

    return self.apply(lambda df: df.drop(_to_drop, axis=1))


def _keep(self, keep):

    columns = self.columns
    if not self.stranded:
        always_keep = "Chromosome Start End".split()
    else:
        always_keep = "Chromosome Start End Strand".split()

    _to_drop = []
    if isinstance(keep, Iterable) or isinstance(keep, list):
        _to_drop = set(columns) - set(keep) - set(always_keep)
    else:
        raise Exception("Column(s) to drop must be in list.")

    self = self.apply(lambda df: df.drop(_to_drop, axis=1))

    columns = self.columns

    keep = [c for c in keep if not c in always_keep]
    new_order = []
    i = 0
    for c in columns:
        if c in always_keep:
            new_order.append(c)
        else:
            new_order.append(keep[i])
            i += 1

    self = self.apply(lambda df: df[new_order])

    return self
