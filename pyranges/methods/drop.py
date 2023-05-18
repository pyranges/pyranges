from collections.abc import Iterable


def _drop(self, drop=None, like=None):
    columns = self.columns

    if isinstance(drop, str):
        drop = [drop]

    if drop is not None:
        for i in "Chromosome Start End".split():
            assert i not in drop, "Cannot drop {}".format(i)

    want_to_drop_strand = isinstance(drop, str) and drop == "Strand" or (isinstance(drop, list) and "Strand" in drop)
    if not self.stranded or want_to_drop_strand:
        always_keep = "Chromosome Start End".split()
    else:
        always_keep = "Chromosome Start End Strand".split()

    _to_drop = []

    if like:
        import re

        r = re.compile(like)
        _to_drop = [c for c in self.columns if r.search(c) is not None]
    elif not drop:
        _to_drop = set(columns) - set(always_keep)
    # elif isinstance(drop, str):
    #     _to_drop = [drop]
    elif isinstance(drop, Iterable) or isinstance(drop, list):
        _to_drop = drop
    else:
        raise Exception("Not valid subsetters!")

    _to_drop = set(_to_drop) - set(always_keep)

    # need to use unstrand to remove "+"/"-" from dict
    if "Strand" in _to_drop and self.stranded:
        self = self.unstrand()
        _to_drop = [d for d in _to_drop if not d == "Strand"]

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
        raise Exception("Column(s) to keep must be in list.")

    self = self.apply(lambda df: df.drop(_to_drop, axis=1))

    columns = self.columns

    keep = [c for c in keep if c not in always_keep]
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
