import re
from collections.abc import Iterable


def _drop(self, drop=None, drop_strand=False):
    columns = self.columns
    if drop_strand:
        always_keep = "Chromosome Start End".split()
    else:
        always_keep = "Chromosome Start End Strand".split()

    _to_drop = []

    if not drop:
        _to_drop = set(columns) - set(always_keep)
    elif isinstance(drop, str):
        r = re.compile(drop)
        _to_drop = [c for c in columns if not r.search(c) == None]
    elif isinstance(drop, Iterable) or isinstance(drop, list):
        _to_drop = drop
    else:
        raise Exception("Not valid subsetters!")

    if not drop_strand:
        assert not set(always_keep).intersection(_to_drop), \
            "Can never drop Chromosome, Start or End. Cannot drop Strand unless drop_strand is True"
    _to_drop = set(_to_drop) - set(always_keep)

    return self.apply(lambda df: df.drop(_to_drop, axis=1))


def _keep(self, keep, drop_strand=False):

    columns = self.columns
    if drop_strand:
        always_keep = "Chromosome Start End".split()
    else:
        always_keep = "Chromosome Start End Strand".split()

    _to_drop = []
    if isinstance(keep, str):
        r = re.compile(keep)
        _to_drop = [c for c in columns if r.search(c) is None]
    elif isinstance(keep, Iterable) or isinstance(keep, list):
        _to_drop = set(columns) - set(keep) - set(always_keep)

    return self.apply(lambda df: df.drop(_to_drop, axis=1))
