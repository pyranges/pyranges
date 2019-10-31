import pandas as pd
import numpy as np

import pyranges as pr
from pyranges.multithreaded import pyrange_apply

from pyranges.methods.statistics import _relative_distance


class StatisticsMethods():

    pr = None

    def __init__(self, pr):

        self.pr = pr

    def jaccard(self, other, **kwargs):

        self = self.pr

        kwargs["sparse"] = {"self": True, "other": True}
        kwargs = pr.pyranges.fill_kwargs(kwargs)
        strand = True if kwargs["strandedness"] else False

        intersection_sum = sum(
            v.sum()
            for v in self.set_intersect(other, **kwargs).lengths(as_dict=True).values())

        union_sum = 0
        for gr in [self, other]:
            union_sum += sum(
                v.sum() for v in gr.merge(strand=strand).lengths(as_dict=True).values())

        denominator = (union_sum - intersection_sum)
        if denominator == 0:
            return 1
        else:
            jc = intersection_sum / denominator

        return jc

    def relative_distance(self, other, **kwargs):

        self = self.pr

        kwargs["sparse"] = {"self": True, "other": True}
        kwargs = pr.pyranges.fill_kwargs(kwargs)

        result = pyrange_apply(_relative_distance, self, other, **kwargs)  # pylint: disable=E1132

        result = pd.Series(np.concatenate(list(result.values())))

        not_nan = ~np.isnan(result)
        result.loc[not_nan] = np.floor(result[not_nan] * 100) / 100
        vc = result.value_counts(dropna=False).to_frame().reset_index()
        vc.columns = "reldist count".split()
        vc.insert(vc.shape[1], "total", len(result))
        vc.insert(vc.shape[1], "fraction", vc["count"] / len(result))
        vc = vc.sort_values("reldist", ascending=True)
        vc = vc.reset_index(drop=True)

        return vc


def mcc(grs, genome, labels=None, strandedness=None):

    from itertools import product

    if labels is None:
        labels = list(range(len(grs)))
        labels = product(labels, labels)
    else:
        assert len(labels) == len(grs)
        labels = product(labels, labels)

    rowdicts = []
    for (lt, lf), (t, f) in zip(labels, product(grs, grs)):
        print(lt, lf)

        if t == f:
            tp = t.length
            fn = 0
            tn = genome.length - tp
            fp = 0
        else:
            f = f.merge()
            c = pr.concat([t, f])
            j = t.join(f, strandedness=strandedness)
            tp_gr = j.new_position("intersection").merge()
            tp = tp_gr.length
            fp = f.subtract(tp_gr, strandedness=strandedness, no_merge=True).merge(strand=False).length
            fn = t.subtract(tp_gr, strandedness=strandedness, no_merge=True).merge(strand=False).length
            tn = genome.subtract(c, strandedness=strandedness, no_merge=True).merge(strand=False).length

        # https://stackoverflow.com/a/56875660/992687
        vals = [tp, fp, fn, tn]
        print(vals)
        print([ type(v) for v in vals ])
        x = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
        print("x", x)
        mcc = ((tp * tn) - (fp * fn)) / np.sqrt(x)

        rowdicts.append({"T": lt, "F": lf, "TP": tp, "FP": fp, "TN": tn, "FN": fn, "MCC": mcc})

    df = pd.DataFrame.from_dict(rowdicts)

    return df


