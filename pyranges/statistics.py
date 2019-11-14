import pandas as pd
import numpy as np

import pyranges as pr
from pyranges.multithreaded import pyrange_apply

from pyranges.methods.statistics import _relative_distance

import numpy as np

def fdr(p_vals):

    from scipy.stats import rankdata
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1

    return fdr


def fisher_exact(n1, d1, n2, d2, **kwargs):
    try:
        from fisher import pvalue_npy
    except:
        import sys
        print("fisher needs to be installed to use fisher exact. pip install fisher or conda install -c bioconda fisher.")
        sys.exit(-1)

    pseudocount = kwargs.get("pseudocount", 0)
    fe_type = kwargs.get("alternative", "twosided")

    n1 = np.array(n1, dtype=np.uint)
    n2 = np.array(n2, dtype=np.uint)
    d1 = np.array(d1, dtype=np.uint)
    d2 = np.array(d2, dtype=np.uint)

    left, right, twosided = pvalue_npy(n1, d1, n2, d2)

    if fe_type == "twosided":
        p_vals = twosided
    elif fe_type == "left":
        p_vals = left
    elif fe_type == "right":
        p_vals = right
    else:
        raise Exception("alternative must be twosided, left or right")


    OR = ((n1 + pseudocount) / (d2 + pseudocount)) / ((n2 + pseudocount) / (d1 + pseudocount))

    df = pd.DataFrame({"OR": OR, "P": p_vals})

    return df


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

from math import sqrt
def _mcc(tp, fp, tn, fn):

    # https://stackoverflow.com/a/56875660/992687
    x = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    return ((tp * tn) - (fp * fn)) / sqrt(x)


def mcc(grs, genome=None, labels=None, strand=False, verbose=False):
    import sys
    from itertools import combinations_with_replacement, chain

    if labels is None:
        _labels = list(range(len(grs)))
        _labels = combinations_with_replacement(_labels, r=2)
    else:
        assert len(labels) == len(grs)
        _labels = combinations_with_replacement(labels, r=2)

    # remove all non-loc columns before computation
    grs = [gr.merge(strand=strand) for gr in grs]

    if genome is not None:
        try:
            genome_length = int(genome)
        except (TypeError, ValueError):
            genome_length = int(genome.End.sum())

        if verbose:
            # check that genome definition does not have many more
            # chromosomes than datafiles
            gr_cs = set(chain(*[gr.chromosomes for gr in grs]))

            g_cs = set(genome.chromosomes)
            surplus = g_cs - gr_cs
            if len(surplus):
                print("The following chromosomes are in the genome, but not the PyRanges:", ", ".join(surplus), file=sys.stderr)


        if strand:
            def make_stranded(df):
                df = df.copy()
                df2 = df.copy()
                df.insert(df.shape[1], "Strand", "+")
                df2.insert(df2.shape[1], "Strand", "-")
                return pd.concat([df, df2])

        genome = genome.apply(make_stranded)
        genome_length = genome.length
    else:
        if len(grs) == 2:
            print("If you do not have a genome and the number of compared pyranges is two, mcc will not give any true negatives.", file=sys.stderr)
        genome_length = pr.concat(grs).merge().length


    strandedness = "same" if strand else None

    rowdicts = []
    for (lt, lf), (t, f) in zip(_labels, combinations_with_replacement(grs, r=2)):
        if verbose:
            print(lt, lf, file=sys.stderr)

        if lt == lf:

            if not strand:
                tp = t.length
                fn = 0
                tn = genome_length - tp
                fp = 0
                rowdicts.append({"T": lt, "F": lf, "TP": tp, "FP": fp, "TN": tn, "FN": fn, "MCC": 1})
            else:
                for strand in "+ -".split():
                    tp = t[strand].length
                    fn = 0
                    tn = genome_length - tp
                    fp = 0
                    rowdicts.append({"T": lt, "F": lf, "Strand": strand, "TP": tp, "FP": fp, "TN": tn, "FN": fn, "MCC": 1})
            continue

        else:
            j = t.join(f, strandedness=strandedness)
            tp_gr = j.new_position("intersection").merge(strand=strand)
            if strand:
                for strand in "+ -".split():
                    tp = tp_gr[strand].length
                    fp = f[strand].length - tp
                    fn = t[strand].length - tp
                    tn = genome_length - (tp + fp + fn)
                    mcc = _mcc(tp, fp, tn, fn)
                    rowdicts.append({"T": lt, "F": lf, "Strand": strand, "TP": tp, "FP": fp, "TN": tn, "FN": fn, "MCC": mcc})
                    rowdicts.append({"T": lf, "F": lt, "Strand": strand, "TP": tp, "FP": fn, "TN": tn, "FN": fp, "MCC": mcc})
            else:
                tp = tp_gr.length
                fp = f.length - tp
                fn = t.length - tp
                tn = genome_length - (tp + fp + fn)
                mcc = _mcc(tp, fp, tn, fn)

                rowdicts.append({"T": lt, "F": lf, "TP": tp, "FP": fp, "TN": tn, "FN": fn, "MCC": mcc})
                rowdicts.append({"T": lf, "F": lt, "TP": tp, "FP": fn, "TN": tn, "FN": fp, "MCC": mcc})

    df = pd.DataFrame.from_dict(rowdicts).sort_values(["T", "F"])

    return df
