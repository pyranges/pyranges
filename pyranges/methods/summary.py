from collections import OrderedDict

import pandas as pd

from tabulate import tabulate


def _summary(self, to_stdout=True, return_df=False):

    lengths = {}
    total_lengths = {}
    lengths["pyrange"] = self.lengths(as_dict=True)
    total_lengths["pyrange"] = [self.length]

    if self.stranded:
        c = self.merge(strand=True)
        lengths["coverage_forward"] = c["+"].lengths(as_dict=True)
        lengths["coverage_reverse"] = c["-"].lengths(as_dict=True)
        total_lengths["coverage_forward"] = [c["+"].length]
        total_lengths["coverage_reverse"] = [c["-"].length]
    else:
        c = self

    c = c.merge(strand=False)
    lengths["coverage_unstranded"] = c.lengths(as_dict=True)
    total_lengths["coverage_unstranded"] = [c.length]

    summaries = OrderedDict()

    # statistics for lengths
    for summary, d in lengths.items():
        if d:
            summaries[summary] = pd.concat(d.values()).describe()

    summary = pd.concat(summaries.values(), axis=1)
    summary.columns = list(summaries)

    df = pd.DataFrame.from_dict(total_lengths)
    df.index = ["sum"]
    summary = pd.concat([summary, df])

    if to_stdout:
        str_repr = tabulate(summary, headers=summary.columns, tablefmt='psql')
        print(str_repr)

    if return_df:
        return summary
