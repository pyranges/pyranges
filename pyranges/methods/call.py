import pandas as pd
import pyranges as pr


def _handle_eval_return(self, result, col, as_pyranges, subset):
    """Handle return from eval.

    If col is set, add/update cols. If subset is True, use return series to subset PyRanges.
    Otherwise return PyRanges or dict of data."""

    if as_pyranges:
        if not result:
            return pr.PyRanges()

        first_hit = list(result.values())[0]

        if isinstance(first_hit, pd.Series):
            if first_hit.dtype == bool and subset:
                return self[result]
            elif col:
                self.__setattr__(col, result)
                return self
            else:
                raise Exception(
                    "Cannot return PyRanges when function returns a Series! Use as_pyranges=False."
                )
        return pr.PyRanges(result)
    else:
        return result


def _call(self, f, strand=None, as_pyranges=True, **kwargs):

    if strand is None:
        strand = self.stranded

    if self.stranded and not strand:
        self = self.unstrand()

    result = self.apply(f, strand=strand, as_pyranges=False, **kwargs)

    # result = _handle_eval_return(self, result, col, as_pyranges, subset)

    return result
