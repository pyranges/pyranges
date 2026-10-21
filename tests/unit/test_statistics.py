import numpy as np
import pandas as pd
import pytest

import pyranges as pr


def test_fisher_exact_matches_known_contingency_tables():
    pytest.importorskip("fisher")
    tables = pd.DataFrame({"TP": [12, 0], "FP": [5, 12], "TN": [29, 10], "FN": [2, 2]})

    result = pr.stats.fisher_exact(tables.TP, tables.FP, tables.TN, tables.FN)

    np.testing.assert_allclose(result.OR, [12 * 2 / (5 * 29), 0])
    np.testing.assert_allclose(result.P, [0.080269, 0.000067], atol=5e-7)
    np.testing.assert_allclose(result.PLeft, [0.044555, 0.000034], atol=5e-7)
    np.testing.assert_allclose(result.PRight, [0.994525, 1.0], atol=5e-7)


@pytest.mark.parametrize("count", [-1, 2**32])
def test_fisher_exact_rejects_counts_outside_uint32_range(count):
    pytest.importorskip("fisher", minversion="0.2.0")
    with pytest.raises(ValueError, match="requires counts between"):
        pr.stats.fisher_exact(pd.Series([count]), [5], [29], [2])
