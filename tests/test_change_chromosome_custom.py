import pandas as pd

import pyranges as pr


def test_change_chromosomes():
    df1 = pd.DataFrame({"Chromosome": ["chr1", "chr2"], "Start": [100, 200], "End": [150, 201]})
    py1 = pr.PyRanges(df1)

    def modify_chrom_series(df):
        df.Chromosome = df.Chromosome.apply(lambda val: val.replace("chr", ""))
        return df

    def fix_chrom(regions):
        return regions.apply(modify_chrom_series)

    print(py1)

    py1 = fix_chrom(py1)

    assert py1.chromosomes == ["1", "2"]
