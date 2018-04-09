# pyranges

Not optimized, not ready for use, only the "in self but not in other" operator
(-) implemented.

The input dataframe needs at least the columns Chromome, Start, End. The rest are arbitrary.

Using all operations on the GRanges object return a new GRanges object, not a
view. (Execpt those that just change the underlying dataframe, which is
accessible with gr.df)

Columns in the dataframe given to the constructor are available as an attribute
of the object. I.e. gr.Score returns the Score series in the example below.

You can add columns to the underlying dataframe just fine and change all columns
except Chromosome, Start and End. If you want to change those columns you should
create a new GRanges object with the new values. Note: Changing the Chromosome, Start
and End columns does not change the datastructures in the GRanges object that
are used to perform quick range operations on the GRange and will therefore lead to
erroneus results.

```
In [1]: from pyranges import GRanges

In [2]: import pandas as pd

In [3]: df = pd.read_table("example_data/lamina.bed", sep="\s+", header=None, names="Chromosome Start End Score".split(), skiprows=1)

In [4]: gr = GRanges(df)

In [5]: gr
Out[5]: <pyranges.pyranges.GRanges at 0x112341f60>

In [6]: print(gr)
+------|--------------|----------|----------|--------------------+
|      | Chromosome   | Start    | End      | Score              |
|------|--------------|----------|----------|--------------------|
| 0    | chr1         | 11323785 | 11617177 | 0.86217008797654   |
| 1    | chr1         | 12645605 | 13926923 | 0.9348914858096831 |
| 2    | chr1         | 14750216 | 15119039 | 0.9459459459459459 |
| ...  | ...          | ...      | ...      | ...                |
| 1341 | chrY         | 13556427 | 13843364 | 0.892655367231638  |
| 1342 | chrY         | 14113371 | 15137286 | 0.9364089775561101 |
| 1343 | chrY         | 15475619 | 19472504 | 0.8138424821002389 |
+------|--------------|----------|----------|--------------------+
GRanges object with 1344 sequences from 24 chromosomes.

In [7]: subset_gr = gr[:111617177]

In [8]: print(subset_gr)
+-----|--------------|-----------|-----------|--------------------+
|     | Chromosome   | Start     | End       | Score              |
|-----|--------------|-----------|-----------|--------------------|
| 0   | chr15        | 50892438  | 53249359  | 0.920298879202989  |
| 1   | chr15        | 23238025  | 23460157  | 0.8214285714285708 |
| 2   | chr15        | 20639479  | 22659952  | 0.8454692556634301 |
| ... | ...          | ...       | ...       | ...                |
| 989 | chr5         | 107749740 | 107997497 | 0.9004975124378108 |
| 990 | chr5         | 108850224 | 109049698 | 1.0                |
| 991 | chr5         | 111540207 | 112050044 | 0.9844054580896691 |
+-----|--------------|-----------|-----------|--------------------+
GRanges object with 992 sequences from 24 chromosomes.

In [9]: print(subset_gr["chr1"])
+-----|--------------|-----------|-----------|--------------------+
|     | Chromosome   | Start     | End       | Score              |
|-----|--------------|-----------|-----------|--------------------|
| 0   | chr1         | 95358151  | 96453308  | 0.952742616033755  |
| 1   | chr1         | 38272060  | 39078902  | 0.940932642487047  |
| 2   | chr1         | 12645605  | 13926923  | 0.9348914858096831 |
| ... | ...          | ...       | ...       | ...                |
| 45  | chr1         | 110470808 | 110667025 | 0.7201646090534979 |
| 46  | chr1         | 108566438 | 108844851 | 1.0                |
| 47  | chr1         | 111580508 | 111778144 | 0.9533678756476679 |
+-----|--------------|-----------|-----------|--------------------+
GRanges object with 48 sequences from 1 chromosomes.

In [10]: gr.Score.head()
Out[10]:
0    0.862170
1    0.934891
2    0.945946
3    0.895175
4    0.892526
Name: Score, dtype: float64
```
