import pandas as pd
import pyranges as pr
import ray

def assert_dfs_equal(gr1, gr2):

    dfs1 = gr1.dfs
    dfs2 = gr2.dfs

    print(list(dfs1.keys()), list(dfs2.keys()))
    assert dfs1.keys() == dfs2.keys()

    for k, v in dfs1.items():

        v2 = dfs2[k]
        try:
            v = ray.get(v)
            v2 = ray.get(v2)
        except:
            pass

        assert_df_equal(v, v2)

def assert_df_equal(df1, df2):

    # if not isinstance(df1, pd.DataFrame):
    #     df1 = ray.get(df1)
    #     df2 = ray.get(df2)

    if "Strand" in df1 and "Strand" in df2:
        sort_on = "Chromosome Start End Strand".split()
        df1.Strand = df1.Strand.astype("object")
        df2.Strand = df2.Strand.astype("object")
    else:
        sort_on = "Chromosome Start End".split()

    if "Strand_b" in df1:
        sort_on += "Start_b End_b Strand_b".split()
        df1.Strand_b = df1.Strand_b.astype("object")
        df2.Strand_b = df2.Strand_b.astype("object")
    elif "Start_b" in df2:
        sort_on += "Start_b End_b".split()

    df1.Chromosome = df1.Chromosome.astype("object")
    df2.Chromosome = df2.Chromosome.astype("object")

    df1 = df1.reset_index(drop=True)
    df2 = df2.reset_index(drop=True)

    df1 = df1.sort_values(sort_on)
    df2 = df2.sort_values(sort_on)

    print("Actual")
    print(df1.to_csv(sep=" "))
    print("Expected")
    print(df2.to_csv(sep=" "))

    print("Actual dtypes")
    print(df1.dtypes)
    print("Expected dtypes")
    print(df2.dtypes)
    # print("dtypes Strand\n", "1",  df1.Strand.dtype, "2", df2.Strand.dtype)
    # print("dtypes Strand\n", df1.Strand.dtype == df2.Strand.dtype)
    # print("dtypes equal\n", df1.dtypes == df2.dtypes)

    print("Actual index")
    print(df1.index)
    print("Expected index")
    print(df2.index)
    print("index equal", df1.index == df2.index)

    pd.testing.assert_frame_equal(df1, df2)

from io import StringIO

def string_to_pyrange(s):

    df = pd.read_table(StringIO(s), sep="\s+", header=0)
    return pr.PyRanges(df)
