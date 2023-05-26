from typing import Tuple, Union, List

import pandas as pd


def get_chromosomes_from_dict(dfs) -> List[str]:
    keys = list(dfs.keys())
    if isinstance(keys[0], tuple):
        chromosomes = [k[0] for k in keys]
    else:
        chromosomes = keys

    return chromosomes


def get_strands_from_dict(dfs) -> Union[List[str], List[Tuple[str, str]]]:
    keys = list(dfs.keys())
    if isinstance(keys[0], tuple):
        strands = [k[1] for k in keys]
    else:
        strands = keys

    return strands


def get_key_from_df(df: pd.DataFrame) -> Union[str, Tuple[str, str]]:
    chromosome = df.Chromosome.head(1).iloc[0]
    if "Strand" in df:
        strand = df.Strand.head(1).iloc[0]
        return chromosome, strand

    return chromosome


def single_value_key(df: pd.DataFrame) -> bool:
    if "Strand" in df:
        return len(df[["Chromosome", "Strand"]].drop_duplicates(["Chromosome", "Strand"])) == 1
    else:
        return len(df.Chromosome.drop_duplicates()) == 1
