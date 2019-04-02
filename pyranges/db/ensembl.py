import pandas as pd
import numpy as np

import pyranges as pr

import MySQLdb


def _find_correct_db_version(conn, genome, version):

    """ensembl dbs have names like danio_rerio_75_11 where the user has no way of knowing the last
magica number. This function connects to the correct database and version number (danio_rerio_75)."""

    query = 'show databases like "{}_core_%"'.format(genome)

    df = pd.read_sql(query, conn)

    if df.empty:
        # query = 'show databases like %_core_%"'
        # df = pd.read_sql(query, conn)
        raise Exception("Genome {} not found.".format(genome))

    df.columns = ["db"]

    df = df[(df.db.str.contains(genome)) & (df.db.str.contains("_core_"))]

    versions = df.db.str.extract(genome + "_core_([0-9]+)_*").astype(int).squeeze()

    if version == "latest":
        idx = versions.idxmax()
        db = df.loc[idx].iloc[0]
    else:
        idx = versions[version == versions].index

        assert len(idx), "No version {} in list of ensembl versions: {}".format(version, versions)
        assert len(idx) == 1, "Multiple rows matching version number {} in {}".format(version, df)

        db = df.loc[idx].iloc[0,0]

    cursor = conn.cursor()
    cursor.execute("use {}".format(db))

    return conn


def ensembl(genome, query, version="latest"):

    host = "ensembldb.ensembl.org"
    user = "anonymous"
    conn = MySQLdb.Connection(host=host, user=user)

    conn = _find_correct_db_version(conn, genome, version)

    df = pd.read_sql(query, conn)

    conn.close()

    return df


def chromosome_lengths(genome, version="latest"):

    query = 'select name,length from seq_region'

    df = ensembl(genome, query, version)

    return pd.Series(data=df.length.values, index=df.name.values)


