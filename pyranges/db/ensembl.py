import os
import pandas as pd

import MySQLdb

from pyranges.db.methods import get_ftp_file


def _find_correct_db_version(conn, genome, version):
    """ensembl dbs have names like danio_rerio_75_11 where the user has no way of knowing the last
magica number. This function connects to the correct database and version number (danio_rerio_75)."""

    query = 'show databases like "{}_core_%"'.format(genome)

    df = pd.read_sql(query, conn)

    if df.empty:
        raise Exception(
            "Genome {} not found. List all possible genomes and versions with pyranges.db.ensembl.genomes()"
            .format(genome))

    df.columns = ["db"]

    df = df[(df.db.str.contains(genome)) & (df.db.str.contains("_core_"))]

    versions = df.db.str.extract(genome +
                                 "_core_([0-9]+)_*").astype(int).squeeze()

    if version == "latest":
        idx = versions.idxmax()
        db = df.loc[idx].iloc[0]
    else:
        idx = versions[version == versions].index

        assert len(
            idx), "No version {} in list of ensembl versions: {}".format(
                version, versions)
        assert len(
            idx) == 1, "Multiple rows matching version number {} in {}".format(
                version, df)

        db = df.loc[idx].iloc[0, 0]

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


def chromosome_sizes(genome, version="latest"):

    query = 'select name,length from seq_region'

    df = ensembl(genome, query, version)

    return pd.Series(data=df.length.values, index=df.name.values)


def genomes():

    host = "ensembldb.ensembl.org"
    user = "anonymous"
    conn = MySQLdb.Connection(host=host, user=user)

    query = "show databases like '%_core_%'"
    df = pd.read_sql(query, conn)
    df.columns = ["db"]

    df = df.db.str.extract("(\w+)_core_(\d+)_*.")
    df.columns = ["Genome", "Version"]

    return df


# ftp://ftp.ensembl.org/pub/release-95/gtf/gopherus_agassizii

from ftplib import FTP


def genes(genome, release="latest", path=None, head=False):

    ftp = FTP("ftp.ensembl.org")
    ftp.login()

    if release == "latest":
        release = max(
            [r.split("-")[1] for r in ftp.nlst("pub/") if "release-" in r],
            key=int)

    _dir = "pub/release-{}/gtf/{}/".format(release, genome)
    dir_listing = pd.Series(ftp.nlst(_dir))

    matches = (dir_listing.str.lower().apply(os.path.basename).str.match(
        "{}.*{}.gtf.gz".format(genome, release)))

    assert matches.sum() == 1, "More than one file matching: {}".format(
        dir_listing[matches])

    if not len(dir_listing):
        raise Exception(
            "No files found for genome {}. Use pyranges.db.ensembl.genomes() to list all available genomes."
            .format(genome))

    binary_loc = dir_listing[matches].iloc[0]

    host = "ftp.ensembl.org"
    return get_ftp_file(host, binary_loc, path, test=head)
