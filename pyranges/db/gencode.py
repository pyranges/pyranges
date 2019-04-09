import pandas as pd

from ftplib import FTP

from natsort import natsorted

from pyranges.db.methods import get_ftp_file


def genomes():

    hr = natsorted(_releases("human"))
    mr = natsorted(_releases("mouse"))

    hg = ["Human"] * len(hr)
    mg = ["Mouse"] * len(mr)

    v = pd.Series(hr + mr).str.extract("release_(\w+)")
    g = pd.Series(hg + mg)

    df = pd.concat([g, v], axis=1, sort=False)
    df.columns = "Genome Version".split()

    return df


def _releases(species="human"):

    assert species in "human mouse".split()

    ftp = FTP("ftp.ebi.ac.uk")
    ftp.login()
    ftp.cwd('pub/databases/gencode/Gencode_' + species)
    dir_listing = ftp.nlst()

    return [l for l in dir_listing if "release" in l]


def genes(species, release="latest", path=None, head=False):

    assert species in "human mouse".split()

    releases = []
    for r in _releases(species):
        r = r.split("_")[1]
        if r.replace("M", "").isdecimal():
            releases.append(r)

    if not release == "latest":
        assert str(
            release) in releases, str(release) + " not in " + str(releases)
    else:
        release = max(releases, key=lambda n: int(n.replace("M", "")))

    # ftp = FTP("ftp.ebi.ac.uk")
    # ftp.login()
    ftp_host = "ftp.ebi.ac.uk"

    binary_loc = 'pub/databases/gencode/Gencode_' + species + "/release_" + release + "/gencode.v{}.annotation.gtf.gz".format(
        release)

    return get_ftp_file(ftp_host, binary_loc, path, test=head)
