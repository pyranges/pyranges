
import os
import pandas as pd

from ftplib import FTP

import tempfile

import pyranges as pr

def releases(species="human"):

    assert species in "human mouse".split()

    ftp = FTP("ftp.ebi.ac.uk")
    ftp.login()
    ftp.cwd('pub/databases/gencode/Gencode_' + species)
    dir_listing = ftp.nlst()

    return [l for l in dir_listing if "release" in l]


def genes(species, release="latest"):

    assert species in "human mouse".split()

    _releases = []
    for r in releases(species):
        r = r.split("_")[1]
        if r.isdecimal():
            _releases.append(r)

    if not release == "latest":
        assert str(release) in _releases, str(release) + " not in " + str(_releases)
    else:
        release = max(_releases, key=int)

    ftp = FTP("ftp.ebi.ac.uk")
    ftp.login()

    binary_loc = 'pub/databases/gencode/Gencode_' + species + "/release_" + release + "/gencode.v{}.annotation.gtf.gz".format(release)

    with tempfile.NamedTemporaryFile(suffix=".gz", delete=False) as t:

        ftp.retrbinary("RETR {}".format(binary_loc), t.write) #tp.retrbinary("RETR {}".format(binary_loc), t.write)
        t.close()

        gr = pr.read_gtf(t.name)

        os.remove(t.name)

    return gr
