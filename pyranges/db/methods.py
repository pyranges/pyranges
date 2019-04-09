import os
import pyranges as pr

import tempfile

import requests


def get_ftp_file(ftp, binary_loc, path=None, test=False):

    url = "http://" + ftp + "/" + binary_loc
    if test:
        # if test only download small part
        from io import BytesIO
        import gzip
        headers = {"Range": "bytes=0-1000000"}
        c = requests.get(url, headers=headers).content
        s = gzip.GzipFile(fileobj=BytesIO(c)).read(1000000)
        f = s.decode().rsplit("\n", 1)[0]
        with tempfile.NamedTemporaryFile(suffix=".gtf", mode="w+") as t:
            t.write(f)
            gr = pr.read_gtf(t.name)
    else:
        if path and os.path.dirname(path):
            os.makedirs(os.path.dirname(path), exist_ok=True)

        c = requests.get(url).content

        if not path:
            with tempfile.NamedTemporaryFile(
                    suffix=".gtf.gz", mode="wb+") as t:
                t.write(c)
                gr = pr.read_gtf(t.name)
        else:
            if os.path.dirname(path):
                os.makedirs(os.path.dirname(path), exist_ok=True)
            fh = open(path, "wb+")
            fh.write(c)
            fh.close()

            gr = pr.read_gtf(path)

    return gr
