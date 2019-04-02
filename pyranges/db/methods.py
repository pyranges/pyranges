
import os
import pyranges as pr
from ftplib import FTP


def get_ftp_file(ftp, binary_loc, path=None):

    if path:
        if os.path.dirname(path):
            os.makedirs(os.path.dirname(path), exist_ok=True)
        fh = open(path, "wb")
        ftp.retrbinary("RETR {}".format(binary_loc), fh.write) 
        fh.close()
        gr = pr.read_gtf(path)

    else:
        with tempfile.NamedTemporaryFile(suffix=".gz", delete=False) as t:

            ftp.retrbinary("RETR {}".format(binary_loc), t.write) 
            t.close()

            gr = pr.read_gtf(t.name)

            os.remove(t.name)

    return gr
