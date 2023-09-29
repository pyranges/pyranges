import doctest
import os
from pathlib import Path


def test_how_to_pages():
    failure_count2, test_count2 = doctest.testfile("../../docs/how_to_pages.rst")
    if not failure_count2:
        print("All tests in the how_to_pages were successful!")
    FILES = ["minigenome.fa", "minigenome.fa.fai", "chipseq.gtf"]
    [Path(f).unlink() for f in FILES]
