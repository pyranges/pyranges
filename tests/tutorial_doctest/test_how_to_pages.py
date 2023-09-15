import doctest
import os


def test_how_to_pages():
    failure_count2, test_count2 = doctest.testfile('../../docs/how_to_pages.rst')
    if not failure_count2:
        print('All tests in the how_to_pages were successful!')
    print()
    os.remove("minigenome.fa")
    os.remove("minigenome.fa.fai")
    os.remove("chipseq.gtf")
