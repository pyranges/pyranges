"""Tests that all example data loads properly."""

from pyranges.data import (f1, f2, chipseq, chipseq_background, aorta, aorta2,
                           ensembl_gtf, cpg, exons, control_bam)


def test_all_data():

    for f in [
            f1,
            f2,
            chipseq,
            chipseq_background,
            cpg,
            exons,
            aorta,
            aorta2,
            ensembl_gtf  #, control_bam
    ]:
        gr = f()
        print(gr)  # to avoid unused variable warning
