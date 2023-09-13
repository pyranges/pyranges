"""Tests that all example data loads properly."""

from pyranges.data import aorta, aorta2, chipseq, chipseq_background, cpg, ensembl_gtf, exons, f1, f2  # control_bam)


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
        ensembl_gtf,  # , control_bam
    ]:
        f()
