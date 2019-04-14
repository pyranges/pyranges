import pkg_resources

import pyranges as pr
import pandas as pd


def get_example_path(basename):

    full_path = pkg_resources.resource_filename(
        "pyranges", "example_data/{}".format(basename))

    if full_path.endswith(".bam"):
        _hack_to_load_idx = pkg_resources.resource_filename(
            "pyranges", "example_data/{}.bai".format(basename))

    return full_path


def f1():

    full_path = get_example_path("f1.bed")

    return pr.read_bed(full_path)


def f2():

    full_path = get_example_path("f2.bed")

    return pr.read_bed(full_path)


def chipseq():

    full_path = get_example_path("chipseq.bed")

    return pr.read_bed(full_path)


def chipseq_background():

    full_path = get_example_path("chipseq_background.bed")

    return pr.read_bed(full_path)


def aorta():

    full_path = get_example_path("aorta.bed")

    return pr.read_bed(full_path)


def aorta2():

    full_path = get_example_path("aorta2.bed")

    return pr.read_bed(full_path)


def ensembl_gtf():

    full_path = get_example_path("ensembl.gtf")

    return pr.read_gtf(full_path)


def control_bam():

    full_path = get_example_path("control.bam")

    return pr.read_bam(full_path)


def cpg():
    full_path = get_example_path("cpg.bed")

    df = pd.read_csv(
        full_path,
        sep="\t",
        header=None,
        names="Chromosome Start End CpG".split())

    return pr.PyRanges(df)


def exons():

    full_path = get_example_path("exons.bed")

    return pr.read_bed(full_path)
