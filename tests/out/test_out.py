def test_write_bed(chip_10, tmp_path):
    outfile = tmp_path / "deleteme.bed"
    chip_10.to_bed(outfile)


def test_write_bed_no_path(chip_10):
    result = chip_10.to_bed()
    assert isinstance(result, str)


def test_write_gtf(chip_10, tmp_path):
    outfile = tmp_path / "deleteme.gtf"
    chip_10.to_gtf(outfile)


def test_write_gff3(chip_10, tmp_path):
    outfile = tmp_path / "deleteme.gff3"
    chip_10.to_gff3(outfile)


def test_write_gtf_no_path(chip_10):
    result = chip_10.to_gtf()
    assert isinstance(result, str)


def test_write_bigwig(chip_10, tmp_path, chromsizes):
    try:
        outfile = tmp_path / "deleteme.bigwig"
        outpath = str(outfile)
        chip_10.to_bigwig(outpath, chromosome_sizes=chromsizes)
    except SystemExit:
        pass
