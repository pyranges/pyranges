# import pytest

# from io import StringIO

# import pyranges as pr

# from pyranges.guessers import (parse_file, guess_columns, guess_chromosome, guess_strand, guess_start_end, guess_delim)


# @pytest.fixture
# def bed_df():

#     return pr.data.chipseq().df


# @pytest.fixture
# def bed_file():

#     return pr.get_example_path("chipseq") + ".bed"

# @pytest.fixture
# def gtf_file():

#     return pr.get_example_path("ensembl.gtf")

# @pytest.fixture
# def number_unique(bed_df):

#     number_unique = {c: len(bed_df[c].unique()) for c in bed_df}

#     return number_unique


# def test_guess_delim(gtf_file):

#     delim = guess_delim(gtf_file).encode()
#     assert delim == " ".encode() # this is the wrong guess


# def test_guess_strand(bed_df, number_unique):

#     result = guess_strand(bed_df, number_unique, "Start", "End")

#     assert result == "Strand"


# def test_guess_start_end(bed_df, number_unique):

#     result = guess_start_end(bed_df, number_unique)

#     assert result == ("Start", "End")


# def test_guess_chromosome(bed_df, number_unique):

#     result = guess_chromosome(bed_df, number_unique, "Strand")

#     assert result == "Chromosome"


# def test_guess_columns(bed_file):

#     result = guess_columns(bed_file)
#     assert result == (0, 1, 2, 5)


# def test_guess_columns(gtf_file):

#     result = guess_columns(gtf_file, sep="\t", header=False)
#     print(result)
#     assert result == (0, 1, 2, 5)
