import doctest, os, sys, subprocess


thepath = "/home/anandia/pyranges/"
sys.path.insert(0, thepath)
from pyranges import *

def test_tutorial():
	subprocess.run(["curl", "-O", "https://mariottigenomicslab.bio.ub.edu/pyranges_data/pyranges_tutorial_data.tar.gz"])
	subprocess.run(["tar", "zxf", "pyranges_tutorial_data.tar.gz"])
	subprocess.run(["pip", "install", "pyfaidx"])

	os.chdir('../tests/')
	failure_count2, test_count2 = doctest.testfile('/home/anandia/pyranges/docs/tutorial.rst')
	if not failure_count2:
	    print('All tests in the tutorial were successful!')

	print()

	os.remove("pyranges_tutorial_data.tar.gz") 
	os.remove("Dgyro_annotation.canonical_CDS.gtf") 
	#os.remove("Dgyro_canonical_CDS.fa") 
	os.remove("Dgyro_canonical_CDS.seq.tsv") 
	os.remove("Dgyro_genome.fa.fai") 
	os.remove("Dgyro_annotation.gff")
	os.remove("Dgyro_genome.fa")

