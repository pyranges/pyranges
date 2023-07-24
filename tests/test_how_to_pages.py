import doctest, os, sys, subprocess
#unittest, 

thepath = "/home/anandia/pyranges/"
sys.path.insert(0, thepath)
from pyranges import *



subprocess.run(["pip", "install", "pyranges_db"])
subprocess.run(["pip", "install", "pyBigWig"])
subprocess.run(["pip", "install", "pyfaidx"])


os.chdir('../tests/')
failure_count2, test_count2 = doctest.testfile('../docs/how_to_pages.rst')
if not failure_count2:
    print('All tests in the how_to_pages were successful!')

print()


os.remove("minigenome.fa") 
os.remove("minigenome.fa.fai") 
os.remove("chipseq.gtf") 


    
        
        
