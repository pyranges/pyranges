import yaml

rmds = [
    "index.Rmd",
    "creating.Rmd",
    "writing.Rmd",
    "subset.Rmd",
    "manipulation.Rmd",
    "concat.Rmd",
    "piping.Rmd",
    "iterate.Rmd",
    "sort.Rmd",
    "summary.Rmd",
    "single_range_methods.Rmd",
    "intersection.Rmd",
    "overlap.Rmd",
    "join.Rmd",
    "nearest.Rmd",
    "jaccard.Rmd",
    "reldist.Rmd",
    "coverage.Rmd",
    "runlengths.Rmd",
    "runlength_dict.Rmd",
    "subsetting_rles.Rmd",
    "subsetting_pyrles.Rmd",
    "multithreading.Rmd",
    "databases.Rmd",
]

from subprocess import call

single_page = "\n".join(open(f).read() for f in rmds[1:])
import re

single_page = re.sub(re.compile("#"), "##", single_page)
single_page = open(rmds[0]).read() + "\n" + single_page
open("one_page.Rmd", "w+").write(single_page)

# cmd = "cat {} > one_page.Rmd".format(" ".join(rmds))
# call("cat {} > one_page.Rmd".format(" ".join(rmds)), shell=True)

call("Rscript compile2.R", shell=True)
