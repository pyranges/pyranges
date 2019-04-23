
mpl = reticulate::import("matplotlib")
mpl$use('TkAgg')

library(reticulate)

use_python("/mnt/work/endrebak/software/anaconda/bin/python")

library(bookdown)

render_book("bookdown::html_document2")
