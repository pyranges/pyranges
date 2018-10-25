
mpl = reticulate::import("matplotlib")
mpl$use('TkAgg')

library(reticulate)

use_python("/mnt/work/endrebak/software/anaconda/bin/python")

library(bookdown)
## gitbook()
render_book("bookdown::gitbook")
