
mpl = reticulate::import("matplotlib")
mpl$use('TkAgg')

library(reticulate)

library(bookdown)

sessionInfo()

render_book("index.Rmd", "bookdown::gitbook")
