repos <- c(
  CRAN     = "https://packagemanager.posit.co/cran/__linux__/noble/latest",
  BioCsoft = "https://packagemanager.posit.co/bioconductor/3.22",
  BioCann  = "https://packagemanager.posit.co/bioconductor/3.22"
)
options(repos = repos)

install.packages(c(
  "BiocManager", "remotes",
  "Rcpp", "yulab.utils", "ggplot2", "ggupset", "patchwork", "DT",
  "jsonlite", "knitr", "rmarkdown", "quarto",
  "dplyr", "tidyr", "tibble", "pillar", "cli", "rlang", "vctrs",
  "lifecycle", "scales", "reshape2", "igraph", "purrr", "plyr",
  "stringr", "ggnewscale", "ggrepel", "ggfun", "ggplotify",
  "ggforce", "gridGraphics", "cowplot", "aplot"
))

install.packages(c(
  "AnnotationDbi", "GO.db", "GOSemSim",
  "BiocParallel", "fgsea", "DOSE", "enrichplot"
))

stopifnot(
  requireNamespace("DOSE"),
  requireNamespace("enrichplot")
)
