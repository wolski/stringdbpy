options(timeout = 300)

# Use PPM for fast binary CRAN installs (amd64 + arm64 on Ubuntu Noble)
ppm <- "https://packagemanager.posit.co/cran/__linux__/noble/latest"
options(repos = c(CRAN = ppm))

# All CRAN deps of enrichplot/DOSE pre-installed as binaries
install.packages(c(
  "BiocManager", "remotes",
  "Rcpp", "yulab.utils", "ggplot2", "ggupset", "patchwork", "DT",
  "jsonlite", "knitr", "rmarkdown", "quarto",
  "dplyr", "tidyr", "tibble", "pillar", "cli", "rlang", "vctrs",
  "lifecycle", "scales", "reshape2", "igraph", "purrr", "plyr",
  "stringr", "ggnewscale", "ggrepel", "ggfun", "ggplotify",
  "ggforce", "gridGraphics", "cowplot", "aplot", "httr",
  "blob", "DBI", "RSQLite", "memoise"
))

# Bioconductor packages — use BiocManager repos, update=FALSE to avoid
# re-downloading the CRAN packages we just installed above
BiocManager::install(
  c("AnnotationDbi", "GO.db", "GOSemSim", "BiocParallel", "fgsea", "DOSE", "enrichplot"),
  ask = FALSE, update = FALSE
)

stopifnot(
  requireNamespace("DOSE"),
  requireNamespace("enrichplot")
)
