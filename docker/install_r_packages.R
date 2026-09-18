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
  "blob", "DBI", "RSQLite", "memoise",
  # protsea's Suggests, minus its dev tools: the image BUILDS protsea's
  # vignette and ships its report templates, so everything they touch is a
  # build-time requirement here. ggridges is not called as `ggridges::` --
  # enrichplot::ridgeplot() loads it internally -- so grepping the sources for
  # namespace calls does not reveal it; keep this list against protsea's
  # DESCRIPTION instead.
  "ggridges", "rprojroot"
))

# Bioconductor packages — use BiocManager repos, update=FALSE to avoid
# re-downloading the CRAN packages we just installed above
BiocManager::install(
  c("AnnotationDbi", "GO.db", "GOSemSim", "BiocParallel", "fgsea", "DOSE", "enrichplot",
    # protsea's round-trip vignette runs clusterProfiler::GSEA(), and the image
    # builds that vignette, so this is a build-time requirement, not optional.
    "clusterProfiler"),
  ask = FALSE, update = FALSE
)

# Fail this layer, rather than a later vignette render, when a package protsea
# needs is absent.
for (pkg in c(
  "DOSE", "enrichplot", "clusterProfiler", "DT", "ggplot2", "ggridges",
  "ggupset", "knitr", "patchwork", "quarto", "rmarkdown", "rprojroot"
)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("required package not installed: ", pkg)
  }
}
