## Keep the vendored Quarto assets used by the report templates in step with the
## installed fgczQuartoTemplate package.
##
## The reports refer to the assets by bare filename and are rendered by the
## `quarto` CLI (from the string_gsea workflow, or by the vignette builder), so
## nothing calls fgczQuartoTemplate::fgcz_render() to stage them at runtime.
## Copy and verify them here before a package or vignette build.
##
##   Rscript data-raw/sync_quarto_assets.R

assets <- c(
  "_metadata.yml",
  "fgcz.scss",
  "fgcz_header_quarto.html",
  "fgcz-plot-finder.html"
)
target_dir <- normalizePath("vignettes", mustWork = TRUE)
source_files <- fgczQuartoTemplate::fgcz_quarto_dir(assets)
target_files <- file.path(target_dir, assets)

fgczQuartoTemplate::fgcz_copy_assets(target_dir)

if (
  !identical(
    unname(tools::md5sum(source_files)),
    unname(tools::md5sum(target_files))
  )
) {
  stop(
    "FGCZ Quarto assets differ after synchronization: ",
    paste(assets, collapse = ", ")
  )
}

message(
  "OK: synchronized FGCZ Quarto assets into vignettes/: ",
  paste(assets, collapse = ", ")
)
