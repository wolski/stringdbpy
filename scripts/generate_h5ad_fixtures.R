#!/usr/bin/env Rscript
# Generate AnnData rank-source fixtures, one per prolfqua modelling facade.
#
# Each fixture is a zip shaped like a B-Fabric DEA result archive, holding the
# two artifacts string_gsea cares about: the `AnnData.h5ad` it reads, and the
# `.rnk` files prolfquapp derived from the same contrasts. The Python tests
# assert the AnnData rank source reproduces those `.rnk` values, which is what
# keeps one definition of "the rank" across the two languages.
#
# Usage:
#   Rscript scripts/generate_h5ad_fixtures.R [output_dir] [model ...]

suppressPackageStartupMessages({
  library(prolfquapp)
})
# prolfquasaint registers the `saint` facade from its .onLoad.
suppressPackageStartupMessages(
  saint_available <- requireNamespace("prolfquasaint", quietly = TRUE)
)

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[[1]] else "tests/data/h5ad"
requested <- if (length(args) >= 2) args[-1] else character()

dataset <- system.file(
  "application/sim_test/dataset_sim.csv",
  package = "prolfquapp"
)
stopifnot(nchar(dataset) > 0)

facades <- prolfqua::list_facades()
models <- if (length(requested) > 0) requested else sort(names(facades))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Keep only the artifacts the fixture is for; a full result directory carries
# parquet, xlsx and fasta copies that would bloat the repository.
fixture_patterns <- c("[.]h5ad$", "[.]rnk$")

build_one <- function(model) {
  workdir <- tempfile(paste0("fixture-", model, "-"))
  dir.create(workdir, recursive = TRUE)
  config <- prolfquapp::make_DEA_config_R6(
    PATH = workdir,
    WORKUNITID = paste0("FIX", model),
    PROJECTID = "3000",
    ORDERID = "6200",
    Normalization = "none",
    # Thresholds do not affect the ranks; permissive values keep the ORA and
    # significance writers from erroring on an empty selection.
    FDR_threshold = 1,
    diff_threshold = 0,
    model = model
  )
  result <- suppressWarnings(prolfquapp::run_dea(
    indir = tempdir(),
    dataset = dataset,
    # A nested facade switches itself to SIM_PEPTIDE.
    software = "prolfquapp.SIM",
    config = config
  ))
  deanalyse <- result$deanalyse
  config$software <- result$software
  reporter <- prolfquapp::DEAReportGenerator$new(deanalyse, config, name = "")
  reporter$write_DEA_all(boxplot = FALSE)
  se <- reporter$make_SummarizedExperiment()
  prolfquapp:::write_summarized_experiment_h5ad(
    se,
    file.path(reporter$resultdir, "AnnData.h5ad")
  )

  keep <- unlist(lapply(
    fixture_patterns,
    function(pattern) {
      list.files(reporter$resultdir, pattern = pattern, full.names = TRUE)
    }
  ))
  # Resolve through the directory: normalizePath() leaves a path to a file that
  # does not exist yet relative, and the zip below runs from the result folder.
  archive <- file.path(
    normalizePath(output_dir, mustWork = TRUE),
    paste0(model, ".zip")
  )
  unlink(archive)
  owd <- setwd(reporter$resultdir)
  on.exit(setwd(owd), add = TRUE)
  status <- utils::zip(archive, basename(keep), flags = "-q")
  if (!identical(status, 0L) || !file.exists(archive)) {
    stop("could not write fixture archive: ", archive)
  }

  cfg <- deanalyse$contrast_results[[deanalyse$default_model]]$get_config()
  contrasts <- deanalyse$contrast_results[[deanalyse$default_model]]
  frame <- contrasts$get_contrasts()
  data.frame(
    model = model,
    status = "ok",
    rnk_files = sum(grepl("[.]rnk$", basename(keep))),
    rows = nrow(frame),
    contrast_col = cfg$contrast_col,
    effect_col = cfg$effect_col,
    score_col = cfg$score_col,
    pvalue_col = if (cfg$has_pvalue()) cfg$pvalue_col else NA_character_,
    fdr_col = cfg$fdr_col,
    directional = cfg$significance_directional,
    model_names = paste(unique(frame[[cfg$model_name_col]]), collapse = "|"),
    estimate_types = if ("estimate_type" %in% colnames(frame)) {
      paste(unique(frame$estimate_type), collapse = "|")
    } else {
      NA_character_
    },
    size_kb = round(file.size(archive) / 1024),
    stringsAsFactors = FALSE
  )
}

rows <- list()
for (model in models) {
  if (identical(model, "saint") && !saint_available) {
    message("skipping saint: prolfquasaint not installed")
    next
  }
  message("=== ", model)
  row <- tryCatch(
    build_one(model),
    error = function(e) {
      data.frame(
        model = model,
        status = paste("FAILED:", conditionMessage(e)),
        rnk_files = NA_integer_,
        rows = NA_integer_,
        contrast_col = NA_character_,
        effect_col = NA_character_,
        score_col = NA_character_,
        pvalue_col = NA_character_,
        fdr_col = NA_character_,
        directional = NA,
        model_names = NA_character_,
        estimate_types = NA_character_,
        size_kb = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  )
  rows[[model]] <- row
}

summary_table <- do.call(rbind, rows)
print(summary_table, row.names = FALSE)
write.csv(
  summary_table,
  file.path(output_dir, "generation_summary.csv"),
  row.names = FALSE
)
