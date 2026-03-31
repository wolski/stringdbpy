get_test_json <- function() {
  system.file("extdata", "WU2848501_gsea_result.json.gz",
              package = "stringGSEAplot")
}

test_that("read_gsea_json returns double-nested list of enrichResult", {
  json_path <- get_test_json()
  skip_if_not(nzchar(json_path), "JSON test data not available")

  results <- read_gsea_json(json_path)

  # Top level: contrasts
  expect_type(results, "list")
  expect_true(length(results) > 0)

  # Second level: categories, each an enrichResult
  for (contrast_name in names(results)) {
    contrast <- results[[contrast_name]]
    expect_type(contrast, "list")
    expect_true(length(contrast) > 0)
    for (cat_name in names(contrast)) {
      expect_s4_class(contrast[[cat_name]], "enrichResult")
    }
  }
})


test_that("enrichResult has expected slots and structure", {
  json_path <- get_test_json()
  skip_if_not(nzchar(json_path), "JSON test data not available")

  results <- read_gsea_json(json_path)
  er <- results[[1]][[1]]

  # @result is a data.frame with required columns
  res_df <- slot(er, "result")
  expect_s3_class(res_df, "data.frame")
  expected_cols <- c("ID", "Description", "GeneRatio", "pvalue", "p.adjust",
                     "qvalue", "geneID", "Count")
  for (col in expected_cols) {
    expect_true(col %in% names(res_df), info = paste("Missing column:", col))
  }

  # GeneRatio has "k/n" format with constant denominator
  expect_true(all(grepl("^\\d+/\\d+$", res_df$GeneRatio)))
  denoms <- sapply(strsplit(res_df$GeneRatio, "/"), `[`, 2)
  expect_equal(length(unique(denoms)), 1)

  # geneID is slash-separated labels (not STRING IDs)
  first_genes <- strsplit(res_df$geneID[1], "/")[[1]]
  expect_true(length(first_genes) > 0)
  expect_false(any(grepl("^\\d+\\.ENSP", first_genes)))

  # @gene is a character vector
  expect_type(slot(er, "gene"), "character")
  expect_true(length(slot(er, "gene")) > 0)

  # @geneSets is a named list
  gs <- slot(er, "geneSets")
  expect_type(gs, "list")
  expect_true(length(gs) > 0)
  expect_true(!is.null(names(gs)))

  # @universe is a character vector
  expect_type(slot(er, "universe"), "character")
  expect_true(length(slot(er, "universe")) > 0)
})


test_that("all contrasts and categories are present", {
  json_path <- get_test_json()
  skip_if_not(nzchar(json_path), "JSON test data not available")

  results <- read_gsea_json(json_path)

  # Should have 2 contrasts for the human RNK data
  expect_true(length(results) >= 2)

  # Each contrast should have multiple categories
  for (contrast_name in names(results)) {
    expect_true(length(results[[contrast_name]]) > 1,
                info = paste("Contrast", contrast_name, "has too few categories"))
  }
})
