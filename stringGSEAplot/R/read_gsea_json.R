#' Read a STRING-GSEA JSON file and return enrichResult objects
#'
#' Reads a JSON file produced by the string-gsea Python pipeline
#' (`GSEAResult.to_json()`) and constructs one
#' [DOSE::enrichResult-class] object per contrast-category combination.
#'
#' @param json_path Path to the JSON file.
#' @return A double-nested named list: `result[["contrast"]][["category"]]`,
#'   where each leaf is an [DOSE::enrichResult-class] object.
#' @export
#' @examples
#' \dontrun{
#' results <- read_gsea_json("WU2848501_gsea_result.json")
#' names(results)  # contrast names
#' names(results[[1]])  # category names
#' enrichplot::dotplot(results[[1]][["SMART"]])
#' }
read_gsea_json <- function(json_path) {
  json_data <- jsonlite::fromJSON(json_path, simplifyVector = FALSE)

  # Validate required top-level keys
  required_keys <- c("data", "rank_lists", "metadata", "links")
  missing <- setdiff(required_keys, names(json_data))
  if (length(missing) > 0) {
    stop(
      "JSON is missing required keys: ", paste(missing, collapse = ", "),
      ". Re-run the pipeline with a current version of string_gsea to ",
      "produce a complete JSON file.",
      call. = FALSE
    )
  }

  contrasts_data <- json_data[["data"]]
  rank_lists <- json_data[["rank_lists"]]

  result_list <- list()
  raw_data    <- list()

  for (contrast_name in names(contrasts_data)) {
    contrast  <- contrasts_data[[contrast_name]]
    gene_pool <- contrast[["gene_pool"]]
    rank_list <- rank_lists[[contrast_name]]
    cats      <- contrast[["categories"]]

    contrast_results <- list()
    for (cat_name in names(cats)) {
      cat_data <- cats[[cat_name]]
      if (length(cat_data[["terms"]]) == 0) next
      contrast_results[[cat_name]] <- build_enrichResult(cat_data, gene_pool, rank_list)
    }
    result_list[[contrast_name]] <- contrast_results
    raw_data[[contrast_name]]    <- list(gene_pool = gene_pool,
                                         rank_list = rank_list,
                                         cats      = cats)
  }

  attr(result_list, "metadata") <- json_data[["metadata"]]
  attr(result_list, "links")    <- json_data[["links"]]
  attr(result_list, "raw_data") <- raw_data

  result_list
}
