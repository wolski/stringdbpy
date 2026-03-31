#' Build a compareClusterResult from multiple contrasts
#'
#' Merges `enrichResult` objects across contrasts (for a single category)
#' into a [DOSE::compareClusterResult-class] object. This enables
#' cross-contrast comparison plots from enrichplot (e.g., `dotplot()`,
#' `emapplot()`).
#'
#' @param result_list Named list of [DOSE::enrichResult-class] objects,
#'   one per contrast. Names are used as cluster labels.
#' @return A [DOSE::compareClusterResult-class] object.
#' @importFrom methods new slot
#' @importClassesFrom DOSE compareClusterResult
#' @export
#' @examples
#' \dontrun{
#' results <- read_gsea_json("gsea_result.json")
#' # Extract GO Process from each contrast
#' go_list <- lapply(results, function(x) x[["GO Process"]])
#' ccr <- build_compareClusterResult(go_list)
#' enrichplot::dotplot(ccr, showCategory = 10)
#' }
build_compareClusterResult <- function(result_list) {
  stopifnot(length(result_list) >= 2)

  # Merge @result data.frames, adding Cluster column
  dfs <- mapply(function(er, nm) {
    df <- methods::slot(er, "result")
    df$Cluster <- nm
    df
  }, result_list, names(result_list), SIMPLIFY = FALSE)
  merged <- do.call(rbind, dfs)
  merged$Cluster <- factor(merged$Cluster, levels = names(result_list))

  # Build geneClusters: named list of contrast -> gene vector
  gene_clusters <- lapply(result_list, function(er) methods::slot(er, "gene"))

  methods::new(
    "compareClusterResult",
    compareClusterResult = merged,
    geneClusters = gene_clusters,
    fun = "stringGSEA",
    .call = match.call(),
    keytype = "STRING",
    readable = TRUE
  )
}
