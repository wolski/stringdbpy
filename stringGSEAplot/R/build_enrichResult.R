#' Resolve STRING protein IDs to human-readable gene labels
#'
#' @param gene_ids Character vector of STRING protein IDs
#'   (e.g., "9606.ENSP00000269305").
#' @param gene_pool Named list from the JSON, keyed by protein_id.
#'   Each element has fields: protein_id, label, input_label, input_value, rank.
#' @return Character vector of gene labels (same length as gene_ids).
#'   IDs not found in the pool are returned as-is.
#' @keywords internal
resolve_gene_ids <- function(gene_ids, gene_pool) {
  vapply(gene_ids, function(gid) {
    entry <- gene_pool[[gid]]
    if (is.null(entry)) gid else entry[["label"]]
  }, character(1), USE.NAMES = FALSE)
}


#' Build an enrichResult S4 object from one category's JSON data
#'
#' @param category_data List from JSON for one category (has fields:
#'   category, contrast, terms).
#' @param gene_pool Named list from the JSON gene_pool for this contrast.
#' @param rank_list Named list from rank_lists for this contrast
#'   (has fields: contrast, entries).
#' @return An [DOSE::enrichResult-class] object.
#' @importFrom methods new
#' @importClassesFrom DOSE enrichResult
#' @export
build_enrichResult <- function(category_data, gene_pool, rank_list) {
  terms <- category_data[["terms"]]
  if (length(terms) == 0) {
    stop("No terms found in category '", category_data[["category"]], "'")
  }

  # GeneRatio denominator must be constant per cluster (= gene pool size)
  # to match clusterProfiler's convention. Term set size goes in BgRatio.
  n_genes <- length(gene_pool)

  # Build @result data.frame (one row per term)
  rows <- lapply(terms, function(term) {
    gene_ids <- unlist(term[["gene_ids"]])
    labels <- resolve_gene_ids(gene_ids, gene_pool)
    k <- term[["genes_mapped"]]
    n_set <- term[["genes_in_set"]]
    data.frame(
      ID = term[["term_id"]],
      Description = term[["description"]],
      GeneRatio = paste0(k, "/", n_genes),
      BgRatio = paste0(n_set, "/", n_set),
      pvalue = term[["fdr"]],
      p.adjust = term[["fdr"]],
      qvalue = term[["fdr"]],
      geneID = paste(labels, collapse = "/"),
      Count = k,
      stringsAsFactors = FALSE
    )
  })
  result_df <- do.call(rbind, rows)
  rownames(result_df) <- result_df$ID

  # Build @geneSets: named list of term_id -> character vector of labels
  gene_sets <- lapply(terms, function(term) {
    resolve_gene_ids(unlist(term[["gene_ids"]]), gene_pool)
  })
  names(gene_sets) <- vapply(terms, function(t) t[["term_id"]], character(1))

  # Build @gene: all gene labels from this contrast's gene_pool
  gene_labels <- vapply(
    gene_pool,
    function(entry) entry[["label"]],
    character(1),
    USE.NAMES = FALSE
  )

  # Build @universe: all ranked gene identifiers
  universe <- names(rank_list[["entries"]])

  methods::new(
    "enrichResult",
    result = result_df,
    pvalueCutoff = 1.0,
    pAdjustMethod = "BH",
    qvalueCutoff = 1.0,
    gene = gene_labels,
    universe = universe,
    geneSets = gene_sets,
    organism = "unknown",
    ontology = category_data[["category"]],
    keytype = "STRING",
    readable = TRUE
  )
}
