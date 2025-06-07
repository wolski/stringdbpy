# --------------------------------------------
# test_enricher_internal_signif.R
# --------------------------------------------

# 1. Install/load DOSE (if not already):
#    BiocManager::install("DOSE")
library(DOSE)


# 2. Create a minimal “GSON‐like” S4 class on the fly
#    (we only need those slots that enricher_internal() checks:
#     @gsid2gene, @gsid2name, @species, and @keytype).
#    The real GSON class is defined in the 'gson' package, but here
#    we mimic it so that enricher_internal() will treat our USER_DATA
#    as "inherits(USER_DATA, 'GSON') == TRUE".
setClass("GSON",
  slots = c(
    gsid2gene = "data.frame",
    gsid2name = "data.frame",
    species   = "character",
    keytype   = "character",
    gsname    = "character"
  )
)

# 3. Build 100 synthetic genes: GENE1 ... GENE100
all_genes <- paste0("GENE", seq_len(100))

# 4. Define 10 overlapping gene sets (GS1–GS10)
#    * GS1 = GENE1:15, GS2 = GENE10:25, ..., GS10 = GENE90:100
PATHID2EXTID <- list(
  GS1  = paste0("GENE", 1:15),
  GS2  = paste0("GENE", 10:25),
  GS3  = paste0("GENE", 20:35),
  GS4  = paste0("GENE", 30:45),
  GS5  = paste0("GENE", 40:55),
  GS6  = paste0("GENE", 50:65),
  GS7  = paste0("GENE", 60:75),
  GS8  = paste0("GENE", 70:85),
  GS9  = paste0("GENE", 80:95),
  GS10 = paste0("GENE", 90:100)
)

# 5. Give each gene set a name (“Description”)
PATHID2NAME <- data.frame(
  gsid = paste0("GS", 1:10),
  name = paste("Pathway", LETTERS[1:10]),
  stringsAsFactors = FALSE
)

# 6. Build an S4 “GSON” object with exactly the slots that
#    enricher_internal() will inspect under `inherits(USER_DATA, "GSON")`.
#
#    - gsid2gene: a data.frame with columns (gsid, gene)
#    - gsid2name: a data.frame with columns (gsid, name)
#    - species:    a character (e.g. "human")
#    - keytype:    any non‐NULL character (we can leave it "ENTREZID")
gsid2gene_df <- do.call(rbind, lapply(names(PATHID2EXTID), function(id) {
  data.frame(
    gsid = id,
    gene = PATHID2EXTID[[id]],
    stringsAsFactors = FALSE
  )
}))
gsid2name_df <- PATHID2NAME

# Now construct the GSON object
USER_DATA <- new("GSON",
  gsid2gene = gsid2gene_df,
  gsid2name = gsid2name_df,
  species   = "human", # <-- this will set @organism
  keytype   = "ENTREZID",
  gsname    = "dummy;myonto"
) # <-- used if enrichResult wants keytype

# 7. Define our query: choose 10 genes that all lie in GS1
#    (so that k = 10, M = 15, N = 100, n = 10 ⇒ highly significant)
gene_list <- paste0("GENE", c(1:10, 100:91, 40:31, 50:41))

# 8. Call enricher_internal() with minGSSize=1 so GS1 (size=15) passes,
#    and allow p‐valueCutoff = 1 so we still see everything before filtering.
result <- DOSE:::enricher_internal(
  gene          = gene_list,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  universe      = NULL, # we let it default to union of all genes in PATHID2EXTID
  minGSSize     = 1, # allow very small sets
  maxGSSize     = 500,
  qvalueCutoff  = 1,
  USER_DATA     = USER_DATA # our fake GSON object
)

# 9. Inspect the enrichResult
if (is.null(result)) {
  stop("No enrichResult returned (something went wrong).")
}

# 10. Print the raw result table:
print(result@result)

# 11. Confirm that organism/ontology slots came from USER_DATA:
message(">>> organism slot: ", result@organism) # should be "human"
message(">>> ontology slot: ", result@ontology) # should be "myonto"?
# But note: enricher_internal()
# expects @gsname to be something like
# "prefix;ontology", so we need to set gsid2name
# accordingly if we want "myonto".
#
# In our case, PATHID2NAME$name are just “Pathway A…J”.
# enricher_internal() does:
#   x@ontology <- gsub(".*;","", USER_DATA@gsname)
# so USER_DATA@gsname must contain a “;”.
#
# To demonstrate, we could instead set:
#     USER_DATA@gsid2name$name <- paste0("fake;", PATHID2NAME$name)
# and then ontology would be "Pathway A", etc.
# For now, since our gsid2name$name had no “;”,
# x@ontology will remain "UNKNOWN".
#
# If you want exactly @ontology = "myonto", then do:
# USER_DATA@gsname <- "something;myonto"
# after constructing USER_DATA above.

# 12. For clarity, let’s explicitly *force* USER_DATA@gsname so that ontology = "myonto":
USER_DATA@gsname <- "dummy;myonto"
# Now rerun enricher_internal() to see the effect:
result2 <- DOSE:::enricher_internal(
  gene          = gene_list,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  universe      = NULL,
  minGSSize     = 1,
  maxGSSize     = 500,
  qvalueCutoff  = 1,
  USER_DATA     = USER_DATA
)

if (!is.null(result2)) {
  print(result2@result)
  message(">>> organism slot: ", result2@organism) # “human”
  message(">>> ontology slot: ", result2@ontology) # “myonto”
} else {
  stop("Second call to enricher_internal() failed.")
}

library(enrichplot)
library(ggplot2)
barplot(result)
barplot(result,
  showCategory = 20,
  title = "Top 20 Enriched Terms",
  x = "Count"
) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dotplot(result)
library(enrichplot)
cnetplot(result)
heatplot(result)
edox2 <- pairwise_termsim(result)
p1 <- treeplot(edox2)
p1
edo <- pairwise_termsim(result)
p1 <- emapplot(edo)
p1
