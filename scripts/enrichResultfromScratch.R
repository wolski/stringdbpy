# ------------------------------------------
# Build enrichResult “By Hand” (including Over)
# ------------------------------------------
rm(list = ls())
# 1. Load required package (defines enrichResult class)
#    If you’re using clusterProfiler instead, load that instead.
library(DOSE)

# 2. Simulate an enrichment result table 'Over'
#    (ID, Description, GeneRatio, BgRatio, RichFactor, FoldEnrichment,
#     zScore, pvalue, p.adjust, qvalue, geneID, Count)
Over <- data.frame(
  ID = c("PATH_A", "PATH_B", "PATH_C"),
  Description = c("Sim Pathway A", "Sim Pathway B", "Sim Pathway C"),
  GeneRatio = c("6/20", "4/20", "5/20"),
  BgRatio = c("10/100", "12/100", "8/100"),
  RichFactor = c(0.30, 0.20, 0.25),
  FoldEnrichment = c(3.0, 1.67, 2.50),
  zScore = c(2.8, 1.5, 2.1),
  pvalue = c(0.005, 0.04, 0.01),
  p.adjust = c(0.0075, 0.06, 0.015),
  qvalue = c(0.01, 0.07, 0.02),
  geneID = c(
    "G1/G2/G3/G4/G5/G6",
    "G10/G11/G12/G13",
    "G20/G21/G22/G23/G24"
  ),
  Count = c(6L, 4L, 5L),
  row.names = c("PATH_A", "PATH_B", "PATH_C"),
  stringsAsFactors = FALSE
)

# 3. Define the original gene list and universe
gene <- c(
  "G1", "G2", "G3", "G4", "G5", "G6",
  "G10", "G11", "G12", "G13",
  "G20", "G21", "G22", "G23", "G24"
)
universe <- paste0("G", 1:100)

# 4. Define geneSets as a named list
geneSets <- list(
  PATH_A = paste0("G", 1:10), # first 10 genes
  PATH_B = paste0("G", 10:20), # genes 10–20
  PATH_C = paste0("G", 20:35) # genes 20–35
)

# 5. Set your cutoffs and methods
pvalueCutoff <- 0.05
pAdjustMethod <- "BH"
qvalueCutoff <- 0.2

# 6. (Optional) Fake a minimal GSON for metadata
setClass("GSON", slots = c(
  gsid2gene     = "data.frame",
  gsid2name     = "data.frame",
  gene2name     = "data.frame",
  species       = "character",
  gsname        = "character",
  version       = "character",
  accessed_date = "character",
  keytype       = "character",
  info          = "character"
))
USER_DATA <- new("GSON",
  gsid2gene     = data.frame(),
  gsid2name     = data.frame(),
  gene2name     = data.frame(),
  species       = "Homo sapiens",
  gsname        = "CUSTOM;myonto",
  version       = "1.0",
  accessed_date = Sys.Date() %>% as.character(),
  keytype       = "ENTREZID",
  info          = ""
)

# 7. Construct the enrichResult object
x <- new("enrichResult",
  result         = Over,
  pvalueCutoff   = pvalueCutoff,
  pAdjustMethod  = pAdjustMethod,
  qvalueCutoff   = qvalueCutoff,
  gene           = as.character(gene),
  universe       = universe,
  geneSets       = geneSets,
  organism       = "UNKNOWN",
  keytype        = "UNKNOWN",
  ontology       = "UNKNOWN",
  readable       = FALSE
)

# 8. Pull metadata from USER_DATA (exactly as enricher_internal does)
if (inherits(USER_DATA, "GSON")) {
  if (!is.null(USER_DATA@keytype)) x@keytype <- USER_DATA@keytype
  if (!is.null(USER_DATA@species)) x@organism <- USER_DATA@species
  if (!is.null(USER_DATA@gsname)) x@ontology <- sub(".*;", "", USER_DATA@gsname)
}

# 9. Inspect
print(x)
print(x@result)
message("Organism: ", x@organism)
message("Ontology: ", x@ontology)
result <- x
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
