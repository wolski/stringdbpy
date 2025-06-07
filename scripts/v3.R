# Install if needed:
# install.packages("gson")
# if (!requireNamespace("DOSE", quietly=TRUE)) BiocManager::install("DOSE")
# if (!requireNamespace("clusterProfiler", quietly=TRUE)) BiocManager::install("clusterProfiler")

library(gson) # for constructing GSON objects :contentReference[oaicite:0]{index=0}
library(DOSE) # defines enrichResult class :contentReference[oaicite:1]{index=1}
library(clusterProfiler) # for downstream plotting (optional)
Over <- data.frame(
  ID = c("TERM_A", "TERM_B", "TERM_C"),
  Description = c(
    "Simulated Pathway A",
    "Simulated Pathway B",
    "Simulated Pathway C"
  ),
  GeneRatio = c("5/20", "3/20", "4/20"),
  BgRatio = c("10/100", "8/100", "12/100"),
  RichFactor = c(0.25, 0.15, 0.20),
  FoldEnrichment = c(2.5, 1.875, 1.667),
  zScore = c(2.3, 1.7, 1.9),
  pvalue = c(0.01, 0.05, 0.02),
  p.adjust = c(0.015, 0.075, 0.03),
  qvalue = c(0.02, 0.08, 0.035),
  geneID = c(
    "G1/G2/G3/G4/G5",
    "G10/G11/G12",
    "G20/G21/G22/G23"
  ),
  Count = c(5L, 3L, 4L),
  row.names = c("TERM_A", "TERM_B", "TERM_C"),
  stringsAsFactors = FALSE
)
# The input gene list (as.character) used to generate Over:
gene <- unique(unlist(strsplit(Over$geneID, "/")))

# A background universe of 100 genes:
universe <- paste0("G", 1:100)

# Fake geneSets matching the IDs in Over:
geneSets <- list(
  TERM_A = paste0("G", 1:10),
  TERM_B = paste0("G", 10:17),
  TERM_C = paste0("G", 20:31)
)
# 1) Prepare the two mapping data.frames:
gsid2gene_df <- data.frame(
  gsid = rep(names(geneSets), lengths(geneSets)),
  gene = unlist(geneSets),
  stringsAsFactors = FALSE
)
gsid2name_df <- data.frame(
  gsid = c("TERM_A", "TERM_B", "TERM_C"),
  name = c("Simulated Pathway A", "Simulated Pathway B", "Simulated Pathway C"),
  stringsAsFactors = FALSE
)

# 2) Call gson() to get a GSON object:
USER_DATA <- gson(
  gsid2gene     = gsid2gene_df,
  gsid2name     = gsid2name_df,
  gene2name     = NULL,
  species       = "Homo sapiens", # → organism slot
  gsname        = "Custom;myonto", # → ontology = "myonto"
  version       = "v1.0",
  accessed_date = Sys.Date() %>% as.character(),
  keytype       = "ENTREZID",
  info          = ""
)
# Your chosen cutoffs & method:
pvalueCutoff <- 0.05
pAdjustMethod <- "BH"
qvalueCutoff <- 0.2

# 1) Instantiate with UNKNOWN metadata:
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

# 2) Mimic the post‐hoc block in enricher_internal():
if (inherits(USER_DATA, "GSON")) {
  if (!is.null(USER_DATA@keytype)) x@keytype <- USER_DATA@keytype
  if (!is.null(USER_DATA@species)) x@organism <- USER_DATA@species
  if (!is.null(USER_DATA@gsname)) x@ontology <- sub(".*;", "", USER_DATA@gsname)
}

# 3) Inspect your new object:
print(x@result)
message("Organism: ", x@organism) # "Homo sapiens"
message("Ontology: ", x@ontology) # "myonto"
