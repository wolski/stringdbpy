Downloaded from:
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip

keeping files nodes.dmp and readme.md from archive

```
ncbiNodes <- read_delim("../../data/mappings/nodes.dmp", delim= "\t|\t", col_names = FALSE)

column_names = c(
  "taxon_id",
  "parent_taxon_id",
  "rank",
  "code",
  "field5",
  "field6",
  "field7",
  "field8",
  "field9",
  "field10",
  "field11",
  "field12",
  "comment"
)
colnames(ncbiNodes) <- column_names
ncbiNodes <- ncbiNodes[,1:3]


readr::write_tsv(ncbiNodes, file = "../../data/mappings/nodes.tsv")


```

both zipped into the NCBI_nodes.zip