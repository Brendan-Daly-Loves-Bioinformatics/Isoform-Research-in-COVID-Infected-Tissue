# Part 5: Gene Ontology Enrichment Analysis
# identify biological processes associated with significant genes

# install required packages if needed
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}

if (!requireNamespace("enrichplot", quietly = TRUE)) {
  BiocManager::install("enrichplot")
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}

# load packages
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

# convert significant Ensembl IDs to Entrez IDs

# remove Ensembl version numbers
sig_ensembl <- sub(
  "\\..*$",
  "",
  sig_genes$ensembl_id
)

# convert Ensembl IDs to Entrez IDs
entrez_ids <- bitr(
  sig_ensembl,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# keep only one-to-one Ensembl to Entrez mappings

# count how many Entrez IDs each Ensembl gene maps to
mapping_counts <- table(entrez_ids$ENSEMBL)

# keep only genes with exactly one Entrez mapping
one_to_one_ensembl <- names(
  mapping_counts[mapping_counts == 1]
)

# create clean Entrez ID table
entrez_ids_clean <- entrez_ids[
  entrez_ids$ENSEMBL %in% one_to_one_ensembl,
]

# remove duplicate Entrez IDs
entrez_ids_clean <- entrez_ids_clean[
  !duplicated(entrez_ids_clean$ENTREZID),
]

# create the GO background universe

# get all genes tested by DESeq2
tested_ensembl <- rownames(dds)

# remove Ensembl version numbers
tested_ensembl <- sub(
  "\\..*$",
  "",
  tested_ensembl
)

# convert tested Ensembl IDs to Entrez IDs
tested_entrez <- bitr(
  tested_ensembl,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# remove duplicate Ensembl IDs
tested_entrez <- tested_entrez[
  !duplicated(tested_entrez$ENSEMBL),
]

# remove duplicate Entrez IDs
tested_entrez <- tested_entrez[
  !duplicated(tested_entrez$ENTREZID),
]

# perform GO Biological Process enrichment analysis

go_results <- enrichGO(
  gene = entrez_ids_clean$ENTREZID,
  universe = tested_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# display results

# view the top enriched biological processes
head(go_results)

# count enriched GO Biological Processes
cat(
  "Number of enriched GO Biological Processes:",
  nrow(as.data.frame(go_results)),
  "
"
)

# number of significant genes with clean mappings
cat(
  "Significant genes with one-to-one Entrez mappings:",
  nrow(entrez_ids_clean),
  "
"
)

# number of genes in the GO background
cat(
  "Genes in GO background universe:",
  nrow(tested_entrez),
  "
"
)

# save GO enrichment results

go_results_table <- as.data.frame(go_results)

write.csv(
  go_results_table,
  "GO_Biological_Process_enrichment.csv",
  row.names = FALSE
)

# plot the top enriched biological processes
png("go_dotplot.png", width = 1400, height = 1000, res = 150)
dotplot(
  go_results,
  showCategory = 15
)
dev.off()