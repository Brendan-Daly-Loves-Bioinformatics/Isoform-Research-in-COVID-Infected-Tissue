# Part 6: Gene Set Enrichment Analysis (GSEA)
# SARS-CoV-2 Infected A549 RNA-seq
#
# GSEA uses the FULL ranked gene list (all ~20,000 genes), not just the 22
# significant DE genes. This detects coordinated but subtle shifts that ORA
# (Part 5 GO enrichment) cannot — critical when few genes pass significance.
#
# PREREQUISITES: Run in the same R session after Part 3 (DE).
# Required objects in memory:
#   res_df - DESeq2 results data frame with columns:
#            gene_id (Ensembl), gene_symbol, log2FoldChange, padj
#   dds    - DESeq2 dataset (optional, used to extract Wald statistic)
#
# If objects are NOT in memory, uncomment the loading section below.

# ---- Setup ----

# load packages
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(DESeq2)

# ---- OPTIONAL: load objects from RDS if not in memory ----
# res_df <- readRDS("res_df.rds")
# dds <- readRDS("dds.rds")

# ---- Step 1: Build the ranked gene list ----

cat("\n=== Building ranked gene list ===\n")

# Ranking metric priority:
#   1. Wald statistic (stat = log2FC / lfcSE) — captures effect size AND
#      precision, bounded, no extreme skew. This is the ideal GSEA ranking.
#   2. Compute Wald from log2FoldChange / lfcSE if both columns exist.
#   3. Extract from results(dds) if the dds object is available.
#   4. Fall back to log2FoldChange (effect size only, bounded, no p-value skew).
#
# We do NOT use sign(log2FC) * -log10(pvalue) because it produces extreme
# skew when p-values span many orders of magnitude (e.g., -229 to +9.66),
# which biases GSEA toward detecting downregulated gene sets.

gene_ids_clean <- sub("\\..*$", "", res_df$gene_id)

if ("stat" %in% colnames(res_df)) {
  rank_metric <- res_df$stat
  cat("Using Wald statistic (stat column) as ranking metric\n")
} else if ("lfcSE" %in% colnames(res_df)) {
  rank_metric <- res_df$log2FoldChange / res_df$lfcSE
  cat("Using computed Wald statistic (log2FC / lfcSE) as ranking metric\n")
} else if (exists("dds", envir = globalenv())) {
  dds_res <- results(dds)
  dds_ids <- sub("\\..*$", "", rownames(dds_res))
  stat_by_id <- setNames(dds_res$stat, dds_ids)
  rank_metric <- stat_by_id[gene_ids_clean]
  cat("Using Wald statistic extracted from results(dds)\n")
} else {
  rank_metric <- res_df$log2FoldChange
  cat("WARNING: Using log2FoldChange as ranking metric (no stat or lfcSE available).\n")
  cat("This does not account for precision. Consider adding stat/lfcSE to res_df.\n")
}

# named vector — names are Ensembl IDs, values are the ranking metric
names(rank_metric) <- gene_ids_clean

# remove NAs, duplicates, and infinite values
rank_metric <- rank_metric[!is.na(rank_metric)]
rank_metric <- rank_metric[is.finite(rank_metric)]
rank_metric <- rank_metric[!duplicated(names(rank_metric))]

# sort descending — GSEA expects this
gene_list <- sort(rank_metric, decreasing = TRUE)

cat("Ranked gene list:", length(gene_list), "genes\n")
cat("Range:", round(min(gene_list), 2), "to", round(max(gene_list), 2), "\n")

# ---- Step 2: GSEA on GO Biological Process ----

cat("\n=== Running GSEA: GO Biological Process ===\n")

gsea_go <- gseGO(
  geneList      = gene_list,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  minGSSize     = 10,
  maxGSSize     = 500,
  verbose       = TRUE
)

n_gsea_go <- nrow(as.data.frame(gsea_go))
cat("GSEA GO BP results:", n_gsea_go, "enriched gene sets (padj < 0.05)\n")

if (n_gsea_go > 0) {
  go_gsea_df <- as.data.frame(gsea_go)
  
  # count upregulated vs downregulated significant sets
  n_down <- sum(go_gsea_df$NES < 0)
  n_up   <- sum(go_gsea_df$NES > 0)
  cat("  Downregulated (NES < 0):", n_down, "gene sets\n")
  cat("  Upregulated   (NES > 0):", n_up,   "gene sets\n")
  
  # save results table (ordered by NES ascending: most negative first)
  go_gsea_df <- go_gsea_df[order(go_gsea_df$NES, decreasing = FALSE), ]
  write.csv(go_gsea_df, "GSEA_GO_Biological_Process.csv", row.names = FALSE)
  cat("Saved: GSEA_GO_Biological_Process.csv\n")
  
  # top 5 downregulated (most negative NES — strongest suppression)
  cat("\nTop 5 downregulated gene sets (negative NES = suppressed in infected):\n")
  down_sets <- head(go_gsea_df[go_gsea_df$NES < 0, ], 5)
  if (nrow(down_sets) > 0) {
    print(down_sets[, c("Description", "NES", "p.adjust", "setSize")])
  } else {
    cat("  None\n")
  }
  
  # top 5 upregulated (most positive NES — strongest activation)
  # Sort positive-NES sets by NES DESCENDING to show the strongest, not weakest
  cat("\nTop 5 upregulated gene sets (positive NES = activated in infected):\n")
  up_sets <- go_gsea_df[go_gsea_df$NES > 0, ]
  up_sets <- up_sets[order(up_sets$NES, decreasing = TRUE), ]
  up_sets <- head(up_sets, 5)
  if (nrow(up_sets) > 0) {
    print(up_sets[, c("Description", "NES", "p.adjust", "setSize")])
  } else {
    cat("  None\n")
  }
}

# ---- Step 3: GSEA on KEGG Pathways ----

cat("\n=== Running GSEA: KEGG Pathways ===\n")

# KEGG requires Entrez IDs — map from Ensembl
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys       = names(gene_list),
  keytype    = "ENSEMBL",
  column     = "ENTREZID",
  multiVals  = "first"
)

gene_list_entrez <- gene_list[!is.na(entrez_ids)]
names(gene_list_entrez) <- entrez_ids[!is.na(entrez_ids)]

# remove any remaining duplicates
gene_list_entrez <- gene_list_entrez[!duplicated(names(gene_list_entrez))]
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

cat("Genes mapped to Entrez IDs:", length(gene_list_entrez), "\n")

gsea_kegg <- gseKEGG(
  geneList      = gene_list_entrez,
  organism      = "hsa",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  minGSSize     = 10,
  maxGSSize     = 500,
  verbose       = TRUE
)

n_gsea_kegg <- nrow(as.data.frame(gsea_kegg))
cat("GSEA KEGG results:", n_gsea_kegg, "enriched pathways (padj < 0.05)\n")

if (n_gsea_kegg > 0) {
  kegg_gsea_df <- as.data.frame(gsea_kegg)
  
  # count upregulated vs downregulated significant pathways
  n_down_k <- sum(kegg_gsea_df$NES < 0)
  n_up_k   <- sum(kegg_gsea_df$NES > 0)
  cat("  Downregulated (NES < 0):", n_down_k, "pathways\n")
  cat("  Upregulated   (NES > 0):", n_up_k,   "pathways\n")
  
  # save results table (ordered by NES ascending)
  kegg_gsea_df <- kegg_gsea_df[order(kegg_gsea_df$NES, decreasing = FALSE), ]
  write.csv(kegg_gsea_df, "GSEA_KEGG_Pathways.csv", row.names = FALSE)
  cat("Saved: GSEA_KEGG_Pathways.csv\n")
  
  # top 5 downregulated
  cat("\nTop 5 downregulated KEGG pathways:\n")
  down_kegg <- head(kegg_gsea_df[kegg_gsea_df$NES < 0, ], 5)
  if (nrow(down_kegg) > 0) {
    print(down_kegg[, c("Description", "NES", "p.adjust", "setSize")])
  } else {
    cat("  None\n")
  }
  
  # top 5 upregulated (sort by NES DESCENDING for strongest)
  cat("\nTop 5 upregulated KEGG pathways:\n")
  up_kegg <- kegg_gsea_df[kegg_gsea_df$NES > 0, ]
  up_kegg <- up_kegg[order(up_kegg$NES, decreasing = TRUE), ]
  up_kegg <- head(up_kegg, 5)
  if (nrow(up_kegg) > 0) {
    print(up_kegg[, c("Description", "NES", "p.adjust", "setSize")])
  } else {
    cat("  None\n")
  }
}

# ---- Step 4: Save objects for visualization script ----

saveRDS(gsea_go, "gsea_go.rds")
saveRDS(gsea_kegg, "gsea_kegg.rds")
cat("\nSaved: gsea_go.rds, gsea_kegg.rds\n")

cat("\n=== GSEA analysis complete ===\n")
cat("Objects created: gsea_go, gsea_kegg\n")