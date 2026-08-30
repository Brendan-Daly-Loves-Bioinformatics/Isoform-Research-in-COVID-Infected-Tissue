# Part 7: Visualization of SARS-CoV-2 Infected A549 RNA-seq Results
# Publication-ready figures for QC, DE, AS (null result), GO ORA, and GSEA
#
# PREREQUISITES: Run in the same R session after Parts 3-6.
# Required objects in memory:
#   vsd         - DESeq2 variance-stabilized transform (Part 3)
#   dds         - DESeq2 dataset (Part 3)
#   res_df      - DESeq2 results data frame (Part 3)
#   sig_table   - significant DE genes table (Part 3)
#   dpsi        - MAJIQ deltapsi results data frame (Part 4)
#   go_results  - enrichGO result object (Part 5, ORA)
#   gsea_go     - gseGO result object (Part 6, GSEA on GO BP)
#   gsea_kegg   - gseKEGG result object (Part 6, GSEA on KEGG)
#
# If objects are NOT in memory, uncomment the loading section below.

# ---- Setup ----

# install missing packages
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
if (!requireNamespace("ggridges", quietly = TRUE)) install.packages("ggridges")

# load packages
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(ggridges)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

# ---- OPTIONAL: load objects from RDS if not in memory ----
# vsd <- readRDS("vsd.rds")
# dds <- readRDS("dds.rds")
# res_df <- readRDS("res_df.rds")
# sig_table <- readRDS("sig_table.rds")
# dpsi <- read.csv("splicing_posterior_all.csv", stringsAsFactors = FALSE)
# go_results <- readRDS("go_results.rds")
# gsea_go <- readRDS("gsea_go.rds")
# gsea_kegg <- readRDS("gsea_kegg.rds")

# ---- Theme and color palette ----

# Okabe-Ito colorblind-friendly palette
COL_CONTROL  <- "#0072B2"   # blue
COL_INFECTED <- "#D55E00"   # vermillion
COL_UP       <- "#D55E00"   # vermillion (upregulated)
COL_DOWN     <- "#0072B2"   # blue (downregulated)
COL_NS       <- "grey80"    # not significant
COL_GREY     <- "grey50"

# condition colors for heatmaps
CONDITION_COLORS <- c(control = COL_CONTROL, infected = COL_INFECTED)

# publication theme
pub_theme <- theme_bw() +
  theme(
    text = element_text(family = "Liberation Sans"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

# helper: save ggplot as PNG
save_png <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  ggsave(filename, plot, width = width, height = height, dpi = dpi,
         units = "in", bg = "white")
  cat("Saved:", filename, "\n")
}

# ---- Figure 1: QC & Sample Summary ----

cat("\n=== Figure 1: QC & Sample Summary ===\n")

# 1A: PCA plot
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pca_pct  <- attr(pca_data, "percentVar") * 100

fig1a <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4, alpha = 0.9) +
  stat_ellipse(level = 0.95, type = "norm", linewidth = 0.8) +
  scale_color_manual(
    values = c(control = COL_CONTROL, infected = COL_INFECTED),
    name = "Condition"
  ) +
  labs(
    title = "PCA: Control vs SARS-CoV-2 Infected",
    x = paste0("PC1 (", round(pca_pct[1], 1), "% variance)"),
    y = paste0("PC2 (", round(pca_pct[2], 1), "% variance)")
  ) +
  pub_theme

save_png(fig1a, "fig1a_pca.png", width = 7, height = 5)

# 1B: Library size barplot
lib_sizes <- data.frame(
  sample = colnames(counts(dds)),
  total_counts = colSums(counts(dds)),
  condition = colData(dds)$condition
)
lib_sizes$sample <- factor(lib_sizes$sample, levels = lib_sizes$sample)

fig1b <- ggplot(lib_sizes, aes(sample, total_counts / 1e6, fill = condition)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(control = COL_CONTROL, infected = COL_INFECTED)) +
  labs(
    title = "Library Sizes",
    x = "Sample",
    y = "Total Counts (millions)"
  ) +
  pub_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_png(fig1b, "fig1b_library_sizes.png", width = 7, height = 5)

# 1C: Sample distance heatmap
vsd_mat <- assay(vsd)
sample_dists <- dist(t(vsd_mat))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- colnames(vsd_mat)
colnames(sample_dist_matrix) <- colnames(vsd_mat)

# annotation for condition
annotation_col <- data.frame(
  condition = colData(dds)$condition,
  row.names = colnames(vsd_mat)
)
ann_colors <- list(condition = CONDITION_COLORS)

png("fig1c_sample_distance_heatmap.png", width = 8, height = 7, units = "in",
    res = 300, bg = "white")
pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
  fontsize = 10,
  main = "Sample Distance Heatmap"
)
dev.off()
cat("Saved: fig1c_sample_distance_heatmap.png\n")

# ---- Figure 2: Differential Expression ----

cat("\n=== Figure 2: Differential Expression ===\n")

# 2A: Volcano plot
# ensure significance column exists
res_df$significance <- "Not significant"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange >= 1]  <- "Upregulated"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange <= -1] <- "Downregulated"
res_df$significance <- factor(res_df$significance,
                              levels = c("Upregulated", "Downregulated", "Not significant"))

# select top genes to label (by padj among significant)
label_genes <- res_df[res_df$significance != "Not significant", ]
label_genes <- label_genes[order(label_genes$padj), ]
label_genes <- head(label_genes, 10)

fig2a <- ggplot(res_df, aes(log2FoldChange, -log10(padj), color = significance)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_point(
    data = res_df[res_df$significance != "Not significant", ],
    size = 2, alpha = 0.9
  ) +
  scale_color_manual(
    values = c(Upregulated = COL_UP, Downregulated = COL_DOWN,
               `Not significant` = COL_NS),
    name = "Significance"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_text_repel(
    data = label_genes,
    aes(label = gene_symbol),
    size = 3, max.overlaps = 20,
    box.padding = 0.4, segment.color = "grey50"
  ) +
  labs(
    title = "Differential Gene Expression: Infected vs Control",
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~(adjusted~p-value))
  ) +
  pub_theme +
  coord_cartesian(ylim = c(0, max(-log10(res_df$padj), na.rm = TRUE) * 1.1))

save_png(fig2a, "fig2a_volcano.png", width = 8, height = 6)

# 2B: Heatmap of significant genes
# get vsd matrix for significant genes
sig_gene_ids <- sig_table$gene_id
heatmap_mat <- vsd_mat[rownames(vsd_mat) %in% sig_gene_ids, , drop = FALSE]

# use gene symbols as row names
rownames(heatmap_mat) <- sapply(rownames(heatmap_mat), function(id) {
  sym <- sig_table$gene_symbol[sig_table$gene_id == id]
  ifelse(length(sym) > 0 && !is.na(sym), sym, id)
})

# z-score rows
heatmap_mat_z <- t(scale(t(heatmap_mat)))

# column annotation
annotation_col_heat <- data.frame(
  condition = colData(dds)$condition,
  row.names = colnames(heatmap_mat_z)
)

png("fig2b_de_heatmap.png", width = 7, height = 8, units = "in",
    res = 300, bg = "white")
pheatmap(
  heatmap_mat_z,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  annotation_col = annotation_col_heat,
  annotation_colors = ann_colors,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(255),
  fontsize = 9,
  fontsize_row = 8,
  main = "Significant DE Genes (padj < 0.05, |log2FC| >= 1)",
  scale = "none"
)
dev.off()
cat("Saved: fig2b_de_heatmap.png\n")

# ---- Figure 3: Alternative Splicing (Null Result) ----

cat("\n=== Figure 3: Alternative Splicing (Null Result) ===\n")

# 3A: Splicing volcano — E(dPSI) vs P(|dPSI| > 0.20)
fig3a <- ggplot(dpsi, aes(dpsi_mean, probability_changing)) +
  geom_point(size = 0.8, alpha = 0.3, color = COL_GREY) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = COL_INFECTED,
             linewidth = 0.7) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "grey50",
             linewidth = 0.5) +
  annotate("text", x = 0, y = 0.97, label = "P = 0.95 threshold",
           color = COL_INFECTED, size = 3, vjust = -0.5) +
  annotate("text", x = 0, y = 0.45,
           label = "No events reach significance\n(max P = 0.532)",
           color = COL_GREY, size = 3.5, hjust = 0.5) +
  labs(
    title = "Differential Splicing: Infected vs Control",
    subtitle = "MAJIQ posterior — 3 biological replicates per group",
    x = expression(E(Delta*PSI)),
    y = expression(P(bgroup("|", Delta*PSI, "|") > 0.20))
  ) +
  pub_theme

save_png(fig3a, "fig3a_splicing_volcano.png", width = 8, height = 6)

# 3B: Histogram of probability_changing
fig3b <- ggplot(dpsi, aes(probability_changing)) +
  geom_histogram(bins = 100, fill = COL_CONTROL, alpha = 0.7,
                 color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = COL_INFECTED,
             linewidth = 0.7) +
  annotate("text", x = 0.95, y = Inf, label = "P = 0.95",
           color = COL_INFECTED, size = 3, vjust = 2, hjust = -0.1) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    title = "Distribution of Posterior Probabilities",
    subtitle = paste0("59,002 events — max P(|dPSI| > 0.20) = ",
                      signif(max(dpsi$probability_changing, na.rm = TRUE), 3)),
    x = expression(P(bgroup("|", Delta*PSI, "|") > 0.20)),
    y = "Number of Events"
  ) +
  pub_theme

save_png(fig3b, "fig3b_probability_histogram.png", width = 7, height = 5)

# 3C: Top 10 candidate events PSI comparison (exploratory)
top_events <- dpsi[order(-dpsi$probability_changing), ]
top_events <- head(top_events, 10)

# map gene symbols
top_events$ensembl_id <- sub("\\..*$", "", top_events$gene_id)
top_events$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = top_events$ensembl_id,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

# build data frame for plotting — make.unique prevents duplicate factor levels
event_labels <- paste0(
  ifelse(is.na(top_events$gene_symbol), top_events$gene_id, top_events$gene_symbol),
  " [", top_events$event_type, "] (P=", signif(top_events$probability_changing, 2), ")"
)
event_labels <- make.unique(event_labels)

psi_compare <- data.frame(
  event = factor(rep(event_labels, each = 2),
                 levels = event_labels),
  group = rep(c("control", "infected"), 10),
  psi = c(rbind(top_events$control_raw_psi_mean,
                top_events$infected_raw_psi_mean))
)

fig3c <- ggplot(psi_compare, aes(group, psi, fill = group)) +
  geom_col(width = 0.6, position = position_dodge(0.7)) +
  facet_wrap(~ event, scales = "free_y", ncol = 2) +
  scale_fill_manual(
    values = c(control = COL_CONTROL, infected = COL_INFECTED),
    name = "Group"
  ) +
  labs(
    title = "Top 10 Candidate Splicing Events (Exploratory)",
    subtitle = "Not statistically significant at n=3 — shown for hypothesis generation",
    x = "Group",
    y = "Mean PSI"
  ) +
  pub_theme +
  theme(
    axis.text.x = element_text(size = 8),
    strip.text = element_text(size = 8),
    panel.spacing = unit(0.8, "lines")
  )

save_png(fig3c, "fig3c_candidate_events.png", width = 9, height = 10)

# ---- Figure 4: GO Enrichment ----

cat("\n=== Figure 4: GO Enrichment ===\n")

# check if there are enriched terms
n_go <- nrow(as.data.frame(go_results))

if (n_go > 0) {
  
  # 4A: Dotplot
  png("fig4a_go_dotplot.png", width = 10, height = 7, units = "in",
      res = 300, bg = "white")
  print(dotplot(go_results, showCategory = min(15, n_go)) +
          ggtitle("GO Biological Process Enrichment") +
          theme(text = element_text(family = "Liberation Sans"),
                plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
  dev.off()
  cat("Saved: fig4a_go_dotplot.png\n")
  
  # 4B: Barplot
  png("fig4b_go_barplot.png", width = 10, height = 7, units = "in",
      res = 300, bg = "white")
  print(barplot(go_results, showCategory = min(15, n_go)) +
          ggtitle("GO Biological Process Enrichment") +
          theme(text = element_text(family = "Liberation Sans"),
                plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
  dev.off()
  cat("Saved: fig4b_go_barplot.png\n")
  
  # 4C: Cnetplot (if >= 3 enriched terms)
  if (n_go >= 3) {
    png("fig4c_go_cnetplot.png", width = 10, height = 8, units = "in",
        res = 300, bg = "white")
    print(cnetplot(go_results, showCategory = min(10, n_go)) +
            ggtitle("Gene-Concept Network") +
            theme(text = element_text(family = "Liberation Sans"),
                  plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
    dev.off()
    cat("Saved: fig4c_go_cnetplot.png\n")
  } else {
    cat("Skipping cnetplot (fewer than 3 enriched terms)\n")
  }
  
} else {
  cat("No enriched GO terms found. Skipping GO figures.\n")
  cat("Consider relaxing pvalueCutoff in Part 5 if appropriate.\n")
}

# ---- Figure 5: GSEA (Gene Set Enrichment Analysis) ----

cat("\n=== Figure 5: GSEA ===\n")

# 5A: GSEA GO dotplot — shows NES (not gene ratio like ORA dotplot)
n_gsea_go <- nrow(as.data.frame(gsea_go))

if (n_gsea_go > 0) {
  
  # tryCatch: split="NES" + facet_grid may not exist in older enrichplot
  p5a <- tryCatch({
    dotplot(gsea_go, showCategory = min(20, n_gsea_go), split = "NES") +
      facet_grid(~ NES > 0, scales = "free", space = "free")
  }, error = function(e) {
    cat("  (split='NES' not supported in this enrichplot version — using plain dotplot)\n")
    dotplot(gsea_go, showCategory = min(20, n_gsea_go))
  })
  p5a <- p5a +
    ggtitle("GSEA: GO Biological Process") +
    theme(text = element_text(family = "Liberation Sans"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  save_png(p5a, "fig5a_gsea_go_dotplot.png", width = 10, height = 8)
  
  # 5B: GSEA ridgeplot — density of fold changes within each gene set
  p5b <- ridgeplot(gsea_go, showCategory = min(15, n_gsea_go)) +
    ggtitle("GSEA Ridgeplot: GO BP (top gene sets)") +
    theme(text = element_text(family = "Liberation Sans"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  save_png(p5b, "fig5b_gsea_ridgeplot.png", width = 10, height = 8)
  
  # 5C: Running enrichment score for the top pathway
  # pick the most significant by adjusted p-value
  top_idx <- which.min(gsea_go@result$p.adjust)
  top_id <- gsea_go@result$ID[top_idx]
  top_desc <- gsea_go@result$Description[top_idx]
  
  p5c <- tryCatch({
    gseaplot2(gsea_go, geneSetID = top_id, pvalue_table = FALSE) +
      ggtitle(paste0("GSEA Running Score: ", top_desc)) +
      theme(text = element_text(family = "Liberation Sans"),
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  }, error = function(e) {
    cat("  (ggtitle not supported on gseaplot2 output — rendering without title)\n")
    gseaplot2(gsea_go, geneSetID = top_id, pvalue_table = FALSE)
  })
  png("fig5c_gsea_running_score.png", width = 10, height = 7, units = "in",
      res = 300, bg = "white")
  print(p5c)
  dev.off()
  cat("Saved: fig5c_gsea_running_score.png\n")
  
} else {
  cat("No significant GSEA GO BP gene sets (padj < 0.05).\n")
  cat("Consider relaxing pvalueCutoff in Part 6 if appropriate.\n")
}

# 5D: GSEA KEGG dotplot (if KEGG has results)
n_gsea_kegg <- nrow(as.data.frame(gsea_kegg))

if (n_gsea_kegg > 0) {
  
  p5d <- tryCatch({
    dotplot(gsea_kegg, showCategory = min(20, n_gsea_kegg), split = "NES") +
      facet_grid(~ NES > 0, scales = "free", space = "free")
  }, error = function(e) {
    cat("  (split='NES' not supported — using plain dotplot)\n")
    dotplot(gsea_kegg, showCategory = min(20, n_gsea_kegg))
  })
  p5d <- p5d +
    ggtitle("GSEA: KEGG Pathways") +
    theme(text = element_text(family = "Liberation Sans"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  save_png(p5d, "fig5d_gsea_kegg_dotplot.png", width = 10, height = 8)
  
} else {
  cat("No significant GSEA KEGG pathways (padj < 0.05).\n")
}

# ---- Summary ----

cat("\n=== Visualization complete ===\n")
cat("Figures saved to current working directory:\n")
cat("  fig1a_pca.png              - PCA plot\n")
cat("  fig1b_library_sizes.png    - Library size barplot\n")
cat("  fig1c_sample_distance_heatmap.png - Sample distance heatmap\n")
cat("  fig2a_volcano.png          - DE volcano plot\n")
cat("  fig2b_de_heatmap.png       - Significant DE genes heatmap\n")
cat("  fig3a_splicing_volcano.png - Splicing volcano (null result)\n")
cat("  fig3b_probability_histogram.png - Posterior probability histogram\n")
cat("  fig3c_candidate_events.png - Top 10 exploratory splicing events\n")
cat("  fig4a_go_dotplot.png       - GO ORA dotplot\n")
cat("  fig4b_go_barplot.png       - GO ORA barplot\n")
cat("  fig4c_go_cnetplot.png      - GO ORA gene-concept network\n")
cat("  fig5a_gsea_go_dotplot.png  - GSEA GO dotplot (NES)\n")
cat("  fig5b_gsea_ridgeplot.png   - GSEA ridgeplot (fold change density)\n")
cat("  fig5c_gsea_running_score.png - GSEA running score (top pathway)\n")
cat("  fig5d_gsea_kegg_dotplot.png  - GSEA KEGG dotplot (if enriched)\n")

