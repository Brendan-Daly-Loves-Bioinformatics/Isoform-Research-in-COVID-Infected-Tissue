# Part 3: Differential Expression Analysis of A549 COVID-19 RNA Sequencing Data
# gene-level expression is compared between control and infected samples

# load Bioconductor package
library(DESeq2)

# check the current working directory
getwd()

# check that the featureCounts output file is available
file.exists("gene_counts.txt")

# load featureCounts output
counts_raw <- read.delim(
  "gene_counts.txt",
  comment.char = "#",
  check.names = FALSE
)

# extract count columns
count_matrix <- counts_raw[, 7:30]

# give samples readable names
colnames(count_matrix) <- paste0("SRR114122", 39:62)

# set gene IDs as row names
rownames(count_matrix) <- counts_raw$Geneid

# merge 4 runs per biological sample by summing counts
bio_samples <- rep(c("Mock_1", "Mock_2", "Mock_3",
                     "CoV2_1", "CoV2_2", "CoV2_3"), each = 4)
conditions  <- rep(c("control", "control", "control",
                     "infected", "infected", "infected"), each = 4)

# sum the 4 run columns within each biological sample
count_collapsed <- t(rowsum(t(count_matrix), group = bio_samples))
colnames(count_collapsed) <- unique(bio_samples)

# 6-sample metadata
sample_info <- data.frame(
  row.names = colnames(count_collapsed),
  condition = factor(conditions[!duplicated(bio_samples)],
                     levels = c("control", "infected"))
)

# verify metadata matches count matrix
stopifnot(identical(
  rownames(sample_info),
  colnames(count_collapsed)
))

# create the DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_collapsed,
  colData   = sample_info,
  design    = ~ condition
)

# filter out genes with very low counts
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# sample-level QC
# library sizes
barplot(colSums(counts(dds)),
        las = 2, cex.names = 0.6,
        main = "Library sizes",
        ylab = "Total counts")

# PCA on variance-stabilized counts
vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = "condition") +
  ggtitle("PCA: control vs infected")

# check the number of genes remaining
dim(dds)

# run DESeq2 differential expression analysis
dds <- DESeq(dds)

# view the experimental design
design(dds)

# extract results comparing infected samples to controls
res <- results(
  dds,
  contrast = c("condition", "infected", "control")
)

# apply LFC shrinkage for accurate fold-change estimates
if (!requireNamespace("ashr", quietly = TRUE)) {
  install.packages("ashr")
}

res_shr <- lfcShrink(
  dds,
  contrast = c("condition", "infected", "control"),
  res = res,
  type = "ashr"
)

# order shrunken results by adjusted p-value
res_shr <- res_shr[order(res_shr$padj), ]

# display the top differential expression results
head(res_shr)

# convert DESeq2 results to a data frame
res_df <- as.data.frame(res_shr)

# add Ensembl gene IDs as a column
res_df$gene_id <- rownames(res_df)

# remove genes without an adjusted p-value
res_df <- res_df[!is.na(res_df$padj), ]

# define significantly differentially expressed genes
sig_genes <- subset(
  res_df,
  padj < 0.05 & abs(log2FoldChange) >= 1
)

# count upregulated and downregulated genes
upregulated <- sum(
  sig_genes$log2FoldChange >= 1
)

downregulated <- sum(
  sig_genes$log2FoldChange <= -1
)

# display summary
nrow(sig_genes)
upregulated
downregulated

# MA plot of differential expression results
png("ma_plot.png", width = 1200, height = 1000, res = 150)
plotMA(
  res_shr,
  ylim = c(-6, 6),
  alpha = 0.05
)
dev.off()

# create a volcano plot of differential expression results

# label genes as significant or not significant
res_df$significance <- "Not significant"

res_df$significance[
  res_df$padj < 0.05 &
    res_df$log2FoldChange >= 1
] <- "Upregulated"

res_df$significance[
  res_df$padj < 0.05 &
    res_df$log2FoldChange <= -1
] <- "Downregulated"

# calculate -log10 adjusted p-value
res_df$neg_log10_padj <- -log10(res_df$padj)

# create volcano plot
png("volcano_plot.png", width = 1200, height = 1000, res = 150)
plot(
  res_df$log2FoldChange,
  res_df$neg_log10_padj,
  pch = 20,
  col = ifelse(res_df$significance == "Upregulated", "#E64B35",
               ifelse(res_df$significance == "Downregulated", "#4DBBD5", "grey70")),
  main = "Differential Gene Expression: Infected vs Control",
  xlab = "Log2 Fold Change",
  ylab = "-Log10 Adjusted P-value"
)

# add significance thresholds
abline(
  v = c(-1, 1),
  lty = 2
)

abline(
  h = -log10(0.05),
  lty = 2
)
dev.off()

# install human gene annotation package if needed
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}

# load the annotation package
library(org.Hs.eg.db)

# remove Ensembl version numbers
res_df$ensembl_id <- sub(
  "\\..*$",
  "",
  res_df$gene_id
)

# convert Ensembl IDs to gene symbols
res_df$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = res_df$ensembl_id,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

# add gene biotype so rRNA/snRNA/pseudogene artifacts can be excluded
res_df$gene_biotype <- mapIds(
  org.Hs.eg.db,
  keys = res_df$ensembl_id,
  keytype = "ENSEMBL",
  column = "GENETYPE",
  multiVals = "first"
)

# view the annotated results
head(res_df[, c(
  "gene_id", "gene_symbol", "gene_biotype",
  "baseMean", "log2FoldChange", "pvalue", "padj"
)])

# recreate significant gene list with gene annotations
sig_genes <- subset(
  res_df,
  padj < 0.05 & abs(log2FoldChange) >= 1 & gene_biotype == "protein-coding"
)

# order by adjusted p-value
sig_genes <- sig_genes[order(sig_genes$padj), ]

# view significant genes
head(
  sig_genes[, c(
    "gene_symbol", "baseMean",
    "log2FoldChange", "pvalue", "padj"
  )],
  20
)

# create a clean table of significant genes
sig_table <- sig_genes[, c(
  "gene_symbol", "gene_id", "gene_biotype",
  "baseMean", "log2FoldChange", "pvalue", "padj"
)]

# remove genes without a gene symbol
sig_table <- sig_table[!is.na(sig_table$gene_symbol), ]

# view the complete significant gene table
sig_table

# save the significant gene table as a CSV file
write.csv(
  sig_table,
  "significant_differentially_expressed_genes.csv",
  row.names = FALSE
)

# display the number of significant genes
cat("Total significant genes:", nrow(sig_table), "
")

# display the number of upregulated genes
cat(
  "Upregulated genes:",
  sum(sig_table$log2FoldChange >= 1),
  "
"
)

# display the number of downregulated genes
cat(
  "Downregulated genes:",
  sum(sig_table$log2FoldChange <= -1),
  "
"
)

# display the 10 most strongly upregulated genes
top_upregulated <- sig_table[
  sig_table$log2FoldChange >= 1,
]

top_upregulated <- top_upregulated[
  order(-top_upregulated$log2FoldChange),
]

head(
  top_upregulated[, c(
    "gene_symbol", "log2FoldChange", "padj"
  )],
  10
)

# display the most strongly downregulated genes
top_downregulated <- sig_table[
  sig_table$log2FoldChange <= -1,
]

top_downregulated <- top_downregulated[
  order(top_downregulated$log2FoldChange),
]

top_downregulated[, c(
  "gene_symbol", "log2FoldChange", "padj"
)]