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

# create sample metadata
sample_info <- data.frame(
  row.names = colnames(count_matrix),
  condition = factor(
    c(rep("control", 12), rep("infected", 12)),
    levels = c("control", "infected")
  )
)

# verify metadata matches count matrix
stopifnot(identical(
  rownames(sample_info),
  colnames(count_matrix)
))

# create the DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ condition
)

# filter out genes with very low counts
keep <- rowSums(counts(dds) >= 10) >= 3

dds <- dds[keep, ]

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

# order results by adjusted p-value
res <- res[order(res$padj), ]

# Display the top differential expression results
head(res)

# convert DESeq2 results to a data frame
res_df <- as.data.frame(res)

# add Ensembl gene IDs as a column
res_df$gene_id <- rownames(res_df)

# Remove genes without an adjusted p-value
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
plotMA(
  res,
  ylim = c(-6, 6),
  alpha = 0.05
)

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
plot(
  res_df$log2FoldChange,
  res_df$neg_log10_padj,
  pch = 20,
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

# view the annotated results
head(res_df[, c(
  "gene_id",
  "gene_symbol",
  "baseMean",
  "log2FoldChange",
  "pvalue",
  "padj"
)])

# recreate significant gene list with gene annotations
sig_genes <- subset(
  res_df,
  padj < 0.05 & abs(log2FoldChange) >= 1
)

# order by adjusted p-value
sig_genes <- sig_genes[order(sig_genes$padj), ]

# view significant genes
head(
  sig_genes[, c(
    "gene_symbol",
    "baseMean",
    "log2FoldChange",
    "pvalue",
    "padj"
  )],
  20
)

# create a clean table of significant genes
sig_table <- sig_genes[, c(
  "gene_symbol",
  "gene_id",
  "baseMean",
  "log2FoldChange",
  "pvalue",
  "padj"
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
cat("Total significant genes:", nrow(sig_table), "\n")

# display the number of upregulated genes
cat(
  "Upregulated genes:",
  sum(sig_table$log2FoldChange >= 1),
  "\n"
)

# display the number of downregulated genes
cat(
  "Downregulated genes:",
  sum(sig_table$log2FoldChange <= -1),
  "\n"
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
    "gene_symbol",
    "log2FoldChange",
    "padj"
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
  "gene_symbol",
  "log2FoldChange",
  "padj"
)]