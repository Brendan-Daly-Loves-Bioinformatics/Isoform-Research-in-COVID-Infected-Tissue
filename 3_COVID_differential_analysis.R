# Load featureCounts output
counts_raw <- read.delim(
  "gene_counts.txt",
  comment.char = "#",
  check.names = FALSE
)

# Extract count columns
count_matrix <- counts_raw[, 7:30]

# Give samples readable names
colnames(count_matrix) <- paste0("SRR114122", 39:62)

# Set gene IDs as row names
rownames(count_matrix) <- counts_raw$Geneid

# Create sample metadata
sample_info <- data.frame(
  row.names = colnames(count_matrix),
  condition = factor(
    c(rep("control", 12), rep("infected", 12)),
    levels = c("control", "infected")
  )
)

# Verify metadata matches count matrix
stopifnot(identical(
  rownames(sample_info),
  colnames(count_matrix)
))

