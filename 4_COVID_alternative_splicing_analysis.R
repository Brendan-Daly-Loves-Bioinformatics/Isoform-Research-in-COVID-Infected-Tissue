# Part 4: SARS-CoV-2 Alternative Splicing Analysis (Posterior-Based)
# Uses MAJIQ deltapsi posterior probability P(|dPSI| > 0.20)
# from 3 mock vs 3 infected biological replicates

# Thresholds
PROB_CHANGING_THRESHOLD <- 0.95   # P(|dPSI| > 0.20) > 0.95 (standard MAJIQ)
DPSI_MEAN_THRESHOLD <- 0.2        # |E(dPSI)| >= 0.2 (effect size)

# Event type decoding
EVENT_TYPE_MAP <- c(
  "s" = "exon_skipping",
  "t" = "alternative_terminal"
)

# Output directory
OUT_DIR <- "."

# 1. Load deltapsi output (skip JSON metadata header)

cat("=== Loading MAJIQ deltapsi output ===
")

deltapsi_file <- "//wsl$/Ubuntu/home/brendan/majiq_quant/deltapsi_bio3.tsv"

stopifnot(file.exists(deltapsi_file))

dpsi <- read.delim(
  deltapsi_file,
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Total events quantified:", nrow(dpsi), "
")
cat("Columns:", paste(names(dpsi), collapse = ", "), "

")

# 2. Diagnostic: distribution of probability_changing

cat("=== Diagnostic: probability_changing distribution ===
")
cat("Quantiles of P(|dPSI| > 0.20):
")
print(quantile(dpsi$probability_changing, na.rm = TRUE,
               probs = c(0, 0.5, 0.9, 0.95, 0.99, 0.999, 1)))

cat("
Events at various thresholds:
")
for (thresh in c(0.5, 0.8, 0.9, 0.95, 0.99)) {
  n <- sum(dpsi$probability_changing > thresh, na.rm = TRUE)
  cat("  P > ", thresh, ":", n, "events
")
}
cat("
")

# 3. Call significant events using posterior probability

cat("=== Significance calling ===
")
cat("Criteria: P(|dPSI| > 0.20) > ", PROB_CHANGING_THRESHOLD,
    " and |E(dPSI)| >= ", DPSI_MEAN_THRESHOLD, "
", sep = "")

significant <- dpsi[
  !is.na(dpsi$probability_changing) &
    dpsi$probability_changing > PROB_CHANGING_THRESHOLD &
    abs(dpsi$dpsi_mean) >= DPSI_MEAN_THRESHOLD,
]

cat("Significant events:", nrow(significant), "

")

if (nrow(significant) == 0) {
  cat("No events reached significance at P > 0.95.
")
  cat("Retrying with relaxed threshold P > 0.90 for exploratory purposes...

")
  
  significant <- dpsi[
    !is.na(dpsi$probability_changing) &
      dpsi$probability_changing > 0.90 &
      abs(dpsi$dpsi_mean) >= DPSI_MEAN_THRESHOLD,
  ]
  
  cat("Events at P > 0.90:", nrow(significant), "

")
  
  if (nrow(significant) == 0) {
    # Save full results for reference and stop
    write.csv(dpsi, file.path(OUT_DIR, "splicing_posterior_all.csv"),
              row.names = FALSE)
    cat("Full posterior results saved to splicing_posterior_all.csv
")
    cat("No significant splicing events found at any explored threshold.
")
    cat("This is an honest null result with n=3 biological replicates.
")
    stop("No significant events found.")
  } else {
    cat("NOTE: Using relaxed threshold P > 0.90 (exploratory, not confirmatory).

")
  }
}

# 4. Decode event types

significant$event_type_label <- EVENT_TYPE_MAP[significant$event_type]

unknown_types <- setdiff(
  unique(significant$event_type),
  names(EVENT_TYPE_MAP)
)
if (length(unknown_types) > 0) {
  warning("Unrecognized event_type codes: ",
          paste(unknown_types, collapse = ", "),
          " -- update EVENT_TYPE_MAP.")
}

# Fix is_intron if character
if ("is_intron" %in% names(significant) &&
    is.character(significant$is_intron)) {
  significant$is_intron <- as.logical(significant$is_intron)
}

# 5. Map Ensembl gene IDs to gene symbols

cat("=== Annotating gene symbols ===
")

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)

# gene_id has version numbers (e.g., ENSG00000293331.1); strip them
significant$ensembl_id <- sub("\\..*$", "", significant$gene_id)

significant$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = significant$ensembl_id,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

significant$gene_biotype <- mapIds(
  org.Hs.eg.db,
  keys = significant$ensembl_id,
  keytype = "ENSEMBL",
  column = "GENETYPE",
  multiVals = "first"
)

cat("Genes with mapped symbols:",
    sum(!is.na(significant$gene_symbol)), "of", nrow(significant), "

")

# 6. Build the full summary table

desired_cols <- c(
  "gene_name", "gene_id", "gene_symbol", "gene_biotype",
  "event_type", "event_type_label",
  "seqid", "strand",
  "ref_exon_start", "ref_exon_end",
  "start", "end",
  "other_exon_start", "other_exon_end",
  "is_intron", "is_denovo", "event_denovo",
  "dpsi_mean", "dpsi_std",
  "probability_changing", "probability_nonchanging",
  "control_raw_psi_mean", "infected_raw_psi_mean",
  "control_raw_coverage", "infected_raw_coverage"
)

available_cols <- intersect(desired_cols, names(significant))
missing_cols <- setdiff(desired_cols, names(significant))
if (length(missing_cols) > 0) {
  warning("Columns not found, omitted: ",
          paste(missing_cols, collapse = ", "))
}

significant_summary <- significant[, available_cols, drop = FALSE]

# Sort by probability_changing (descending), then |dpsi_mean| (descending)
significant_summary <- significant_summary[
  order(-significant_summary$probability_changing,
        -abs(significant_summary$dpsi_mean)),
]
rownames(significant_summary) <- NULL

# 7. Collapse reciprocal event pairs

pair_key <- paste(
  significant_summary$gene_name,
  significant_summary$event_type,
  significant_summary$ref_exon_start,
  significant_summary$ref_exon_end,
  sep = "_"
)

significant_summary$pair_key <- pair_key
significant_summary$n_paired <- as.integer(ave(pair_key, pair_key, FUN = length))

collapsed_list <- lapply(
  split(significant_summary, pair_key),
  function(df) {
    if (nrow(df) == 1) return(df)
    # Prefer positive dpsi_mean; tie-break by higher probability_changing
    pos <- df[df$dpsi_mean > 0, , drop = FALSE]
    if (nrow(pos) >= 1) {
      df <- pos[which.max(pos$probability_changing), , drop = FALSE]
    } else {
      df <- df[which.max(df$probability_changing), , drop = FALSE]
    }
    df
  }
)

collapsed_summary <- do.call(rbind, collapsed_list)
collapsed_summary$pair_key <- NULL

# Re-sort by probability_changing descending
collapsed_summary <- collapsed_summary[
  order(-collapsed_summary$probability_changing),
]
rownames(collapsed_summary) <- NULL

# 8. Display results

cat("
=== Full significant events (junction-level) ===
")
print(significant_summary)

cat("
=== Collapsed events (one row per splice decision) ===
")
print(collapsed_summary)

# 9. Summary statistics

cat("
=== Summary statistics ===
")
cat("Total quantified events:", nrow(dpsi), "
")
cat("Junction-level significant:", nrow(significant_summary), "
")
cat("Collapsed splice decisions:", nrow(collapsed_summary), "
")
cat("Unique genes:",
    length(unique(significant_summary$gene_symbol[!is.na(significant_summary$gene_symbol)])),
    "
")

cat("
Events by type (junction-level):
")
print(table(significant_summary$event_type_label))

cat("
Events by type (collapsed):
")
print(table(collapsed_summary$event_type_label))

cat("
Genes with multiple splice decisions:
")
gene_counts <- table(collapsed_summary$gene_symbol[!is.na(collapsed_summary$gene_symbol)])
print(gene_counts[gene_counts > 1])

cat("
E(dPSI) range:
")
cat("  min:", round(min(significant_summary$dpsi_mean), 3), "
")
cat("  max:", round(max(significant_summary$dpsi_mean), 3), "
")

cat("
Probability_changing range:
")
cat("  min:", signif(min(significant_summary$probability_changing), 3), "
")
cat("  max:", signif(max(significant_summary$probability_changing), 3), "
")

# 10. Main results table (display copy with rounding)

main_cols <- intersect(
  c("gene_symbol", "event_type_label",
    "control_raw_psi_mean", "infected_raw_psi_mean",
    "dpsi_mean", "dpsi_std",
    "probability_changing", "probability_nonchanging"),
  names(collapsed_summary)
)

main_results <- collapsed_summary[, main_cols, drop = FALSE]

main_results$dpsi_mean <- round(main_results$dpsi_mean, 3)
main_results$dpsi_std <- round(main_results$dpsi_std, 4)
main_results$control_raw_psi_mean <- round(main_results$control_raw_psi_mean, 3)
main_results$infected_raw_psi_mean <- round(main_results$infected_raw_psi_mean, 3)
main_results$probability_changing <- signif(main_results$probability_changing, 3)
main_results$probability_nonchanging <- signif(main_results$probability_nonchanging, 3)

cat("
=== Main results (collapsed, rounded for display) ===
")
print(main_results)

# 11. Save results

# Full posterior results (all events, for reference and visualization)
write.csv(
  dpsi,
  file.path(OUT_DIR, "splicing_posterior_all.csv"),
  row.names = FALSE
)

# Junction-level significant events
write.csv(
  significant_summary,
  file.path(OUT_DIR, "significant_splicing_summary.csv"),
  row.names = FALSE
)

# Collapsed (one row per splice decision)
write.csv(
  collapsed_summary,
  file.path(OUT_DIR, "significant_splicing_collapsed.csv"),
  row.names = FALSE
)

# Display table
write.csv(
  main_results,
  file.path(OUT_DIR, "significant_splicing_main_results.csv"),
  row.names = FALSE
)

cat("
=== Results saved to:", OUT_DIR, "===
")
cat("  splicing_posterior_all.csv (all events with posterior probabilities)
")
cat("  significant_splicing_summary.csv (significant junction-level)
")
cat("  significant_splicing_collapsed.csv (one row per splice decision)
")
cat("  significant_splicing_main_results.csv (display table)
")