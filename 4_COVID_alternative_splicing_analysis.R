# SARS-CoV-2 Alternative Splicing Analysis

# IMPORTANT LIMITATIONS:
#   - PSI values are bounded [0,1] and often bimodal; the t-test
#     assumes normality.  With 11 mock + 12 infected samples the
#     CLT helps, but bimodal PSI distributions may still violate
#     assumptions.  Check residuals if interpreting borderline
#     calls.
#   - Many events have missing (NA) PSI values due to insufficient
#     junction coverage.  The MIN_N_PER_GROUP filter controls this
#     but reduces the testable universe.  Events with low n per
#     group have less reliable tests even after filtering.

SAMPLE_GROUPS <- c(
  SRR11412240 = "mock",
  SRR11412241 = "mock",
  SRR11412242 = "mock",
  SRR11412243 = "mock",
  SRR11412244 = "mock",
  SRR11412245 = "mock",
  SRR11412246 = "mock",
  SRR11412247 = "mock",
  SRR11412248 = "mock",
  SRR11412249 = "mock",
  SRR11412250 = "mock",
  SRR11412251 = "infected",
  SRR11412252 = "infected",
  SRR11412253 = "infected",
  SRR11412254 = "infected",
  SRR11412255 = "infected",
  SRR11412256 = "infected",
  SRR11412257 = "infected",
  SRR11412258 = "infected",
  SRR11412259 = "infected",
  SRR11412260 = "infected",
  SRR11412261 = "infected",
  SRR11412262 = "infected"
)

# Thresholds
MIN_N_PER_GROUP    <- 3       # minimum non-NA samples per group
DELTA_PSI_THRESHOLD <- 0.2    # minimum |delta_psi| for significance
FDR_THRESHOLD       <- 0.05   # BH-adjusted p-value cutoff

# Event type decoding
# Adjust to match your splicing tool's encoding.
EVENT_TYPE_MAP <- c(
  "s" = "exon_skipping",
  "t" = "alternative_terminal"
)

# Output directory
OUT_DIR <- "."

# 1. Validate inputs


cat("=== Input data dimensions ===\n")
cat("psi:           ", dim(psi), "\n")
cat("annotation:    ", dim(annotation), "\n")
cat("psi_results:   ", dim(psi_results), "\n\n")

# Identify sample columns (SRR accessions)
sample_cols <- names(SAMPLE_GROUPS)

# Check all sample columns exist in psi_results
missing_samples <- setdiff(sample_cols, names(psi_results))
if (length(missing_samples) > 0) {
  stop("These sample IDs in SAMPLE_GROUPS are not in psi_results: ",
       paste(missing_samples, collapse = ", "))
}

# Check ec_idx exists for the join
stopifnot("ec_idx" %in% names(psi_results))
stopifnot("ec_idx" %in% names(psi))

# annotation doesn't carry ec_idx -- it shares row order with psi.
# Bridge the key by copying ec_idx from psi (same 127191 rows).
if (!"ec_idx" %in% names(annotation)) {
  stopifnot(nrow(annotation) == nrow(psi))
  annotation$ec_idx <- psi$ec_idx
  cat("Added ec_idx to annotation from psi (positional match).\n")
}

# Report group sizes
n_mock     <- sum(SAMPLE_GROUPS == "mock")
n_infected <- sum(SAMPLE_GROUPS == "infected")
cat("Group sizes: mock =", n_mock, ", infected =", n_infected, "\n\n")

# Diagnostic: non-NA PSI count per sample (helps spot coverage issues)
cat("=== Non-NA PSI counts per sample ===\n")
non_na_counts <- colSums(!is.na(psi_results[, sample_cols]))
print(round(non_na_counts / nrow(psi_results) * 100, 1))
cat("(values are % of events with estimable PSI per sample)\n\n")

# 2. Prepare data for testing

mock_cols     <- sample_cols[SAMPLE_GROUPS == "mock"]
infected_cols <- sample_cols[SAMPLE_GROUPS == "infected"]

# Work on a copy so we don't modify the original
test_data <- psi_results

# Recompute mock_n and infected_n from the data to verify
# against the existing columns (validates group assignment).
test_data$mock_n_calc     <- rowSums(!is.na(test_data[, mock_cols]))
test_data$infected_n_calc <- rowSums(!is.na(test_data[, infected_cols]))

# computed n should match existing n columns

n_mismatch_mock <- sum(test_data$mock_n_calc != test_data$mock_n, na.rm = TRUE)
n_mismatch_inf  <- sum(test_data$infected_n_calc != test_data$infected_n, na.rm = TRUE)

if (n_mismatch_mock > 0 || n_mismatch_inf > 0) {
  warning("Computed n differs from existing n columns for ",
          n_mismatch_mock, " mock and ", n_mismatch_inf,
          " infected entries. Check your SAMPLE_GROUPS assignment.")
} else {
  cat("Group assignment verified: computed n matches existing n columns.\n\n")
}

# 3. Filter to testable events

# Keep only events with enough non-NA samples in both groups.
# This defines the universe of hypotheses we will test.
testable <- test_data[
  test_data$mock_n_calc >= MIN_N_PER_GROUP &
    test_data$infected_n_calc >= MIN_N_PER_GROUP,
]

cat("=== Filtering to testable events ===\n")
cat("Total events in psi_results:  ", nrow(test_data), "\n")
cat("Testable (n >= ", MIN_N_PER_GROUP, " per group): ", nrow(testable), "\n\n", sep = "")

if (nrow(testable) == 0) {
  stop("No events meet the minimum-n criterion. Lower MIN_N_PER_GROUP.")
}

# 4. Statistical testing (Welch's t-test + BH correction)

cat("Running Welch's t-test for ", nrow(testable), " events...\n", sep = "")

# run t-test on two vectors, handling edge cases
run_welch_ttest <- function(mock_vals, inf_vals) {
  m <- mock_vals[!is.na(mock_vals)]
  i <- inf_vals[!is.na(inf_vals)]
  
  # Need at least 2 per group for variance estimation
  if (length(m) < 2 || length(i) < 2) return(NA_real_)
  
  # If both groups have zero variance and identical means, p is undefined
  if (var(m) == 0 && var(i) == 0 && mean(m) == mean(i)) return(NA_real_)
  
  tryCatch(
    t.test(m, i, var.equal = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

# Extract sample matrices for speed
mock_mat     <- as.matrix(testable[, mock_cols])
infected_mat <- as.matrix(testable[, infected_cols])

# Run row-wise t-test
p_values <- vapply(
  seq_len(nrow(testable)),
  function(j) run_welch_ttest(mock_mat[j, ], infected_mat[j, ]),
  numeric(1)
)

testable$p_value <- p_values

# BH FDR correction over all tested events
testable$adj_p_value <- p.adjust(testable$p_value, method = "BH")

n_na_p <- sum(is.na(testable$p_value))
cat("Completed. ", n_na_p, " events had insufficient data for a test.\n\n", sep = "")

# 5. Define significant events

# Apply both criteria: FDR threshold AND effect-size threshold
significant <- testable[
  !is.na(testable$adj_p_value) &
    testable$adj_p_value <= FDR_THRESHOLD &
    abs(testable$delta_psi) >= DELTA_PSI_THRESHOLD,
]

cat("=== Significance calling ===\n")
cat("Criteria: FDR <= ", FDR_THRESHOLD, " and |delta_psi| >= ", DELTA_PSI_THRESHOLD, "\n", sep = "")
cat("Significant events: ", nrow(significant), "\n\n", sep = "")

if (nrow(significant) == 0) {
  cat("No events reached significance.\n")
  cat("Consider relaxing FDR_THRESHOLD or DELTA_PSI_THRESHOLD.\n")
  cat("Diagnostic: distribution of adj_p_value among testable events:\n")
  print(quantile(testable$adj_p_value, na.rm = TRUE, probs = c(0, .01, .05, .1, .25, .5)))
  
  # Still save the full test results for inspection
  write.csv(
    testable[, c("ec_idx", "mock_n_calc", "infected_n_calc",
                 "mock_mean", "infected_mean", "delta_psi",
                 "p_value", "adj_p_value")],
    file.path(OUT_DIR, "splicing_test_results_all.csv"),
    row.names = FALSE
  )
  cat("\nFull test results saved to splicing_test_results_all.csv\n")
  stop("No significant events found.")
}

# 6. Join with annotation

significant_annotated <- merge(
  significant,
  annotation,
  by = "ec_idx",
  suffixes = c("", ".anno")
)

# Verify every significant event found a match
if (nrow(significant_annotated) != nrow(significant)) {
  n_missing <- nrow(significant) - nrow(significant_annotated)
  stop(n_missing, " significant ec_idx values could not be matched to annotation.")
}

# 7. Decode event_type and fix is_intron

significant_annotated$event_type_label <-
  EVENT_TYPE_MAP[significant_annotated$event_type]

unknown_types <- setdiff(
  unique(significant_annotated$event_type),
  names(EVENT_TYPE_MAP)
)
if (length(unknown_types) > 0) {
  warning("Unrecognized event_type codes: ",
          paste(unknown_types, collapse = ", "),
          " -- update EVENT_TYPE_MAP.")
}

if ("is_intron" %in% names(significant_annotated) &&
    is.character(significant_annotated$is_intron)) {
  significant_annotated$is_intron <-
    as.logical(significant_annotated$is_intron)
}

# 8. Build the full summary table

desired_cols <- c(
  "ec_idx",
  "gene_name",
  "gene_id",
  "event_type",
  "event_type_label",
  "seqid",
  "strand",
  "ref_exon_start",
  "ref_exon_end",
  "start",
  "end",
  "other_exon_start",
  "other_exon_end",
  "is_intron",
  "delta_psi",
  "p_value",
  "adj_p_value",
  "mock_n_calc",
  "infected_n_calc",
  "mock_mean",
  "infected_mean"
)

available_cols <- intersect(desired_cols, names(significant_annotated))
missing_cols   <- setdiff(desired_cols, names(significant_annotated))

if (length(missing_cols) > 0) {
  warning("Columns not found, omitted: ",
          paste(missing_cols, collapse = ", "))
}

significant_summary <- significant_annotated[, available_cols, drop = FALSE]

# Sort by magnitude of splicing change (descending)
significant_summary <- significant_summary[
  order(-abs(significant_summary$delta_psi)),
]
rownames(significant_summary) <- NULL

# 9. Collapse reciprocal event pairs

pair_key <- paste(
  significant_summary$gene_name,
  significant_summary$event_type,
  significant_summary$ref_exon_start,
  significant_summary$ref_exon_end,
  sep = "_"
)

significant_summary$pair_key <- pair_key
significant_summary$n_paired  <- as.integer(ave(pair_key, pair_key, FUN = length))

# Split by pair key, pick representative, recombine
collapsed_list <- lapply(
  split(significant_summary, pair_key),
  function(df) {
    if (nrow(df) == 1) return(df)
    # Prefer positive delta_psi; tie-break by larger |delta_psi|
    pos <- df[df$delta_psi > 0, , drop = FALSE]
    if (nrow(pos) >= 1) {
      df <- pos[which.max(abs(pos$delta_psi)), , drop = FALSE]
    } else {
      df <- df[which.max(abs(df$delta_psi)), , drop = FALSE]
    }
    df
  }
)

collapsed_summary <- do.call(rbind, collapsed_list)
collapsed_summary$pair_key <- NULL

# Re-sort
collapsed_summary <- collapsed_summary[
  order(-abs(collapsed_summary$delta_psi)),
]
rownames(collapsed_summary) <- NULL

# 10. Display results

cat("\n=== Full significant events (junction-level) ===\n")
print(significant_summary)

cat("\n=== Collapsed events (one row per splice decision) ===\n")
print(collapsed_summary)

# 11. Summary statistics

cat("\n=== Summary statistics ===\n")
cat("Testable events:              ", nrow(testable), "\n")
cat("Junction-level significant:   ", nrow(significant_summary), "\n")
cat("Collapsed splice decisions:   ", nrow(collapsed_summary), "\n")
cat("Unique genes:                 ",
    length(unique(significant_summary$gene_name)), "\n")

cat("\nEvents by type (junction-level):\n")
print(table(significant_summary$event_type_label))

cat("\nEvents by type (collapsed):\n")
print(table(collapsed_summary$event_type_label))

cat("\nGenes with multiple splice decisions:\n")
gene_counts <- table(collapsed_summary$gene_name)
print(gene_counts[gene_counts > 1])

cat("\nDelta PSI range:\n")
cat("  min:", round(min(significant_summary$delta_psi), 3), "\n")
cat("  max:", round(max(significant_summary$delta_psi), 3), "\n")

cat("\nP-value range (significant events):\n")
cat("  raw p min:  ", signif(min(significant_summary$p_value), 3), "\n")
cat("  adj p min:  ", signif(min(significant_summary$adj_p_value), 3), "\n")
cat("  adj p max:  ", signif(max(significant_summary$adj_p_value), 3), "\n")

# 12. Main results table (display copy with rounding)

main_cols <- intersect(
  c("gene_name", "event_type_label",
    "mock_mean", "infected_mean", "delta_psi",
    "p_value", "adj_p_value",
    "mock_n_calc", "infected_n_calc"),
  names(collapsed_summary)
)

main_results <- collapsed_summary[, main_cols, drop = FALSE]

# Round only this display copy
main_results$mock_mean      <- round(main_results$mock_mean, 3)
main_results$infected_mean  <- round(main_results$infected_mean, 3)
main_results$delta_psi      <- round(main_results$delta_psi, 3)
if ("p_value" %in% names(main_results)) {
  main_results$p_value <- signif(main_results$p_value, 3)
}
if ("adj_p_value" %in% names(main_results)) {
  main_results$adj_p_value <- signif(main_results$adj_p_value, 3)
}

cat("\n=== Main results (collapsed, rounded for display) ===\n")
print(main_results)

# 13. Save results

# Full test results (all testable events, for reference)
write.csv(
  testable[, c("ec_idx", "mock_n_calc", "infected_n_calc",
               "mock_mean", "infected_mean", "delta_psi",
               "p_value", "adj_p_value")],
  file.path(OUT_DIR, "splicing_test_results_all.csv"),
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

cat("\n=== Results saved to:", OUT_DIR, "===\n")
cat("  splicing_test_results_all.csv        (all testable events with stats)\n")
cat("  significant_splicing_summary.csv      (significant junction-level)\n")
cat("  significant_splicing_collapsed.csv    (one row per splice decision)\n")
cat("  significant_splicing_main_results.csv (display table)\n")
