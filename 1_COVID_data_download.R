# Part 1: Download the A549 Series2 RNA-seq data from GSE147507
# 24 SRR runs = 6 biological samples (3 Mock + 3 SARS-CoV-2), 4 runs each

data_dir <- "//wsl$/Ubuntu/home/brendan/data/RNAseq_data/GSE147507/A549_COVID"
fastq_dir <- file.path(data_dir, "FASTQ")
dir.create(fastq_dir, recursive = TRUE, showWarnings = FALSE)

# 24 run accessions (4 runs per biological sample)
run_ids <- c(
  "SRR11412239", "SRR11412240", "SRR11412241", "SRR11412242",  # Mock_1
  "SRR11412243", "SRR11412244", "SRR11412245", "SRR11412246",  # Mock_2
  "SRR11412247", "SRR11412248", "SRR11412249", "SRR11412250",  # Mock_3
  "SRR11412251", "SRR11412252", "SRR11412253", "SRR11412254",  # SARS-CoV-2_1
  "SRR11412255", "SRR11412256", "SRR11412257", "SRR11412258",  # SARS-CoV-2_2
  "SRR11412259", "SRR11412260", "SRR11412261", "SRR11412262"   # SARS-CoV-2_3
)

# map runs to biological samples for downstream merging
sample_map <- data.frame(
  run = run_ids,
  sample = rep(c("Mock_1","Mock_2","Mock_3","CoV2_1","CoV2_2","CoV2_3"), each = 4),
  condition = rep(c("Mock","Mock","Mock","SARS_CoV_2","SARS_CoV_2","SARS_CoV_2"), each = 4),
  stringsAsFactors = FALSE
)

# uncomment to download
# for (run in run_ids) {
#   status <- system2("fasterq-dump", args = c(
#     run, "--outdir", fastq_dir, "--threads", "8", "--gzip", "--temp", fastq_dir
#   ))
#   if (status != 0) warning("fasterq-dump failed for ", run)
# }

# verify downloads
fastq_files <- file.path(fastq_dir, paste0(run_ids, ".fastq.gz"))
sum(file.exists(fastq_files))  # should be 24

# verify downloads — show which are missing
missing <- run_ids[!file.exists(fastq_files)]
if (length(missing) > 0) {
  cat("Missing", length(missing), "files:
")
  print(missing)
} else {
  cat("All 24 FASTQ files present.
")
}