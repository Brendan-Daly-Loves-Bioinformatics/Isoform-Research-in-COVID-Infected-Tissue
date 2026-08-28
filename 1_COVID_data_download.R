# Download the A549 COVID-19 RNA-seq data used for this research
# (Files are very large, use with caution)

# file directories
data_dir <- "//wsl$/Ubuntu/home/brendan/data/RNAseq_data/GSE147507/A549_COVID"
fastq_dir <- file.path(data_dir, "FASTQ")

# create new directories
dir.create(fastq_dir, recursive = TRUE, showWarnings = FALSE)

list.files(fastq_dir)

# define the list of 24 sample IDs used in analysis
sample_ids <- c(
  "SRR11412239", "SRR11412240", "SRR11412241", "SRR11412242",
  "SRR11412243", "SRR11412244", "SRR11412245", "SRR11412246",
  "SRR11412247", "SRR11412248", "SRR11412249", "SRR11412250",
  "SRR11412251", "SRR11412252", "SRR11412253", "SRR11412254",
  "SRR11412255", "SRR11412256", "SRR11412257", "SRR11412258",
  "SRR11412259", "SRR11412260", "SRR11412261", "SRR11412262"
)

# confirms the 24 samples are present for use in R
length(sample_ids)

# uncomment the chunk below to download FASTQ files on a new system.

# for (sample in sample_ids) {
#   system2(
#     "prefetch",
#     args = sample
#   )
#
#   system2(
#     "fasterq-dump",
#     args = c(
#       sample,
#       "--outdir", fastq_dir,
#       "--threads", "8"
#     )
#   )
# }

# verify that all FASTQ files were downloaded correctly for use
fastq_files <- file.path(
  fastq_dir,
  paste0(sample_ids, ".fastq")
)

file.exists(fastq_files)

sum(file.exists(fastq_files))