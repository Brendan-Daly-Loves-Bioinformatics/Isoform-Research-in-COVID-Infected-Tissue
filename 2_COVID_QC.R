# Part 2: Quality Control of A549 COVID-19 RNA-seq Short Reads
# raw and trimmed FASTQ files are assessed before downstream alignment and analysis

# location of the RNA-seq data
data_dir <- "//wsl$/Ubuntu/home/brendan/data/RNAseq_data/GSE147507/A549_COVID"
fastq_dir <- file.path(data_dir, "FASTQ")

# directory for FastQC output
qc_dir <- file.path(fastq_dir, "fastqc_reports")

# create the QC output directory
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# identify raw FASTQ files
fastq_files <- list.files(
  fastq_dir,
  pattern = "\\.fastq$",
  full.names = TRUE
)

# extract sample IDs from FASTQ filenames
sample_ids <- sub(
  "\\.fastq$",
  "",
  basename(fastq_files)
)

length(sample_ids)

# confirm that all expected FASTQ files are present
length(sample_ids) == 24

# confirm that all FASTQ files have non-zero file size
file_sizes <- file.info(fastq_files, extra_cols = FALSE)$size

all(file_sizes > 0)

# location of previously trimmed FASTQ files
trimmed_dir <- file.path(fastq_dir, "fastp_trimmed")

# use the compressed trimmed FASTQ files for post-trimming QC
trimmed_files <- list.files(
  trimmed_dir,
  pattern = "_trimmed\\.fastq\\.gz$",
  full.names = TRUE
)

# confirm that one trimmed file is available for each sample
length(trimmed_files)

# extract sample IDs from the trimmed FASTQ files
trimmed_ids <- sub(
  "_trimmed\\.fastq\\.gz$",
  "",
  basename(trimmed_files)
)

# confirm that all original samples have a trimmed FASTQ file
all(sample_ids %in% trimmed_ids)

# directory for FastQC results from trimmed reads
trimmed_qc_dir <- file.path(trimmed_dir, "fastqc_reports")

dir.create(
  trimmed_qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# verify that FastQC generated one HTML and one ZIP report per sample
trimmed_qc_files <- list.files(
  trimmed_qc_dir,
  pattern = "_fastqc\\.(html|zip)$",
  full.names = TRUE
)

length(trimmed_qc_files)

# separate HTML and ZIP reports
trimmed_html <- list.files(
  trimmed_qc_dir,
  pattern = "_fastqc\\.html$",
  full.names = TRUE
)

trimmed_zip <- list.files(
  trimmed_qc_dir,
  pattern = "_fastqc\\.zip$",
  full.names = TRUE
)

# confirm that every trimmed sample has both report types
length(trimmed_html) == length(trimmed_ids)
length(trimmed_zip) == length(trimmed_ids)