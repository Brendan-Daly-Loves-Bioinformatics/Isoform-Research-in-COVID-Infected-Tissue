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

# verify raw FastQC reports exist (one HTML + one ZIP per sample)
raw_html <- list.files(qc_dir, pattern = "_fastqc\\.html$", full.names = TRUE)
raw_zip  <- list.files(qc_dir, pattern = "_fastqc\\.zip$",  full.names = TRUE)

length(raw_html)
length(raw_zip)

# confirm that every raw sample has both report types
length(raw_html) == length(sample_ids)
length(raw_zip)  == length(sample_ids)

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

# confirm no extra or unexpected trimmed files
setequal(sample_ids, trimmed_ids)

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
length(trimmed_zip)  == length(trimmed_ids)

# aggregate raw and trimmed FastQC reports with MultiQC (run via WSL)
multiqc_dir <- file.path(data_dir, "multiqc")
dir.create(multiqc_dir, recursive = TRUE, showWarnings = FALSE)

# convert WSL Windows paths to native Linux paths for MultiQC
linux_qc_dir      <- sub("//wsl\\$/Ubuntu", "", qc_dir)
linux_trimmed_qc  <- sub("//wsl\\$/Ubuntu", "", trimmed_qc_dir)
linux_multiqc_dir <- sub("//wsl\\$/Ubuntu", "", multiqc_dir)

system(paste("wsl ~/.local/bin/multiqc", linux_qc_dir, linux_trimmed_qc, "-o", linux_multiqc_dir, "--dirs"))

list.files(multiqc_dir, pattern = "\\.html$")

# verify that the MultiQC report was created
list.files(multiqc_dir, pattern = "\\.html$")