# Alternative-Splicing-Research-in-COVID-Infected-Tissue
Research of SARS-CoV-2 infection of A549 ACE2-expressing cells in human lung tissue
I am using R Bioconductor and its vast tools to analyze short-read RNAseq data, and research further into how COVID-19 can effect alternative splicing in A549 ACE2-expressing cells. Changed alternative splicing can cause many issues in host cells, such as inactive variants of immune response genes, enhanced viral replication, altered host proteins, and abnormal protein isoforms that can cause cancer.

Transcriptomics Pipeline

1. Use SRA Toolkit to download SRA files and fasterq-dump to convert them to fastq files
2. Run FastQC for quality control reports on the data, analyze and fix any issues, and trim adapters off the reads using fastp
3. Use STAR (short read) or minimap2 (long read) to map RNAseq to a reference human genome, hg38. (may take hours/days)
4. Simultaneously as step 3, samtools should be used to sort the minimap2 .bam file outputs into BAM and BAI files.
5. Use STAR to make TAB files with the splice junction that minimap2 won't produce if needed.
6. I ran MAJIQ on the BAM files to make gene-specific splice graphs.
   
Data Download

Quality Control

Differential Expression Analysis

Alternative Splicing Analysis

Gene Ontology Enrichment

Gene Set Enrichment Analysis

Visualizations

Sources Cited
