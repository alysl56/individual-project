# DMSO RNA-seq Project

This repository contains scripts and instructions for processing RNA-seq data
from human cell lines (A549, U937, Calu-3, HepG2) under DMSO-treated and untreated conditions.
All computations were performed on the Ada HPC cluster (University of Nottingham) using Slurm.

---

## 1. Data Sources

- RNA-seq data were downloaded from **NCBI GEO/SRA**.
- For each cell line, 5 GEO series (GSE) were selected.
- In each GSE, 2 samples (GSM) were selected.
- Each condition (with DMSO / without DMSO) therefore contains 10 biological replicates, i.e. 20 samples per cell line.
- Download links were obtained via [SRA Explorer](https://sra-explorer.info/).

Scripts for downloading raw FASTQ files are provided in:
projects/scripts/download_fastq/

Raw FASTQ files were downloaded using **wget** (GNU Wget 1.20+).

---

## 2. Reference Transcriptome

We used the **GENCODE v38** human transcriptome (GRCh38) as the reference:

- Transcript FASTA: `gencode.v38.transcripts.fa`
- Annotation GTF: `gencode.v38.annotation.gtf`

Due to file size, these files and the Salmon index are **not included** in this repository.  
They can be downloaded from [GENCODE](https://www.gencodegenes.org/human/release_38.html).

---

## 3. Building the Salmon Index

The Salmon index was built using **Salmon v1.9.0**:

```bash
salmon index \
  -t gencode.v38.transcripts.fa \
  -i gencode_v38_index \
  -k 31
```

- `-t` = transcriptome FASTA  
- `-i` = output directory for index  
- `-k` = k-mer size (31 is default and recommended)  
- Auxiliary logs from Salmon  

---
## 4. Preprocessing (QC + Trimming)

Preprocessing of raw FASTQ files was performed using **fastp v0.22.0**:

- **Input**: raw paired-end FASTQ files  
- **Process**: quality control and trimming  
  - Adapter removal  
  - Filtering low-quality reads  
  - Trimming short reads  
- **Output**:  
  - Trimmed FASTQ files (`*_trimmed.fastq.gz`)  
  - QC reports (`*.html`, `*.json`)  

Scripts are named following the convention:
qc_fastp_<CellLine>_<Condition>.sh
Example: `qc_fastp_U937_withDMSO.sh`

---

## 5. Quantification (Salmon)

Transcript quantification was performed using **Salmon v1.9.0** with the **GENCODE v38** reference transcriptome.

### Reference preparation
Before running `salmon quant`, the reference transcriptome and annotation were downloaded and indexed:

```bash
# Download reference GTF and transcript FASTA
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz

# Unzip
gunzip gencode.v38.annotation.gtf.gz
gunzip gencode.v38.transcripts.fa.gz

# Build Salmon index
salmon index -t gencode.v38.transcripts.fa -i gencode_v38_index --gencode
gencode.v38.annotation.gtf: gene annotation file (used later for tx2gene mapping)

gencode.v38.transcripts.fa: reference transcriptome

gencode_v38_index: Salmon index used for quantification

Input: trimmed FASTQ files (*_1.trimmed.fastq.gz, *_2.trimmed.fastq.gz)

Process: salmon quant with --validateMappings enabled

Output:

Quantification directory for each sample (<SampleID>_quant/)

Expression estimates in quant.sf

Scripts are named following the convention:
quant_salmon_<CellLine>.sh
Example: quant_salmon_U937.sh

**Output structure:**
projects/DMSO_<CellLine>_RNAseq/salmon_quant/with_DMSO/<SampleID>quant/
projects/DMSO<CellLine>_RNAseq/salmon_quant/without_DMSO/<SampleID>_quant/

Each `<SampleID>_quant/` directory contains:
- `quant.sf` (expression estimates)
- Auxiliary logs from Salmon
