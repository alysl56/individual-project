# DMSO RNA-seq Project

This repository contains scripts and instructions for processing RNA-seq data
from human cell lines (A549, U937, Calu-3, HepG2) under DMSO-treated and untreated conditions.
All computations were performed on the Ada HPC cluster (University of Nottingham) using Slurm.

---

## 1. Data Sources

- Raw RNA-seq data were downloaded from **NCBI GEO/SRA**.  
- Four human cell lines were analyzed: **A549, Calu-3, HepG2, U937**.  
- For each cell line, samples were collected under two conditions (**with DMSO / without DMSO**).  
- Each condition contains **10 biological replicates**, i.e. **20 samples per cell line**.  
- Download links were obtained via [SRA Explorer](https://sra-explorer.info/).  

Raw data links are provided in:  
[`projects/raw_data/data_links/`](projects/raw_data/data_links/)  
- Example: [A549_withDMSO_links.txt](projects/raw_data/data_links/A549_withDMSO_links.txt), [Calu3_withoutDMSO_links.txt](projects/raw_data/data_links/Calu3_withoutDMSO_links.txt)

Download scripts are provided in:  
[`projects/scripts/download_scripts/`](projects/scripts/download_scripts/)  
- Example: [dl_wget_A549_withDMSO.sh](projects/scripts/download_scripts/dl_wget_A549_withDMSO.sh), [dl_wget_U937_withoutDMSO.sh](projects/scripts/download_scripts/dl_wget_U937_withoutDMSO.sh)

Raw FASTQ files were downloaded using **wget** (GNU Wget 1.20+).


---

---

## 2. Reference Transcriptome

We used the **GENCODE v38** human transcriptome (GRCh38) as the reference.  

- Transcript FASTA: [gencode.v38.transcripts.fa.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz)  
- Annotation GTF: [gencode.v38.annotation.gtf.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz)  

Due to file size, these files and the Salmon index are **not included** in this repository.  
They can be downloaded directly from [GENCODE](https://www.gencodegenes.org/human/release_38.html).  

**Download example:**  
```bash
# Create reference folder
mkdir -p references
cd references

# Download GENCODE v38 transcriptome and annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz

# Unzip
gunzip gencode.v38.transcripts.fa.gz
gunzip gencode.v38.annotation.gtf.gz
```
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

QC and trimming scripts are provided in:  
`projects/scripts/QC_trim_scripts/`  

- Example: [qc_fastp_A549_withDMSO.sh](projects/scripts/QC_trim_scripts/qc_fastp_A549_withDMSO.sh), [qc_fastp_U937_withoutDMSO.sh](projects/scripts/QC_trim_scripts/qc_fastp_U937_withoutDMSO.sh)  

Raw QC reports (`*.html`, `*.json`) are stored in:  
`projects/qc_reports/`  


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
```

-gencode.v38.annotation.gtf: gene annotation file (used later for tx2gene mapping)

-gencode.v38.transcripts.fa: reference transcriptome

-gencode_v38_index: Salmon index used for quantification

- **Input**: trimmed FASTQ files (*_1.trimmed.fastq.gz, *_2.trimmed.fastq.gz)

- **Process**: salmon quant with --validateMappings enabled

- **Output**:
Quantification directory for each sample (<SampleID>_quant/)

Expression estimates in quant.sf

Scripts are named following the convention:
quant_salmon_<CellLine>.sh
Example: quant_salmon_U937.sh


**Output structure:**
projects/DMSO_<CellLine>_RNAseq/salmon_quant/with_DMSO/<SampleID>_quant/
projects/DMSO_<CellLine>_RNAseq/salmon_quant/without_DMSO/<SampleID>_quant/

Each `<SampleID>_quant/` directory contains:
- `quant.sf` (expression estimates)
- Auxiliary logs from Salmon
## 6. Sample Metadata (samplesheet.csv)

To facilitate downstream transcript-level summarization with **tximport**,  
a metadata table (`samplesheet.csv`) was generated for each cell line.

- **Input**: Quantification directories from Salmon (`<SampleID>_quant/`)
- **Process**: Generate `samplesheet.csv` listing each sample, its condition, and quantification path
- **Output**: One `samplesheet.csv` per cell line, stored in the project root directory

The CSV has the following format:

sample,condition,quant_dir
SRRxxxxxxx,with_DMSO,/path/to/salmon_quant/with_DMSO/SRRxxxxxxx_quant
SRRyyyyyyy,without_DMSO,/path/to/salmon_quant/without_DMSO/SRRyyyyyyy_quant

Scripts are named following the convention:

`make_samplesheet_<CellLine>.sbatch`

Examples:

- `make_samplesheet_U937.sbatch`
- `make_samplesheet_HepG2.sbatch`
- `make_samplesheet_A549.sbatch`
- `make_samplesheet_Calu3.sbatch`

Each script was executed with `sbatch`, and produced the corresponding `samplesheet.csv` in:

projects/DMSO_<CellLine>_RNAseq/samplesheet.csv

---
## 7. Import with tximport (R)

Transcript-level quantifications from Salmon were summarized to gene-level counts using **tximport v1.30+** in R.

### Reference mapping (tx2gene)

A transcript-to-gene mapping file (`tx2gene.csv`) was generated from the **GENCODE v38 annotation** (`gencode.v38.annotation.gtf`).  
This file maps transcript IDs to gene IDs and is required by `tximport`.

- **Input**: `gencode.v38.annotation.gtf`  
- **Process**: Extract transcript–gene pairs  
- **Output**: `tx2gene.csv` (stored in `/gpfs01/home/alysl56/references/`)

The script used to generate this file is included in the repository:

`projects/scripts/references_scripts/make_tx2gene_AWK.sbatch`

### tximport execution

R was run in an **interactive Slurm session** on the Ada HPC cluster:

```bash
srun -p defq -c 4 --mem=24G -t 02:00:00 --pty bash
module load R-uoneasy/4.4.1-gfbf-2023b-rstudio
R
Within R, tximport was executed with the generated samplesheet.csv and tx2gene.csv:
