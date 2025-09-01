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
[`raw_data/data_links/`](raw_data/data_links/)  
- Example: [A549_withDMSO_links.txt](raw_data/data_links/A549_withDMSO_links.txt), [Calu3_withoutDMSO_links.txt](raw_data/data_links/Calu3_withoutDMSO_links.txt)

Download scripts are provided in:  
[`scripts/download_scripts/`](scripts/download_scripts/)  
- Example: [dl_wget_A549_withDMSO.sh](scripts/download_scripts/dl_wget_A549_withDMSO.sh), [dl_wget_U937_withoutDMSO.sh](scripts/download_scripts/dl_wget_U937_withoutDMSO.sh)

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
  - Trimmed FASTQ files (`*_trimmed.fastq.gz`, not uploaded due to large size; can be regenerated using scripts)  
  - QC reports (`*.html`, `*.json`)  
QC and trimming scripts are provided in:  
[`scripts/QC_trim_scripts/`](scripts/QC_trim_scripts/)  

- Example: [qc_fastp_A549_withDMSO.sh](scripts/QC_trim_scripts/qc_fastp_A549_withDMSO.sh), [qc_fastp_U937_withoutDMSO.sh](scripts/QC_trim_scripts/qc_fastp_U937_withoutDMSO.sh)  

QC reports are organized by cell line and condition in:  
[`qc_reports/`](qc_reports/)

- Example: [A549_with_DMSO](qc_reports/A549_with_DMSO), [HepG2_without_DMSO](qc_reports/HepG2_without_DMSO)  

---

## 5. Quantification (Salmon)

Transcript quantification was performed using **Salmon v1.9.0** with the **GENCODE v38** reference transcriptome.
Transcript quantification was performed using **Salmon v1.9.0** with the **GENCODE v38** reference transcriptome.

- **Input**: trimmed FASTQ files (`*_trimmed.fastq.gz`)  
- **Process**: quantification with `salmon quant` (`--validateMappings` enabled)  
- **Output**:  
  - One `quant.sf` file per sample (not uploaded due to size; can be regenerated using scripts)  
  - Auxiliary Salmon logs  

Quantification scripts are provided in:  
[`scripts/salmon_quant_scripts/`](scripts/salmon_quant_scripts/)

- Example: [quant_salmon_A549.sh](scripts/salmon_quant_scripts/quant_salmon_A549.sh), [quant_salmon_U937.sh](scripts/salmon_quant_scripts/quant_salmon_U937.sh)  

Each sample produces a quantification directory `<SampleID>_quant/` containing:  
- `quant.sf` (transcript-level expression estimates)  
- Auxiliary logs from Salmon  

**Output structure:**
projects/DMSO_<CellLine>_RNAseq/salmon_quant/with_DMSO/<SampleID>_quant/
projects/DMSO_<CellLine>_RNAseq/salmon_quant/without_DMSO/<SampleID>_quant/

## 6. Sample Metadata (samplesheet.csv)

To facilitate downstream transcript-level summarization with **tximport**,  
a metadata table (`samplesheet.csv`) was generated for each cell line.

- **Input**: Quantification directories from Salmon (`<SampleID>_quant/`)  
- **Process**: Generate `samplesheet.csv` listing each sample, its condition, and quantification path  
- **Output**: One `samplesheet.csv` per cell line (uploaded to repository for reproducibility)  

The CSV has the following format:

sample,condition,quant_dir
SRRxxxxxxx,with_DMSO,/path/to/salmon_quant/with_DMSO/SRRxxxxxxx_quant
SRRyyyyyyy,without_DMSO,/path/to/salmon_quant/without_DMSO/SRRyyyyyyy_quant


Scripts for generating the metadata are stored in:  
[`scripts/make_samplesheet_scripts/`](scripts/make_samplesheet_scripts)  

- Example: [`make_samplesheet_A549.sbatch`](scripts/make_samplesheet_scripts/make_samplesheet_A549.sbatch),  
  [`make_samplesheet_U937.sbatch`](scripts/make_samplesheet_scripts/make_samplesheet_U937.sbatch)

The generated metadata files are stored in:  
[`samplesheet/`](samplesheet)  

- Example: [`A549_samplesheet.csv`](samplesheet/A549_samplesheet.csv),  
  [`U937_samplesheet.csv`](samplesheet/U937_samplesheet.csv)


---
## 7. Import with tximport (R) + Differential Expression (DESeq2)

Transcript-level quantifications from Salmon were summarized to gene-level counts using **tximport (v1.34.0;)**, followed by differential expression analysis with **DESeq2 (v1.46.0)**.

- **Environment**:  
  Analyses were performed on the **Ada HPC cluster** in two ways:  
  - Interactive R session (for testing/exploration)  
  - Batch mode (Slurm scripts for full analysis)

- **Required R packages**:  
  `DESeq2 (v1.46.0)`, `tximport (v1.34.0)`, `clusterProfiler (v4.14.6)`,  
  `org.Hs.eg.db (v3.20.0)`, `ggplot2 (v3.5.2)`, `pheatmap (v1.0.13)`, `VennDiagram (v1.7.3)`

- **Input**:  
  - Salmon quantification directories (`<SampleID>_quant/`)  
  - Transcript-to-gene mapping file: [`tx2gene.csv`](intermediate/tx2gene.csv)  

- **Process**:  
  1. Generate `tx2gene.csv` using  
     Scripts: [`make_tx2gene_AWK.sbatch`](scripts/references_scripts/make_tx2gene_AWK.sbatch)
  2. Run **tximport** for gene-level summarization  
     Scripts: [`run_tximport_*.R`](scripts/tximport_DESeq2_scripts/)  
  3. Run **DESeq2** differential expression analysis  
     Scripts: [`run_DESeq2_*.R`](scripts/tximport_DESeq2_scripts/)  

- **Output**:  
  - Normalized counts per cell line: [`normalized_counts.csv`](intermediate/)  
    - Example: [`A549_normalized_counts.csv`](intermediate/A549_normalized_counts.csv)  
  - Differential expression results per cell line: [`*_DESeq2_results.csv`](DE_results/)  
    - Example: [`A549_DESeq2_results.csv`](DE_results/A549_DESeq2_results.csv)  

