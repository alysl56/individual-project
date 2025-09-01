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

The Salmon index was built using **Salmon v1.9.0** (Patro et al., 2017):
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

Preprocessing of raw FASTQ files was performed using **fastp v0.22.0** (Chen et al., 2018):

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
  Analyses were performed on the **Ada HPC cluster** in two modes:  

  - **Interactive R (for testing/exploration):**
    ```bash
    srun -p defq -c 4 --mem=24G -t 02:00:00 --pty bash
    module load R-uoneasy/4.4.1-gfbf-2023b-rstudio
    R
    ```

  - **Batch (Slurm, for full analysis):**  
    Example: [`run_DESeq2_A549.sbatch`](projects/scripts/tximport_DESeq2_scripts/run_DESeq2_A549.sbatch)

- **Required R packages**:  
  ```r
  library(DESeq2)           # v1.46.0 (Love et al., 2014)
  library(tximport)         # v1.34.0 (Soneson et al., 2015)
  library(clusterProfiler)  # v4.14.6 (Yu et al., 2012)
  library(org.Hs.eg.db)     # v3.20.0 (Carlson, 2019)
  library(ggplot2)          # v3.5.2 (Wickham, 2016)
  library(pheatmap)         # v1.0.13 (Kolde, 2019)
  library(VennDiagram)      # v1.7.3 (Chen and Boutros, 2011)

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


---
## 8. Single Cell-Line DEG Analysis

For each cell line, differential expression results (`*_DESeq2_results.csv`) were further summarized:

- **Input**: DESeq2 results (`*_DESeq2_results.csv`)
- **Process**: 
  - Count total genes, significant DEGs (padj < 0.05, |log2FC| > 1)
  - Extract top 10 up-regulated and down-regulated genes
  - Generate volcano plots for visualization
- **Output**: 
  - DEG summary (`*_DEG_summary.csv`)
  - Top 10 up/down gene lists (`*_top10_up.csv`, `*_top10_down.csv`)
  - Volcano plots (`*.pdf`)

Scripts are provided in:  
[`scripts/degstats_scripts/`](/scripts/degstats_scripts/)

- Examples:  
  - [`run_degstats_A549.R`](/scripts/degstats_scripts/run_degstats_A549.R), [`run_degstats_A549.sbatch`](/scripts/degstats_scripts/run_degstats_A549.sbatch) 
  - [`run_volcano_A549.R`](/scripts/degstats_scripts/run_volcano_A549.R),[`run_volcano_A549.sbatch`](/scripts/degstats_scripts/run_volcano_A549.sbatch`) 

Results are stored in:  
[DE_results/`](DE_results/)

- Example files:  
  - [`A549_DESeq2_results.csv`](DE_results/A549_DESeq2_results.csv) 
  - [`A549_DEG_summary.csv`](DE_results/A549_DEG_summary.csv) 
  - [`A549_top10_up.csv`](DE_results/A549_top10_up.csv) 


---
## 9. Single Cell-Line Gene Detection Totals

To evaluate overall gene detection per sample,  
we calculated **gene totals** (number of detected genes) for each cell line.  

- **Input**:  
  - `samplesheet.csv` (sample metadata)  
  - `gene_counts.csv` (raw counts per gene)  
  - `gene_tpm.csv` (TPM-normalized values)  

- **Process**:  
  - Count total number of genes per sample  
  - Calculate number of genes with counts > 0  
  - Calculate number of genes with TPM ≥ 0.1  
  - Summarize mean and SD of detected genes by condition (with/without DMSO)  

- **Output**:  
  - Per-sample gene totals (`*_gene_totals_per_sample.csv`)  
  - Condition-level summaries (`*_gene_totals_by_condition.csv`)  

Scripts are provided in:  
[`projects/scripts/gene_totals/`](scripts/gene_totals)

Examples:  
- [`run_gene_totals_A549.R`](scripts/gene_totals/run_gene_totals_A549.R)  
- [`run_gene_totals_A549.sbatch`](scripts/gene_totals/run_gene_totals_A549.sbatch)

Results are available in:  
[`coverage_summary/`](coverage_summary)  

- Example output: [`A549_gene_totals_per_sample.csv`](coverage_summary/A549_gene_totals_per_sample.csv), [`A549_gene_totals_by_condition.csv`](coverage_summary/A549_gene_totals_per_sample.csv)


---
## 10. Cross-cell line integration (Venn → Core genes → Bar plots)

- **Input**:  
  Per–cell-line DE tables from DESeq2:  
  [`DE_results/`](DE_results/)  

- **Scripts**:  
  - 4-line / 3-line Venn:  
    [`scripts/deg_venn/run_venn_crossline.R`](scripts/deg_venn/run_venn_crossline.R)  
    Example: [`run_venn_crossline.sbatch`](scripts/deg_venn/run_venn_crossline.sbatch)  
  - Core gene tables:  
    [`scripts/deg_venn/run_core_tables.R`](scripts/deg_venn/run_core_tables.R)  
    Example: [`run_core_tables.sbatch`](scripts/deg_venn/run_core_tables.sbatch)    
  - Core-6 bar plots:  
    - Generate table: [`scripts/deg_venn/gene_bars/make_six_gene_table_from_dds.R`](scripts/deg_venn/gene_bars/make_six_gene_table_from_dds.R)  
      Example: [`make_six_gene_table_from_dds.sbatch`](scripts/deg_venn/gene_bars/make_six_gene_table_from_dds.sbatch)  
    - Plot bars: [`scripts/deg_venn/gene_bars/run_gene_bars.R`](scripts/deg_venn/gene_bars/run_gene_bars.R)  
      Example: [`run_gene_bars.sbatch`](scripts/deg_venn/gene_bars/run_gene_bars.sbatch)  

- **Process**:  
  1. Construct **four-line all-up Venn** (A549, Calu-3, HepG2, U937).  
  2. Construct **three-epithelial Venn (A549/Calu-3/HepG2)** for both all-up and all-down.  
  3. Export intersection tables and identify shared sets.  
  4. Extract **Core-6 genes** (RPPH1, RNU1-1, RNU4-1, RNU4-2, SNORA74A, RNU2-1).  
  5. Summarize expression of Core-6 across with/without DMSO conditions.  
  6. Generate bar plots (mean ± s.e.m.) for visualization.  

- **Output**:  
  - Venn plots and intersection tables: [`DE_results/`](DE_results/)  
    - Example: [`Venn4_All_padj0.05_LFC1_regions.tsv`](DE_results/Venn4_All_padj0.05_LFC1_regions.tsv)  , [`common3_All_table.tsv`](DE_results/common3_All_table.tsv`) 
    - Core gene list: [`core6_selected.tsv`](DE_results/core6_selected.tsv`)   
  - Core-6 expression and bar plots: [`DE_results`](DE_results)  
    - Example: [`six_genes_expression_long_table.csv`](DE_results/six_genes_expression_long_table.csv)  

- **Execution examples**:  
  ```bash
  Rscript projects/scripts/deg_venn/run_venn_crossline.R
  sbatch projects/scripts/deg_venn/run_core_tables.sbatch
  Rscript projects/scripts/deg_venn/gene_bars/make_six_gene_table_from_dds.R
  sbatch projects/scripts/deg_venn/gene_bars/run_gene_bars.sbatch


---
## 11. Functional Enrichment Analysis

After identifying DEGs, functional enrichment (GO:BP) analyses were performed for each cell line (A549, Calu-3, HepG2, U937) using **clusterProfiler**. Two complementary approaches were applied:

- **Overrepresentation Analysis (ORA)**
- **Gene Set Enrichment Analysis (GSEA)**


- **Scripts:** [`scripts/enrich_scripts/`](scripts/enrich_scripts/)
- **Results:** [`enrichment_results/`](enrichment_results/)

---

### Input
- Differential expression results from DESeq2 ( e.g. [`A549_DESeq2_results.csv`](DE_results/A549_DESeq2_results.csv))
- DEG summary tables [`*_DEG_summary.csv`](DE_results/))
- Top DEG subsets [`*_top10_up.csv`](DE_results/))

---

### Scripts and Usage

- [`run_enrich_4lines.R`](scripts/enrich_scripts/run_enrich_4lines.R)
  R script to run ORA and GSEA for all four cell lines.  
  **Batch:** [`run_enrich_4lines.sbatch`](scripts/enrich_scripts/run_enrich_4lines.sbatch)

- Example outputs for A549:  
  - [`A549_ORA_GO_BP_up.tsv`](enrichment_results/A549_ORA_GO_BP_up.tsv) (ORA results, up-regulated DEGs)  
  - [`A549_ORA_GO_BP_down.tsv`](enrichment_results/A549_ORA_GO_BP_down.tsv) (ORA results, down-regulated DEGs)  
  - [`A549_GSEA_GO_BP.tsv`](enrichment_results/A549_GSEA_GO_BP.tsv)(GSEA results, all ranked DEGs)

Same format was used for **Calu3**, **HepG2**, and **U937**.

---

### Output
Results stored in:  
[`enrichment_results/`](enrichment_results)


---

 

---

## 12. Cross-line Enrichment Integration

To assess shared functional patterns across cell lines, enrichment results were integrated and summarized using custom R scripts:

[`scripts/enrich_scripts/`](scripts/enrich_scripts/)

### Example scripts
- [`make_crossline_top10.R`](scripts/enrich_scripts/make_crossline_top10.R) : Extracts the top-10 representative GO:BP terms across ≥2 cell lines  
- [`make_crossline_heatmap.R`](scripts/enrich_scripts/make_crossline_heatmap.R) : Generates heatmap of representative terms  
- [`make_crossline_minitable.R`](scripts/enrich_scripts/make_crossline_minitable.R): Produces compact summary tables of shared terms  
- [`make_crossline_representative.R`](scripts/enrich_scripts/make_crossline_representative.R): Clusters redundant GO terms and selects representatives  

### Example outputs
- Representative enrichment tables [`crossline_ORA_common_rep_table.tsv`](enrichment_results/crossline_ORA_common_rep_table.tsv)  


---

**Notes**:  
- In A549 and Calu-3, DMSO upregulated immune/stress pathways and downregulated chromatin/nucleosome assembly and sensory perception.  
- In HepG2, enrichment was biased towards metal ion homeostasis and developmental processes.  
- U937 showed very few DEGs, but ORA/GSEA results are still available in the enrichment_results/ folder.  

---

## Acknowledgements

This project was completed as part of my MSc dissertation at the University of Nottingham.  
I would like to thank **Prof. Patrick Tighe** **Dr. Hannah Jackson**and **Dr. Laura Dean** for their supervision, constructive feedback, and continuous support throughout the project.

---

## References


