# individual-project
# U937 RNA-Seq Data Download (with and without DMSO)

This project documents the download and management of **RNA-Seq data for the U937 cell line under DMSO treatment and untreated control**.  
The data are obtained from the **NCBI GEO/SRA databases** and processed using the **University of Nottingham HPC (Ada cluster)**.

---

## 📂 Directory Structure

projects/
├── raw_data/
│ └── data_links/
│ ├── u937_withDMSO_links.txt # Download links for with DMSO samples
│ └── u937_withoutDMSO_links.txt # Download links for without DMSO samples
│
└── scripts/
└── download_fastq/
├── download_u937_withDMSO.sh # SLURM batch script for with DMSO
└── download_u937_withoutDMSO.sh # SLURM batch script for without DMSO

---

## 📄 File Descriptions

- **u937_withDMSO_links.txt**  
  List of FTP links to paired-end FASTQ files (with DMSO condition).  

- **u937_withoutDMSO_links.txt**  
  List of FTP links to paired-end FASTQ files (without DMSO condition).  

- **download_u937_withDMSO.sh**  
  SLURM batch script to download FASTQ files for the with DMSO condition.  
  - Partition: `defq`  
  - Resources: 1 node, 1 task, 1 CPU, 8 GB memory  
  - Max runtime: 48 hours  
  - Logs saved in the same directory  

- **download_u937_withoutDMSO.sh**  
  SLURM batch script to download FASTQ files for the without DMSO condition.  
  - Partition: `defq`  
  - Resources: 1 node, 1 task, 1 CPU, 8 GB memory  
  - Max runtime: 48 hours  
  - Logs saved in the same directory  

---

## ⚙️ Usage

### 1. Navigate to working directory
```bash
cd /gpfs01/home/alysl56/projects/DMSO_U937_RNAseq/raw_data/with_DMSO
cd /gpfs01/home/alysl56/projects/DMSO_U937_RNAseq/raw_data/without_DMSO

### 2. Submit the download job
For with DMSO:
sbatch /gpfs01/home/alysl56/projects/DMSO_U937_RNAseq/scripts/download_fastq/download_u937_withDMSO.sh
For without DMSO:
sbatch /gpfs01/home/alysl56/projects/DMSO_U937_RNAseq/scripts/download_fastq/download_u937_withoutDMSO.sh

