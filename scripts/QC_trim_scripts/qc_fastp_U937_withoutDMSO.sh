#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --time=48:00:00
#SBATCH --job-name=fastp_u937_withoutDMSO
#SBATCH --output=/gpfs01/home/alysl56/projects/DMSO_U937_RNAseq/qc_reports/without_DMSO/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/alysl56/projects/DMSO_U937_RNAseq/qc_reports/without_DMSO/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alysl56@exmail.nottingham.ac.uk

source $HOME/.bash_profile
conda activate rna_seq_env

RAW_DIR=/gpfs01/home/alysl56/projects/DMSO_U937_RNAseq/raw_data/without_DMSO
TRIM_DIR=/gpfs01/home/alysl56/projects/DMSO_U937_RNAseq/raw_data/trimmed_data/without_DMSO
QC_DIR=/gpfs01/home/alysl56/projects/DMSO_U937_RNAseq/qc_reports/without_DMSO

mkdir -p "$TRIM_DIR" "$QC_DIR"
cd "$RAW_DIR"

for fq1 in *_1.fastq.gz; do
    fq2=${fq1/_1.fastq.gz/_2.fastq.gz}
    base=${fq1%%_1.fastq.gz}
    fastp -i "$fq1" -I "$fq2" \
          -o "${TRIM_DIR}/${base}_1.trimmed.fastq.gz" \
          -O "${TRIM_DIR}/${base}_2.trimmed.fastq.gz" \
          -q 20 -u 30 -n 5 -l 50 -w 1 \
          -h "${QC_DIR}/${base}.html" \
          -j "${QC_DIR}/${base}.json"
done

conda deactivate
