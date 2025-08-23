#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --time=48:00:00
#SBATCH --job-name=quant_salmon_Calu3
#SBATCH --output=/gpfs01/home/alysl56/projects/DMSO_Calu-3_RNAseq/salmon_quant/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/alysl56/projects/DMSO_Calu-3_RNAseq/salmon_quant/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alysl56@exmail.nottingham.ac.uk

source $HOME/.bash_profile
conda activate salmon_env

INDEX=/gpfs01/home/alysl56/references/gencode_v38_index
TRIM_DIR=/gpfs01/home/alysl56/projects/DMSO_Calu-3_RNAseq/raw_data/trimmed_data
OUT_DIR=/gpfs01/home/alysl56/projects/DMSO_Calu-3_RNAseq/salmon_quant

mkdir -p "${OUT_DIR}/with_DMSO" "${OUT_DIR}/without_DMSO"

for COND in with_DMSO without_DMSO; do
    cd ${TRIM_DIR}/${COND}
    for fq1 in *_1.trimmed.fastq.gz; do
        sample=${fq1%_1.trimmed.fastq.gz}
        fq2=${sample}_2.trimmed.fastq.gz
        salmon quant \
            -i ${INDEX} \
            -l A \
            -1 ${fq1} \
            -2 ${fq2} \
            -p 1 \
            --validateMappings \
            -o ${OUT_DIR}/${COND}/${sample}_quant
    done
done

conda deactivate
