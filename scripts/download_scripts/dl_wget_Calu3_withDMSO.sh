#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --time=48:00:00
#SBATCH --job-name=dl_Calu3_withDMSO
#SBATCH --output=/gpfs01/home/alysl56/projects/DMSO_Calu3_RNAseq/raw_data/with_DMSO/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/alysl56/projects/DMSO_Calu3_RNAseq/raw_data/with_DMSO/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alysl56@exmail.nottingham.ac.uk

source $HOME/.bash_profile
conda activate rna_seq_env

cd /gpfs01/home/alysl56/projects/DMSO_Calu3_RNAseq/raw_data/with_DMSO
wget -c --tries=20 --timeout=60 --no-verbose -i Calu3_withDMSO_links.txt

conda deactivate
