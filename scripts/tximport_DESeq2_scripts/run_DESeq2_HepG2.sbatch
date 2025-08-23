#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12g
#SBATCH --time=04:00:00
#SBATCH --job-name=deseq2_HepG2
#SBATCH --output=/gpfs01/home/alysl56/projects/DMSO_HepG2_RNAseq/DE_analysis/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/alysl56/projects/DMSO_HepG2_RNAseq/DE_analysis/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alysl56@exmail.nottingham.ac.uk
module purge
module load R-uoneasy/4.4.1-gfbf-2023b-rstudio
Rscript /gpfs01/home/alysl56/projects/DMSO_HepG2_RNAseq/DE_analysis/run_DESeq2_HepG2.R
