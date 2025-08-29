#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --time=00:30:00
#SBATCH --job-name=venn_crossline
#SBATCH --output=/gpfs01/home/alysl56/projects/DMSO_A549_RNAseq/DE_analysis/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/alysl56/projects/DMSO_A549_RNAseq/DE_analysis/slurm-%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=alysl56@exmail.nottingham.ac.uk

module purge
module load R-uoneasy/4.4.1-gfbf-2023b-rstudio

Rscript /gpfs01/home/alysl56/projects/DMSO_crossline/DE_analysis/run_venn_crossline.R
