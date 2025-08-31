#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2g
#SBATCH --time=00:05:00
#SBATCH --job-name=mini_table
#SBATCH --output=/gpfs01/home/alysl56/projects/DE_analysis/logs/slurm-mini_table-%j.out
#SBATCH --error=/gpfs01/home/alysl56/projects/DE_analysis/logs/slurm-mini_table-%j.err

module purge
module load R-uoneasy/4.4.1-gfbf-2023b-rstudio

Rscript /gpfs01/home/alysl56/projects/DE_analysis/enrich/make_crossline_minitable.R
