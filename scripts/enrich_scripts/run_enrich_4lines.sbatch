#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16g
#SBATCH --time=04:00:00
#SBATCH --job-name=enrich_4lines
#SBATCH --output=/gpfs01/home/alysl56/projects/DE_analysis/logs/slurm-enrich_4lines-%j.out
#SBATCH --error=/gpfs01/home/alysl56/projects/DE_analysis/logs/slurm-enrich_4lines-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=alysl56@email.nottingham.ac.uk

module purge
module load R-uoneasy/4.4.1-gfbf-2023b-rstudio

Rscript /gpfs01/home/alysl56/projects/DE_analysis/run_enrich_4lines.R
