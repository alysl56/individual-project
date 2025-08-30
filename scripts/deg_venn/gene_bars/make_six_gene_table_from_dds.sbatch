#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --time=00:15:00
#SBATCH --job-name=six_gene_table
#SBATCH --output=/gpfs01/home/alysl56/projects/DE_analysis/logs/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/alysl56/projects/DE_analysis/logs/slurm-%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=alysl56@exmail.nottingham.ac.uk

set -euo pipefail
module --force purge || true
module load R-uoneasy/4.4.1-gfbf-2023b-rstudio || true
Rscript --vanilla "$SLURM_SUBMIT_DIR/make_six_gene_table_from_dds.R"
