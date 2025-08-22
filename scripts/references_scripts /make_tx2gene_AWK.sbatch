#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=00:20:00
#SBATCH --job-name=make_tx2gene_awk
#SBATCH --output=/gpfs01/home/alysl56/references/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/alysl56/references/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alysl56@exmail.nottingham.ac.uk

cd /gpfs01/home/alysl56/references
{ echo "transcript_id,gene_id"; awk '$3=="transcript"{match($0,/transcript_id "([^"]+)"/,a);match($0,/gene_id "([^"]+)"/,b);if(a[1]&&b[1])print a[1]","b[1]}' gencode.v38.annotation.gtf | sort -u; } > tx2gene.csv
