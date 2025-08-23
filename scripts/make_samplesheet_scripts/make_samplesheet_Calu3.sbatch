#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=01:00:00
#SBATCH --job-name=make_samplesheet_Calu3
#SBATCH --output=/gpfs01/home/alysl56/projects/DMSO_Calu-3_RNAseq/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/alysl56/projects/DMSO_Calu-3_RNAseq/slurm-%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=alysl56@exmail.nottingham.ac.uk

BASE="/gpfs01/home/alysl56/projects/DMSO_Calu-3_RNAseq"
OUT="$BASE/samplesheet.csv"

echo "sample,condition,quant_dir" > "$OUT"

for COND in with_DMSO without_DMSO; do
  for D in "$BASE/salmon_quant/$COND"/*_quant; do
    [ -d "$D" ] || continue
    S=$(basename "$D"); S=${S%_quant}
    echo "$S,$COND,$D" >> "$OUT"
  done
done
