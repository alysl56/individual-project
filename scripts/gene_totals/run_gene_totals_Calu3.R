library(readr)
base <- "/gpfs01/home/alysl56/projects/DMSO_Calu-3_RNAseq"
samples <- read.csv(file.path(base,"samplesheet.csv"), colClasses=c("character","character","character"))
counts <- read.csv(file.path(base,"DE_analysis","gene_counts.csv"), row.names=1, check.names=FALSE)
tpm <- read.csv(file.path(base,"DE_analysis","gene_tpm.csv"), row.names=1, check.names=FALSE)
samples$quant_dir <- trimws(samples$quant_dir)
id2cond <- setNames(samples$condition, samples$sample)
det_counts_gt0 <- colSums(as.matrix(counts) > 0)
det_tpm_ge0.1 <- colSums(as.matrix(tpm) >= 0.1)
total_genes <- nrow(counts)
out_sample <- data.frame(
  sample = colnames(counts),
  condition = unname(id2cond[colnames(counts)]),
  genes_total = total_genes,
  detected_counts_gt0 = as.integer(det_counts_gt0),
  detected_tpm_ge0.1 = as.integer(det_tpm_ge0.1),
  stringsAsFactors = FALSE
)
conds <- unique(out_sample$condition)
summ <- do.call(rbind, lapply(conds, function(cn){
  ix <- out_sample$condition == cn
  data.frame(
    condition = cn,
    n_samples = sum(ix),
    mean_counts_gt0 = mean(out_sample$detected_counts_gt0[ix]),
    sd_counts_gt0 = ifelse(sum(ix)>1, sd(out_sample$detected_counts_gt0[ix]), NA),
    mean_tpm_ge0.1 = mean(out_sample$detected_tpm_ge0.1[ix]),
    sd_tpm_ge0.1 = ifelse(sum(ix)>1, sd(out_sample$detected_tpm_ge0.1[ix]), NA),
    stringsAsFactors = FALSE
  )
}))
dir.create(file.path(base,"DE_analysis"), showWarnings=FALSE, recursive=TRUE)
write.csv(out_sample, file.path(base,"DE_analysis","Calu3_gene_totals_per_sample.csv"), row.names=FALSE)
write.csv(summ, file.path(base,"DE_analysis","Calu3_gene_totals_by_condition.csv"), row.names=FALSE)
