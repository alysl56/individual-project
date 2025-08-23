library(tximport)
library(readr)
library(DESeq2)

base <- "/gpfs01/home/alysl56/projects/DMSO_Calu-3_RNAseq"
txi <- readRDS(file.path(base,"DE_analysis","txi_counts.rds"))
samples <- read.csv(file.path(base,"samplesheet.csv"), colClasses=c("character","character","character"))
samples$quant_dir <- trimws(samples$quant_dir)

cts <- txi$counts
cts <- cts[rowSums(cts)>0,,drop=FALSE]
cts[is.na(cts)] <- 0
samples <- samples[match(colnames(cts), samples$sample), , drop=FALSE]
mode(cts) <- "integer"

samples$condition <- factor(samples$condition, levels=c("without_DMSO","with_DMSO"))
dds <- DESeqDataSetFromMatrix(countData=cts, colData=samples, design=~condition)
ok <- TRUE
tryCatch({
  dds <<- estimateSizeFactors(dds, type="poscounts")
  if (any(is.na(sizeFactors(dds)))) ok <<- FALSE
}, error=function(e){ ok <<- FALSE })
if (!ok) {
  sf <- DESeq2::estimateSizeFactorsForMatrix(cts+1)
  sizeFactors(dds) <- sf
}
dds <- DESeq(dds, sfType="poscounts")
res <- results(dds, contrast=c("condition","with_DMSO","without_DMSO"))
normCounts <- counts(dds, normalized=TRUE)

out_dir <- file.path(base,"DE_analysis")
dir.create(out_dir, showWarnings=FALSE)
write.csv(as.data.frame(res), file.path(out_dir,"Calu3_DESeq2_results.csv"))
write.csv(as.data.frame(normCounts), file.path(out_dir,"Calu3_normalized_counts.csv"))
saveRDS(dds, file.path(out_dir,"Calu3_dds.rds"))
