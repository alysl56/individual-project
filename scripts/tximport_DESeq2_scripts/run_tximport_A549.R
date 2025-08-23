library(tximport)
library(readr)

base <- "/gpfs01/home/alysl56/projects/DMSO_A549_RNAseq"

samples <- read.csv(file.path(base, "samplesheet.csv"), colClasses = c("character","character","character"))
samples$quant_dir <- trimws(samples$quant_dir)

tx2gene <- read.csv("/gpfs01/home/alysl56/references/tx2gene.csv", colClasses = c("character","character"))

files <- file.path(samples$quant_dir, "quant.sf")
names(files) <- samples$sample
stopifnot(all(file.exists(files)))

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

dir.create(file.path(base, "DE_analysis"), showWarnings = FALSE, recursive = TRUE)
write.csv(txi$counts,    file.path(base, "DE_analysis", "gene_counts.csv"))
write.csv(txi$abundance, file.path(base, "DE_analysis", "gene_tpm.csv"))
saveRDS(txi,             file.path(base, "DE_analysis", "txi_counts.rds"))
