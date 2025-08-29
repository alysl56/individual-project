library(readr)
library(AnnotationDbi)
library(org.Hs.eg.db)
base <- "/gpfs01/home/alysl56/projects/DMSO_A549_RNAseq"
f <- file.path(base,"DE_analysis","A549_DESeq2_results.csv")
res <- read.csv(f, row.names=1, check.names=FALSE)
res <- res[!is.na(res$padj), , drop=FALSE]
top_n <- 10
total_genes <- nrow(res)
sig <- res[res$padj < 0.05 & abs(res$log2FoldChange) > 1, , drop=FALSE]
up <- sig[sig$log2FoldChange > 0, , drop=FALSE]
down <- sig[sig$log2FoldChange < 0, , drop=FALSE]
o_up <- up[order(up$padj, -up$log2FoldChange), , drop=FALSE]
o_down <- down[order(down$padj, down$log2FoldChange), , drop=FALSE]
top_up <- head(o_up, top_n)
top_down <- head(o_down, top_n)
annotate_tbl <- function(x){
  if(nrow(x)==0) return(data.frame(Gene=character(), Symbol=character(), log2FoldChange=numeric(), padj=numeric()))
  ens <- sub("\\.\\d+$","", rownames(x))
  sym <- mapIds(org.Hs.eg.db, keys=ens, keytype="ENSEMBL", column="SYMBOL", multiVals="first")
  data.frame(Gene=rownames(x), Symbol=unname(sym[ens]), log2FoldChange=x$log2FoldChange, padj=x$padj, check.names=FALSE)
}
summ <- data.frame(metric=c("genes_total","deg_total","up","down"), value=c(total_genes,nrow(sig),nrow(up),nrow(down)))
dir.create(file.path(base,"DE_analysis"), showWarnings=FALSE, recursive=TRUE)
write.csv(summ, file.path(base,"DE_analysis","A549_DEG_summary.csv"), row.names=FALSE)
write.csv(annotate_tbl(top_up), file.path(base,"DE_analysis","A549_top10_up.csv"), row.names=FALSE)
write.csv(annotate_tbl(top_down), file.path(base,"DE_analysis","A549_top10_down.csv"), row.names=FALSE)
