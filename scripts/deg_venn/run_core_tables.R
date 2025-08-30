#!/usr/bin/env Rscript
.libPaths('~/Rlibs')
options(stringsAsFactors=FALSE, bitmapType="cairo")

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble)
  library(AnnotationDbi); library(org.Hs.eg.db)
  library(pheatmap)
})

f_A549  <- "/gpfs01/home/alysl56/projects/DMSO_A549_RNAseq/DE_analysis/A549_DESeq2_results.csv"
f_Calu3 <- "/gpfs01/home/alysl56/projects/DMSO_Calu-3_RNAseq/DE_analysis/Calu3_DESeq2_results.csv"
f_HepG2 <- "/gpfs01/home/alysl56/projects/DMSO_HepG2_RNAseq/DE_analysis/HepG2_DESeq2_results.csv"
f_U937  <- "/gpfs01/home/alysl56/projects/DMSO_U937_RNAseq/DE_analysis/U937_DESeq2_results.csv"

out_dir <- "/gpfs01/home/alysl56/projects/DE_analysis/venn_table"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

thr_p <- 0.05; thr_fc <- 1

read_try <- function(f){
  x <- try(read.csv(f, check.names=FALSE), silent=TRUE)
  if(inherits(x,"try-error")) x <- read.csv(f, row.names=1, check.names=FALSE)
  x
}
read_deg <- function(f){
  x <- read_try(f); nms <- tolower(names(x))
  if(!("symbol" %in% nms)){
    if(!("gene" %in% nms)){ y<-try(read.csv(f,row.names=1,check.names=FALSE),silent=TRUE); if(!inherits(y,"try-error")){x<-y; x$gene<-rownames(x); nms<-tolower(names(x))}}
    id_col <- if("gene"%in%nms) "gene" else if("ensembl"%in%nms) "ensembl" else if("ensembl_gene_id"%in%nms) "ensembl_gene_id" else if("gene_id"%in%nms) "gene_id" else names(x)[1]
    ids <- gsub("\\.\\d+$","", x[[id_col]])
    ann <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=ids, keytype="ENSEMBL", columns="SYMBOL"))
    ann <- ann[!duplicated(ann$ENSEMBL),]
    x$ENSEMBL_clean <- gsub("\\.\\d+$","", x[[id_col]])
    x <- merge(x, ann, by.x="ENSEMBL_clean", by.y="ENSEMBL", all.x=TRUE, sort=FALSE)
    x$SYMBOL <- ifelse(is.na(x$SYMBOL)|x$SYMBOL=="", NA, x$SYMBOL)
  } else {
    x$SYMBOL <- x[[names(x)[which(nms=="symbol")[1]]]]
  }
  if(!("log2FoldChange" %in% names(x))){ m<-match("log2foldchange",tolower(names(x))); if(!is.na(m)) names(x)[m]<-"log2FoldChange" }
  if(!("padj" %in% names(x))){ m<-match("padj",tolower(names(x))); if(!is.na(m)) names(x)[m]<-"padj" }
  keep <- intersect(c("SYMBOL","log2FoldChange","padj","baseMean"), names(x))
  x <- x[, keep, drop=FALSE]
  x <- x[!is.na(x$SYMBOL) & !is.na(x$log2FoldChange) & !is.na(x$padj), ]
  x <- x[!duplicated(x$SYMBOL), ]; rownames(x) <- x$SYMBOL; x
}
A549  <- read_deg(f_A549)
Calu3 <- read_deg(f_Calu3)
HepG2 <- read_deg(f_HepG2)
U937  <- read_deg(f_U937)

sel <- function(df, type=c("all","up","down")){
  type <- match.arg(type)
  if(type=="all") df %>% filter(padj<thr_p, abs(log2FoldChange)>thr_fc) %>% pull(SYMBOL) %>% unique()
  else if(type=="up") df %>% filter(padj<thr_p, log2FoldChange>thr_fc) %>% pull(SYMBOL) %>% unique()
  else df %>% filter(padj<thr_p, log2FoldChange< -thr_fc) %>% pull(SYMBOL) %>% unique()
}

all_A <- sel(A549,"all");  up_A <- sel(A549,"up");  dn_A <- sel(A549,"down")
all_C <- sel(Calu3,"all"); up_C <- sel(Calu3,"up"); dn_C <- sel(Calu3,"down")
all_H <- sel(HepG2,"all"); up_H <- sel(HepG2,"up"); dn_H <- sel(HepG2,"down")

common3_all  <- Reduce(intersect, list(all_A, all_C, all_H))
common3_up   <- Reduce(intersect, list(up_A,  up_C,  up_H))
common3_down <- Reduce(intersect, list(dn_A,  dn_C,  dn_H))

write_tsv(tibble(gene=common3_all),  file.path(out_dir,"common3_All.tsv"))
write_tsv(tibble(gene=common3_up),   file.path(out_dir,"common3_Up.tsv"))
write_tsv(tibble(gene=common3_down), file.path(out_dir,"common3_Down.tsv"))

pull_cols <- function(df, genes, line){
  cols <- intersect(c("log2FoldChange","padj","baseMean"), names(df))
  out <- df[genes, cols, drop=FALSE]
  colnames(out) <- paste(line, colnames(out), sep="_"); out$SYMBOL <- rownames(out); out
}

genes_all  <- common3_all
genes_up   <- common3_up
genes_down <- common3_down

tab_all <- Reduce(function(x,y) dplyr::full_join(x,y, by="SYMBOL"),
                  list(pull_cols(A549, genes_all,"A549"),
                       pull_cols(Calu3,genes_all,"Calu3"),
                       pull_cols(HepG2,genes_all,"HepG2"))) %>% dplyr::relocate(SYMBOL)

tab_up <- Reduce(function(x,y) dplyr::full_join(x,y, by="SYMBOL"),
                 list(pull_cols(A549, genes_up,"A549"),
                      pull_cols(Calu3,genes_up,"Calu3"),
                      pull_cols(HepG2,genes_up,"HepG2"))) %>% dplyr::relocate(SYMBOL)

tab_down <- Reduce(function(x,y) dplyr::full_join(x,y, by="SYMBOL"),
                   list(pull_cols(A549, genes_down,"A549"),
                        pull_cols(Calu3,genes_down,"Calu3"),
                        pull_cols(HepG2,genes_down,"HepG2"))) %>% dplyr::relocate(SYMBOL)

write_tsv(tab_all,  file.path(out_dir,"common3_All_table.tsv"))
write_tsv(tab_up,   file.path(out_dir,"common3_Up_table.tsv"))
write_tsv(tab_down, file.path(out_dir,"common3_Down_table.tsv"))

rank_top_anydir <- function(tab, k=6){
  lfc_cols <- grep("log2FoldChange$", colnames(tab), value=TRUE)
  padj_cols <- grep("padj$", colnames(tab), value=TRUE)
  m <- as.matrix(tab[, lfc_cols])
  tab$mean_abs_log2FC <- rowMeans(abs(m))
  tab$max_padj <- apply(as.matrix(tab[, padj_cols]), 1, max, na.rm=TRUE)
  tab <- tab[order(-tab$mean_abs_log2FC, tab$max_padj), ]
  head(tab, min(k, nrow(tab)))
}

core6 <- rank_top_anydir(tab_all, k=6)
write_tsv(core6, file.path(out_dir,"core6_selected.tsv"))

lfc_mat <- function(genes){
  a <- A549[genes,"log2FoldChange"];  c <- Calu3[genes,"log2FoldChange"]
  h <- HepG2[genes,"log2FoldChange"]; u <- if(all(genes %in% rownames(U937))) U937[genes,"log2FoldChange"] else setNames(rep(NA_real_, length(genes)), genes)
  m <- cbind(A549=a, Calu3=c, HepG2=h, U937=u); rownames(m) <- genes; m
}

if(nrow(core6) > 0){
  mat6 <- lfc_mat(core6$SYMBOL)
  pngf <- file.path(out_dir,"core6_log2FC_heatmap.png")
  if (requireNamespace("ragg", quietly=TRUE)) {
    ragg::agg_png(pngf, width=1200, height=800, units="px", res=150)
  } else {
    grDevices::png(pngf, width=1200, height=800, res=150, type="cairo")
  }
  pheatmap(mat6, cluster_rows=FALSE, cluster_cols=FALSE, display_numbers=TRUE, number_format="%.2f", main="Core-6 log2FC (Top-6 by mean |log2FC|)")
  dev.off()
}
