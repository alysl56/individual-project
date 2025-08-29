#!/usr/bin/env Rscript
.libPaths('~/Rlibs')
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(VennDiagram)
  library(grid)
  library(grDevices)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(tibble)
  library(readr)
  library(dplyr)
})

files <- c(
  A549   = "/gpfs01/home/alysl56/projects/DMSO_A549_RNAseq/DE_analysis/A549_DESeq2_results.csv",
  "Calu-3" = "/gpfs01/home/alysl56/projects/DMSO_Calu-3_RNAseq/DE_analysis/Calu3_DESeq2_results.csv",
  U937   = "/gpfs01/home/alysl56/projects/DMSO_U937_RNAseq/DE_analysis/U937_DESeq2_results.csv",
  HepG2  = "/gpfs01/home/alysl56/projects/DMSO_HepG2_RNAseq/DE_analysis/HepG2_DESeq2_results.csv"
)

thr_fc <- 1
thr_p  <- 0.05
out_dir <- "/gpfs01/home/alysl56/projects/DE_analysis/venn_plots"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_try <- function(f){
  x <- try(read.csv(f, check.names = FALSE), silent = TRUE)
  if(inherits(x, "try-error")) x <- read.csv(f, row.names = 1, check.names = FALSE)
  x
}

read_deg <- function(f){
  x <- read_try(f)
  nms <- tolower(names(x))
  if(!("symbol" %in% nms)){
    if(!("gene" %in% nms)){
      y <- try(read.csv(f, row.names = 1, check.names = FALSE), silent = TRUE)
      if(!inherits(y,"try-error")){ x <- y; x$gene <- rownames(x); nms <- tolower(names(x)) }
    }
    id_col <- if("gene" %in% nms) "gene" else if("ensembl" %in% nms) "ensembl" else if("ensembl_gene_id" %in% nms) "ensembl_gene_id" else if("gene_id" %in% nms) "gene_id" else names(x)[1]
    ids <- gsub("\\.\\d+$","", x[[id_col]])
    ann <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = ids, keytype = "ENSEMBL", columns = "SYMBOL"))
    ann <- ann[!duplicated(ann$ENSEMBL), ]
    x$ENSEMBL_clean <- gsub("\\.\\d+$","", x[[id_col]])
    x <- merge(x, ann, by.x = "ENSEMBL_clean", by.y = "ENSEMBL", all.x = TRUE, sort = FALSE)
    x$SYMBOL <- ifelse(is.na(x$SYMBOL) | x$SYMBOL=="", NA, x$SYMBOL)
  } else {
    sycol <- names(x)[which(nms=="symbol")[1]]
    x$SYMBOL <- x[[sycol]]
  }
  if(!("log2FoldChange" %in% names(x))) { m <- match("log2foldchange", tolower(names(x))); if(!is.na(m)) names(x)[m] <- "log2FoldChange" }
  if(!("padj" %in% names(x)))         { m <- match("padj",           tolower(names(x))); if(!is.na(m)) names(x)[m] <- "padj" }
  x <- x[!is.na(x$SYMBOL) & !is.na(x$padj) & !is.na(x$log2FoldChange), c("SYMBOL","padj","log2FoldChange")]
  x <- x[!duplicated(x$SYMBOL), ]
  x
}

lst <- lapply(files, read_deg)

get_sets <- function(L, thr_p, thr_fc){
  all  <- lapply(L, function(d) unique(d$SYMBOL[d$padj < thr_p & abs(d$log2FoldChange) >  thr_fc]))
  up   <- lapply(L, function(d) unique(d$SYMBOL[d$padj < thr_p & d$log2FoldChange >   thr_fc]))
  down <- lapply(L, function(d) unique(d$SYMBOL[d$padj < thr_p & d$log2FoldChange <  -thr_fc]))
  list(all=all, up=up, down=down)
}
sets <- get_sets(lst, thr_p, thr_fc)

col_map <- c(A549="#F8766D","Calu-3"="#7CAE00",U937="#C77CFF",HepG2="#00BFC4")
plot_vd <- function(xlist, title_main, stub){
  fill_cols <- unname(col_map[names(xlist)])
  pngf <- file.path(out_dir, paste0(stub, ".png"))
  pdff <- file.path(out_dir, paste0(stub, ".pdf"))
  vd <- VennDiagram::venn.diagram(x = xlist, category.names = names(xlist), filename = NULL, cex = 1.1, cat.cex = 1.1, fill = fill_cols, alpha = 0.6, main = title_main)
  if (requireNamespace("ragg", quietly = TRUE)) { ragg::agg_png(pngf, width = 7.5, height = 6.5, units = "in", res = 300); grid::grid.newpage(); grid::grid.draw(vd); grDevices::dev.off()
  } else { grDevices::png(pngf, width = 7.5, height = 6.5, units = "in", res = 300, type = "cairo"); grid::grid.newpage(); grid::grid.draw(vd); grDevices::dev.off() }
  grDevices::pdf(pdff, width = 7.5, height = 6.5); grid::grid.newpage(); grid::grid.draw(vd); grDevices::dev.off()
  labs <- names(xlist)
  combs <- unlist(lapply(1:length(labs), function(k) combn(labs, k, simplify = FALSE)), recursive = FALSE)
  rows <- lapply(combs, function(cs){
    s <- Reduce(intersect, xlist[cs])
    others <- setdiff(labs, cs)
    if(length(others)>0) s <- setdiff(s, unique(unlist(xlist[others])))
    data.frame(region=paste(cs, collapse="&"), n=length(s), genes=paste(s, collapse=";"), stringsAsFactors=FALSE)
  })
  reg_df <- do.call(rbind, rows)
  write.table(reg_df, file.path(out_dir, paste0(stub, "_regions.tsv")), sep="\t", row.names=FALSE, quote=FALSE)
  write.csv(reg_df,   file.path(out_dir, paste0(stub, "_regions.csv")), row.names=FALSE)
}

x4_all <- sets$all
x4_up  <- sets$up

plot_vd(x4_all, sprintf("DEG Venn Diagram: A549, HepG2, Calu-3, U937 (All, padj<%.2f, |log2FC|>%s)", thr_p, thr_fc), sprintf("Venn4_All_padj%.2f_LFC%s", thr_p, thr_fc))
plot_vd(x4_up,  sprintf("DEG Venn Diagram: A549, HepG2, Calu-3, U937 (Up, padj<%.2f, |log2FC|>%s)",  thr_p, thr_fc), sprintf("Venn4_Up_padj%.2f_LFC%s",  thr_p, thr_fc))

drop_u937 <- function(x) { y <- x; y$U937 <- NULL; y }
x3 <- lapply(sets, drop_u937)

plot_vd(x3$all,  sprintf("DEG Venn Diagram: A549, HepG2, Calu-3 (All, padj<%.2f, |log2FC|>%s)", thr_p, thr_fc), sprintf("Venn3noU937_All_padj%.2f_LFC%s", thr_p, thr_fc))
plot_vd(x3$up,   sprintf("DEG Venn Diagram: A549, HepG2, Calu-3 (Up, padj<%.2f, |log2FC|>%s)",  thr_p, thr_fc), sprintf("Venn3noU937_Up_padj%.2f_LFC%s",  thr_p, thr_fc))
plot_vd(x3$down, sprintf("DEG Venn Diagram: A549, HepG2, Calu-3 (Down, padj<%.2f, |log2FC|>%s)",thr_p, thr_fc), sprintf("Venn3noU937_Down_padj%.2f_LFC%s",thr_p, thr_fc))
