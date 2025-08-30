#!/usr/bin/env Rscript
.libPaths('~/Rlibs')
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tidyr)
  library(org.Hs.eg.db); library(AnnotationDbi)
  library(clusterProfiler)
  library(ggplot2)
})

thr_p <- 0.05
thr_fc <- 1

in_csv <- c(
  "A549"   = "/gpfs01/home/alysl56/projects/DMSO_A549_RNAseq/DE_analysis/A549_DESeq2_results.csv",
  "Calu-3" = "/gpfs01/home/alysl56/projects/DMSO_Calu-3_RNAseq/DE_analysis/Calu3_DESeq2_results.csv",
  "HepG2"  = "/gpfs01/home/alysl56/projects/DMSO_HepG2_RNAseq/DE_analysis/HepG2_DESeq2_results.csv",
  "U937"   = "/gpfs01/home/alysl56/projects/DMSO_U937_RNAseq/DE_analysis/U937_DESeq2_results.csv"
)

out_root <- "/gpfs01/home/alysl56/projects/DE_analysis/enrich"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

has_reactome <- requireNamespace("ReactomePA", quietly = TRUE)

read_res <- function(f){
  x <- utils::read.csv(f, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  x <- tibble::rownames_to_column(x, "gene")
  x$gene_clean <- sub("\\.\\d+$","", x$gene)
  ann <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = x$gene_clean, keytype = "ENSEMBL", columns = c("SYMBOL")))
  ann <- ann %>% distinct(ENSEMBL, .keep_all = TRUE)
  x <- x %>% left_join(ann, by = c("gene_clean"="ENSEMBL"))
  x$SYMBOL <- ifelse(is.na(x$SYMBOL) | x$SYMBOL=="", x$gene_clean, x$SYMBOL)
  x
}

plot_dot_bar <- function(df, file_prefix, topn = 20){
  if(is.null(df)) return(invisible(NULL))
  dd <- as.data.frame(df)
  utils::write.table(dd, paste0(file_prefix, ".tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  if(nrow(dd) == 0) return(invisible(NULL))
  top <- dd %>% mutate(neglog10FDR = -log10(p.adjust)) %>% arrange(p.adjust) %>% head(min(topn, nrow(dd)))
  p_dot <- ggplot(top, aes(x = neglog10FDR, y = reorder(Description, neglog10FDR))) +
    geom_point(aes(size = Count, color = neglog10FDR)) +
    scale_color_gradient(low="#6baed6", high="#08306b") +
    labs(x = "-log10(FDR)", y = NULL, title = "Dotplot") +
    theme_minimal(base_size=11)
  p_bar <- ggplot(top, aes(x = neglog10FDR, y = reorder(Description, neglog10FDR))) +
    geom_col() +
    labs(x = "-log10(FDR)", y = NULL, title = "Barplot") +
    theme_minimal(base_size=11)
  ggsave(paste0(file_prefix, "_dotplot.png"), p_dot, width=7.5, height=5.8, dpi=300)
  ggsave(paste0(file_prefix, "_barplot.png"), p_bar, width=7.5, height=5.8, dpi=300)
  ggsave(paste0(file_prefix, "_dotplot.pdf"), p_dot, width=7.5, height=5.8)
  ggsave(paste0(file_prefix, "_barplot.pdf"), p_bar, width=7.5, height=5.8)
  invisible(NULL)
}

save_gsea_tbl <- function(gse, file_prefix){
  if(is.null(gse)) return(invisible(NULL))
  df <- as.data.frame(gse) %>% arrange(p.adjust, -abs(NES))
  utils::write.table(df, paste0(file_prefix, ".tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  invisible(NULL)
}

for(cl in names(in_csv)){
  out_dir <- file.path(out_root, cl); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  res <- read_res(in_csv[[cl]])
  df  <- res %>% filter(!is.na(padj) & !is.na(log2FoldChange))
  df$stat2 <- if("stat" %in% names(df)) df$stat else sign(df$log2FoldChange) * -log10(pmax(df$pvalue, 1e-300))
  up   <- df %>% filter(padj < thr_p, log2FoldChange >  thr_fc) %>% pull(SYMBOL) %>% unique()
  down <- df %>% filter(padj < thr_p, log2FoldChange < -thr_fc) %>% pull(SYMBOL) %>% unique()
  univ <- df %>% pull(SYMBOL) %>% unique()
  ranks <- df$stat2; names(ranks) <- df$SYMBOL
  ranks <- ranks[is.finite(ranks) & !is.na(names(ranks))]
  ranks <- ranks[!duplicated(names(ranks))]

  ora_bp_up <- if(length(up) >= 5)
    enrichGO(gene = up, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
             universe = univ, pAdjustMethod = "BH", qvalueCutoff = 0.05,
             minGSSize = 10, maxGSSize = 5000, readable = TRUE) else NULL
  ora_bp_dn <- if(length(down) >= 5)
    enrichGO(gene = down, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
             universe = univ, pAdjustMethod = "BH", qvalueCutoff = 0.05,
             minGSSize = 10, maxGSSize = 5000, readable = TRUE) else NULL
  plot_dot_bar(ora_bp_up, file.path(out_dir, "ORA_GO_BP_up"))
  plot_dot_bar(ora_bp_dn, file.path(out_dir, "ORA_GO_BP_down"))

  gsea_bp <- if(length(ranks) >= 100)
    gseGO(geneList = sort(ranks, decreasing = TRUE),
          OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
          minGSSize = 10, maxGSSize = 5000, pAdjustMethod = "BH") else NULL
  save_gsea_tbl(gsea_bp, file.path(out_dir, "GSEA_GO_BP"))

  if(has_reactome){
    ReactomePA <- asNamespace("ReactomePA")
    ora_re_up <- if(length(up) >= 5)
      ReactomePA$enrichPathway(gene = up, organism = "human", pAdjustMethod = "BH",
                               universe = univ, qvalueCutoff = 0.05,
                               minGSSize = 10, maxGSSize = 5000, readable = TRUE) else NULL
    ora_re_dn <- if(length(down) >= 5)
      ReactomePA$enrichPathway(gene = down, organism = "human", pAdjustMethod = "BH",
                               universe = univ, qvalueCutoff = 0.05,
                               minGSSize = 10, maxGSSize = 5000, readable = TRUE) else NULL
    plot_dot_bar(ora_re_up, file.path(out_dir, "ORA_Reactome_up"))
    plot_dot_bar(ora_re_dn, file.path(out_dir, "ORA_Reactome_down"))

    gsea_re <- if(length(ranks) >= 100)
      ReactomePA$gsePathway(geneList = sort(ranks, decreasing = TRUE),
                            organism = "human", pAdjustMethod = "BH",
                            minGSSize = 10, maxGSSize = 5000) else NULL
    save_gsea_tbl(gsea_re, file.path(out_dir, "GSEA_Reactome"))
  }
}
cat("done\n")
