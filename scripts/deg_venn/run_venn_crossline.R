suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(ggVennDiagram)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(readr)
})

files <- c(
  A549   = "/gpfs01/home/alysl56/projects/DMSO_A549_RNAseq/DE_analysis/A549_DESeq2_results.csv",
  `Calu-3` = "/gpfs01/home/alysl56/projects/DMSO_Calu-3_RNAseq/DE_analysis/Calu3_DESeq2_results.csv",
  U937   = "/gpfs01/home/alysl56/projects/DMSO_U937_RNAseq/DE_analysis/U937_DESeq2_results.csv",
  HepG2  = "/gpfs01/home/alysl56/projects/DMSO_HepG2_RNAseq/DE_analysis/HepG2_DESeq2_results.csv"
)

thr_fc <- 1
thr_p  <- 0.05

out_dir <- "/gpfs01/home/alysl56/projects/DE_analysis/venn_plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

read_deg <- function(f){
  x <- read.csv(f, check.names = FALSE)
  nms <- tolower(names(x))
  sym_col <- which(nms %in% c("symbol","gene_symbol","hgnc_symbol"))
  if(length(sym_col)){
    x$SYMBOL <- x[[sym_col[1]]]
  } else {
    if(!("gene" %in% nms)) {
      x <- tibble::rownames_to_column(x, "gene")
      nms <- tolower(names(x))
    }
    id_col <- which(nms %in% c("gene","ensembl","ensembl_gene_id","gene_id"))[1]
    ids <- gsub("\\.\\d+$","", x[[id_col]])
    ann <- AnnotationDbi::select(org.Hs.eg.db, keys = ids, keytype = "ENSEMBL", columns = "SYMBOL")
    ann <- ann %>% distinct(ENSEMBL, .keep_all = TRUE)
    x <- x %>% mutate(ENSEMBL_clean = gsub("\\.\\d+$","", .data[[id_col]])) %>%
      left_join(ann, by = c("ENSEMBL_clean" = "ENSEMBL"))
  }
  names(x)[match("log2foldchange", tolower(names(x)))] <- "log2FoldChange"
  names(x)[match("padj", tolower(names(x)))] <- "padj"
  x %>%
    mutate(SYMBOL = ifelse(is.na(SYMBOL) | SYMBOL=="", NA, SYMBOL)) %>%
    filter(!is.na(SYMBOL), !is.na(padj), !is.na(log2FoldChange)) %>%
    distinct(SYMBOL, .keep_all = TRUE)
}

degs <- lapply(files, read_deg)

make_sets <- function(lst, thr_p, thr_fc){
  all  <- lapply(lst, \(d) d %>% filter(padj < thr_p, abs(log2FoldChange) > thr_fc) %>% pull(SYMBOL) %>% unique())
  up   <- lapply(lst, \(d) d %>% filter(padj < thr_p, log2FoldChange >  thr_fc) %>% pull(SYMBOL) %>% unique())
  down <- lapply(lst, \(d) d %>% filter(padj < thr_p, log2FoldChange < -thr_fc) %>% pull(SYMBOL) %>% unique())
  list(all=all, up=up, down=down)
}
sets <- make_sets(degs, thr_p, thr_fc)

plot_venn4 <- function(xlist, title_main, fn_stub){
  p <- ggVennDiagram(xlist, label = "count", label_size = 3) +
    ggtitle(title_main) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(out_dir, paste0(fn_stub, ".png")), p, width = 7.5, height = 6.5, dpi = 300)
  ggsave(file.path(out_dir, paste0(fn_stub, ".pdf")), p, width = 7.5, height = 6.5)
  pd <- ggVennDiagram::process_data(xlist)
  reg_df <- tibble::tibble(
    region = names(pd@region$item),
    n      = vapply(pd@region$item, length, integer(1)),
    genes  = vapply(pd@region$item, function(v) paste(v, collapse=";"), character(1))
  )
  readr::write_tsv(reg_df, file.path(out_dir, paste0(fn_stub, "_regions.tsv")))
  readr::write_csv(reg_df, file.path(out_dir, paste0(fn_stub, "_regions.csv")))
  invisible(p)
}

ttl_all  <- sprintf("All DEGs (padj < %.02f & |log2FC| > %s)", thr_p, thr_fc)
ttl_down <- sprintf("Down-regulated DEGs (padj < %.02f & |log2FC| > %s)", thr_p, thr_fc)
ttl_up   <- sprintf("Up-regulated DEGs (padj < %.02f & |log2FC| > %s)", thr_p, thr_fc)

plot_venn4(sets$all,  ttl_all,  sprintf("Venn_All_padj%.2f_LFC%s",  thr_p, thr_fc))
plot_venn4(sets$down, ttl_down, sprintf("Venn_Down_padj%.2f_LFC%s", thr_p, thr_fc))
plot_venn4(sets$up,   ttl_up,   sprintf("Venn_Up_padj%.2f_LFC%s",   thr_p, thr_fc))

up_list   <- sets$up;  names(up_list) <- names(files)
down_list <- sets$down; names(down_list) <- names(files)
all_genes <- unique(unlist(c(up_list, down_list)))
mat_up   <- sapply(up_list,   function(s) all_genes %in% s) * 1
mat_down <- sapply(down_list, function(s) all_genes %in% s) * 1
core_up   <- all_genes[rowSums(mat_up)   >= 3]
core_down <- all_genes[rowSums(mat_down) >= 3]
readr::write_tsv(tibble(type="core_up", gene=core_up),   file.path(out_dir, "core_up_genes.tsv"))
readr::write_tsv(tibble(type="core_down", gene=core_down), file.path(out_dir, "core_down_genes.tsv"))
spec_up_U937 <- setdiff(up_list$U937, unique(unlist(up_list[names(up_list)!="U937"])))
readr::write_tsv(tibble(line="U937", direction="UP", gene=spec_up_U937),
                 file.path(out_dir, "specific_U937_UP.tsv"))
