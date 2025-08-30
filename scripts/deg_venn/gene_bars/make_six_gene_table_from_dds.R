#!/usr/bin/env Rscript
.libPaths('~/Rlibs')
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(DESeq2); library(SummarizedExperiment); library(AnnotationDbi); library(org.Hs.eg.db) })

core6_tsv <- "/gpfs01/home/alysl56/projects/DE_analysis/venn_table/core6_selected.tsv"
out_csv   <- "/gpfs01/home/alysl56/projects/DE_analysis/venn_table/six_genes_expression_long_table.csv"
out_sum   <- "/gpfs01/home/alysl56/projects/DE_analysis/venn_table/six_genes_condition_summary.tsv"

dds_paths <- c(
  "A549"   = "/gpfs01/home/alysl56/projects/DMSO_A549_RNAseq/DE_analysis/A549_dds.rds",
  "Calu-3" = "/gpfs01/home/alysl56/projects/DMSO_Calu-3_RNAseq/DE_analysis/Calu3_dds.rds",
  "HepG2"  = "/gpfs01/home/alysl56/projects/DMSO_HepG2_RNAseq/DE_analysis/HepG2_dds.rds"
)

genes <- utils::read.delim(core6_tsv, check.names = FALSE, stringsAsFactors = FALSE)$SYMBOL |> unique()

map_symbols <- function(ids){
  ids2 <- gsub("\\.\\d+$","", ids)
  ann <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=ids2, keytype="ENSEMBL", columns="SYMBOL"))
  ann <- ann[!is.na(ann$ENSEMBL),]
  ann <- ann[!duplicated(ann$ENSEMBL),]
  dplyr::tibble(gene_id = ids, ENS = ids2) |>
    dplyr::left_join(dplyr::tibble(ENS = ann$ENSEMBL, SYMBOL = ann$SYMBOL), by="ENS") |>
    dplyr::mutate(SYMBOL = ifelse(is.na(SYMBOL) | SYMBOL=="", ENS, SYMBOL))
}

build_cond <- function(cd, samples){
  txt <- apply(cd, 1, function(r) paste(tolower(trimws(as.character(r))), collapse="|"))
  sname <- tolower(samples)
  neg <- grepl("without|no[_-]?dmso|\\bctrl\\b|control|vehicle|untreat|untreated|minus|\\bneg\\b|mock|blank|naive|^0$|\\bfalse\\b", txt) |
         grepl("without|no[_-]?dmso|\\bctrl\\b|control|vehicle|untreat|untreated|minus|\\bneg\\b|mock|blank|naive|^0$|\\bfalse\\b", sname)
  pos <- grepl("with[_-]?dmso|\\bdmso\\b|treated|treat|plus|\\bpos\\b|^1$|\\btrue\\b", txt) |
         grepl("with[_-]?dmso|\\bdmso\\b|treated|treat|plus|\\bpos\\b|^1$|\\btrue\\b", sname)
  out <- ifelse(neg, "−", ifelse(pos, "+", NA_character_))
  names(out) <- samples
  out
}

one_line <- function(line, path){
  dds <- readRDS(path)
  nc  <- as.data.frame(counts(dds, normalized = TRUE))
  nc$gene_id <- rownames(nc)
  sym <- rownames(nc)
  if(all(grepl("^ENSG", sym))){
    m <- map_symbols(sym)
    nc <- m |> dplyr::select(gene_id, SYMBOL) |> dplyr::left_join(nc, by=c("gene_id"="gene_id"))
  } else {
    nc$SYMBOL <- sym
  }
  cd <- as.data.frame(SummarizedExperiment::colData(dds))
  smp <- rownames(cd)
  cond <- build_cond(cd, smp)
  nc |>
    dplyr::select(SYMBOL, dplyr::all_of(smp)) |>
    tidyr::pivot_longer(-SYMBOL, names_to="sample", values_to="norm_count") |>
    dplyr::mutate(condition = cond[sample], cell_line = line) |>
    dplyr::filter(SYMBOL %in% genes) |>
    dplyr::transmute(symbol = SYMBOL, cell_line, sample, condition, norm_count)
}

res <- dplyr::bind_rows(mapply(one_line, names(dds_paths), dds_paths, SIMPLIFY = FALSE))
res$cell_line <- factor(res$cell_line, levels = c("A549","Calu-3","HepG2"))
res <- dplyr::arrange(res, symbol, cell_line, condition, sample)

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, out_csv, row.names = FALSE)

sumtab <- res |>
  dplyr::mutate(condition = dplyr::if_else(is.na(condition), "NA", condition)) |>
  dplyr::count(cell_line, condition) |>
  dplyr::arrange(cell_line, condition)
utils::write.table(sumtab, out_sum, sep="\t", row.names = FALSE, quote = FALSE)

cat("written:", out_csv, "\n")
cat("summary:", out_sum, "\n")
