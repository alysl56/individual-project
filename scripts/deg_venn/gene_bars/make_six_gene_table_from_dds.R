#!/usr/bin/env Rscript
.libPaths('~/Rlibs')
suppressPackageStartupMessages({
  library(dplyr); library(tidyr)
  library(DESeq2); library(SummarizedExperiment)
  library(AnnotationDbi); library(org.Hs.eg.db)
})

core6_tsv <- "/gpfs01/home/alysl56/projects/DE_analysis/venn_table/core6_selected.tsv"
out_csv   <- "/gpfs01/home/alysl56/projects/DE_analysis/venn_table/six_genes_expression_long_table.csv"

dds_paths <- c(
  "A549"   = "/gpfs01/home/alysl56/projects/DMSO_A549_RNAseq/DE_analysis/A549_dds.rds",
  "Calu-3" = "/gpfs01/home/alysl56/projects/DMSO_Calu-3_RNAseq/DE_analysis/Calu3_dds.rds",
  "HepG2"  = "/gpfs01/home/alysl56/projects/DMSO_HepG2_RNAseq/DE_analysis/HepG2_dds.rds",
  "U937"   = "/gpfs01/home/alysl56/projects/DMSO_U937_RNAseq/DE_analysis/U937_dds.rds"
)

genes <- utils::read.delim(core6_tsv, check.names = FALSE, stringsAsFactors = FALSE)$SYMBOL |> unique()

map_symbols <- function(ids){
  ids2 <- gsub("\\.\\d+$","", ids)
  ann <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=ids2, keytype="ENSEMBL", columns="SYMBOL"))
  ann <- ann[!duplicated(ann$ENSEMBL),]
  dplyr::tibble(gene_id = ids, ENS = ids2) |>
    dplyr::left_join(dplyr::tibble(ENS = ann$ENSEMBL, SYMBOL = ann$SYMBOL), by="ENS") |>
    dplyr::mutate(SYMBOL = ifelse(is.na(SYMBOL) | SYMBOL=="", ENS, SYMBOL))
}

pick_condition <- function(cd){
  nm <- tolower(colnames(cd))
  i  <- which(nm %in% c("condition","group","treatment","status"))
  if(length(i)==0) i <- grep("cond|treat|group|status", nm)
  stopifnot(length(i)>=1)
  x <- as.character(cd[[ i[1] ]])
  y <- tolower(trimws(x))
  z <- ifelse(grepl("with|dmso|treat|plus|pos", y), "+",
              ifelse(grepl("without|ctrl|control|vehicle|untreat|minus|neg|no", y), "−", y))
  z[z %in% c("1","treated")] <- "+"
  z[z %in% c("0","control","ctrl")] <- "−"
  z
}

one_line <- function(line, path){
  dds <- readRDS(path)
  nc  <- as.data.frame(counts(dds, normalized = TRUE))
  nc$gene_id <- rownames(nc)
  sym <- rownames(nc)
  if(all(grepl("^ENSG", sym))) {
    m <- map_symbols(sym)
    nc <- m |> dplyr::select(gene_id, SYMBOL) |> dplyr::left_join(nc, by=c("gene_id"="gene_id"))
  } else {
    nc$SYMBOL <- sym
  }
  cd <- as.data.frame(SummarizedExperiment::colData(dds))
  stopifnot(all(colnames(nc)[!(colnames(nc) %in% c("gene_id","SYMBOL"))] %in% rownames(cd)))
  cond <- pick_condition(cd); names(cond) <- rownames(cd)
  long <- nc |>
    dplyr::select(SYMBOL, dplyr::all_of(rownames(cd))) |>
    tidyr::pivot_longer(-SYMBOL, names_to="sample", values_to="norm_count") |>
    dplyr::mutate(condition = cond[sample], cell_line = line) |>
    dplyr::filter(SYMBOL %in% genes) |>
    dplyr::transmute(symbol = SYMBOL, cell_line, sample, condition, norm_count)
  long
}

res <- dplyr::bind_rows(mapply(one_line, names(dds_paths), dds_paths, SIMPLIFY = FALSE))
res$cell_line <- factor(res$cell_line, levels = c("A549","Calu-3","HepG2","U937"))
res$condition <- ifelse(res$condition %in% c("+","−"), res$condition,
                        ifelse(tolower(res$condition) %in% c("with","with_dmso","withdmso","dmso","treated"), "+",
                               ifelse(tolower(res$condition) %in% c("without","without_dmso","withoutdmso","control","ctrl","untreated","vehicle"), "−", res$condition)))
res <- dplyr::arrange(res, symbol, cell_line, condition, sample)
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, out_csv, row.names = FALSE)
cat("written:", out_csv, "\n")
