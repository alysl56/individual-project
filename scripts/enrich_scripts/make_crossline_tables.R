#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(stringr)})

base <- "/gpfs01/home/alysl56/projects/DE_analysis/enrich"
lines <- c("A549","Calu-3","HepG2","U937")

read_tsv <- function(fp){
  if(!file.exists(fp)) return(NULL)
  x <- try(utils::read.delim(fp, check.names=FALSE, stringsAsFactors=FALSE), silent=TRUE)
  if(inherits(x,"try-error")) return(NULL)
  x
}

collect_ora <- function(direction=c("up","down")){
  direction <- match.arg(direction)
  dfs <- list()
  for(cl in lines){
    fp <- file.path(base, cl, paste0("ORA_GO_BP_", direction, ".tsv"))
    tb <- read_tsv(fp)
    if(is.null(tb) || !all(c("Description","p.adjust") %in% colnames(tb))) next
    tb <- tb %>% select(Description, p.adjust) %>% mutate(line = cl)
    dfs[[cl]] <- tb
  }
  if(length(dfs)==0) return(NULL)
  bind_rows(dfs) %>%
    mutate(Description = as.character(Description)) %>%
    mutate(is_sig = ifelse(!is.na(p.adjust) & p.adjust <= 0.05, 1L, 0L)) %>%
    select(Description, line, p.adjust, is_sig) %>%
    tidyr::pivot_wider(names_from=line, values_from=c(p.adjust, is_sig))
}

ora_up   <- collect_ora("up")
ora_down <- collect_ora("down")

write_out <- function(df, fn){
  if(is.null(df)){cat("no data:", fn, "\n"); return(invisible(NULL))}
  padj_cols <- grep("^p.adjust_", names(df), value=TRUE)
  df$n_lines_sig <- df %>% select(starts_with("is_sig_")) %>% rowSums(na.rm=TRUE)
  out <- df %>% arrange(desc(n_lines_sig), !!!rlang::syms(padj_cols))
  utils::write.table(out, file.path(base, fn), sep="\t", quote=FALSE, row.names=FALSE)
}

write_out(ora_up,   "crossline_GO_ORA_up.tsv")
write_out(ora_down, "crossline_GO_ORA_down.tsv")

collect_gsea <- function(){
  dfs <- list()
  for(cl in lines){
    fp <- file.path(base, cl, "GSEA_GO_BP.tsv")
    tb <- read_tsv(fp)
    if(is.null(tb)) next
    have <- intersect(c("ID","Description","NES","p.adjust"), colnames(tb))
    if(!all(c("Description","NES","p.adjust") %in% have)) next
    tb <- tb %>% select(Description, NES, p.adjust) %>% mutate(line=cl)
    dfs[[cl]] <- tb
  }
  if(length(dfs)==0) return(NULL)
  bind_rows(dfs) %>%
    tidyr::pivot_wider(names_from=line,
                       values_from=c(NES, p.adjust),
                       names_sep="_")
}

gsea_all <- collect_gsea()
if(!is.null(gsea_all)){
  utils::write.table(gsea_all, file.path(base, "crossline_GO_GSEA.tsv"),
                     sep="\t", quote=FALSE, row.names=FALSE)
}
cat("done\n")
