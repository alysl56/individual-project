#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(dplyr); library(stringr)})

base <- "/gpfs01/home/alysl56/projects/DE_analysis/enrich"
infile <- file.path(base, "crossline_ORA_common_table.tsv")
outfile_rep <- file.path(base, "crossline_ORA_common_rep_table.tsv")
outfile_clusters <- file.path(base, "crossline_ORA_clusters.tsv")

df <- utils::read.delim(infile, check.names = FALSE, stringsAsFactors = FALSE)

numify <- function(x){y<-suppressWarnings(as.numeric(x)); y}
fdr_cols <- c("A549_FDR","Calu3_FDR","HepG2_FDR","U937_FDR")
for(cn in fdr_cols){if(cn %in% names(df)) df[[cn]] <- numify(df[[cn]])}

df$minFDR <- apply(df[, fdr_cols], 1, function(v){
  v <- suppressWarnings(as.numeric(v)); if(all(is.na(v))) Inf else min(v, na.rm=TRUE)
})
df$sumFDR <- apply(df[, fdr_cols], 1, function(v){
  v <- suppressWarnings(as.numeric(v)); if(all(is.na(v))) Inf else sum(v, na.rm=TRUE)
})

stopwords <- c("of","the","and","to","in","by","for","with","via","on","into","within",
               "process","processes","pathway","pathways","activity","activities",
               "biological","cellular","regulation","positive","negative","regulations")
tokenize <- function(s){
  s <- tolower(s)
  s <- gsub("[^a-z0-9]+"," ", s)
  tt <- unique(strsplit(trimws(s), "\\s+")[[1]])
  tt <- tt[!(tt %in% stopwords)]
  tt[nchar(tt)>0]
}
jac <- function(a,b){
  if(length(a)==0 || length(b)==0) return(0)
  length(intersect(a,b))/length(union(a,b))
}

thr <- 0.60
df <- df %>% mutate(idx = row_number())
df$tokens <- lapply(df$Term, tokenize)

clusters <- list()
assigned <- rep(FALSE, nrow(df))

make_cluster <- function(sub){
  ord <- order(-sub$Lines_sig, sub$minFDR, sub$sumFDR, nchar(sub$Term))
  sub[ord, , drop=FALSE]
}

res_rep <- list()
res_map <- list()

for(dir_level in unique(df$Direction)){
  sub_all <- df %>% filter(Direction==dir_level)
  ord_seed <- order(-sub_all$Lines_sig, sub_all$minFDR, sub_all$sumFDR, nchar(sub_all$Term))
  sub_all <- sub_all[ord_seed, , drop=FALSE]
  used <- rep(FALSE, nrow(sub_all))
  for(i in seq_len(nrow(sub_all))){
    if(used[i]) next
    seed <- sub_all[i, , drop=FALSE]
    mem_idx <- i
    sim_keep <- 1
    for(j in (i+1):nrow(sub_all)){
      if(is.na(j) || j>nrow(sub_all)) break
      if(used[j]) next
      s <- jac(seed$tokens[[1]], sub_all$tokens[[j]])
      if(s >= thr){ mem_idx <- c(mem_idx, j); sim_keep <- c(sim_keep, s); used[j] <- TRUE }
    }
    cl <- sub_all[mem_idx, , drop=FALSE]
    cl <- make_cluster(cl)
    rep_row <- cl[1, , drop=FALSE]
    res_rep[[length(res_rep)+1]] <- rep_row %>%
      select(Term, Direction, all_of(fdr_cols)) %>%
      mutate(ClusterSize = nrow(cl))
    res_map[[length(res_map)+1]] <- data.frame(
      Representative = rep_row$Term,
      Direction = rep_row$Direction,
      Member = cl$Term,
      stringsAsFactors = FALSE
    )
  }
}

rep_tab <- bind_rows(res_rep) %>%
  arrange(Direction, desc(ClusterSize), A549_FDR, Calu3_FDR, HepG2_FDR, U937_FDR)

fmt <- function(x){ ifelse(is.na(x),"", formatC(x, format="e", digits=2)) }
rep_tab_out <- rep_tab %>%
  mutate(across(all_of(fdr_cols), fmt))

utils::write.table(rep_tab_out, outfile_rep, sep="\t", quote=FALSE, row.names=FALSE)
utils::write.table(bind_rows(res_map), outfile_clusters, sep="\t", quote=FALSE, row.names=FALSE)
cat("done\n")
