#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(dplyr); library(stringr)})

base <- "/gpfs01/home/alysl56/projects/DE_analysis/enrich"
infile <- file.path(base, "crossline_ORA_common_rep_table.tsv")
out_top <- file.path(base, "crossline_ORA_top10_table.tsv")
out_log <- file.path(base, "crossline_ORA_top10_log.tsv")

df <- utils::read.delim(infile, check.names = FALSE, stringsAsFactors = FALSE)

numify <- function(x){suppressWarnings(as.numeric(x))}
fdr_cols <- c("A549_FDR","Calu3_FDR","HepG2_FDR","U937_FDR")
for(cn in fdr_cols){ if(cn %in% names(df)) df[[cn]] <- numify(df[[cn]]) }

df$minFDR <- apply(df[, fdr_cols], 1, function(v){
  v <- suppressWarnings(as.numeric(v)); if(all(is.na(v))) Inf else min(v, na.rm=TRUE)
})
df$sumFDR <- apply(df[, fdr_cols], 1, function(v){
  v <- suppressWarnings(as.numeric(v)); if(all(is.na(v))) Inf else sum(v, na.rm=TRUE)
})

class_theme <- function(term){
  s <- tolower(term)
  if(grepl("immune|antigen|mhc|defense|interferon|inflammatory|leukocyte|lymph", s)) return("Immune/Antigen")
  if(grepl("sensory|smell|olfactory|taste|chemical stimulus", s)) return("Sensory/Chemical")
  if(grepl("nucleosome|chromatin|cenp|histone|protein-?dna", s)) return("Chromatin/Nucleosome")
  if(grepl("rna|splic|mrna|ribonucle|transcript|processing", s)) return("RNA/Splicing")
  if(grepl("metabolic|lipid|oxid|redox|fatty acid|cholesterol", s)) return("Metabolic/Redox")
  "Other"
}
df$Theme <- vapply(df$Term, class_theme, character(1))

rank_df <- function(x){
  x %>% arrange(desc(Lines_sig), minFDR, sumFDR, nchar(Term))
}

target_total <- 10
per_dir_cap  <- 5
per_theme_cap <- 3

sel_log <- list()
pick_dir <- function(dd){
  res <- dd[0,]
  themes_used <- setNames(integer(0), character(0))
  cand <- rank_df(dd)
  for(i in seq_len(nrow(cand))){
    r <- cand[i,]
    th <- r$Theme
    themes_used[th] <- ifelse(!is.na(themes_used[th]), themes_used[th]+1L, 1L)
    if(themes_used[th] > per_theme_cap){
      themes_used[th] <- themes_used[th]-1L
      next
    }
    res <- bind_rows(res, r)
    sel_log[[length(sel_log)+1]] <<- data.frame(
      Term=r$Term, Direction=r$Direction, Theme=th,
      Reason="ranked", Lines_sig=r$Lines_sig, minFDR=r$minFDR, stringsAsFactors=FALSE)
    if(nrow(res) >= per_dir_cap) break
  }
  res
}

up_sel   <- pick_dir(df %>% filter(Direction=="Up"))
down_sel <- pick_dir(df %>% filter(Direction=="Down"))
top <- bind_rows(up_sel, down_sel)
if(nrow(top) < target_total){
  pool <- anti_join(rank_df(df), top, by="Term")
  for(i in seq_len(nrow(pool))){
    if(nrow(top) >= target_total) break
    r <- pool[i,]; top <- bind_rows(top, r)
    sel_log[[length(sel_log)+1]] <- data.frame(
      Term=r$Term, Direction=r$Direction, Theme=r$Theme,
      Reason="fill_to_target", Lines_sig=r$Lines_sig, minFDR=r$minFDR, stringsAsFactors=FALSE)
  }
}
top <- rank_df(top) %>% slice_head(n=target_total)

fmt <- function(x){ ifelse(is.na(x),"", formatC(x, format="e", digits=2)) }
top_out <- top %>% select(Term, Direction, all_of(fdr_cols), Lines_sig, Theme) %>%
  mutate(across(all_of(fdr_cols), fmt))

utils::write.table(top_out, out_top, sep="\t", quote=FALSE, row.names=FALSE)
utils::write.table(bind_rows(sel_log), out_log, sep="\t", quote=FALSE, row.names=FALSE)
cat("done\n")
