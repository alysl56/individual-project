#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(dplyr); library(stringr)})

base <- "/gpfs01/home/alysl56/projects/DE_analysis/enrich"
upf   <- file.path(base, "crossline_GO_ORA_up.tsv")
dnf   <- file.path(base, "crossline_GO_ORA_down.tsv")

rd <- function(fp){
  if(!file.exists(fp)) stop("file not found: ", fp)
  utils::read.delim(fp, check.names = TRUE, stringsAsFactors = FALSE)
}

find_col <- function(df, line){
  nm <- names(df)
  cand <- grep(paste0("^p\\.adjust_", gsub("[^A-Za-z0-9]", "\\\\W+", line), "$"),
               nm, perl=TRUE, value=TRUE)
  if(length(cand)==0){
    key <- paste0("p.adjust_", gsub("[^A-Za-z0-9]+","", line))
    cand <- nm[tolower(gsub("[^A-Za-z0-9]+","", nm)) == tolower(key)]
  }
  if(length(cand)) cand[1] else NA_character_
}
getcol <- function(df, cname){
  if(!is.na(cname) && cname %in% names(df)) df[[cname]] else rep(NA_real_, nrow(df))
}
fmt <- function(x){
  z <- suppressWarnings(as.numeric(x))
  ifelse(is.na(z), "", formatC(z, format="e", digits=2))
}

mk <- function(df, direction){
  if(nrow(df)==0) return(df[0,])
  ca <- find_col(df, "A549")
  cc <- find_col(df, "Calu-3"); if(is.na(cc)) cc <- find_col(df, "Calu.3")
  ch <- find_col(df, "HepG2")
  cu <- find_col(df, "U937")
  tibble(
    Term = df$Description,
    Direction = direction,
    A549_FDR  = getcol(df, ca),
    Calu3_FDR = getcol(df, cc),
    HepG2_FDR = getcol(df, ch),
    U937_FDR  = getcol(df, cu),
    Lines_sig = df$n_lines_sig
  )
}

up <- mk(rd(upf), "Up")
dn <- mk(rd(dnf), "Down")

mini_all <- bind_rows(up, dn) %>%
  arrange(desc(Lines_sig), A549_FDR, Calu3_FDR, HepG2_FDR, U937_FDR)

mini_common <- mini_all %>%
  filter(Lines_sig >= 2) %>%
  mutate(across(c(A549_FDR, Calu3_FDR, HepG2_FDR, U937_FDR), fmt))

utils::write.table(mini_common,
  file.path(base, "crossline_ORA_common_table.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE)

utils::write.table(mini_all %>%
    mutate(across(c(A549_FDR, Calu3_FDR, HepG2_FDR, U937_FDR), fmt)),
  file.path(base, "crossline_ORA_all_table.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE)

cat("done\n")
