#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(stringr)})

base <- "/gpfs01/home/alysl56/projects/DE_analysis/enrich"
f_top10 <- file.path(base, "crossline_ORA_top10_table.tsv")
f_rep   <- file.path(base, "crossline_ORA_common_rep_table.tsv")
infile  <- if (file.exists(f_top10)) f_top10 else f_rep
df0 <- utils::read.delim(infile, check.names = FALSE, stringsAsFactors = FALSE)

numify <- function(x){suppressWarnings(as.numeric(x))}
cols <- c("A549_FDR","Calu3_FDR","HepG2_FDR","U937_FDR")
for (cn in cols) if (cn %in% names(df0)) df0[[cn]] <- numify(df0[[cn]])
if (!"Theme" %in% names(df0)) df0$Theme <- "All"

mk_long <- function(df){
  b1 <- df %>% transmute(Term, Direction, Theme, Lines_sig, Cell="A549",  FDR=A549_FDR)
  b2 <- df %>% transmute(Term, Direction, Theme, Lines_sig, Cell="Calu-3", FDR=Calu3_FDR)
  b3 <- df %>% transmute(Term, Direction, Theme, Lines_sig, Cell="HepG2",  FDR=HepG2_FDR)
  b4 <- df %>% transmute(Term, Direction, Theme, Lines_sig, Cell="U937",   FDR=U937_FDR)
  bind_rows(b1,b2,b3,b4)
}

dd <- mk_long(df0) %>%
  mutate(score = ifelse(is.na(FDR) | FDR>0.05, NA_real_, -log10(pmax(FDR, 1e-300))),
         score = ifelse(Direction=="Up",  score, -score),
         Cell  = factor(Cell, levels=c("A549","Calu-3","HepG2","U937")))

cap <- 10
dd$score <- pmin(pmax(dd$score, -cap), cap)

ord_terms <- df0 %>% arrange(Direction, desc(Lines_sig), if_any(all_of(cols), ~.x)) %>% pull(Term)
dd$Term <- factor(dd$Term, levels=rev(unique(ord_terms)))

p <- ggplot(dd, aes(x=Cell, y=Term, fill=score)) +
  geom_tile(color="white", linewidth=0.2) +
  scale_fill_gradient2(low="#2C7BE5", mid="grey92", high="#D62728",
                       midpoint=0, limits=c(-cap,cap), oob=scales::squish,
                       name="signed -log10(FDR)") +
  labs(x=NULL, y=NULL,
       title="Cross-line ORA summary heatmap",
       subtitle=ifelse(basename(infile)=="crossline_ORA_top10_table.tsv","Top-10 representative terms","Representative terms")) +
  theme_minimal(base_size=11) +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(angle=0, hjust=0.5),
        plot.title=element_text(face="bold"))
if (length(unique(df0$Theme))>1) {
  p <- p + facet_grid(Theme ~ Direction, scales="free_y", space="free_y")
} else {
  p <- p + facet_grid(. ~ Direction, scales="free_y", space="free_y")
}

ggsave(file.path(base,"crossline_ORA_heatmap.pdf"), p, width=8.2, height=6.5)
ggsave(file.path(base,"crossline_ORA_heatmap.png"), p, width=8.2, height=6.5, dpi=300)
cat("done\n")
