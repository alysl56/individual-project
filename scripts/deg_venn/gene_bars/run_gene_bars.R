#!/usr/bin/env Rscript
.libPaths('~/Rlibs')
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(cowplot) })

expr_csv <- "/gpfs01/home/alysl56/projects/DE_analysis/venn_table/six_genes_expression_long_table.csv"
core6_tsv <- "/gpfs01/home/alysl56/projects/DE_analysis/venn_table/core6_selected.tsv"
out_dir   <- "/gpfs01/home/alysl56/projects/DE_analysis/gene_plots_bar"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

x <- utils::read.csv(expr_csv, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("symbol","cell_line","condition","norm_count") %in% names(x)))
core <- utils::read.delim(core6_tsv, check.names = FALSE, stringsAsFactors = FALSE)
gene_col <- if ("SYMBOL" %in% names(core)) "SYMBOL" else names(core)[1]
genes <- unique(core[[gene_col]])

map_col <- function(prefix, nm){
  c1 <- paste0(prefix, "_", nm)
  c2 <- paste0(sub("Calu3","Calu-3", prefix), "_", nm)
  if (c1 %in% names(core)) c1 else if (c2 %in% names(core)) c2 else NA_character_
}
lfc_cols  <- c(A549 = map_col("A549","log2FoldChange"),
               `Calu-3` = map_col("Calu3","log2FoldChange"),
               HepG2  = map_col("HepG2","log2FoldChange"))
padj_cols <- c(A549 = map_col("A549","padj"),
               `Calu-3` = map_col("Calu3","padj"),
               HepG2  = map_col("HepG2","padj"))

mk_lab <- function(line){
  data.frame(
    SYMBOL = core[[gene_col]],
    cell_line = line,
    lfc  = suppressWarnings(as.numeric(core[[lfc_cols[[line]]]])),
    padj = suppressWarnings(as.numeric(core[[padj_cols[[line]]]])),
    stringsAsFactors = FALSE
  )
}
lab_long <- do.call(rbind, lapply(names(lfc_cols), mk_lab))
lab_long$lab_lfc <- paste0("log2FC=", sprintf("%.2f", lab_long$lfc))
lab_long$lab_p   <- paste0("padj=",   format(lab_long$padj, scientific = TRUE, digits = 2))

canon <- function(v){
  vv <- tolower(trimws(v))
  if (vv %in% c("-", "−", "minus", "neg", "0", "false")) return("minus")
  if (vv %in% c("+", "plus", "pos", "1", "true"))         return("plus")
  if (grepl("\\b(without|no[_-]?dmso|\\bctrl\\b|control|vehicle|untreat|untreated|minus|\\bneg\\b|mock|blank|naive)\\b", vv)) return("minus")
  if (grepl("\\b(with[_-]?dmso|\\bdmso\\b|treated|treat|plus|\\bpos\\b)\\b", vv)) return("plus")
  NA_character_
}

dat <- x %>%
  mutate(cond = vapply(condition, canon, character(1))) %>%
  filter(cell_line %in% c("A549","Calu-3","HepG2"), !is.na(cond)) %>%
  mutate(cond = factor(cond, levels = c("plus","minus")),
         cell_line = factor(cell_line, levels = c("A549","Calu-3","HepG2"))) %>%
  filter(symbol %in% genes)

sumdat <- dat %>%
  group_by(symbol, cell_line, cond) %>%
  summarise(mean = mean(as.numeric(norm_count), na.rm = TRUE),
            sd   = sd(as.numeric(norm_count),   na.rm = TRUE),
            n    = dplyr::n(),
            se   = sd / pmax(n,1),
            .groups = "drop")

SHIFT_LFC <- 1.14
SHIFT_PADJ <- 1.26

plot_one <- function(g){
  df <- sumdat %>% filter(symbol == g); if(nrow(df) == 0) return(NULL)
  ann <- lab_long %>% filter(SYMBOL == g)
  y_top <- max(df$mean + df$se, na.rm = TRUE)
  pads <- data.frame(cell_line = levels(df$cell_line)) %>%
    left_join(ann, by = "cell_line") %>%
    mutate(y_lfc = y_top * SHIFT_LFC, y_p = y_top * SHIFT_PADJ)
  pd <- position_dodge(width = 0.6)
  p <- ggplot(df, aes(x = cell_line, y = mean, fill = cond)) +
    geom_col(position = pd, width = 0.6) +
    geom_errorbar(aes(ymin = pmax(mean - se, 0), ymax = mean + se),
                  position = pd, width = 0.25, linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.0f", mean), y = mean),
              position = pd, vjust = -0.9, size = 3.0, color = "black") +
    scale_fill_manual(values = c(plus = "#2C7BE5", minus = "#6c757d"),
                      breaks = c("plus","minus"),
                      limits = c("plus","minus"),
                      labels = c("+ with DMSO","- without DMSO"),
                      drop = FALSE) +
    guides(fill = guide_legend(override.aes = list(alpha = 1))) +
    geom_text(data = pads, aes(x = cell_line, y = y_lfc, label = lab_lfc),
              inherit.aes = FALSE, size = 3.3, vjust = 0, lineheight = 1.0) +
    geom_text(data = pads, aes(x = cell_line, y = y_p, label = lab_p),
              inherit.aes = FALSE, size = 3.3, vjust = 0, lineheight = 1.0) +
    labs(title = paste0("Expression of ", g),
         x = "Cell line",
         y = "Normalized counts (DESeq2)",
         fill = "Condition") +
    coord_cartesian(ylim = c(0, y_top * 1.34), clip = "off") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold"),
          plot.margin = margin(10, 34, 10, 10),
          legend.key = element_rect(fill = NA, colour = NA))
  ggsave(file.path(out_dir, paste0(g, "_barplot_3lines.png")), p, width = 8.0, height = 5.2, dpi = 300)
  ggsave(file.path(out_dir, paste0(g, "_barplot_3lines.pdf")),  p, width = 8.0, height = 5.2)
  p
}

plots <- lapply(genes, plot_one)
plots <- plots[!vapply(plots, is.null, logical(1))]
if(length(plots) > 0){
  comb <- cowplot::plot_grid(plotlist = plots, ncol = 2, align = "hv")
  ggsave(file.path(out_dir, "Core6_barplots_3lines.png"), comb, width = 10.8, height = 8.6, dpi = 300)
  ggsave(file.path(out_dir, "Core6_barplots_3lines.pdf"), comb, width = 10.8, height = 8.6)
}
