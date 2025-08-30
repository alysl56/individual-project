#!/usr/bin/env Rscript
.libPaths('~/Rlibs')
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(cowplot) })

expr_csv <- "/gpfs01/home/alysl56/projects/DE_analysis/venn_table/six_genes_expression_long_table.csv"
core6_tsv <- "/gpfs01/home/alysl56/projects/DE_analysis/venn_table/core6_selected.tsv"
out_dir   <- "/gpfs01/home/alysl56/projects/DE_analysis/gene_plots_bar"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

x <- utils::read.csv(expr_csv, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("symbol","cell_line","condition","norm_count") %in% names(x)))
core <- utils::read.delim(core6_tsv, check.names = FALSE, stringsAsFactors = FALSE)
gene_col <- if ("SYMBOL" %in% names(core)) "SYMBOL" else names(core)[1]
genes <- unique(core[[gene_col]])

map_col <- function(prefix, nm){ c1 <- paste0(prefix, "_", nm); c2 <- paste0(sub("Calu3","Calu-3", prefix), "_", nm); if(c1 %in% names(core)) c1 else if(c2 %in% names(core)) c2 else NA_character_ }
lfc_cols  <- c(A549 = map_col("A549","log2FoldChange"), `Calu-3` = map_col("Calu3","log2FoldChange"), HepG2 = map_col("HepG2","log2FoldChange"))
padj_cols <- c(A549 = map_col("A549","padj"),            `Calu-3` = map_col("Calu3","padj"),            HepG2 = map_col("HepG2","padj"))

mk_lab <- function(line){ data.frame(SYMBOL = core[[gene_col]], cell_line = line, lfc = as.numeric(core[[lfc_cols[[line]]]]), padj = as.numeric(core[[padj_cols[[line]]]]), stringsAsFactors = FALSE) }
lab_long <- do.call(rbind, lapply(names(lfc_cols), mk_lab))
lab_long$label <- paste0("log2FC=", sprintf("%.2f", lab_long$lfc), "\npadj=", format(lab_long$padj, scientific = TRUE, digits = 2))

canon <- function(v){
  vv <- tolower(trimws(v))
  ifelse(grepl("\\b(without|no[_-]?dmso|ctrl|control|vehicle|untreat|untreated|minus|neg|no)\\b", vv), "−",
  ifelse(grepl("\\b(with|with[_-]?dmso|dmso|treated|treat|plus|pos)\\b", vv), "+", v))
}

dat <- x %>%
  mutate(condition = canon(condition),
         cell_line = factor(cell_line, levels = c("A549","HepG2","Calu-3","U937"))) %>%
  filter(symbol %in% genes)

sumdat <- dat %>%
  group_by(symbol, cell_line, condition) %>%
  summarise(mean = mean(as.numeric(norm_count), na.rm = TRUE),
            sd   = sd(as.numeric(norm_count),   na.rm = TRUE),
            n    = dplyr::n(),
            se   = sd / pmax(n,1),
            .groups = "drop")

plot_one <- function(g){
  df <- sumdat %>% filter(symbol == g); if(nrow(df)==0) return(NULL)
  ann <- lab_long %>% filter(SYMBOL == g)
  pad <- df %>% group_by(cell_line) %>% summarise(ypos = max(mean + se, na.rm = TRUE) * 1.12, .groups = "drop") %>% left_join(ann, by = "cell_line")
  pd <- position_dodge(width = 0.6)
  p <- ggplot(df, aes(x = cell_line, y = mean, fill = condition)) +
    geom_col(position = pd, width = 0.6) +
    geom_errorbar(aes(ymin = pmax(mean - se, 0), ymax = mean + se), position = pd, width = 0.25, linewidth = 0.4) +
    scale_fill_manual(values = c("+" = "#2C7BE5", "−" = "#E55353"), name = "Condition", labels = c("+ with DMSO","− without DMSO")) +
    geom_text(data = pad, aes(x = cell_line, y = ypos, label = label), inherit.aes = FALSE, size = 3.2, vjust = 0, lineheight = 1.0) +
    labs(title = paste0("Expression of ", g), x = "Cell line", y = "Normalized counts (DESeq2)") +
    coord_cartesian(clip = "off") + theme_minimal(base_size = 12) +
    theme(legend.position = "right", plot.title = element_text(face = "bold"), plot.margin = margin(10, 20, 10, 10))
  ggsave(file.path(out_dir, paste0(g, "_barplot.png")), p, width = 7.5, height = 4.6, dpi = 300)
  ggsave(file.path(out_dir, paste0(g, "_barplot.pdf")),  p, width = 7.5, height = 4.6)
  p
}

plots <- lapply(genes, plot_one); plots <- plots[!vapply(plots, is.null, logical(1))]
if(length(plots)>0){ comb <- cowplot::plot_grid(plotlist = plots, ncol = 2, align = "hv")
  ggsave(file.path(out_dir, "Core6_barplots_combined.png"), comb, width = 10, height = 8, dpi = 300)
  ggsave(file.path(out_dir, "Core6_barplots_combined.pdf"), comb, width = 10, height = 8) }
utils::write.csv(sumdat, file.path(out_dir, "Core6_expression_summary.csv"), row.names = FALSE)

