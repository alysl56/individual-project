
library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(cowplot))

base <- "/gpfs01/home/alysl56/projects/DMSO_U937_RNAseq"
f <- file.path(base, "DE_analysis", "U937_DESeq2_results.csv")

res <- read.csv(f, row.names = 1, check.names = FALSE) %>% tibble::rownames_to_column("gene")
res$gene_clean <- str_replace(res$gene, "\\.\\d+$", "")
ann <- AnnotationDbi::select(org.Hs.eg.db, keys = res$gene_clean, keytype = "ENSEMBL", columns = c("SYMBOL"))
ann <- ann %>% distinct(ENSEMBL, .keep_all = TRUE)
res <- res %>% left_join(ann, by = c("gene_clean" = "ENSEMBL"))
res$SYMBOL <- ifelse(is.na(res$SYMBOL) | res$SYMBOL == "", res$gene_clean, res$SYMBOL)
df <- res %>% filter(!is.na(padj) & !is.na(log2FoldChange))

thr_fc <- 1
thr_p  <- 0.05
df$neglog10padj <- -log10(df$padj)
df$status <- "NS"
df$status[df$padj < thr_p & df$log2FoldChange >  thr_fc] <- "UP"
df$status[df$padj < thr_p & df$log2FoldChange < -thr_fc] <- "DOWN"
df$status <- factor(df$status, levels = c("DOWN","NS","UP"))

n_up <- sum(df$status == "UP")
n_dn <- sum(df$status == "DOWN")
n_ns <- sum(df$status == "NS")

top_up <- df %>% filter(status == "UP")   %>% arrange(padj, desc(abs(log2FoldChange))) %>% slice_head(n = 10)
top_dn <- df %>% filter(status == "DOWN") %>% arrange(padj, desc(abs(log2FoldChange))) %>% slice_head(n = 10)
labs_df <- bind_rows(top_up, top_dn) %>% distinct(gene, .keep_all = TRUE)

x_min <- min(df$log2FoldChange, na.rm = TRUE)
x_max <- max(df$log2FoldChange, na.rm = TRUE)
x_rng <- x_max - x_min
x_right <- x_max + 0.33 * x_rng

y_top <- max(df$neglog10padj, na.rm = TRUE)
lab_y_max <- if (nrow(labs_df) > 0) max(labs_df$neglog10padj, na.rm = TRUE) else y_top
y0 <- max(0.20 * y_top, lab_y_max - 1.2)
dy <- 0.55

size_header <- 3.2
size_text   <- 3.0
size_bullet <- 3.2

p <- ggplot(df, aes(x = log2FoldChange, y = neglog10padj, color = status)) +
  geom_hline(yintercept = -log10(thr_p), linetype = "dashed", linewidth = 0.3, color = "grey55") +
  geom_vline(xintercept = c(-thr_fc, thr_fc), linetype = "dashed", linewidth = 0.3, color = "grey55") +
  geom_point(size = 1.2, alpha = 0.7, stroke = 0) +
  scale_color_manual(values = c(DOWN = "#2ca02c", NS = "grey75", UP = "#d62728"), name = NULL) +
  geom_text_repel(
    data = labs_df,
    aes(label = SYMBOL),
    size = 3.0, max.overlaps = 100, min.segment.length = 0,
    box.padding = 0.35, point.padding = 0.25, seed = 42,
    segment.color = "grey60",
    xlim = c(x_min, x_max)
  ) +
  scale_x_continuous(limits = c(x_min, x_right)) +
  labs(
    x = expression(log[2]*"Fold Change"),
    y = expression(-log[10]*"(adjusted p-value)"),
    title    = "Volcano plot",
    subtitle = "U937: with DMSO vs without DMSO"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey25"),
    plot.margin = margin(t = 6, r = 22, b = 6, l = 6)
  )

x_mark <- x_max + 0.02 * x_rng
x_text <- x_max + 0.06 * x_rng

p <- p +
  annotate("text", x = x_text, y = y0 + 2*dy, label = "padj < 0.05",  hjust = 0, vjust = 0.5, size = size_header) +
  annotate("text", x = x_text, y = y0 + 1*dy, label = "|log2FC| > 1", hjust = 0, vjust = 0.5, size = size_header) +
  annotate("text", x = x_mark, y = y0 - 0*dy, label = "a", hjust = 0, vjust = 0.5, size = size_bullet, color = "#2ca02c") +
  annotate("text", x = x_text, y = y0 - 0*dy, label = paste0("DOWN ", n_dn), hjust = 0, vjust = 0.5, size = size_text) +
  annotate("text", x = x_mark, y = y0 - 1*dy, label = " ^` ", hjust = 0, vjust = 0.5, size = size_bullet, color = "grey70") +
  annotate("text", x = x_text, y = y0 - 1*dy, label = paste0("NS ", n_ns),   hjust = 0, vjust = 0.5, size = size_text) +
  annotate("text", x = x_mark, y = y0 - 2*dy, label = "a", hjust = 0, vjust = 0.5, size = size_bullet, color = "#d62728") +
  annotate("text", x = x_text, y = y0 - 2*dy, label = paste0("UP ", n_up),   hjust = 0, vjust = 0.5, size = size_text)

dir.create(file.path(base, "DE_analysis"), showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(base, "DE_analysis", "U937_volcano.pdf"), p, width = 6.6, height = 6.2)
ggsave(file.path(base, "DE_analysis", "U937_volcano.png"), p, width = 6.6, height = 6.2, dpi = 300)

write.csv(dplyr::select(labs_df, gene, SYMBOL, log2FoldChange, padj),
          file.path(base, "DE_analysis", "U937_volcano_labels.csv"),
          row.names = FALSE)
