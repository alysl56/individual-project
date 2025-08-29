library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(cowplot))

base <- "/gpfs01/home/alysl56/projects/DMSO_U937_RNAseq"
f <- file.path(base,"DE_analysis","U937_DESeq2_results.csv")

res <- read.csv(f, row.names = 1, check.names = FALSE) %>% tibble::rownames_to_column("gene")
res$gene_clean <- str_replace(res$gene, "\\.\\d+$", "")
ann <- AnnotationDbi::select(org.Hs.eg.db, keys = res$gene_clean, keytype = "ENSEMBL", columns = c("SYMBOL"))
ann <- ann %>% distinct(ENSEMBL, .keep_all = TRUE)
res <- res %>% left_join(ann, by = c("gene_clean" = "ENSEMBL"))
res$SYMBOL <- ifelse(is.na(res$SYMBOL) | res$SYMBOL=="", res$gene_clean, res$SYMBOL)
df <- res %>% filter(!is.na(padj) & !is.na(log2FoldChange))

thr_fc <- 1
thr_p  <- 0.05
df$neglog10padj <- -log10(df$padj)
df$status <- "NS"
df$status[df$padj < thr_p & df$log2FoldChange >  thr_fc] <- "UP"
df$status[df$padj < thr_p & df$log2FoldChange < -thr_fc] <- "DOWN"
df$status <- factor(df$status, levels = c("DOWN","NS","UP"))

n_up <- sum(df$status=="UP")
n_dn <- sum(df$status=="DOWN")
n_ns <- sum(df$status=="NS")

top_up <- df %>% filter(status=="UP")   %>% arrange(padj, desc(abs(log2FoldChange))) %>% slice_head(n=10)
top_dn <- df %>% filter(status=="DOWN") %>% arrange(padj, desc(abs(log2FoldChange))) %>% slice_head(n=10)
labs_df <- bind_rows(top_up, top_dn) %>% distinct(gene, .keep_all = TRUE)

p_core <- ggplot(df, aes(x = log2FoldChange, y = neglog10padj, color = status)) +
  geom_hline(yintercept = -log10(thr_p), linetype = "dashed", linewidth = 0.3, color = "grey55") +
  geom_vline(xintercept = c(-thr_fc, thr_fc), linetype = "dashed", linewidth = 0.3, color = "grey55") +
  geom_point(size = 1.2, alpha = 0.7, stroke = 0) +
  scale_color_manual(
    values = c(DOWN = "#2ca02c", NS = "grey75", UP = "#d62728"),
    breaks = c("DOWN","NS","UP"),
    labels = c(paste0("DOWN ", n_dn), paste0("NS ", n_ns), paste0("UP ", n_up)),
    drop   = FALSE,
    name   = NULL
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  geom_text_repel(
    data = labs_df,
    aes(label = SYMBOL),
    size = 3.2, max.overlaps = 100, min.segment.length = 0,
    box.padding = 0.35, point.padding = 0.25, seed = 42,
    segment.color = "grey60"
  ) +
  labs(
    x = expression(log[2]*"Fold Change"),
    y = expression(-log[10]*"(adjusted p-value)"),
    title    = "Volcano plot",
    subtitle = "U937: with DMSO vs without DMSO"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position      = c(0.18, 0.80),
    legend.justification = c(0,0),
    legend.background    = element_rect(fill = alpha("white",0.75), color = NA),
    plot.title           = element_text(face="bold"),
    plot.subtitle        = element_text(color="grey25"),
    plot.margin          = margin(t=6,r=20,b=6,l=6)
  ) +
  coord_cartesian(clip = "off")

p_final <- cowplot::ggdraw(p_core) +
  cowplot::draw_label("padj < 0.05\n|log2FC| > 1",
                      x = 0.33, y = 0.89, hjust = 0, vjust = 0,
                      fontface = "plain", size = 10, color = "black",
                      lineheight = 1.05)

dir.create(file.path(base,"DE_analysis"), showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(base,"DE_analysis","U937_volcano.pdf"), p_final, width = 6.6, height = 6.2)
ggsave(file.path(base,"DE_analysis","U937_volcano.png"), p_final, width = 6.6, height = 6.2, dpi = 300)

write.csv(dplyr::select(labs_df, gene, SYMBOL, log2FoldChange, padj),
          file.path(base,"DE_analysis","U937_volcano_labels.csv"), row.names = FALSE)
