# ============================================================
# Volcano Plots for Four GEO Datasets
# GSE278723, GSE173955, GSE48350, GSE5281
# Threshold: adj.P < 0.05, |logFC| > 0.58
# ============================================================

library(ggplot2)
library(ggrepel)
library(readxl)
library(dplyr)
library(patchwork)

# ---- 参数 ----
FC_CUTOFF    <- 0.58
ADJ_P_CUTOFF <- 0.05

# ---- 文件路径 ----
file_path <- "D:/My Desktop/返修/补充文件/Supplementary Table S7.xlsx"
out_dir   <- "D:/My Desktop/返修/补充文件"

# ---- 读取数据 ----
gse278723 <- read_excel(file_path, sheet = "GSE278723")
gse173955 <- read_excel(file_path, sheet = "GSE173955")
gse48350  <- read_excel(file_path, sheet = "GSE48350")
gse5281   <- read_excel(file_path, sheet = "GSE5281")

# ---- 统一列名函数 ----
# gene_col  : 基因名列
# logfc_col : logFC列
# adjp_col  : adjusted p-value列
prepare_data <- function(df, gene_col, logfc_col, adjp_col) {
  df %>%
    rename(gene      = all_of(gene_col),
           logFC     = all_of(logfc_col),
           adj.P.Val = all_of(adjp_col)) %>%
    filter(!is.na(logFC), !is.na(adj.P.Val)) %>%
    filter(!grepl("^AFFX", as.character(gene), ignore.case = TRUE)) %>%
    filter(!is.na(gene) & gene != "" & gene != "NA") %>%
    mutate(
      logFC       = as.numeric(logFC),
      adj.P.Val   = as.numeric(adj.P.Val),
      neg_log10_p = -log10(adj.P.Val),
      status = case_when(
        logFC >  FC_CUTOFF & adj.P.Val < ADJ_P_CUTOFF ~ "Upregulated",
        logFC < -FC_CUTOFF & adj.P.Val < ADJ_P_CUTOFF ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
}

# ---- 准备各数据集（注意每个数据集列名不同）----
df1 <- prepare_data(gse278723,
                    gene_col  = "gene",
                    logfc_col = "logFC",
                    adjp_col  = "adj.P.Val")

df2 <- prepare_data(gse173955,
                    gene_col  = "Symbol",
                    logfc_col = "log2FoldChange",
                    adjp_col  = "padj")

df3 <- prepare_data(gse48350,
                    gene_col  = "Gene.symbol",
                    logfc_col = "logFC",
                    adjp_col  = "adj.P.Val")

df4 <- prepare_data(gse5281,
                    gene_col  = "Gene.symbol",
                    logfc_col = "logFC",
                    adjp_col  = "adj.P.Val")

# ---- 火山图函数 ----
make_volcano <- function(df, title, top_n = 10) {

  colors <- c("Upregulated"     = "#E41A1C",
              "Downregulated"   = "#377EB8",
              "Not significant" = "grey70")

  n_up   <- sum(df$status == "Upregulated",   na.rm = TRUE)
  n_down <- sum(df$status == "Downregulated", na.rm = TRUE)
  n_ns   <- sum(df$status == "Not significant", na.rm = TRUE)

  top_genes <- df %>%
    filter(status != "Not significant") %>%
    arrange(adj.P.Val) %>%
    slice_head(n = top_n)

  y_max <- min(max(df$neg_log10_p, na.rm = TRUE) * 1.1, 300)

  ggplot(df, aes(x = logFC, y = neg_log10_p, color = status)) +
    geom_point(alpha = 0.5, size = 1.0) +
    scale_color_manual(
      values = colors,
      labels = c(
        paste0("Upregulated (n=", n_up, ")"),
        paste0("Downregulated (n=", n_down, ")"),
        paste0("Not significant (n=", n_ns, ")")
      )
    ) +
    geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF),
               linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_hline(yintercept = -log10(ADJ_P_CUTOFF),
               linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_text_repel(
      data          = top_genes,
      aes(label     = gene),
      size          = 2.5,
      max.overlaps  = 20,
      segment.color = "grey50",
      segment.size  = 0.3,
      color         = "black"
    ) +
    coord_cartesian(ylim = c(0, y_max)) +
    labs(
      title = title,
      x     = expression(log[2]~"Fold Change"),
      y     = expression(-log[10]~"(adjusted p-value)"),
      color = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.position  = "bottom",
      legend.text      = element_text(size = 8),
      axis.text        = element_text(color = "black"),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
    )
}

# ---- 生成四张图 ----
p1 <- make_volcano(df1, "GSE278723")
p2 <- make_volcano(df2, "GSE173955")
p3 <- make_volcano(df3, "GSE48350")
p4 <- make_volcano(df4, "GSE5281")

# ---- 拼合 ----
combined <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "Volcano Plots of Differentially Expressed Genes in Four AD Hippocampal Datasets",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
    )
  )

# ---- 保存 ----
ggsave(file.path(out_dir, "volcano_plots_combined.pdf"),
       plot   = combined,
       width  = 12,
       height = 10)

ggsave(file.path(out_dir, "volcano_plots_combined.tiff"),
       plot        = combined,
       width       = 12,
       height      = 10,
       dpi         = 300,
       compression = "lzw")

message("✅ 完成！文件已保存至：", out_dir)
