# ============================================================
# Venn Diagram of DEG Overlap Across Four AD Hippocampal Datasets
# GSE278723, GSE173955, GSE48350, GSE5281
# Threshold: adj.P < 0.05, |logFC| > 0.58
# ============================================================

library(readxl)
library(dplyr)
library(ggVennDiagram)
library(ggplot2)

# ---- 文件路径 ----
file_path <- "D:/My Desktop/返修/补充文件/Supplementary Table S7.xlsx"
out_dir   <- "D:/My Desktop/返修/补充文件"

# ---- 读取数据 ----
gse278723 <- read_excel(file_path, sheet = "GSE278723")
gse173955 <- read_excel(file_path, sheet = "GSE173955")
gse48350  <- read_excel(file_path, sheet = "GSE48350")
gse5281   <- read_excel(file_path, sheet = "GSE5281")

# ---- 筛选显著DEG并提取基因名 ----
FC_CUTOFF    <- 0.58
ADJ_P_CUTOFF <- 0.05

get_sig_genes <- function(df, gene_col, logfc_col, adjp_col) {
  df %>%
    rename(gene      = all_of(gene_col),
           logFC     = all_of(logfc_col),
           adj.P.Val = all_of(adjp_col)) %>%
    mutate(logFC     = as.numeric(logFC),
           adj.P.Val = as.numeric(adj.P.Val)) %>%
    filter(!is.na(logFC), !is.na(adj.P.Val)) %>%
    filter(!grepl("^AFFX", as.character(gene), ignore.case = TRUE)) %>%
    filter(!is.na(gene) & gene != "" & gene != "NA") %>%
    filter(abs(logFC) > FC_CUTOFF & adj.P.Val < ADJ_P_CUTOFF) %>%
    pull(gene) %>%
    unique()
}

genes_278723 <- get_sig_genes(gse278723, "gene",        "logFC",         "adj.P.Val")
genes_173955 <- get_sig_genes(gse173955, "Symbol",      "log2FoldChange","padj")
genes_48350  <- get_sig_genes(gse48350,  "Gene.symbol", "logFC",         "adj.P.Val")
genes_5281   <- get_sig_genes(gse5281,   "Gene.symbol", "logFC",         "adj.P.Val")

# ---- 打印各数据集DEG数量 ----
cat("GSE278723 significant DEGs:", length(genes_278723), "\n")
cat("GSE173955 significant DEGs:", length(genes_173955), "\n")
cat("GSE48350  significant DEGs:", length(genes_48350),  "\n")
cat("GSE5281   significant DEGs:", length(genes_5281),   "\n")

# ---- 构建列表 ----
gene_list <- list(
  GSE278723 = genes_278723,
  GSE173955 = genes_173955,
  GSE48350  = genes_48350,
  GSE5281   = genes_5281
)

# ---- 绘制Venn图 ----
p <- ggVennDiagram(
  gene_list,
  label_alpha = 0,
  edge_size   = 0.5
) +
  scale_fill_gradient(low = "#F0F8FF", high = "#2166AC") +
  scale_color_manual(values = c(
    "GSE278723" = "#E41A1C",
    "GSE173955" = "#377EB8",
    "GSE48350"  = "#4DAF4A",
    "GSE5281"   = "#FF7F00"
  )) +
  labs(
    title = "Overlap of Significant DEGs Across Four AD Hippocampal Datasets",
    fill  = "Gene Count"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "right"
  )

# ---- 保存 ----
ggsave(file.path(out_dir, "venn_diagram_DEG_overlap.pdf"),
       plot   = p,
       width  = 8,
       height = 7)

ggsave(file.path(out_dir, "venn_diagram_DEG_overlap.tiff"),
       plot        = p,
       width       = 8,
       height      = 7,
       dpi         = 300,
       compression = "lzw")

message("✅ Venn图已保存至：", out_dir)