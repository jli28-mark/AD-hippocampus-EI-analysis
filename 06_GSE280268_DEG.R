## =============================================================
## GSE280268 Differential Expression Analysis
## Platform: GPL21697 NextSeq 550 (Homo sapiens)
## Comparison: AD vs HealthyControl (No dementia)
## Brain region: Hippocampus (HI) only
## Data: Pre-normalized library counts (All_libnorml_AD_HI.counts.txt)
## Note: This script was reconstructed based on the analysis
##       described in the STAR Methods section.
##       The original script was not retained.
## =============================================================

## ---- 0. 环境 ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("limma", "edgeR", "org.Hs.eg.db", "AnnotationDbi")) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)
}
library(limma)
library(edgeR)
library(dplyr)
library(tibble)
library(AnnotationDbi)
library(org.Hs.eg.db)

## ---- 1. 路径设置 ----
# 请将 data_dir 改为你本地存放 GSE280268 数据的路径
data_dir <- "path/to/GSE280268"
counts_f <- file.path(data_dir, "GSE280268_All_libnorml_AD_HI.counts.txt")
out_csv  <- file.path(data_dir, "DEG_AD_vs_HC_HI_all.csv")

## ---- 2. 读取数据 ----
counts <- read.delim(counts_f, check.names = FALSE, row.names = 1)

# AD samples: AD1H, AD2H, AD3H, AD4H, AD5H, AD6H
# HC samples: HC1H1, HC2H1, HC3H1, HC1H2, HC2H2, HC3H2
ad_cols <- grep("^AD[0-9]+H$", colnames(counts), value = TRUE)
hc_cols <- grep("^HC[0-9]+H[0-9]+$", colnames(counts), value = TRUE)

cat("AD样本数：", length(ad_cols), "\n")
cat("HC样本数：", length(hc_cols), "\n")

counts_hi <- counts[, c(ad_cols, hc_cols), drop = FALSE]

## ---- 3. 构建样本信息 ----
meta <- data.frame(
  SampleID = c(ad_cols, hc_cols),
  Group    = factor(c(rep("AD", length(ad_cols)),
                      rep("HC", length(hc_cols))),
                    levels = c("HC", "AD")),
  stringsAsFactors = FALSE
)

## ---- 4. 设计矩阵 ----
design <- model.matrix(~ Group, data = meta)
colnames(design) <- make.names(colnames(design))

## ---- 5. limma差异分析 ----
# 数据已library normalized，取log2后直接用limma
expr_log2 <- log2(counts_hi + 1)
fit  <- lmFit(as.matrix(expr_log2), design)
fit  <- eBayes(fit, trend = TRUE, robust = TRUE)

## ---- 6. 提取结果 ----
res <- topTable(fit, coef = "GroupAD", number = Inf, sort.by = "P")
res <- res %>% rownames_to_column("ENSEMBL")

## ---- 7. 基因注释（ENSEMBL -> SYMBOL）----
res$ENSEMBL_clean <- sub("\\..*$", "", res$ENSEMBL)

anno <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(res$ENSEMBL_clean),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID", "GENENAME")
) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

res <- res %>%
  left_join(anno, by = c("ENSEMBL_clean" = "ENSEMBL"))

## ---- 8. 导出全部结果 ----
write.csv(res, out_csv, row.names = FALSE)
cat("已保存全部结果：", nrow(res), "个基因\n")

## ---- 9. 筛选显著DEGs ----
deg_sig <- res %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 0.58)

write.csv(deg_sig,
          file.path(data_dir, "DEG_AD_vs_HC_HI_significant.csv"),
          row.names = FALSE)

cat("显著DEG数（adj.P<0.05, |logFC|>0.58）：", nrow(deg_sig), "\n")
