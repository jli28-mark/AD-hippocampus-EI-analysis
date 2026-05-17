## ---- 0. 环境 ----
library(edgeR)
library(limma)
library(tidyverse)

data_dir <- "D:/My Desktop/5.19AD/raw data/2.GSE278723"
counts_f <- file.path(data_dir, "GSE278723_Counts_filtered.txt")
meta_f   <- file.path(data_dir, "GSE278723_sample_information.txt")
out_csv  <- file.path(data_dir, "DEG_GSE278723_AD_vs_NC_blockBySubject_regionAdj.csv")

## ---- 1. 读取数据 ----
# 计数矩阵：第一列基因，其他列为样本
counts <- read.delim(counts_f, check.names = FALSE)
gene_col <- 1
genes <- counts[[gene_col]]
mat <- counts[ , -(gene_col), drop = FALSE]
rownames(mat) <- make.unique(genes)

# 元数据
meta <- read.delim(meta_f, sep = "\t", check.names = FALSE)
colnames(meta) <- c("SampleID","Group")  # 保证列名一致
meta$SampleID <- as.character(meta$SampleID)
meta$Group    <- factor(meta$Group, levels = c("NC","AD"))

# 从列名解析 subject 与 region
# 列名形如 "01.CA1"；subject=01，region=CA1
parse_id <- function(x){
  tibble(
    SampleID = x,
    subject  = sub("\\..*$", "", x),
    region   = sub("^.*\\.", "", x)
  )
}
design_df <- parse_id(colnames(mat)) %>%
  left_join(meta, by = "SampleID")

# 确保列顺序与设计一致
stopifnot(all(design_df$SampleID == colnames(mat)))

design_df$subject <- factor(design_df$subject)    # 01–28
design_df$region  <- factor(design_df$region)     # CA1–CA4

## ---- 2. 低表达过滤 + 归一化 ----
y <- DGEList(counts = as.matrix(mat))
# 以主对比 Group 作为分组参考来做过滤
keep <- filterByExpr(y, group = design_df$Group)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")

## ---- 3. 设计矩阵（主效应：Group；协变量：region；重复测量：subject）----
design <- model.matrix(~ Group + region, data = design_df)

# voom + 重复测量相关性估计
v  <- voom(y, design, plot = FALSE)
dupcor1 <- duplicateCorrelation(v, design, block = design_df$subject)
v  <- voom(y, design, plot = FALSE, block = design_df$subject, correlation = dupcor1$consensus.correlation)
fit <- lmFit(v, design, block = design_df$subject, correlation = dupcor1$consensus.correlation)
fit <- eBayes(fit)

## ---- 4. 提取 AD vs NC 主效应 ----
# 设计中 GroupAD 系数即 AD- NC
res <- topTable(fit, coef = "GroupAD", number = Inf, sort.by = "P")
res$gene <- rownames(res)
res <- res %>% relocate(gene)
write.csv(res, out_csv, row.names = FALSE)

## 可选：阈值筛选
deg <- res %>% filter(adj.P.Val < 0.05, abs(logFC) > 0.58)
write.csv(deg, sub(".csv$", "_sig.csv", out_csv), row.names = FALSE)

## ---- 5.（可选）每个亚区单独分析作为敏感性 ----
do_per_region <- FALSE
if (do_per_region) {
  reg_levels <- levels(design_df$region)
  for (rg in reg_levels) {
    idx <- design_df$region == rg
    y2 <- DGEList(counts = as.matrix(mat[, idx, drop=FALSE]))
    keep2 <- filterByExpr(y2, group = droplevels(design_df$Group[idx]))
    y2 <- y2[keep2,, keep.lib.sizes = FALSE]
    y2 <- calcNormFactors(y2)

    des2 <- model.matrix(~ Group, data = droplevels(design_df[idx,]))
    v2 <- voom(y2, des2, plot = FALSE)
    fit2 <- eBayes(lmFit(v2, des2))

    res2 <- topTable(fit2, coef = "GroupAD", number = Inf)
    res2$gene <- rownames(res2)
    write.csv(res2 %>% relocate(gene),
              file.path(data_dir, paste0("DEG_GSE278723_", rg, "_AD_vs_NC.csv")),
              row.names = FALSE)
  }
}
