## ===== 0) 安装与加载 =====
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("oligo", "Biobase", "tibble", "readr")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE)

## HTA-2_0 平台建议安装的包（用于正确识别平台与后续注释）
# 设计包（pd）：pd.hta.2.0；注释包：hta20transcriptcluster.db（差异分析后用）
for (p in c("pd.hta.2.0")) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE)

library(oligo)
library(Biobase)
library(tibble)
library(readr)

## ===== 1) 目录与文件 =====
# 原始 CEL（或 CEL.gz）所在目录（不要使用 setwd）
cel_dir   <- "affymetrix/GSE173954/GSE173954_RAW"
# 可选：你的样本信息（上一条你给的 Targets.txt）
targets_fp <- "affymetrix/GSE173954/Targets.txt"

# 输出目录与文件名
out_dir <- "affymetrix/GSE173954/processed"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rdata_fp <- file.path(out_dir, "GSE173954_after_rma_core.Rdata")
expr_fp  <- file.path(out_dir, "GSE173954_RMA_coreexpr.tsv.gz")

# 替代 list.celfiles()
cel_files <- list.files(
  path = cel_dir,
  pattern = "\\.CEL(\\.gz)?$",   # 匹配 .CEL 或 .CEL.gz
  full.names = TRUE,
  ignore.case = TRUE
)
length(cel_files)
head(cel_files)


## ===== 3) 读取 CEL 并 RMA（HTA-2_0 用 target='core' 得到 gene level）=====
# 读入时会用到 pd.hta.2.0（若没装会提示）
raw <- oligo::read.celfiles(cel_files)

# RMA（gene/transcript cluster level）
# 可选 target：'core'（gene层，常用）或 'probeset'（exon层）
eset <- oligo::rma(raw, target = "core")
expr <- Biobase::exprs(eset)  # 已经是 log2
message("表达矩阵维度：", paste(dim(expr), collapse = " x "))

## ===== 4) （可选）对齐 Targets，并把 pheno 写入 pData =====
if (file.exists(targets_fp)) {
  pheno <- readr::read_tsv(targets_fp, show_col_types = FALSE)
  # 期望至少包含 FileName, sample_id, group（你的表有：FileName、sample_id、tissue_type、group、age、gender）
  # 用文件名对齐（去除路径，仅保留文件名）
  cel_base <- basename(sampleNames(eset))
  # 如果 sampleNames 是解压后的 .CEL，需要同时考虑去掉扩展名：
  cel_base2 <- sub("\\.gz$", "", cel_base, ignore.case = TRUE)
  
  # 与 Targets 的 FileName 对齐（也去掉可能的路径差异）
  tf_base <- basename(pheno$FileName)
  tf_base2 <- sub("\\.CEL$", ".CEL.gz", tf_base, ignore.case = TRUE)  # 容错：有的 Targets 记录 .CEL
  # 构造匹配索引
  idx <- match(cel_base, tf_base)
  if (anyNA(idx)) idx <- match(cel_base2, tf_base)      # 再试一次（去掉.gz）
  if (anyNA(idx)) idx <- match(cel_base, tf_base2)      # 再试一次（把.CEL映射成.CEL.gz）
  if (anyNA(idx)) idx <- match(cel_base2, tf_base2)     # 最后试一次

  if (anyNA(idx)) {
    warning("有样本无法在 Targets 中找到匹配，将跳过写入 pData。")
  } else {
    # 按表达矩阵列顺序重排 targets，并写入 pData
    ph <- as.data.frame(pheno[idx, , drop = FALSE])
    rownames(ph) <- sampleNames(eset)
    Biobase::pData(eset) <- ph
  }
}

## ===== 5) 保存对象与矩阵 =====
save(eset, expr, cel_files, file = rdata_fp)
# 也可以导出表达矩阵（制表符分隔 + gzip）
# 带行名写出（探针/cluster ID）
write.table(expr, file = gzfile(expr_fp), sep = "\t", quote = FALSE, col.names = NA)

message("已保存：\n  R对象: ", rdata_fp, "\n  矩阵:  ", expr_fp)
