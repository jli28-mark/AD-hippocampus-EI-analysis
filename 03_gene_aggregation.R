library(AnnotationDbi)
library(hta20transcriptcluster.db)
library(limma)
library(dplyr)
library(tibble)

# 1) 取表达矩阵 & 对应的探针ID
expr <- exprs(eset)                     # 70523 x 18，行为 PROBEID
probe_ids <- rownames(expr)

# 2) 映射：PROBEID -> SYMBOL
anno <- AnnotationDbi::select(
  hta20transcriptcluster.db,
  keys    = probe_ids,
  keytype = "PROBEID",
  columns = c("SYMBOL")
) %>% 
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  distinct(PROBEID, SYMBOL, .keep_all = TRUE)

# 3) 用映射表对齐表达矩阵的行顺序
expr2 <- expr[anno$PROBEID, , drop = FALSE]
stopifnot(identical(rownames(expr2), anno$PROBEID))

# 4) 聚合到基因层（均值；可改 median）
expr_gene <- limma::avereps(expr2, ID = anno$SYMBOL, FUN = mean)

# 5) 保存
dir.create(file.path(out_dir, "genelevel"), showWarnings = FALSE)
save(expr_gene, file = file.path(out_dir, "genelevel/GSE173954_expr_geneMean.Rdata"))
write.table(expr_gene, file = file.path(out_dir, "genelevel/GSE173954_expr_geneMean.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
dim(expr_gene); expr_gene[1:3,1:3]
