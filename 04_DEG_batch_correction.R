suppressPackageStartupMessages({
  library(Biobase); library(limma); library(sva)
  library(AnnotationDbi); library(hta20transcriptcluster.db)
  library(dplyr); library(tibble)
})

## 0) 拿到表达与表型
stopifnot(exists("eset"))
expr  <- exprs(eset)          # RMA 输出（log2）
pheno <- pData(eset)

## 1) 统一命名：都用纯 GSM ID（去掉 .CEL / .CEL.gz）
colnames(expr)    <- sub("\\.CEL(\\.gz)?$", "", colnames(expr),    ignore.case = TRUE)
rownames(pheno)   <- sub("\\.CEL(\\.gz)?$", "", rownames(pheno),   ignore.case = TRUE)

## 2) 如果 Targets 里有 sample_id（推荐用它作为标准名），也顺便覆盖一下
if ("sample_id" %in% colnames(pheno)) {
  rownames(pheno) <- pheno$sample_id
}
# 确保表达矩阵列顺序与 pheno 行顺序一致
expr <- expr[, rownames(pheno), drop = FALSE]
stopifnot(identical(colnames(expr), rownames(pheno)))

## 3) 分组（仅 AD vs Control；不纳入性别/年龄）
pheno$group <- factor(pheno$group, levels = c("Control","AD"))
stopifnot(!any(is.na(pheno$group)))

## 4) 批次：从 protocolData 推断（YYYY-MM），若只有一个批次将自动跳过 ComBat
pd <- protocolData(eset)@data
if ("dates" %in% names(pd)) {
  pheno$batch <- factor(substr(pd$dates, 1, 7))
} else if ("ScanDate" %in% names(pd)) {
  pheno$batch <- factor(substr(pd$ScanDate, 1, 7))
} else {
  pheno$batch <- factor(rep("B1", ncol(expr)))
}
cat("批次频数：\n"); print(table(pheno$batch))

table(pheno$batch, pheno$group)


## 5) 去批次（仅保留组效应）
if (nlevels(pheno$batch) > 1) {
  mod      <- model.matrix(~ pheno$group)   # 只含 group
  expr_bc  <- ComBat(dat = expr, batch = pheno$batch, mod = mod,
                     par.prior = TRUE, prior.plots = FALSE)
} else {
  message("只检测到一个批次：跳过 ComBat。")
  expr_bc <- expr
}

## 6) 可选过滤：保留变异度前 50% 提升信噪比
sdv    <- apply(expr_bc, 1, sd)
expr_f <- expr_bc[sdv >= quantile(sdv, 0.5, na.rm = TRUE), , drop = FALSE]

## 7) 仅 AD vs Control 的 limma 差异分析
design <- model.matrix(~ 0 + pheno$group)
colnames(design) <- c("Control","AD")
fit  <- lmFit(expr_f, design)
cont <- makeContrasts(ADvsCTRL = AD - Control, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cont), trend = TRUE, robust = TRUE)

tt <- topTable(fit2, coef = "ADvsCTRL", number = Inf, sort.by = "P")  # 探针层

if (!exists("tt_anno")) {
  library(AnnotationDbi)
  library(hta20transcriptcluster.db)
  tt2 <- tt %>% tibble::rownames_to_column("PROBEID")
  anno <- AnnotationDbi::select(
    hta20transcriptcluster.db,
    keys     = tt2$PROBEID,
    keytype  = "PROBEID",
    columns  = c("SYMBOL","GENENAME","ENTREZID","ENSEMBL")
  )
  tt_anno <- tt2 %>% dplyr::left_join(anno, by="PROBEID")
  message("✅ 已完成新的基因注释")
} else {
  message("💡 已检测到 tt_anno 对象，跳过重新注释")
}

library(dplyr)

tt_gene <- tt_anno %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%       # 去掉未注释探针
  group_by(SYMBOL) %>%                            # 按基因分组
  slice_min(order_by = P.Value, n = 1, with_ties = FALSE) %>%  # 每个基因取最显著探针
  ungroup() %>%
  arrange(adj.P.Val)                              # 按FDR排序

# 检查一下结果
head(tt_gene)
cat("基因层总数：", nrow(tt_gene), "\n")


## 9) 导出
out_dir <- "affymetrix/GSE173954/results_onlyGroup"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(tt_anno, file.path(out_dir, "DEG_probelevel_annotated.csv"), row.names = FALSE)
write.csv(tt_gene, file.path(out_dir, "DEG_genelevel.csv"),            row.names = FALSE)

cat("\n显著基因数（adj.P<0.05）：", sum(tt_gene$adj.P.Val < 0.05, na.rm=TRUE), "\n")
cat("较宽松（adj.P<0.1）：",           sum(tt_gene$adj.P.Val < 0.1,  na.rm=TRUE), "\n")
