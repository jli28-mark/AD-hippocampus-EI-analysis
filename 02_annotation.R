library(AnnotationDbi)
library(hta20transcriptcluster.db)
library(dplyr)
library(tibble)

# ① 取探针/cluster ID（任选其一来源）
probe_ids <- rownames(expr)           # 若在做全表注释
# probe_ids <- rownames(tt)          # 若只给 DEG 结果注释
# probe_ids <- featureNames(eset)    # 也可以从 ExpressionSet 直接取

# ② 用官方注释包做映射
anno <- AnnotationDbi::select(
  hta20transcriptcluster.db,
  keys    = probe_ids,
  keytype = "PROBEID",
  columns = c("SYMBOL","GENENAME","ENTREZID","ENSEMBL")
)

# ③（可选）看一眼一对多的情况
# 每个 probe 映射到多少基因
map_count <- table(anno$PROBEID)
plot(table(map_count), xlim=c(1,50),
     main="HTA-2_0 transcript cluster → #genes",
     xlab="#genes per cluster", ylab="#clusters")

# ④ 保存映射表
save(anno, file = file.path(out_dir, "probe2gene_hta20.Rdata"))
write.csv(anno, file.path(out_dir, "probe2gene_hta20.csv"), row.names = FALSE)
