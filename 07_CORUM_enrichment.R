###############################################
###  CORUM enrichment for mitochondrial genes ###
###############################################

library(dplyr)
library(tidyr)
library(readxl)

## 1. 文件路径（只需要改这两行） ---------------------------
corum_file <- "D:/My Desktop/5.19AD/PPI11.15/CORUM/CORUM.xlsx"
gene_file  <- "D:/My Desktop/5.19AD/PPI11.15/CORUM/mitochondria38.txt"

## 2. 读入线粒体基因 ---------------------------------------
genes <- read.table(gene_file, stringsAsFactors = FALSE)[,1]
genes <- unique(toupper(genes))
length(genes)   # 检查是不是 38 个

## 3. 读入 CORUM --------------------------------------------
corum <- read_xlsx(corum_file)

## 如果有 organism 字段，保留 Human
if ("organism" %in% colnames(corum)) {
  corum_h <- corum %>% filter(organism == "Human")
} else {
  corum_h <- corum
}

## 4. 拆分子单元基因名（subunits_gene_name） ----------------
corum_long <- corum_h %>%
  select(complex_id, complex_name, subunits_gene_name) %>%
  filter(!is.na(subunits_gene_name), subunits_gene_name != "") %>%
  separate_rows(subunits_gene_name, sep = ";") %>%
  rename(member_gene = subunits_gene_name) %>%
  mutate(member_gene = toupper(member_gene))

## CORUM 所有基因作为背景
all_genes_corum <- unique(corum_long$member_gene)
total_genes <- length(all_genes_corum)

## 5. 计算每个复合体的 overlap --------------------------------
res <- corum_long %>%
  group_by(complex_id, complex_name) %>%
  summarise(
    complex_size = n(),
    overlap = sum(member_gene %in% genes),
    overlap_genes = paste(member_gene[member_gene %in% genes], collapse = ";"),
    .groups = "drop"
  ) %>%
  filter(overlap > 0)

## 6. 超几何检验 --------------------------------------------
gene_in_set <- length(genes)

res$pvalue <- phyper(
  q = res$overlap - 1,
  m = res$complex_size,
  n = total_genes - res$complex_size,
  k = gene_in_set,
  lower.tail = FALSE
)

res$FDR <- p.adjust(res$pvalue, method = "BH")

## 7. 输出结果 ------------------------------------------------
out_file <- "D:/My Desktop/5.19AD/PPI11.15/CORUM/CORUM_enrichment_mito38.csv"
write.csv(res, out_file, row.names = FALSE)

print("Mitochondrial CORUM analysis finished. Output:")
print(out_file)
