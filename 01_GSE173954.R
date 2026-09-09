
###############################################################

# GSE173954 COMPLETE DIFFERENTIAL EXPRESSION ANALYSIS
#
# Tissue: Human hippocampus
# Platform: GPL17586
#           Affymetrix Human Transcriptome Array 2.0
#
# Groups:
#   AD      = GSM5283431-GSM5283438 (n = 8)
#   Control = GSM5283439-GSM5283448 (n = 10)
#
# Primary analysis:
#   AD vs Control
#
# Primary model:
#   Expression ~ Group
#
# Sensitivity analysis:
#   Expression ~ age + gender + Group
#
# DEG threshold:
#   BH-adjusted P < 0.05
#   |log2FC| >= 0.58
#
# Input:
#   CEL.gz files
#
# Output:
#   QC plots
#   normalized expression matrix
#   complete limma results
#   annotated gene-level results
#   DEG table
#   sensitivity analysis
#   DEG summary
###############################################################


###############################################################
# 0. CLEAN ENVIRONMENT
###############################################################

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  timeout = 600
)


###############################################################
# 1. FILE PATHS
###############################################################

# Raw CEL folder
raw_dir <- "D:/My Desktop/09.09返修/14_data/GSE173954/GSE173954_RAW"

# Metadata
metadata_file <- "D:/My Desktop/09.09返修/14_data/GSE173954/GSE173954.metadata.xlsx"

# Output folder
out_dir <- "D:/My Desktop/09.09返修/14_data/GSE173954/GSE173954_analysis"

if (!dir.exists(out_dir)) {
  dir.create(
    out_dir,
    recursive = TRUE
  )
}


###############################################################
# 2. INSTALL REQUIRED PACKAGES
###############################################################

cran_packages <- c(
  "readxl",
  "ggplot2",
  "pheatmap"
)

bioc_packages <- c(
  "Biobase",
  "oligo",
  "limma",
  "AnnotationDbi",
  "pd.hta.2.0",
  "hta20transcriptcluster.db"
)


# Install CRAN packages
for (pkg in cran_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    install.packages(
      pkg,
      dependencies = TRUE
    )
  }
}


# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {

  install.packages("BiocManager")
}


# Install Bioconductor packages
for (pkg in bioc_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
  }
}


###############################################################
# 3. LOAD PACKAGES
###############################################################

suppressPackageStartupMessages({

  library(readxl)

  library(Biobase)

  library(oligo)

  library(limma)

  library(AnnotationDbi)

  library(pd.hta.2.0)

  library(hta20transcriptcluster.db)

  library(ggplot2)

  library(pheatmap)

})


###############################################################
# 4. READ METADATA
###############################################################

metadata <- read_excel(metadata_file)

metadata <- as.data.frame(metadata)


# Clean column names
colnames(metadata) <- trimws(colnames(metadata))


cat("\n========================================\n")
cat("METADATA COLUMN NAMES\n")
cat("========================================\n")

print(colnames(metadata))


###############################################################
# 5. STANDARDIZE METADATA COLUMN NAMES
###############################################################

# Search sample column
sample_col <- grep(
  "^sample$|gsm",
  colnames(metadata),
  ignore.case = TRUE,
  value = TRUE
)[1]

# Search group column
group_col <- grep(
  "^group$|diagnosis|disease",
  colnames(metadata),
  ignore.case = TRUE,
  value = TRUE
)[1]

# Search age
age_col <- grep(
  "^age$",
  colnames(metadata),
  ignore.case = TRUE,
  value = TRUE
)[1]

# Search sex/gender
gender_col <- grep(
  "^gender$|^sex$",
  colnames(metadata),
  ignore.case = TRUE,
  value = TRUE
)[1]


if (is.na(sample_col)) {
  stop("Could not identify Sample/GSM column in metadata.")
}

if (is.na(group_col)) {
  stop("Could not identify Group column in metadata.")
}


# Rename
colnames(metadata)[
  colnames(metadata) == sample_col
] <- "Sample"

colnames(metadata)[
  colnames(metadata) == group_col
] <- "Group"


if (!is.na(age_col)) {

  colnames(metadata)[
    colnames(metadata) == age_col
  ] <- "age"

}


if (!is.na(gender_col)) {

  colnames(metadata)[
    colnames(metadata) == gender_col
  ] <- "gender"

}


###############################################################
# 6. CLEAN METADATA VALUES
###############################################################

metadata$Sample <- trimws(
  as.character(metadata$Sample)
)

metadata$Group <- trimws(
  as.character(metadata$Group)
)


# Normalize group names
metadata$Group[
  tolower(metadata$Group) %in%
    c(
      "control",
      "ctrl",
      "normal",
      "non-ad",
      "non ad",
      "nonad"
    )
] <- "Control"


metadata$Group[
  tolower(metadata$Group) %in%
    c(
      "ad",
      "alzheimer",
      "alzheimer's disease",
      "alzheimers disease"
    )
] <- "AD"


metadata$Group <- factor(
  metadata$Group,
  levels = c(
    "Control",
    "AD"
  )
)


###############################################################
# 7. CLEAN AGE / GENDER
###############################################################

if ("age" %in% colnames(metadata)) {

  metadata$age <- as.numeric(
    gsub(
      "[^0-9.]",
      "",
      as.character(metadata$age)
    )
  )

}


if ("gender" %in% colnames(metadata)) {

  metadata$gender <- trimws(
    as.character(metadata$gender)
  )

  metadata$gender <- tolower(
    metadata$gender
  )

  metadata$gender[
    metadata$gender %in% c("m", "male")
  ] <- "Male"

  metadata$gender[
    metadata$gender %in% c("f", "female")
  ] <- "Female"

  metadata$gender <- factor(metadata$gender)

}


###############################################################
# 8. CHECK METADATA
###############################################################

cat("\n========================================\n")
cat("METADATA\n")
cat("========================================\n")

print(metadata)


cat("\nGroup counts:\n")

print(
  table(
    metadata$Group,
    useNA = "ifany"
  )
)


if (any(is.na(metadata$Group))) {

  stop(
    "Some samples cannot be assigned to AD or Control."
  )

}


if (nrow(metadata) != 18) {

  warning(
    paste0(
      "Metadata contains ",
      nrow(metadata),
      " rows; expected 18."
    )
  )

}


###############################################################
# 9. CHECK EXPECTED GSM IDs
###############################################################

expected_AD <- paste0(
  "GSM",
  5283431:5283438
)

expected_Control <- paste0(
  "GSM",
  5283439:5283448
)

expected_samples <- c(
  expected_AD,
  expected_Control
)


missing_samples <- setdiff(
  expected_samples,
  metadata$Sample
)


if (length(missing_samples) > 0) {

  stop(
    paste(
      "Missing expected samples:",
      paste(
        missing_samples,
        collapse = ", "
      )
    )
  )

}


###############################################################
# 10. FIND CEL FILES
###############################################################

cel_files <- list.files(
  raw_dir,
  pattern = "\\.CEL(\\.gz)?$",
  full.names = TRUE,
  ignore.case = TRUE
)


cat("\n========================================\n")
cat("CEL FILES\n")
cat("========================================\n")

cat(
  "Number of CEL files found:",
  length(cel_files),
  "\n"
)


if (length(cel_files) != 18) {

  stop(
    paste0(
      "Expected 18 CEL files but found ",
      length(cel_files),
      "."
    )
  )

}


###############################################################
# 11. EXTRACT GSM FROM CEL FILENAMES
###############################################################

gsm_from_file <- sub(
  "^(GSM[0-9]+).*",
  "\\1",
  basename(cel_files)
)


cel_information <- data.frame(
  GSM = gsm_from_file,
  CELfile = basename(cel_files)
)


print(cel_information)


###############################################################
# 12. MATCH CEL FILES WITH METADATA
###############################################################

match_index <- match(
  metadata$Sample,
  gsm_from_file
)


if (any(is.na(match_index))) {

  stop(
    "Some metadata GSM IDs cannot be matched to CEL files."
  )

}


# reorder CEL files according to metadata
cel_files <- cel_files[
  match_index
]


# final GSM check
gsm_check <- sub(
  "^(GSM[0-9]+).*",
  "\\1",
  basename(cel_files)
)


if (!all(
  gsm_check == metadata$Sample
)) {

  stop(
    "CEL files and metadata order do not match."
  )

}


cat(
  "\nCEL files and metadata matched successfully.\n"
)


###############################################################
# 13. SAVE METADATA USED
###############################################################

write.csv(
  metadata,
  file.path(
    out_dir,
    "GSE173954_metadata_used.csv"
  ),
  row.names = FALSE
)


###############################################################
# 14. DESCRIPTIVE AGE/SEX INFORMATION
###############################################################

if ("age" %in% colnames(metadata)) {

  cat("\n========================================\n")
  cat("AGE SUMMARY\n")
  cat("========================================\n")

  print(
    aggregate(
      age ~ Group,
      data = metadata,
      FUN = function(x) {
        c(
          n = length(x),
          mean = mean(x, na.rm = TRUE),
          SD = sd(x, na.rm = TRUE),
          median = median(x, na.rm = TRUE),
          min = min(x, na.rm = TRUE),
          max = max(x, na.rm = TRUE)
        )
      }
    )
  )

}


if ("gender" %in% colnames(metadata)) {

  cat("\n========================================\n")
  cat("SEX DISTRIBUTION\n")
  cat("========================================\n")

  print(
    table(
      metadata$Group,
      metadata$gender
    )
  )

}


###############################################################
# 15. READ CEL FILES
###############################################################

cat("\n========================================\n")
cat("READING CEL FILES\n")
cat("========================================\n")


raw_data <- read.celfiles(
  cel_files
)


cat(
  "CEL files successfully loaded.\n"
)


###############################################################
# 16. RAW DATA QC
###############################################################

tiff(
  file.path(
    out_dir,
    "QC_01_Raw_Boxplot.tiff"
  ),
  width = 3400,
  height = 2400,
  res = 300,
  compression = "lzw"
)


boxplot(
  raw_data,
  target = "core",
  las = 2,
  main = "GSE173954 raw CEL intensity",
  ylab = "Raw intensity"
)

dev.off()


###############################################################
# 17. RMA NORMALIZATION
###############################################################

cat("\n========================================\n")
cat("RMA NORMALIZATION\n")
cat("========================================\n")


eset <- rma(
  raw_data,
  target = "core"
)


expr <- exprs(eset)


# Rename samples using GSM numbers
colnames(expr) <- metadata$Sample


cat(
  "Expression matrix:",
  nrow(expr),
  "features x",
  ncol(expr),
  "samples\n"
)


###############################################################
# 18. SAVE NORMALIZED EXPRESSION MATRIX
###############################################################

write.csv(
  expr,
  file.path(
    out_dir,
    "GSE173954_RMA_expression_matrix_transcript_cluster.csv"
  )
)


###############################################################
# 19. NORMALIZED DATA BOXPLOT
###############################################################

tiff(
  file.path(
    out_dir,
    "QC_02_RMA_Boxplot.tiff"
  ),
  width = 3400,
  height = 2400,
  res = 300,
  compression = "lzw"
)


boxplot(
  expr,
  las = 2,
  outline = FALSE,
  main = "GSE173954 RMA-normalized expression",
  ylab = "log2 expression"
)

dev.off()


###############################################################
# 20. EXPRESSION DENSITY
###############################################################

tiff(
  file.path(
    out_dir,
    "QC_03_RMA_Density.tiff"
  ),
  width = 3000,
  height = 2200,
  res = 300,
  compression = "lzw"
)


limma::plotDensities(
  expr,
  main = "GSE173954 expression density after RMA",
  xlab = "log2 expression"
)

dev.off()


###############################################################
# 21. PCA
###############################################################

pca <- prcomp(
  t(expr),
  center = TRUE,
  scale. = FALSE
)


pca_variance <- (
  pca$sdev^2 /
    sum(pca$sdev^2)
) * 100


pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = metadata$Group
)


if ("gender" %in% colnames(metadata)) {

  pca_df$gender <- metadata$gender

}


if ("age" %in% colnames(metadata)) {

  pca_df$age <- metadata$age

}


p_pca <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    shape = Group
  )
) +

  geom_point(
    size = 4
  ) +

  geom_text(
    aes(
      label = Sample
    ),
    size = 2.5,
    vjust = -0.8
  ) +

  xlab(
    paste0(
      "PC1 (",
      round(
        pca_variance[1],
        1
      ),
      "%)"
    )
  ) +

  ylab(
    paste0(
      "PC2 (",
      round(
        pca_variance[2],
        1
      ),
      "%)"
    )
  ) +

  ggtitle(
    "GSE173954 PCA"
  ) +

  theme_bw(
    base_size = 13
  )


ggsave(
  file.path(
    out_dir,
    "QC_04_PCA.tiff"
  ),
  p_pca,
  width = 8,
  height = 6,
  dpi = 300,
  compression = "lzw"
)


write.csv(
  pca_df,
  file.path(
    out_dir,
    "GSE173954_PCA_coordinates.csv"
  ),
  row.names = FALSE
)


###############################################################
# 22. SAMPLE CORRELATION
###############################################################

cor_matrix <- cor(
  expr,
  method = "pearson"
)


annotation_sample <- data.frame(
  Group = metadata$Group
)


if ("gender" %in% colnames(metadata)) {

  annotation_sample$gender <- metadata$gender

}


rownames(annotation_sample) <- metadata$Sample


tiff(
  file.path(
    out_dir,
    "QC_05_Sample_Correlation.tiff"
  ),
  width = 3200,
  height = 3000,
  res = 300,
  compression = "lzw"
)


pheatmap(
  cor_matrix,
  annotation_col = annotation_sample,
  annotation_row = annotation_sample,
  main = "GSE173954 sample correlation"
)


dev.off()


###############################################################
# 23. HIERARCHICAL CLUSTERING
###############################################################

sample_distance <- dist(
  t(expr)
)


hc <- hclust(
  sample_distance,
  method = "complete"
)


tiff(
  file.path(
    out_dir,
    "QC_06_Hierarchical_Clustering.tiff"
  ),
  width = 3400,
  height = 2400,
  res = 300,
  compression = "lzw"
)


plot(
  hc,
  labels = paste0(
    metadata$Sample,
    "_",
    metadata$Group
  ),
  main = "GSE173954 hierarchical clustering",
  xlab = "",
  sub = "",
  cex = 0.75
)


dev.off()


###############################################################
# 24. PRIMARY DIFFERENTIAL EXPRESSION MODEL
#
# Expression ~ Group
###############################################################

cat("\n========================================\n")
cat("PRIMARY LIMMA ANALYSIS\n")
cat("AD vs CONTROL\n")
cat("========================================\n")


design <- model.matrix(
  ~ 0 + Group,
  data = metadata
)


colnames(design) <- c(
  "Control",
  "AD"
)


rownames(design) <- metadata$Sample


print(design)


###############################################################
# 25. LIMMA FIT
###############################################################

fit <- lmFit(
  expr,
  design
)


contrast_matrix <- makeContrasts(
  AD_vs_Control = AD - Control,
  levels = design
)


fit2 <- contrasts.fit(
  fit,
  contrast_matrix
)


fit2 <- eBayes(
  fit2,
  trend = TRUE
)


###############################################################
# 26. EXTRACT ALL FEATURES
###############################################################

results_feature <- topTable(
  fit2,
  coef = "AD_vs_Control",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)


results_feature$ProbeID <- rownames(
  results_feature
)


rownames(results_feature) <- NULL


# Put ProbeID first
results_feature <- results_feature[
  ,
  c(
    "ProbeID",
    setdiff(
      colnames(results_feature),
      "ProbeID"
    )
  )
]


write.csv(
  results_feature,
  file.path(
    out_dir,
    "GSE173954_01_limma_all_transcript_clusters.csv"
  ),
  row.names = FALSE
)


###############################################################
# 27. HTA 2.0 ANNOTATION
###############################################################

cat("\n========================================\n")
cat("ANNOTATION\n")
cat("========================================\n")


probe_ids <- as.character(
  results_feature$ProbeID
)


# Check available annotation fields
available_columns <- AnnotationDbi::columns(
  hta20transcriptcluster.db
)


cat(
  "Available HTA annotation columns:\n"
)

print(
  available_columns
)


# Annotation
anno <- AnnotationDbi::select(
  hta20transcriptcluster.db,
  keys = probe_ids,
  keytype = "PROBEID",
  columns = c(
    "SYMBOL",
    "ENTREZID",
    "GENENAME"
  )
)


###############################################################
# 28. CLEAN ANNOTATION
###############################################################

anno$SYMBOL <- trimws(
  as.character(
    anno$SYMBOL
  )
)


anno$SYMBOL[
  anno$SYMBOL == ""
] <- NA


# Remove exact duplicated mapping rows
anno <- unique(
  anno
)


###############################################################
# 29. MERGE EXPRESSION RESULTS + ANNOTATION
###############################################################

results_annotated <- merge(
  results_feature,
  anno,
  by.x = "ProbeID",
  by.y = "PROBEID",
  all.x = TRUE,
  sort = FALSE
)


write.csv(
  results_annotated,
  file.path(
    out_dir,
    "GSE173954_02_all_features_annotated.csv"
  ),
  row.names = FALSE
)


###############################################################
# 30. REMOVE FEATURES WITHOUT GENE SYMBOL
###############################################################

results_gene <- results_annotated[
  !is.na(results_annotated$SYMBOL),
]


###############################################################
# 31. RESOLVE DUPLICATED GENE SYMBOLS
#
# If multiple transcript clusters map to one gene:
# retain the transcript cluster with highest AveExpr
###############################################################

results_gene <- results_gene[
  order(
    results_gene$SYMBOL,
    -results_gene$AveExpr
  ),
]


results_gene_unique <- results_gene[
  !duplicated(
    results_gene$SYMBOL
  ),
]


###############################################################
# 32. ADD DIRECTION
###############################################################

results_gene_unique$Direction <- "Not_significant"


results_gene_unique$Direction[
  results_gene_unique$adj.P.Val < 0.05 &
    results_gene_unique$logFC >= 0.58
] <- "Up"


results_gene_unique$Direction[
  results_gene_unique$adj.P.Val < 0.05 &
    results_gene_unique$logFC <= -0.58
] <- "Down"


###############################################################
# 33. SAVE ALL GENE-LEVEL RESULTS
###############################################################

results_gene_unique <- results_gene_unique[
  order(
    results_gene_unique$adj.P.Val
  ),
]


write.csv(
  results_gene_unique,
  file.path(
    out_dir,
    "GSE173954_03_all_genes_limma.csv"
  ),
  row.names = FALSE
)


###############################################################
# 34. EXTRACT SIGNIFICANT DEGs
###############################################################

DEG <- results_gene_unique[
  !is.na(results_gene_unique$adj.P.Val) &
    results_gene_unique$adj.P.Val < 0.05 &
    abs(results_gene_unique$logFC) >= 0.58,
]


DEG <- DEG[
  order(
    DEG$adj.P.Val
  ),
]


write.csv(
  DEG,
  file.path(
    out_dir,
    "GSE173954_04_DEG_adjP0.05_logFC0.58.csv"
  ),
  row.names = FALSE
)


###############################################################
# 35. SEPARATE UP/DOWN
###############################################################

DEG_up <- DEG[
  DEG$Direction == "Up",
]


DEG_down <- DEG[
  DEG$Direction == "Down",
]


write.csv(
  DEG_up,
  file.path(
    out_dir,
    "GSE173954_05_DEG_UP.csv"
  ),
  row.names = FALSE
)


write.csv(
  DEG_down,
  file.path(
    out_dir,
    "GSE173954_06_DEG_DOWN.csv"
  ),
  row.names = FALSE
)


###############################################################
# 36. DEG SUMMARY
###############################################################

n_total <- nrow(DEG)

n_up <- nrow(DEG_up)

n_down <- nrow(DEG_down)


summary_primary <- data.frame(

  Dataset = "GSE173954",

  Platform = "GPL17586",

  Technology =
    "Affymetrix Human Transcriptome Array 2.0",

  Tissue =
    "Human hippocampus",

  AD_n =
    sum(metadata$Group == "AD"),

  Control_n =
    sum(metadata$Group == "Control"),

  Model =
    "AD vs Control",

  Adjusted_P_cutoff =
    0.05,

  Abs_log2FC_cutoff =
    0.58,

  Total_DEGs =
    n_total,

  Upregulated =
    n_up,

  Downregulated =
    n_down

)


write.csv(
  summary_primary,
  file.path(
    out_dir,
    "GSE173954_07_DEG_summary.csv"
  ),
  row.names = FALSE
)


cat("\n========================================\n")
cat("PRIMARY ANALYSIS SUMMARY\n")
cat("========================================\n")

print(
  summary_primary
)


###############################################################
# 37. VOLCANO PLOT
###############################################################

volcano_df <- results_gene_unique


volcano_df$minusLog10FDR <- -log10(
  pmax(
    volcano_df$adj.P.Val,
    .Machine$double.xmin
  )
)


volcano_plot <- ggplot(
  volcano_df,
  aes(
    x = logFC,
    y = minusLog10FDR,
    shape = Direction
  )
) +

  geom_point(
    alpha = 0.55,
    size = 1.4
  ) +

  geom_vline(
    xintercept = c(
      -0.58,
      0.58
    ),
    linetype = "dashed"
  ) +

  geom_hline(
    yintercept =
      -log10(0.05),
    linetype = "dashed"
  ) +

  labs(
    title = paste0(
      "GSE173954: ",
      n_total,
      " DEGs"
    ),
    subtitle = paste0(
      n_up,
      " upregulated; ",
      n_down,
      " downregulated"
    ),
    x = "log2 fold change",
    y = "-log10 adjusted P value"
  ) +

  theme_bw(
    base_size = 13
  )


ggsave(
  file.path(
    out_dir,
    "GSE173954_08_Volcano.tiff"
  ),
  volcano_plot,
  width = 7,
  height = 6,
  dpi = 300,
  compression = "lzw"
)


###############################################################
# 38. OPTIONAL TOP-50 DEG HEATMAP
###############################################################

if (nrow(DEG) >= 2) {

  top_n <- min(
    50,
    nrow(DEG)
  )


  top_genes <- DEG$ProbeID[
    1:top_n
  ]


  heat_expr <- expr[
    rownames(expr) %in% top_genes,
    ,
    drop = FALSE
  ]


  # Map ProbeID -> SYMBOL
  symbol_map <- DEG$SYMBOL[
    match(
      rownames(heat_expr),
      DEG$ProbeID
    )
  ]


  # Make row names unique
  rownames(heat_expr) <- make.unique(
    ifelse(
      is.na(symbol_map),
      rownames(heat_expr),
      symbol_map
    )
  )


  # z-score by gene
  heat_z <- t(
    scale(
      t(heat_expr)
    )
  )


  tiff(
    file.path(
      out_dir,
      "GSE173954_09_Top50_DEG_heatmap.tiff"
    ),
    width = 2800,
    height = 3600,
    res = 300,
    compression = "lzw"
  )


  pheatmap(
    heat_z,
    annotation_col = annotation_sample,
    show_colnames = TRUE,
    show_rownames = TRUE,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    main = "Top DEGs in GSE173954"
  )


  dev.off()

}


###############################################################
# 39. SENSITIVITY ANALYSIS
#
# age + gender + Group
#
# This does NOT replace primary analysis.
###############################################################

run_sensitivity <- all(
  c(
    "age",
    "gender"
  ) %in%
    colnames(metadata)
)


if (run_sensitivity) {


  complete_covariates <- complete.cases(
    metadata[
      ,
      c(
        "age",
        "gender",
        "Group"
      )
    ]
  )


  if (
    all(complete_covariates) &&
      length(unique(metadata$gender)) >= 2
  ) {


    cat("\n========================================\n")
    cat("SENSITIVITY ANALYSIS\n")
    cat("AGE + SEX ADJUSTMENT\n")
    cat("========================================\n")


    design_adjusted <- model.matrix(
      ~ age + gender + Group,
      data = metadata
    )


    print(
      design_adjusted
    )


    fit_adjusted <- lmFit(
      expr,
      design_adjusted
    )


    fit_adjusted <- eBayes(
      fit_adjusted,
      trend = TRUE
    )


    if (!"GroupAD" %in%
        colnames(
          design_adjusted
        )) {

      stop(
        paste0(
          "Cannot find GroupAD coefficient. ",
          "Design columns are: ",
          paste(
            colnames(
              design_adjusted
            ),
            collapse = ", "
          )
        )
      )

    }


    adjusted_feature <- topTable(
      fit_adjusted,
      coef = "GroupAD",
      number = Inf,
      adjust.method = "BH",
      sort.by = "P"
    )


    adjusted_feature$ProbeID <-
      rownames(
        adjusted_feature
      )


    rownames(
      adjusted_feature
    ) <- NULL


    ###########################################################
    # Add annotation
    ###########################################################

    adjusted_annotated <- merge(
      adjusted_feature,
      anno,
      by.x = "ProbeID",
      by.y = "PROBEID",
      all.x = TRUE,
      sort = FALSE
    )


    adjusted_annotated <- adjusted_annotated[
      !is.na(
        adjusted_annotated$SYMBOL
      ),
    ]


    adjusted_annotated <- adjusted_annotated[
      order(
        adjusted_annotated$SYMBOL,
        -adjusted_annotated$AveExpr
      ),
    ]


    adjusted_gene <- adjusted_annotated[
      !duplicated(
        adjusted_annotated$SYMBOL
      ),
    ]


    adjusted_gene$Direction <-
      "Not_significant"


    adjusted_gene$Direction[
      adjusted_gene$adj.P.Val < 0.05 &
        adjusted_gene$logFC >= 0.58
    ] <- "Up"


    adjusted_gene$Direction[
      adjusted_gene$adj.P.Val < 0.05 &
        adjusted_gene$logFC <= -0.58
    ] <- "Down"


    adjusted_DEG <- adjusted_gene[
      adjusted_gene$adj.P.Val < 0.05 &
        abs(
          adjusted_gene$logFC
        ) >= 0.58,
    ]


    adjusted_DEG <- adjusted_DEG[
      order(
        adjusted_DEG$adj.P.Val
      ),
    ]


    write.csv(
      adjusted_gene,
      file.path(
        out_dir,
        "GSE173954_10_age_sex_adjusted_all_genes.csv"
      ),
      row.names = FALSE
    )


    write.csv(
      adjusted_DEG,
      file.path(
        out_dir,
        "GSE173954_11_age_sex_adjusted_DEG.csv"
      ),
      row.names = FALSE
    )


    ###########################################################
    # Sensitivity summary
    ###########################################################

    sensitivity_summary <- data.frame(

      Dataset =
        "GSE173954",

      Model =
        "Group + age + gender",

      Total_DEGs =
        nrow(
          adjusted_DEG
        ),

      Upregulated =
        sum(
          adjusted_DEG$Direction ==
            "Up"
        ),

      Downregulated =
        sum(
          adjusted_DEG$Direction ==
            "Down"
        )

    )


    write.csv(
      sensitivity_summary,
      file.path(
        out_dir,
        "GSE173954_12_age_sex_adjusted_summary.csv"
      ),
      row.names = FALSE
    )


    ###########################################################
    # Compare primary vs adjusted results
    ###########################################################

    comparison <- merge(

      results_gene_unique[
        ,
        c(
          "SYMBOL",
          "logFC",
          "adj.P.Val",
          "Direction"
        )
      ],

      adjusted_gene[
        ,
        c(
          "SYMBOL",
          "logFC",
          "adj.P.Val",
          "Direction"
        )
      ],

      by = "SYMBOL",

      suffixes = c(
        "_Primary",
        "_Adjusted"
      )

    )


    write.csv(
      comparison,
      file.path(
        out_dir,
        "GSE173954_13_primary_vs_age_sex_adjusted.csv"
      ),
      row.names = FALSE
    )


    ###########################################################
    # DEG overlap
    ###########################################################

    primary_symbols <-
      unique(
        DEG$SYMBOL
      )


    adjusted_symbols <-
      unique(
        adjusted_DEG$SYMBOL
      )


    overlap_symbols <-
      intersect(
        primary_symbols,
        adjusted_symbols
      )


    overlap_summary <- data.frame(

      Primary_DEGs =
        length(
          primary_symbols
        ),

      Adjusted_DEGs =
        length(
          adjusted_symbols
        ),

      Shared_DEGs =
        length(
          overlap_symbols
        ),

      Primary_overlap_percent =
        ifelse(
          length(primary_symbols) > 0,
          length(overlap_symbols) /
            length(primary_symbols) *
            100,
          NA
        )

    )


    write.csv(
      overlap_summary,
      file.path(
        out_dir,
        "GSE173954_14_sensitivity_overlap_summary.csv"
      ),
      row.names = FALSE
    )


    write.csv(
      data.frame(
        SYMBOL =
          overlap_symbols
      ),
      file.path(
        out_dir,
        "GSE173954_15_shared_DEGs_primary_adjusted.csv"
      ),
      row.names = FALSE
    )


    cat("\nSensitivity analysis summary:\n")

    print(
      sensitivity_summary
    )


    cat("\nPrimary-adjusted overlap:\n")

    print(
      overlap_summary
    )

  } else {

    warning(
      "Age/sex sensitivity analysis skipped because covariate data are incomplete."
    )

  }

}


###############################################################
# 40. SAVE R SESSION INFORMATION
###############################################################

sink(
  file.path(
    out_dir,
    "GSE173954_sessionInfo.txt"
  )
)


sessionInfo()


sink()


###############################################################
# 41. SAVE WORKSPACE
###############################################################

save(
  metadata,
  expr,
  results_gene_unique,
  DEG,
  summary_primary,
  file = file.path(
    out_dir,
    "GSE173954_analysis_objects.RData"
  )
)


###############################################################
# 42. FINAL OUTPUT
###############################################################

cat("\n\n")
cat("=============================================\n")
cat("GSE173954 ANALYSIS COMPLETED\n")
cat("=============================================\n")

cat(
  "\nAD:",
  sum(
    metadata$Group ==
      "AD"
  )
)

cat(
  "\nControl:",
  sum(
    metadata$Group ==
      "Control"
  )
)

cat(
  "\n\nPrimary model: AD vs Control"
)

cat(
  "\nDEG threshold: adj.P < 0.05 and |log2FC| >= 0.58"
)

cat(
  "\n\nTotal DEGs:",
  n_total
)

cat(
  "\nUpregulated:",
  n_up
)

cat(
  "\nDownregulated:",
  n_down
)

cat(
  "\n\nResults saved to:\n"
)

cat(
  out_dir
)

cat("\n=============================================\n")



###############################################################
# GSE173954 OUTLIER SENSITIVITY ANALYSIS
#
# Compare:
#   Main analysis: 8 AD + 10 Control = 18 samples
#   Sensitivity:   remove GSM5283438
#                  7 AD + 10 Control = 17 samples
#
# REQUIREMENTS:
# This code should be run AFTER the previous complete script.
#
# Required objects already present:
#   expr
#   metadata
#   results_gene_unique
#   DEG
#   anno
#   out_dir
###############################################################


###############################################################
# 1. DEFINE OUTLIER
###############################################################

outlier_sample <- "GSM5283438"

cat("\n\n")
cat("=====================================================\n")
cat("OUTLIER SENSITIVITY ANALYSIS\n")
cat("Removing:", outlier_sample, "\n")
cat("=====================================================\n")


###############################################################
# 2. CHECK OBJECTS
###############################################################

required_objects <- c(
  "expr",
  "metadata",
  "results_gene_unique",
  "DEG",
  "anno",
  "out_dir"
)

missing_objects <- required_objects[
  !sapply(
    required_objects,
    exists
  )
]

if (length(missing_objects) > 0) {

  stop(
    paste0(
      "Missing required objects: ",
      paste(
        missing_objects,
        collapse = ", "
      ),
      ". Run the main GSE173954 analysis first."
    )
  )

}


###############################################################
# 3. CHECK OUTLIER EXISTS
###############################################################

if (!outlier_sample %in% colnames(expr)) {

  stop(
    paste0(
      outlier_sample,
      " is not present in expression matrix."
    )
  )

}

if (!outlier_sample %in% metadata$Sample) {

  stop(
    paste0(
      outlier_sample,
      " is not present in metadata."
    )
  )

}


###############################################################
# 4. CREATE SENSITIVITY OUTPUT FOLDER
###############################################################

sens_dir <- file.path(
  out_dir,
  "Sensitivity_remove_GSM5283438"
)

if (!dir.exists(sens_dir)) {

  dir.create(
    sens_dir,
    recursive = TRUE
  )

}


###############################################################
# 5. REMOVE GSM5283438
###############################################################

keep_samples <- setdiff(
  colnames(expr),
  outlier_sample
)


expr_sens <- expr[
  ,
  keep_samples,
  drop = FALSE
]


metadata_sens <- metadata[
  metadata$Sample %in% keep_samples,
  ,
  drop = FALSE
]


# Reorder metadata exactly as expression matrix
metadata_sens <- metadata_sens[
  match(
    colnames(expr_sens),
    metadata_sens$Sample
  ),
  ,
  drop = FALSE
]


stopifnot(
  all(
    metadata_sens$Sample ==
      colnames(expr_sens)
  )
)


cat(
  "\nSamples after removal:",
  ncol(expr_sens),
  "\n"
)

cat(
  "AD:",
  sum(
    metadata_sens$Group == "AD"
  ),
  "\n"
)

cat(
  "Control:",
  sum(
    metadata_sens$Group == "Control"
  ),
  "\n"
)


write.csv(
  metadata_sens,
  file.path(
    sens_dir,
    "metadata_without_GSM5283438.csv"
  ),
  row.names = FALSE
)


###############################################################
# 6. PCA AFTER OUTLIER REMOVAL
###############################################################

pca_sens <- prcomp(
  t(expr_sens),
  center = TRUE,
  scale. = FALSE
)


pca_var_sens <- (
  pca_sens$sdev^2 /
    sum(pca_sens$sdev^2)
) * 100


pca_df_sens <- data.frame(

  Sample = rownames(
    pca_sens$x
  ),

  PC1 = pca_sens$x[, 1],

  PC2 = pca_sens$x[, 2],

  Group = metadata_sens$Group

)


p_pca_sens <- ggplot(
  pca_df_sens,
  aes(
    x = PC1,
    y = PC2,
    shape = Group
  )
) +

  geom_point(
    size = 4
  ) +

  geom_text(
    aes(
      label = Sample
    ),
    size = 2.5,
    vjust = -0.8
  ) +

  xlab(
    paste0(
      "PC1 (",
      round(
        pca_var_sens[1],
        1
      ),
      "%)"
    )
  ) +

  ylab(
    paste0(
      "PC2 (",
      round(
        pca_var_sens[2],
        1
      ),
      "%)"
    )
  ) +

  ggtitle(
    "GSE173954 PCA after removal of GSM5283438"
  ) +

  theme_bw(
    base_size = 13
  )


ggsave(
  file.path(
    sens_dir,
    "01_PCA_without_GSM5283438.tiff"
  ),
  p_pca_sens,
  width = 8,
  height = 6,
  dpi = 300,
  compression = "lzw"
)


write.csv(
  pca_df_sens,
  file.path(
    sens_dir,
    "PCA_coordinates_without_GSM5283438.csv"
  ),
  row.names = FALSE
)


###############################################################
# 7. SAMPLE CORRELATION AFTER REMOVAL
###############################################################

cor_sens <- cor(
  expr_sens,
  method = "pearson"
)


annotation_sens <- data.frame(
  Group = metadata_sens$Group
)

rownames(
  annotation_sens
) <- metadata_sens$Sample


tiff(
  file.path(
    sens_dir,
    "02_Sample_Correlation_without_GSM5283438.tiff"
  ),
  width = 3000,
  height = 2800,
  res = 300,
  compression = "lzw"
)


pheatmap(
  cor_sens,
  annotation_col = annotation_sens,
  annotation_row = annotation_sens,
  main = "GSE173954 sample correlation after GSM5283438 removal"
)


dev.off()


###############################################################
# 8. HIERARCHICAL CLUSTERING AFTER REMOVAL
###############################################################

sample_dist_sens <- dist(
  t(expr_sens)
)


hc_sens <- hclust(
  sample_dist_sens,
  method = "complete"
)


tiff(
  file.path(
    sens_dir,
    "03_Hierarchical_Clustering_without_GSM5283438.tiff"
  ),
  width = 3400,
  height = 2400,
  res = 300,
  compression = "lzw"
)


plot(
  hc_sens,
  labels = paste0(
    metadata_sens$Sample,
    "_",
    metadata_sens$Group
  ),
  main = "GSE173954 clustering after GSM5283438 removal",
  xlab = "",
  sub = "",
  cex = 0.75
)


dev.off()


###############################################################
# 9. PRIMARY LIMMA MODEL AFTER OUTLIER REMOVAL
#
# Still use:
# Expression ~ Group
###############################################################

design_sens <- model.matrix(
  ~ 0 + Group,
  data = metadata_sens
)


colnames(
  design_sens
) <- c(
  "Control",
  "AD"
)


rownames(
  design_sens
) <- metadata_sens$Sample


fit_sens <- lmFit(
  expr_sens,
  design_sens
)


contrast_sens <- makeContrasts(
  AD_vs_Control = AD - Control,
  levels = design_sens
)


fit_sens2 <- contrasts.fit(
  fit_sens,
  contrast_sens
)


fit_sens2 <- eBayes(
  fit_sens2,
  trend = TRUE
)


###############################################################
# 10. EXTRACT FEATURE-LEVEL RESULTS
###############################################################

res_sens_feature <- topTable(
  fit_sens2,
  coef = "AD_vs_Control",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)


res_sens_feature$ProbeID <- rownames(
  res_sens_feature
)


rownames(
  res_sens_feature
) <- NULL


###############################################################
# 11. ANNOTATION
###############################################################

res_sens_annotated <- merge(
  res_sens_feature,
  anno,
  by.x = "ProbeID",
  by.y = "PROBEID",
  all.x = TRUE,
  sort = FALSE
)


res_sens_annotated <- res_sens_annotated[
  !is.na(
    res_sens_annotated$SYMBOL
  ),
  ,
  drop = FALSE
]


###############################################################
# 12. ONE FEATURE PER GENE
#
# Same rule as main analysis:
# retain feature with highest AveExpr
###############################################################

res_sens_annotated <- res_sens_annotated[
  order(
    res_sens_annotated$SYMBOL,
    -res_sens_annotated$AveExpr
  ),
]


res_sens_gene <- res_sens_annotated[
  !duplicated(
    res_sens_annotated$SYMBOL
  ),
]


###############################################################
# 13. DEG CLASSIFICATION
###############################################################

res_sens_gene$Direction <- "Not_significant"


res_sens_gene$Direction[
  res_sens_gene$adj.P.Val < 0.05 &
    res_sens_gene$logFC >= 0.58
] <- "Up"


res_sens_gene$Direction[
  res_sens_gene$adj.P.Val < 0.05 &
    res_sens_gene$logFC <= -0.58
] <- "Down"


DEG_sens <- res_sens_gene[
  !is.na(
    res_sens_gene$adj.P.Val
  ) &
    res_sens_gene$adj.P.Val < 0.05 &
    abs(
      res_sens_gene$logFC
    ) >= 0.58,
]


DEG_sens <- DEG_sens[
  order(
    DEG_sens$adj.P.Val
  ),
]


###############################################################
# 14. SAVE SENSITIVITY DEG RESULTS
###############################################################

write.csv(
  res_sens_gene,
  file.path(
    sens_dir,
    "04_all_genes_without_GSM5283438.csv"
  ),
  row.names = FALSE
)


write.csv(
  DEG_sens,
  file.path(
    sens_dir,
    "05_DEG_without_GSM5283438.csv"
  ),
  row.names = FALSE
)


###############################################################
# 15. DEG COUNTS
###############################################################

n_sens <- nrow(
  DEG_sens
)

n_sens_up <- sum(
  DEG_sens$Direction == "Up"
)

n_sens_down <- sum(
  DEG_sens$Direction == "Down"
)


cat(
  "\n========================================\n"
)

cat(
  "DEG RESULTS WITHOUT GSM5283438\n"
)

cat(
  "========================================\n"
)

cat(
  "Total DEG:",
  n_sens,
  "\n"
)

cat(
  "Up:",
  n_sens_up,
  "\n"
)

cat(
  "Down:",
  n_sens_down,
  "\n"
)


###############################################################
# 16. COMPARE ORIGINAL vs OUTLIER-REMOVED DEGs
###############################################################

primary_symbols <- unique(
  DEG$SYMBOL
)

sens_symbols <- unique(
  DEG_sens$SYMBOL
)


shared_symbols <- intersect(
  primary_symbols,
  sens_symbols
)


lost_symbols <- setdiff(
  primary_symbols,
  sens_symbols
)


new_symbols <- setdiff(
  sens_symbols,
  primary_symbols
)


union_symbols <- union(
  primary_symbols,
  sens_symbols
)


jaccard_index <- ifelse(
  length(union_symbols) > 0,
  length(shared_symbols) /
    length(union_symbols),
  NA
)


primary_overlap_pct <- ifelse(
  length(primary_symbols) > 0,
  length(shared_symbols) /
    length(primary_symbols) *
    100,
  NA
)


sens_overlap_pct <- ifelse(
  length(sens_symbols) > 0,
  length(shared_symbols) /
    length(sens_symbols) *
    100,
  NA
)


###############################################################
# 17. SUMMARY TABLE
###############################################################

comparison_summary <- data.frame(

  Analysis = c(
    "18_samples_primary",
    "17_samples_without_GSM5283438"
  ),

  AD_n = c(
    sum(
      metadata$Group == "AD"
    ),
    sum(
      metadata_sens$Group == "AD"
    )
  ),

  Control_n = c(
    sum(
      metadata$Group == "Control"
    ),
    sum(
      metadata_sens$Group == "Control"
    )
  ),

  Total_DEG = c(
    nrow(DEG),
    nrow(DEG_sens)
  ),

  Upregulated = c(
    sum(
      DEG$Direction == "Up"
    ),
    n_sens_up
  ),

  Downregulated = c(
    sum(
      DEG$Direction == "Down"
    ),
    n_sens_down
  )

)


write.csv(
  comparison_summary,
  file.path(
    sens_dir,
    "06_DEG_count_comparison.csv"
  ),
  row.names = FALSE
)


###############################################################
# 18. OVERLAP SUMMARY
###############################################################

overlap_summary <- data.frame(

  Primary_DEGs =
    length(
      primary_symbols
    ),

  Without_GSM5283438_DEGs =
    length(
      sens_symbols
    ),

  Shared_DEGs =
    length(
      shared_symbols
    ),

  Lost_after_removal =
    length(
      lost_symbols
    ),

  Newly_detected_after_removal =
    length(
      new_symbols
    ),

  Primary_overlap_percent =
    primary_overlap_pct,

  Sensitivity_overlap_percent =
    sens_overlap_pct,

  Jaccard_index =
    jaccard_index

)


write.csv(
  overlap_summary,
  file.path(
    sens_dir,
    "07_DEG_overlap_summary.csv"
  ),
  row.names = FALSE
)


###############################################################
# 19. SAVE SHARED / LOST / NEW DEGs
###############################################################

write.csv(
  data.frame(
    SYMBOL = shared_symbols
  ),
  file.path(
    sens_dir,
    "08_shared_DEGs.csv"
  ),
  row.names = FALSE
)


write.csv(
  data.frame(
    SYMBOL = lost_symbols
  ),
  file.path(
    sens_dir,
    "09_DEGs_lost_after_removal.csv"
  ),
  row.names = FALSE
)


write.csv(
  data.frame(
    SYMBOL = new_symbols
  ),
  file.path(
    sens_dir,
    "10_new_DEGs_after_removal.csv"
  ),
  row.names = FALSE
)


###############################################################
# 20. COMPARE ALL GENE logFC VALUES
###############################################################

logfc_compare <- merge(

  results_gene_unique[
    ,
    c(
      "SYMBOL",
      "logFC",
      "adj.P.Val"
    )
  ],

  res_sens_gene[
    ,
    c(
      "SYMBOL",
      "logFC",
      "adj.P.Val"
    )
  ],

  by = "SYMBOL",

  suffixes = c(
    "_18samples",
    "_17samples"
  )

)


###############################################################
# 21. CALCULATE logFC CORRELATION
###############################################################

logfc_cor_pearson <- cor(
  logfc_compare$logFC_18samples,
  logfc_compare$logFC_17samples,
  method = "pearson",
  use = "complete.obs"
)


logfc_cor_spearman <- cor(
  logfc_compare$logFC_18samples,
  logfc_compare$logFC_17samples,
  method = "spearman",
  use = "complete.obs"
)


cat(
  "\nlogFC Pearson correlation:",
  round(
    logfc_cor_pearson,
    4
  ),
  "\n"
)


cat(
  "logFC Spearman correlation:",
  round(
    logfc_cor_spearman,
    4
  ),
  "\n"
)


write.csv(
  logfc_compare,
  file.path(
    sens_dir,
    "11_gene_logFC_comparison.csv"
  ),
  row.names = FALSE
)


###############################################################
# 22. LOGFC CORRELATION PLOT
###############################################################

logfc_plot <- ggplot(
  logfc_compare,
  aes(
    x = logFC_18samples,
    y = logFC_17samples
  )
) +

  geom_point(
    alpha = 0.35,
    size = 1.2
  ) +

  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +

  labs(
    title =
      "GSE173954 effect-size sensitivity analysis",

    subtitle =
      paste0(
        "Pearson r = ",
        round(
          logfc_cor_pearson,
          3
        ),
        "; Spearman rho = ",
        round(
          logfc_cor_spearman,
          3
        )
      ),

    x =
      "log2FC: 18 samples",

    y =
      "log2FC: without GSM5283438"
  ) +

  theme_bw(
    base_size = 13
  )


ggsave(
  file.path(
    sens_dir,
    "12_logFC_correlation_18_vs_17_samples.tiff"
  ),
  logfc_plot,
  width = 7,
  height = 6,
  dpi = 300,
  compression = "lzw"
)


###############################################################
# 23. DIRECTION CONCORDANCE AMONG SHARED DEGs
###############################################################

shared_compare <- logfc_compare[
  logfc_compare$SYMBOL %in%
    shared_symbols,
  ,
  drop = FALSE
]


if (nrow(shared_compare) > 0) {

  shared_compare$Same_direction <-
    sign(
      shared_compare$logFC_18samples
    ) ==
    sign(
      shared_compare$logFC_17samples
    )


  direction_concordance <- mean(
    shared_compare$Same_direction
  ) * 100

} else {

  direction_concordance <- NA

}


###############################################################
# 24. COMPLETE SENSITIVITY SUMMARY
###############################################################

sensitivity_metrics <- data.frame(

  Original_sample_n =
    ncol(expr),

  Sensitivity_sample_n =
    ncol(expr_sens),

  Removed_sample =
    outlier_sample,

  Original_DEG_n =
    nrow(DEG),

  Sensitivity_DEG_n =
    nrow(DEG_sens),

  Shared_DEG_n =
    length(shared_symbols),

  Lost_DEG_n =
    length(lost_symbols),

  New_DEG_n =
    length(new_symbols),

  Jaccard_index =
    jaccard_index,

  Primary_DEG_overlap_percent =
    primary_overlap_pct,

  Sensitivity_DEG_overlap_percent =
    sens_overlap_pct,

  Shared_DEG_direction_concordance_percent =
    direction_concordance,

  All_gene_logFC_Pearson =
    logfc_cor_pearson,

  All_gene_logFC_Spearman =
    logfc_cor_spearman

)


write.csv(
  sensitivity_metrics,
  file.path(
    sens_dir,
    "13_FINAL_sensitivity_metrics.csv"
  ),
  row.names = FALSE
)


###############################################################
# 25. FINAL CONSOLE OUTPUT
###############################################################

cat("\n\n")
cat("=====================================================\n")
cat("SENSITIVITY ANALYSIS COMPLETED\n")
cat("=====================================================\n")

cat(
  "\n18-sample analysis DEGs:",
  nrow(DEG)
)

cat(
  "\n17-sample analysis DEGs:",
  nrow(DEG_sens)
)

cat(
  "\nShared DEGs:",
  length(shared_symbols)
)

cat(
  "\nLost after removal:",
  length(lost_symbols)
)

cat(
  "\nNew after removal:",
  length(new_symbols)
)

cat(
  "\nJaccard index:",
  round(
    jaccard_index,
    3
  )
)

cat(
  "\nPrimary DEG overlap:",
  round(
    primary_overlap_pct,
    1
  ),
  "%"
)

cat(
  "\nSensitivity DEG overlap:",
  round(
    sens_overlap_pct,
    1
  ),
  "%"
)

cat(
  "\nShared DEG direction concordance:",
  round(
    direction_concordance,
    1
  ),
  "%"
)

cat(
  "\nAll-gene logFC Pearson r:",
  round(
    logfc_cor_pearson,
    4
  )
)

cat(
  "\nAll-gene logFC Spearman rho:",
  round(
    logfc_cor_spearman,
    4
  )
)

cat(
  "\n\nResults saved to:\n",
  sens_dir,
  "\n"
)

cat("=====================================================\n")
