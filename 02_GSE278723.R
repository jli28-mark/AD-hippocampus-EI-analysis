#####################################################################
# GSE278723
# Human hippocampal subfields: CA1, CA2, CA3, CA4
#
# Comparison:
#   AD vs NC
#
# Excluded:
#   CAA
#   AD+CAA
#   EC / NX
#
# Each hippocampal subfield is analysed independently:
#   CA1: AD vs NC
#   CA2: AD vs NC
#   CA3: AD vs NC
#   CA4: AD vs NC
#
# Method:
#   DESeq2
#
# DEG threshold:
#   BH-adjusted P < 0.05
#   |log2FC| >= 0.58
#####################################################################


#####################################################################
# 0. CLEAN ENVIRONMENT
#####################################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)


#####################################################################
# 1. PATHS
#####################################################################

base_dir <- "D:/My Desktop/09.09返修/14_data/GSE278723"

count_file <- file.path(
  base_dir,
  "GSE278723_Counts.txt"
)

metadata_file <- file.path(
  base_dir,
  "GSE278723.metadata.xlsx"
)

out_dir <- file.path(
  base_dir,
  "GSE278723_analysis"
)

if (!dir.exists(out_dir)) {
  dir.create(
    out_dir,
    recursive = TRUE
  )
}


#####################################################################
# 2. INSTALL PACKAGES
#####################################################################

cran_packages <- c(
  "readxl",
  "data.table",
  "ggplot2",
  "pheatmap"
)

bioc_packages <- c(
  "DESeq2",
  "AnnotationDbi",
  "org.Hs.eg.db"
)


for (pkg in cran_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    install.packages(
      pkg,
      dependencies = TRUE
    )
  }
}


if (!requireNamespace("BiocManager", quietly = TRUE)) {

  install.packages("BiocManager")
}


for (pkg in bioc_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
  }
}


#####################################################################
# 3. LOAD PACKAGES
#####################################################################

suppressPackageStartupMessages({

  library(readxl)
  library(data.table)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(AnnotationDbi)
  library(org.Hs.eg.db)

})


#####################################################################
# 4. READ METADATA
#####################################################################

metadata <- read_excel(
  metadata_file
)

metadata <- as.data.frame(
  metadata
)

colnames(metadata) <- trimws(
  colnames(metadata)
)


cat("\n========================================\n")
cat("METADATA INFORMATION\n")
cat("========================================\n")

print(
  head(metadata)
)

print(
  colnames(metadata)
)


required_cols <- c(
  "title",
  "data",
  "Sample",
  "group",
  "Region"
)


if (!all(required_cols %in% colnames(metadata))) {

  stop(
    paste0(
      "Metadata must contain: ",
      paste(
        required_cols,
        collapse = ", "
      )
    )
  )
}


#####################################################################
# 5. CLEAN METADATA
#####################################################################

metadata$title <- trimws(
  as.character(metadata$title)
)

metadata$group <- trimws(
  as.character(metadata$group)
)

metadata$Region <- trimws(
  as.character(metadata$Region)
)


cat("\nDisease groups:\n")

print(
  table(metadata$group)
)


cat("\nBrain regions:\n")

print(
  table(metadata$Region)
)


cat("\nGroup x Region:\n")

print(
  table(
    metadata$group,
    metadata$Region
  )
)


#####################################################################
# 6. RETAIN ONLY AD AND NC
#
# Exclude:
# CAA
# AD+CAA
#####################################################################

metadata_use <- metadata[
  metadata$group %in% c(
    "AD",
    "NC"
  ),
  ,
  drop = FALSE
]


#####################################################################
# 7. RETAIN ONLY HIPPOCAMPAL CA SUBFIELDS
#
# Exclude EC
#####################################################################

regions_keep <- c(
  "CA1",
  "CA2",
  "CA3",
  "CA4"
)


metadata_use <- metadata_use[
  metadata_use$Region %in%
    regions_keep,
  ,
  drop = FALSE
]


cat("\n========================================\n")
cat("SAMPLES RETAINED FOR ANALYSIS\n")
cat("========================================\n")

print(
  table(
    metadata_use$group,
    metadata_use$Region
  )
)


#####################################################################
# 8. READ COUNT MATRIX
#####################################################################

counts_raw <- fread(
  count_file,
  data.table = FALSE,
  check.names = FALSE
)


cat("\n========================================\n")
cat("COUNT MATRIX\n")
cat("========================================\n")

cat(
  "Rows:",
  nrow(counts_raw),
  "\n"
)

cat(
  "Columns:",
  ncol(counts_raw),
  "\n"
)


cat("\nFirst 10 column names:\n")

print(
  colnames(counts_raw)[1:min(
    10,
    ncol(counts_raw)
  )]
)


#####################################################################
# 9. IDENTIFY GENE-ID COLUMN
#
# The first column is assumed to contain gene identifiers.
#####################################################################

gene_col <- colnames(
  counts_raw
)[1]


cat(
  "\nGene identifier column:",
  gene_col,
  "\n"
)


gene_ids <- as.character(
  counts_raw[[gene_col]]
)


#####################################################################
# 10. BUILD COUNT MATRIX
#####################################################################

count_matrix <- counts_raw[
  ,
  -1,
  drop = FALSE
]


rownames(count_matrix) <- gene_ids


#####################################################################
# 11. REMOVE EMPTY / DUPLICATED GENE IDS
#####################################################################

keep_gene <- !is.na(
  rownames(count_matrix)
) &
  rownames(count_matrix) != ""


count_matrix <- count_matrix[
  keep_gene,
  ,
  drop = FALSE
]


# If duplicated gene identifiers exist, sum them
if (anyDuplicated(
  rownames(count_matrix)
) > 0) {

  cat(
    "\nDuplicated gene IDs found; summing counts.\n"
  )

  temp <- data.frame(
    GeneID = rownames(count_matrix),
    count_matrix,
    check.names = FALSE
  )

  temp <- aggregate(
    . ~ GeneID,
    data = temp,
    FUN = sum
  )

  rownames(temp) <- temp$GeneID

  temp$GeneID <- NULL

  count_matrix <- temp
}


#####################################################################
# 12. CHECK WHETHER COUNTS ARE NUMERIC
#####################################################################

count_matrix[] <- lapply(
  count_matrix,
  as.numeric
)


if (anyNA(count_matrix)) {

  stop(
    paste0(
      "NA values were generated while converting counts to numeric. ",
      "Please inspect GSE278723_Counts.txt."
    )
  )
}


#####################################################################
# 13. CHECK WHETHER DATA LOOK LIKE RAW COUNTS
#####################################################################

non_integer_fraction <- mean(
  abs(
    as.matrix(count_matrix) -
      round(as.matrix(count_matrix))
  ) > 1e-8
)


cat(
  "\nFraction of non-integer entries:",
  non_integer_fraction,
  "\n"
)


if (non_integer_fraction > 0.001) {

  stop(
    paste0(
      "The matrix contains substantial non-integer values. ",
      "DESeq2 requires raw integer counts. ",
      "Please confirm that GSE278723_Counts.txt is a raw count matrix."
    )
  )
}


count_matrix <- round(
  as.matrix(count_matrix)
)


storage.mode(
  count_matrix
) <- "integer"


#####################################################################
# 14. NORMALIZE COUNT SAMPLE NAMES
#
# Count names:
#   01.CA1
#   01.CA2
#
# Metadata title:
#   01CA1
#   01CA2
#
# Therefore convert count names:
#   01.CA1 -> 01CA1
#
# Also:
#   .NX corresponds to EC and is excluded anyway.
#####################################################################

original_count_names <- colnames(
  count_matrix
)


normalized_count_names <- original_count_names


# NX -> EC for consistency
normalized_count_names <- sub(
  "\\.NX$",
  "EC",
  normalized_count_names,
  ignore.case = TRUE
)


# Remove dot before CA
normalized_count_names <- gsub(
  "\\.",
  "",
  normalized_count_names
)


colnames(
  count_matrix
) <- normalized_count_names


#####################################################################
# 15. SAMPLE MATCHING CHECK
#####################################################################

metadata_samples <- metadata_use$title


missing_in_counts <- setdiff(
  metadata_samples,
  colnames(count_matrix)
)


if (length(missing_in_counts) > 0) {

  cat("\n========================================\n")
  cat("WARNING: METADATA SAMPLES NOT FOUND\n")
  cat("========================================\n")

  print(
    missing_in_counts
  )

  stop(
    paste0(
      "Some required AD/NC CA1-CA4 samples cannot be matched ",
      "to the count matrix. Do not proceed until this is resolved."
    )
  )
}


#####################################################################
# 16. REMOVE ALL NON-REQUIRED COUNT COLUMNS
#####################################################################

count_matrix_use <- count_matrix[
  ,
  metadata_samples,
  drop = FALSE
]


metadata_use <- metadata_use[
  match(
    colnames(count_matrix_use),
    metadata_use$title
  ),
  ,
  drop = FALSE
]


stopifnot(
  all(
    colnames(count_matrix_use) ==
      metadata_use$title
  )
)


#####################################################################
# 17. SAVE FINAL SAMPLE SELECTION
#####################################################################

write.csv(
  metadata_use,
  file.path(
    out_dir,
    "GSE278723_samples_used_AD_vs_NC_CA1_CA4.csv"
  ),
  row.names = FALSE
)


sample_summary <- as.data.frame.matrix(
  table(
    metadata_use$Region,
    metadata_use$group
  )
)


write.csv(
  sample_summary,
  file.path(
    out_dir,
    "GSE278723_sample_number_summary.csv"
  )
)


cat("\nFinal sample numbers:\n")

print(
  table(
    metadata_use$group,
    metadata_use$Region
  )
)


#####################################################################
# 18. OPTIONAL GENE ANNOTATION FUNCTION
#
# If rownames are Ensembl gene IDs:
# ENSG00000123456.1 -> ENSG00000123456
#
# If the file already contains gene symbols, the original ID is retained.
#####################################################################

annotate_gene_ids <- function(result_df) {

  ids <- as.character(
    result_df$GeneID
  )


  ids_clean <- sub(
    "\\..*$",
    "",
    ids
  )


  is_ensembl <- mean(
    grepl(
      "^ENSG[0-9]+$",
      ids_clean
    )
  ) > 0.5


  if (is_ensembl) {

    cat(
      "\nEnsembl gene IDs detected. Mapping to gene symbols...\n"
    )


    anno <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = unique(ids_clean),
      keytype = "ENSEMBL",
      columns = c(
        "SYMBOL",
        "ENTREZID",
        "GENENAME"
      )
    )


    anno <- anno[
      !duplicated(
        anno$ENSEMBL
      ),
    ]


    map_symbol <- anno$SYMBOL[
      match(
        ids_clean,
        anno$ENSEMBL
      )
    ]


    map_entrez <- anno$ENTREZID[
      match(
        ids_clean,
        anno$ENSEMBL
      )
    ]


    map_name <- anno$GENENAME[
      match(
        ids_clean,
        anno$ENSEMBL
      )
    ]


    result_df$ENSEMBL <- ids_clean
    result_df$SYMBOL <- map_symbol
    result_df$ENTREZID <- map_entrez
    result_df$GENENAME <- map_name

  } else {

    result_df$SYMBOL <- ids

  }


  return(
    result_df
  )
}


#####################################################################
# 19. FUNCTION TO ANALYSE ONE HIPPOCAMPAL SUBFIELD
#####################################################################

run_region_DESeq2 <- function(region_name) {


  cat("\n\n")
  cat("====================================================\n")
  cat("ANALYSING ", region_name, "\n", sep = "")
  cat("AD vs NC\n")
  cat("====================================================\n")


  ###############################################################
  # 19.1 CREATE REGION FOLDER
  ###############################################################

  region_dir <- file.path(
    out_dir,
    region_name
  )


  if (!dir.exists(region_dir)) {

    dir.create(
      region_dir,
      recursive = TRUE
    )
  }


  ###############################################################
  # 19.2 SELECT METADATA
  ###############################################################

  meta_region <- metadata_use[
    metadata_use$Region ==
      region_name,
    ,
    drop = FALSE
  ]


  # Explicit reference group:
  # NC = reference
  # AD = comparison
  meta_region$group <- factor(
    meta_region$group,
    levels = c(
      "NC",
      "AD"
    )
  )


  cat("\nSample numbers:\n")

  print(
    table(
      meta_region$group
    )
  )


  ###############################################################
  # 19.3 SELECT COUNT MATRIX
  ###############################################################

  counts_region <- count_matrix_use[
    ,
    meta_region$title,
    drop = FALSE
  ]


  stopifnot(
    all(
      colnames(counts_region) ==
        meta_region$title
    )
  )


  rownames(
    meta_region
  ) <- meta_region$title


  ###############################################################
  # 19.4 REMOVE GENES WITH ZERO COUNTS ACROSS ALL SAMPLES
  #
  # This is not a DEG-based selection.
  ###############################################################

  counts_region <- counts_region[
    rowSums(counts_region) > 0,
    ,
    drop = FALSE
  ]


  cat(
    "Genes with at least one read:",
    nrow(counts_region),
    "\n"
  )


  ###############################################################
  # 19.5 CREATE DESEQ2 OBJECT
  ###############################################################

  dds <- DESeqDataSetFromMatrix(

    countData =
      counts_region,

    colData =
      meta_region,

    design =
      ~ group

  )


  ###############################################################
  # 19.6 DESEQ2
  ###############################################################

  dds <- DESeq(
    dds
  )


  ###############################################################
  # 19.7 NORMALIZED COUNTS
  ###############################################################

  normalized_counts <- counts(
    dds,
    normalized = TRUE
  )


  write.csv(
    normalized_counts,
    file.path(
      region_dir,
      paste0(
        "GSE278723_",
        region_name,
        "_normalized_counts.csv"
      )
    )
  )


  ###############################################################
  # 19.8 DIFFERENTIAL EXPRESSION
  #
  # AD vs NC
  ###############################################################

  res <- results(
    dds,
    contrast = c(
      "group",
      "AD",
      "NC"
    ),
    alpha = 0.05
  )


  res_df <- as.data.frame(
    res
  )


  res_df$GeneID <- rownames(
    res_df
  )


  rownames(
    res_df
  ) <- NULL


  ###############################################################
  # 19.9 ANNOTATION
  ###############################################################

  res_df <- annotate_gene_ids(
    res_df
  )


  ###############################################################
  # 19.10 DIRECTION
  ###############################################################

  res_df$Direction <- "Not_significant"


  res_df$Direction[
    !is.na(res_df$padj) &
      res_df$padj < 0.05 &
      res_df$log2FoldChange >= 0.58
  ] <- "Up"


  res_df$Direction[
    !is.na(res_df$padj) &
      res_df$padj < 0.05 &
      res_df$log2FoldChange <= -0.58
  ] <- "Down"


  ###############################################################
  # 19.11 DEG TABLE
  ###############################################################

  DEG_region <- res_df[
    !is.na(res_df$padj) &
      res_df$padj < 0.05 &
      abs(
        res_df$log2FoldChange
      ) >= 0.58,
    ,
    drop = FALSE
  ]


  DEG_region <- DEG_region[
    order(
      DEG_region$padj
    ),
    ,
    drop = FALSE
  ]


  ###############################################################
  # 19.12 SAVE ALL RESULTS
  ###############################################################

  res_df <- res_df[
    order(
      res_df$padj
    ),
    ,
    drop = FALSE
  ]


  write.csv(
    res_df,
    file.path(
      region_dir,
      paste0(
        "GSE278723_",
        region_name,
        "_all_genes_AD_vs_NC.csv"
      )
    ),
    row.names = FALSE
  )


  write.csv(
    DEG_region,
    file.path(
      region_dir,
      paste0(
        "GSE278723_",
        region_name,
        "_DEG_AD_vs_NC.csv"
      )
    ),
    row.names = FALSE
  )


  write.csv(
    DEG_region[
      DEG_region$Direction == "Up",
    ],
    file.path(
      region_dir,
      paste0(
        "GSE278723_",
        region_name,
        "_DEG_UP.csv"
      )
    ),
    row.names = FALSE
  )


  write.csv(
    DEG_region[
      DEG_region$Direction == "Down",
    ],
    file.path(
      region_dir,
      paste0(
        "GSE278723_",
        region_name,
        "_DEG_DOWN.csv"
      )
    ),
    row.names = FALSE
  )


  ###############################################################
  # 19.13 PCA USING VST
  ###############################################################

  vst_data <- vst(
    dds,
    blind = TRUE
  )


  pca_data <- plotPCA(
    vst_data,
    intgroup = "group",
    returnData = TRUE
  )


  percentVar <- round(
    100 *
      attr(
        pca_data,
        "percentVar"
      )
  )


  pca_data$Sample <- rownames(
    pca_data
  )


  p_pca <- ggplot(
    pca_data,
    aes(
      x = PC1,
      y = PC2,
      shape = group
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
      vjust = -0.7
    ) +

    xlab(
      paste0(
        "PC1: ",
        percentVar[1],
        "% variance"
      )
    ) +

    ylab(
      paste0(
        "PC2: ",
        percentVar[2],
        "% variance"
      )
    ) +

    ggtitle(
      paste0(
        "GSE278723 ",
        region_name,
        ": AD vs NC"
      )
    ) +

    theme_bw(
      base_size = 13
    )


  ggsave(
    file.path(
      region_dir,
      paste0(
        "GSE278723_",
        region_name,
        "_PCA.tiff"
      )
    ),
    p_pca,
    width = 8,
    height = 6,
    dpi = 300,
    compression = "lzw"
  )


  ###############################################################
  # 19.14 SAMPLE CORRELATION
  ###############################################################

  vst_matrix <- assay(
    vst_data
  )


  cor_matrix <- cor(
    vst_matrix,
    method = "pearson"
  )


  annotation_col <- data.frame(
    Group = meta_region$group
  )


  rownames(
    annotation_col
  ) <- meta_region$title


  tiff(
    file.path(
      region_dir,
      paste0(
        "GSE278723_",
        region_name,
        "_SampleCorrelation.tiff"
      )
    ),
    width = 3000,
    height = 2800,
    res = 300,
    compression = "lzw"
  )


  pheatmap(
    cor_matrix,
    annotation_col =
      annotation_col,
    annotation_row =
      annotation_col,
    main = paste0(
      "GSE278723 ",
      region_name,
      " sample correlation"
    )
  )


  dev.off()


  ###############################################################
  # 19.15 HIERARCHICAL CLUSTERING
  ###############################################################

  sample_dist <- dist(
    t(vst_matrix)
  )


  hc <- hclust(
    sample_dist
  )


  tiff(
    file.path(
      region_dir,
      paste0(
        "GSE278723_",
        region_name,
        "_HierarchicalClustering.tiff"
      )
    ),
    width = 3000,
    height = 2200,
    res = 300,
    compression = "lzw"
  )


  plot(
    hc,
    labels = paste0(
      meta_region$title,
      "_",
      meta_region$group
    ),
    main = paste0(
      "GSE278723 ",
      region_name
    ),
    xlab = "",
    sub = "",
    cex = 0.8
  )


  dev.off()


  ###############################################################
  # 19.16 REGION SUMMARY
  ###############################################################

  summary_region <- data.frame(

    Dataset =
      "GSE278723",

    Region =
      region_name,

    AD_n =
      sum(
        meta_region$group == "AD"
      ),

    NC_n =
      sum(
        meta_region$group == "NC"
      ),

    Total_DEG =
      nrow(
        DEG_region
      ),

    Up =
      sum(
        DEG_region$Direction == "Up"
      ),

    Down =
      sum(
        DEG_region$Direction == "Down"
      ),

    padj_cutoff =
      0.05,

    abs_log2FC_cutoff =
      0.58

  )


  write.csv(
    summary_region,
    file.path(
      region_dir,
      paste0(
        "GSE278723_",
        region_name,
        "_summary.csv"
      )
    ),
    row.names = FALSE
  )


  cat("\n")
  print(
    summary_region
  )


  ###############################################################
  # RETURN RESULTS
  ###############################################################

  return(
    list(

      all =
        res_df,

      DEG =
        DEG_region,

      summary =
        summary_region,

      dds =
        dds

    )
  )
}


#####################################################################
# 20. RUN FOUR HIPPOCAMPAL SUBFIELDS
#####################################################################

result_CA1 <- run_region_DESeq2(
  "CA1"
)

result_CA2 <- run_region_DESeq2(
  "CA2"
)

result_CA3 <- run_region_DESeq2(
  "CA3"
)

result_CA4 <- run_region_DESeq2(
  "CA4"
)


#####################################################################
# 21. COMBINE REGION SUMMARIES
#####################################################################

summary_all <- rbind(

  result_CA1$summary,

  result_CA2$summary,

  result_CA3$summary,

  result_CA4$summary

)


write.csv(
  summary_all,
  file.path(
    out_dir,
    "GSE278723_FINAL_region_DEG_summary.csv"
  ),
  row.names = FALSE
)


cat("\n\n")
cat("====================================================\n")
cat("FINAL REGION-SPECIFIC DEG SUMMARY\n")
cat("====================================================\n")

print(
  summary_all
)


#####################################################################
# 22. CREATE FOUR-REGION GENE-LEVEL EFFECT SUMMARY
#
# IMPORTANT:
# This does NOT merge samples or recompute DEGs.
#
# It simply places CA1-CA4 effect estimates side-by-side.
#####################################################################

extract_region <- function(
    result_object,
    region_name) {


  x <- result_object$all[
    ,
    c(
      "GeneID",
      "SYMBOL",
      "log2FoldChange",
      "pvalue",
      "padj",
      "Direction"
    )
  ]


  colnames(x) <- c(

    "GeneID",

    paste0(
      region_name,
      "_SYMBOL"
    ),

    paste0(
      region_name,
      "_log2FC"
    ),

    paste0(
      region_name,
      "_P"
    ),

    paste0(
      region_name,
      "_FDR"
    ),

    paste0(
      region_name,
      "_Direction"
    )

  )


  return(x)
}


CA1_table <- extract_region(
  result_CA1,
  "CA1"
)

CA2_table <- extract_region(
  result_CA2,
  "CA2"
)

CA3_table <- extract_region(
  result_CA3,
  "CA3"
)

CA4_table <- extract_region(
  result_CA4,
  "CA4"
)


region_effects <- Reduce(
  function(x, y) {
    merge(
      x,
      y,
      by = "GeneID",
      all = TRUE
    )
  },
  list(
    CA1_table,
    CA2_table,
    CA3_table,
    CA4_table
  )
)


#####################################################################
# 23. CREATE ONE CONSENSUS SYMBOL COLUMN
#####################################################################

symbol_cols <- grep(
  "_SYMBOL$",
  colnames(region_effects),
  value = TRUE
)


region_effects$SYMBOL <- apply(
  region_effects[
    ,
    symbol_cols,
    drop = FALSE
  ],
  1,
  function(x) {

    x <- x[
      !is.na(x) &
        x != ""
    ]

    if (length(x) == 0) {
      return(NA)
    }

    return(x[1])
  }
)


#####################################################################
#####################################################################
# 24. COUNT NUMBER OF SIGNIFICANT REGIONS FOR EACH GENE
#####################################################################

region_names <- c(
  "CA1",
  "CA2",
  "CA3",
  "CA4"
)

sig_matrix <- sapply(
  region_names,
  function(r) {

    fdr <- region_effects[[paste0(r, "_FDR")]]
    fc  <- region_effects[[paste0(r, "_log2FC")]]

    !is.na(fdr) &
      fdr < 0.05 &
      !is.na(fc) &
      abs(fc) >= 0.58
  }
)

region_effects$N_significant_regions <- rowSums(
  sig_matrix
)


#####################################################################
# 25. COUNT UP / DOWN REGIONS
#####################################################################

up_matrix <- sapply(
  region_names,
  function(r) {

    fdr <- region_effects[[paste0(r, "_FDR")]]
    fc  <- region_effects[[paste0(r, "_log2FC")]]

    !is.na(fdr) &
      fdr < 0.05 &
      !is.na(fc) &
      fc >= 0.58
  }
)

down_matrix <- sapply(
  region_names,
  function(r) {

    fdr <- region_effects[[paste0(r, "_FDR")]]
    fc  <- region_effects[[paste0(r, "_log2FC")]]

    !is.na(fdr) &
      fdr < 0.05 &
      !is.na(fc) &
      fc <= -0.58
  }
)

region_effects$N_up_regions <- rowSums(
  up_matrix
)

region_effects$N_down_regions <- rowSums(
  down_matrix
)


#####################################################################
# 26. DEFINE REGIONAL CONSISTENCY
#####################################################################

region_effects$Regional_pattern <- "Not_consensus"

# Significant in >=2 regions and consistently UP
region_effects$Regional_pattern[
  region_effects$N_up_regions >= 2 &
    region_effects$N_down_regions == 0
] <- "Consistent_Up"

# Significant in >=2 regions and consistently DOWN
region_effects$Regional_pattern[
  region_effects$N_down_regions >= 2 &
    region_effects$N_up_regions == 0
] <- "Consistent_Down"

# Significant in opposite directions across regions
region_effects$Regional_pattern[
  region_effects$N_up_regions >= 1 &
    region_effects$N_down_regions >= 1
] <- "Discordant"

#####################################################################
# 27. SAVE COMPLETE FOUR-REGION MATRIX
#####################################################################

write.csv(
  region_effects,
  file.path(
    out_dir,
    "GSE278723_CA1_CA4_gene_effect_summary.csv"
  ),
  row.names = FALSE
)


#####################################################################
# 28. OPTIONAL CONSENSUS DEG SET
#
# Definition:
# Significant in >=2 hippocampal subfields
# AND same direction.
#
# IMPORTANT:
# Keep this separate from the four primary regional DEG analyses.
#####################################################################

consensus_DEG <- region_effects[
  region_effects$Regional_pattern %in%
    c(
      "Consistent_Up",
      "Consistent_Down"
    ),
  ,
  drop = FALSE
]


write.csv(
  consensus_DEG,
  file.path(
    out_dir,
    "GSE278723_consensus_DEG_at_least_2_regions_same_direction.csv"
  ),
  row.names = FALSE
)


#####################################################################
# 29. SAVE DISCORDANT GENES
#####################################################################

discordant_DEG <- region_effects[
  region_effects$Regional_pattern ==
    "Discordant",
  ,
  drop = FALSE
]


write.csv(
  discordant_DEG,
  file.path(
    out_dir,
    "GSE278723_regionally_discordant_DEGs.csv"
  ),
  row.names = FALSE
)


#####################################################################
# 30. CONSENSUS SUMMARY
#####################################################################

consensus_summary <- data.frame(

  Dataset =
    "GSE278723",

  CA1_DEG =
    nrow(
      result_CA1$DEG
    ),

  CA2_DEG =
    nrow(
      result_CA2$DEG
    ),

  CA3_DEG =
    nrow(
      result_CA3$DEG
    ),

  CA4_DEG =
    nrow(
      result_CA4$DEG
    ),

  Consensus_DEG_2plus_regions =
    nrow(
      consensus_DEG
    ),

  Consensus_UP =
    sum(
      consensus_DEG$Regional_pattern ==
        "Consistent_Up"
    ),

  Consensus_DOWN =
    sum(
      consensus_DEG$Regional_pattern ==
        "Consistent_Down"
    ),

  Discordant_genes =
    nrow(
      discordant_DEG
    )

)


write.csv(
  consensus_summary,
  file.path(
    out_dir,
    "GSE278723_FINAL_consensus_summary.csv"
  ),
  row.names = FALSE
)


#####################################################################
# 31. SAVE SESSION INFORMATION
#####################################################################

sink(
  file.path(
    out_dir,
    "GSE278723_sessionInfo.txt"
  )
)

sessionInfo()

sink()


#####################################################################
# 32. SAVE IMPORTANT R OBJECTS
#####################################################################

save(
  result_CA1,
  result_CA2,
  result_CA3,
  result_CA4,
  summary_all,
  region_effects,
  consensus_DEG,
  file = file.path(
    out_dir,
    "GSE278723_analysis_objects.RData"
  )
)


#####################################################################
# 33. FINAL OUTPUT
#####################################################################

cat("\n\n")
cat("====================================================\n")
cat("GSE278723 ANALYSIS COMPLETED\n")
cat("====================================================\n")

print(
  summary_all
)

cat("\nConsensus summary:\n")

print(
  consensus_summary
)

cat(
  "\nResults saved to:\n",
  out_dir,
  "\n"
)

cat("====================================================\n")
