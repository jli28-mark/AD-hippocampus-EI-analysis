############################################################
# ==========================================================
# Leave-one-study-out dataset contribution analysis
# Existing 10 KEGG pathways
#
# IMPORTANT:
# This analysis DOES NOT rerun KEGG enrichment.
#
# Aim:
# For each of the 10 predefined pathways,
# quantify how many pathway-associated DEGs are lost
# when each study is removed from the four-study union.
# ==========================================================
############################################################


rm(list = ls())
gc()

options(stringsAsFactors = FALSE)


############################################################
# 1. Packages
############################################################

cran_pkgs <- c(
  "readxl",
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2",
  "openxlsx"
)

for (pkg in cran_pkgs) {

  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }

}


suppressPackageStartupMessages({

  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(openxlsx)

})


############################################################
# 2. File paths
############################################################

deg_file <- paste0(
  "D:/My Desktop/09.09返修/Union gene/",
  "Four_studies_DEG_filtered_cleaned_union.xlsx"
)


kegg_file <- paste0(
  "D:/My Desktop/09.09返修/KEGG/Full_union_KEGG/",
  "Full_union_KEGG_pathway_gene_long_table.csv"
)


output_dir <- paste0(
  "D:/My Desktop/09.09返修/KEGG/",
  "Leave_one_study_out_KEGG"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 3. Check files
############################################################

if (!file.exists(deg_file)) {

  stop(
    paste0(
      "DEG file not found:\n",
      deg_file
    )
  )

}


if (!file.exists(kegg_file)) {

  stop(
    paste0(
      "KEGG pathway-gene file not found:\n",
      kegg_file
    )
  )

}


############################################################
# 4. Four studies
############################################################

studies <- c(
  "GSE278723",
  "GSE173955",
  "GSE48350",
  "GSE5281"
)


############################################################
# 5. Fixed 10 pathways
############################################################

target_pathways <- c(

  "Retrograde endocannabinoid signaling",

  "Synaptic vesicle cycle",

  "Dopaminergic synapse",

  "Long-term potentiation",

  "Calcium signaling pathway",

  "Glutamatergic synapse",

  "Cholinergic synapse",

  "GABAergic synapse",

  "Long-term depression",

  "Neuroactive ligand-receptor interaction"

)


############################################################
# 6. Read KEGG pathway-gene long table
############################################################

kegg_long <- readr::read_csv(
  kegg_file,
  show_col_types = FALSE
)


cat(
  "\n========================================\n",
  "KEGG TABLE\n",
  "========================================\n"
)

print(
  names(kegg_long)
)


############################################################
# Check required columns
############################################################

required_kegg_cols <- c(
  "Description",
  "SYMBOL"
)


missing_kegg_cols <- setdiff(
  required_kegg_cols,
  names(kegg_long)
)


if (length(missing_kegg_cols) > 0) {

  stop(
    paste0(
      "Missing KEGG columns: ",
      paste(
        missing_kegg_cols,
        collapse = ", "
      )
    )
  )

}


############################################################
# 7. Extract only the 10 predefined pathways
############################################################

target_genes <- kegg_long %>%

  dplyr::filter(
    Description %in% target_pathways
  ) %>%

  dplyr::filter(
    !is.na(SYMBOL),
    SYMBOL != ""
  ) %>%

  dplyr::transmute(
    Description = as.character(Description),
    Gene = as.character(SYMBOL)
  ) %>%

  dplyr::distinct(
    Description,
    Gene
  )


############################################################
# Check pathways found
############################################################

pathway_check <- data.frame(
  Description = target_pathways,
  stringsAsFactors = FALSE
) %>%

  dplyr::mutate(

    Found = ifelse(
      Description %in%
        unique(
          target_genes$Description
        ),
      "Yes",
      "No"
    )

  )


cat(
  "\n========================================\n",
  "TARGET PATHWAYS FOUND\n",
  "========================================\n"
)

print(
  pathway_check
)


cat(
  "\nNumber of pathways found:",
  dplyr::n_distinct(
    target_genes$Description
  ),
  "\n"
)


cat(
  "Unique genes across 10 pathways:",
  dplyr::n_distinct(
    target_genes$Gene
  ),
  "\n"
)


############################################################
# 8. Read DEG workbook
############################################################

sheet_names <- readxl::excel_sheets(
  deg_file
)


cat(
  "\n========================================\n",
  "EXCEL SHEETS\n",
  "========================================\n"
)

print(
  sheet_names
)


if (!"All_DEGs_long" %in% sheet_names) {

  stop(
    "Sheet 'All_DEGs_long' was not found."
  )

}


deg_raw <- readxl::read_excel(
  deg_file,
  sheet = "All_DEGs_long"
)


deg_raw <- as.data.frame(
  deg_raw
)


cat(
  "\nColumns in All_DEGs_long:\n"
)

print(
  names(deg_raw)
)


############################################################
# 9. Auto-detect DEG columns
############################################################

find_col <- function(
    df,
    candidates
) {

  hit <- candidates[
    candidates %in%
      names(df)
  ]


  if (length(hit) == 0) {
    return(NA_character_)
  }


  hit[1]

}


gene_col <- find_col(
  deg_raw,
  c(
    "Gene",
    "gene",
    "SYMBOL",
    "Symbol",
    "symbol",
    "Gene.symbol"
  )
)


dataset_col <- find_col(
  deg_raw,
  c(
    "Dataset",
    "dataset",
    "Study",
    "study"
  )
)


logfc_col <- find_col(
  deg_raw,
  c(
    "logFC",
    "log2FoldChange",
    "log2FC",
    "Log2FC"
  )
)


padj_col <- find_col(
  deg_raw,
  c(
    "padj",
    "adj.P.Val",
    "FDR",
    "adjusted_p"
  )
)


cat(
  "\nDetected DEG columns:\n"
)

cat(
  "Gene    :",
  gene_col,
  "\n"
)

cat(
  "Dataset :",
  dataset_col,
  "\n"
)

cat(
  "logFC   :",
  logfc_col,
  "\n"
)

cat(
  "padj    :",
  padj_col,
  "\n"
)


if (
  is.na(gene_col) ||
    is.na(dataset_col)
) {

  stop(
    paste0(
      "Could not identify Gene or Dataset columns.\n",
      "Available columns:\n",
      paste(
        names(deg_raw),
        collapse = ", "
      )
    )
  )

}


############################################################
# 10. Standardize DEG data
############################################################

deg <- data.frame(

  Gene = trimws(
    as.character(
      deg_raw[[gene_col]]
    )
  ),

  Dataset = trimws(
    as.character(
      deg_raw[[dataset_col]]
    )
  ),

  stringsAsFactors = FALSE

)


############################################################
# Add logFC if available
############################################################

if (!is.na(logfc_col)) {

  deg$logFC <- suppressWarnings(
    as.numeric(
      deg_raw[[logfc_col]]
    )
  )

} else {

  deg$logFC <- NA_real_

}


############################################################
# Add padj if available
############################################################

if (!is.na(padj_col)) {

  deg$padj <- suppressWarnings(
    as.numeric(
      deg_raw[[padj_col]]
    )
  )

} else {

  deg$padj <- NA_real_

}


############################################################
# 11. Correct known Excel-converted gene symbols
############################################################

gene_fix <- c(

  "1-Sep" = "SEPTIN1",

  "2-Mar" = "MARCHF2",

  "6-Mar" = "MARCHF6",

  "7-Sep" = "SEPTIN7",

  "9-Sep" = "SEPTIN9",

  "15-Sep" = "SELENOF"

)


deg$Gene <- ifelse(

  deg$Gene %in%
    names(gene_fix),

  unname(
    gene_fix[
      deg$Gene
    ]
  ),

  deg$Gene

)


############################################################
# 12. Basic cleaning
############################################################

deg <- deg %>%

  dplyr::filter(

    !is.na(Gene),

    Gene != "",

    !is.na(Dataset),

    Dataset != "",

    Dataset %in% studies

  )


############################################################
# IMPORTANT:
#
# All_DEGs_long should already contain only genes that
# met the DEG criterion in each individual dataset.
#
# If logFC/padj exist, apply the unified DEG threshold again.
############################################################

if (
  !all(
    is.na(
      deg$logFC
    )
  ) &&
    !all(
      is.na(
        deg$padj
      )
    )
) {

  deg <- deg %>%

    dplyr::filter(

      !is.na(logFC),

      !is.na(padj),

      padj < 0.05,

      abs(logFC) >= 0.58

    )

}


############################################################
# Remove multiple-gene annotations if still present
############################################################

deg <- deg %>%

  dplyr::filter(
    !grepl(
      "///",
      Gene,
      fixed = TRUE
    )
  )


############################################################
# 13. Deduplicate within each dataset
############################################################

if (
  !all(
    is.na(
      deg$padj
    )
  )
) {

  deg <- deg %>%

    dplyr::arrange(
      Dataset,
      Gene,
      padj,
      dplyr::desc(
        abs(logFC)
      )
    ) %>%

    dplyr::distinct(
      Dataset,
      Gene,
      .keep_all = TRUE
    )

} else {

  deg <- deg %>%

    dplyr::distinct(
      Dataset,
      Gene,
      .keep_all = TRUE
    )

}


############################################################
# 14. DEG count QC
############################################################

deg_count <- deg %>%

  dplyr::count(
    Dataset,
    name = "N_DEGs"
  )


cat(
  "\n========================================\n",
  "DEG COUNT PER DATASET\n",
  "========================================\n"
)

print(
  deg_count
)


############################################################
# 15. Create gene × study DEG membership matrix
#
# 1 = gene is a DEG in that dataset
# 0 = gene is not a DEG in that dataset
############################################################

gene_study <- deg %>%

  dplyr::mutate(
    Present = 1L
  ) %>%

  dplyr::select(
    Gene,
    Dataset,
    Present
  ) %>%

  tidyr::pivot_wider(

    names_from = Dataset,

    values_from = Present,

    values_fill = 0

  )


############################################################
# Ensure all four study columns exist
############################################################

for (s in studies) {

  if (!s %in% names(gene_study)) {

    gene_study[[s]] <- 0L

  }

}


############################################################
# 16. Merge pathway membership with dataset DEG membership
############################################################

pathway_gene_study <- target_genes %>%

  dplyr::left_join(
    gene_study,
    by = "Gene"
  )


############################################################
# Replace missing study membership with 0
############################################################

for (s in studies) {

  pathway_gene_study[[s]][
    is.na(
      pathway_gene_study[[s]]
    )
  ] <- 0L

}


############################################################
# 17. QC:
# Every KEGG pathway-associated union DEG should ideally
# appear as a DEG in at least one of the four datasets.
############################################################

pathway_gene_study$N_studies_DEG <- rowSums(

  pathway_gene_study[
    ,
    studies,
    drop = FALSE
  ]

)


n_zero <- sum(
  pathway_gene_study$N_studies_DEG == 0
)


cat(
  "\nPathway genes not found as DEG in any dataset:",
  n_zero,
  "\n"
)


############################################################
# Save any problematic genes
############################################################

zero_deg_genes <- pathway_gene_study %>%

  dplyr::filter(
    N_studies_DEG == 0
  )


write.csv(

  zero_deg_genes,

  file.path(
    output_dir,
    "QC_pathway_genes_not_found_in_DEG_tables.csv"
  ),

  row.names = FALSE

)


############################################################
# 18. Full pathway gene count
############################################################

full_counts <- pathway_gene_study %>%

  dplyr::group_by(
    Description
  ) %>%

  dplyr::summarise(

    Full_N_genes =
      dplyr::n_distinct(
        Gene
      ),

    .groups = "drop"

  )


############################################################
# 19. Leave-one-study-out contribution analysis
############################################################

loso_list <- list()


for (leaveout in studies) {


  remaining_studies <- setdiff(
    studies,
    leaveout
  )


  temp <- pathway_gene_study


  ##########################################################
  # Is the gene present in at least one remaining dataset?
  ##########################################################

  temp$Remaining_after_leaveout <- rowSums(

    temp[
      ,
      remaining_studies,
      drop = FALSE
    ]

  ) > 0


  ##########################################################
  # Is the gene present in the omitted dataset?
  ##########################################################

  temp$Present_in_leftout <- (
    temp[[leaveout]] == 1
  )


  ##########################################################
  # Unique-to-leftout:
  # present in omitted dataset
  # AND absent from all remaining datasets
  ##########################################################

  temp$Unique_to_leftout <- (

    temp$Present_in_leftout &

      !temp$Remaining_after_leaveout

  )


  ##########################################################
  # Number of studies containing each gene
  ##########################################################

  temp$N_studies_before_leaveout <- rowSums(

    temp[
      ,
      studies,
      drop = FALSE
    ]

  )


  ##########################################################
  # Summarise per pathway
  ##########################################################

  temp_summary <- temp %>%

    dplyr::group_by(
      Description
    ) %>%

    dplyr::summarise(

      Genes_present_in_leftout =
        sum(
          Present_in_leftout
        ),

      Remaining_N_genes =
        sum(
          Remaining_after_leaveout
        ),

      Unique_contribution =
        sum(
          Unique_to_leftout
        ),

      Genes_shared_with_other_studies =
        sum(
          Present_in_leftout &
            Remaining_after_leaveout
        ),

      .groups = "drop"

    ) %>%

    dplyr::left_join(
      full_counts,
      by = "Description"
    ) %>%

    dplyr::mutate(

      Dataset_removed = leaveout,

      Lost_N_genes =
        Full_N_genes -
        Remaining_N_genes,

      Retention_percent =
        100 *
        Remaining_N_genes /
        Full_N_genes,

      Loss_percent =
        100 *
        Lost_N_genes /
        Full_N_genes,

      Unique_contribution_percent =
        100 *
        Unique_contribution /
        Full_N_genes

    )


  loso_list[[leaveout]] <- temp_summary

}


############################################################
# 20. Combine LOSO results
############################################################

loso_summary <- dplyr::bind_rows(
  loso_list
)


############################################################
# Set factor order
############################################################

loso_summary$Description <- factor(

  loso_summary$Description,

  levels = target_pathways

)


loso_summary$Dataset_removed <- factor(

  loso_summary$Dataset_removed,

  levels = studies

)


loso_summary <- loso_summary %>%

  dplyr::arrange(
    Description,
    Dataset_removed
  )


############################################################
# 21. Determine most influential dataset for each pathway
#
# Primary criterion:
# largest number of genes lost
#
# Tie:
# larger unique contribution
############################################################

most_influential <- loso_summary %>%

  dplyr::group_by(
    Description
  ) %>%

  dplyr::arrange(

    dplyr::desc(
      Lost_N_genes
    ),

    dplyr::desc(
      Unique_contribution
    ),

    .by_group = TRUE

  ) %>%

  dplyr::slice(
    1
  ) %>%

  dplyr::ungroup()


############################################################
# 22. Create matrices
############################################################

retention_matrix <- loso_summary %>%

  dplyr::select(
    Description,
    Dataset_removed,
    Retention_percent
  ) %>%

  tidyr::pivot_wider(
    names_from = Dataset_removed,
    values_from = Retention_percent
  )


loss_matrix <- loso_summary %>%

  dplyr::select(
    Description,
    Dataset_removed,
    Loss_percent
  ) %>%

  tidyr::pivot_wider(
    names_from = Dataset_removed,
    values_from = Loss_percent
  )


lost_gene_matrix <- loso_summary %>%

  dplyr::select(
    Description,
    Dataset_removed,
    Lost_N_genes
  ) %>%

  tidyr::pivot_wider(
    names_from = Dataset_removed,
    values_from = Lost_N_genes
  )


unique_matrix <- loso_summary %>%

  dplyr::select(
    Description,
    Dataset_removed,
    Unique_contribution
  ) %>%

  tidyr::pivot_wider(
    names_from = Dataset_removed,
    values_from = Unique_contribution
  )


############################################################
# 23. Dataset overall contribution across 10 pathways
############################################################

dataset_summary <- loso_summary %>%

  dplyr::group_by(
    Dataset_removed
  ) %>%

  dplyr::summarise(

    Total_pathway_gene_occurrences_lost =
      sum(
        Lost_N_genes
      ),

    Mean_loss_percent =
      mean(
        Loss_percent
      ),

    Median_loss_percent =
      median(
        Loss_percent
      ),

    Total_unique_contribution =
      sum(
        Unique_contribution
      ),

    Mean_retention_percent =
      mean(
        Retention_percent
      ),

    .groups = "drop"

  ) %>%

  dplyr::arrange(
    dplyr::desc(
      Mean_loss_percent
    )
  )


############################################################
# 24. Save CSV files
############################################################

write.csv(

  loso_summary,

  file.path(
    output_dir,
    "LOSO_10_pathways_dataset_contribution.csv"
  ),

  row.names = FALSE

)


write.csv(

  most_influential,

  file.path(
    output_dir,
    "LOSO_most_influential_dataset_per_pathway.csv"
  ),

  row.names = FALSE

)


write.csv(

  dataset_summary,

  file.path(
    output_dir,
    "LOSO_overall_dataset_contribution.csv"
  ),

  row.names = FALSE

)


write.csv(

  pathway_gene_study,

  file.path(
    output_dir,
    "LOSO_pathway_gene_study_membership.csv"
  ),

  row.names = FALSE

)


############################################################
# 25. Excel workbook
############################################################

excel_file <- file.path(
  output_dir,
  "LOSO_10_pathways_dataset_contribution.xlsx"
)


wb <- openxlsx::createWorkbook()


openxlsx::addWorksheet(
  wb,
  "LOSO_summary"
)

openxlsx::writeData(
  wb,
  "LOSO_summary",
  loso_summary
)


openxlsx::addWorksheet(
  wb,
  "Most_influential"
)

openxlsx::writeData(
  wb,
  "Most_influential",
  most_influential
)


openxlsx::addWorksheet(
  wb,
  "Dataset_summary"
)

openxlsx::writeData(
  wb,
  "Dataset_summary",
  dataset_summary
)


openxlsx::addWorksheet(
  wb,
  "Retention_matrix"
)

openxlsx::writeData(
  wb,
  "Retention_matrix",
  retention_matrix
)


openxlsx::addWorksheet(
  wb,
  "Loss_percent_matrix"
)

openxlsx::writeData(
  wb,
  "Loss_percent_matrix",
  loss_matrix
)


openxlsx::addWorksheet(
  wb,
  "Lost_gene_matrix"
)

openxlsx::writeData(
  wb,
  "Lost_gene_matrix",
  lost_gene_matrix
)


openxlsx::addWorksheet(
  wb,
  "Unique_contribution"
)

openxlsx::writeData(
  wb,
  "Unique_contribution",
  unique_matrix
)


openxlsx::addWorksheet(
  wb,
  "Pathway_gene_study"
)

openxlsx::writeData(
  wb,
  "Pathway_gene_study",
  pathway_gene_study
)


openxlsx::addWorksheet(
  wb,
  "Pathway_check"
)

openxlsx::writeData(
  wb,
  "Pathway_check",
  pathway_check
)


for (sh in names(wb)) {

  openxlsx::freezePane(
    wb,
    sh,
    firstRow = TRUE
  )

}


openxlsx::saveWorkbook(

  wb,

  excel_file,

  overwrite = TRUE

)


############################################################
# ==========================================================
# 26. Publication-quality figures
# ==========================================================
############################################################


############################################################
# Plotting dataframe
############################################################

plot_df <- loso_summary %>%

  dplyr::mutate(

    Description = factor(
      Description,
      levels = rev(
        target_pathways
      )
    ),

    Dataset_removed = factor(
      Dataset_removed,
      levels = studies
    )

  )


############################################################
# 27. Figure 1
# Heatmap: percentage of pathway genes LOST
#
# This is the clearest figure for contribution:
# higher value = stronger dependence on that dataset.
############################################################

p_loss <- ggplot(

  plot_df,

  aes(
    x = Dataset_removed,
    y = Description,
    fill = Loss_percent
  )

) +

  geom_tile(
    colour = "white",
    linewidth = 0.6
  ) +

  geom_text(

    aes(
      label = sprintf(
        "%.1f%%",
        Loss_percent
      )
    ),

    size = 3.3

  ) +

  scale_fill_viridis_c(

    option = "C",

    name = "Genes lost (%)"

  ) +

  scale_y_discrete(

    labels = function(x) {

      stringr::str_wrap(
        x,
        width = 38
      )

    }

  ) +

  labs(
    x = "Dataset excluded",
    y = NULL
  ) +

  theme_classic(
    base_size = 11
  ) +

  theme(

    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      size = 10,
      colour = "black"
    ),

    axis.text.y = element_text(
      size = 9.5,
      colour = "black"
    ),

    axis.title.x = element_text(
      size = 11
    ),

    axis.ticks = element_blank(),

    axis.line = element_blank(),

    legend.position = "right",

    legend.title = element_text(
      size = 10
    ),

    legend.text = element_text(
      size = 9
    ),

    plot.margin = margin(
      10,
      15,
      10,
      10
    )

  )


############################################################
# Save heatmap
############################################################

ggsave(

  filename = file.path(
    output_dir,
    "Figure_LOSO_pathway_gene_loss_heatmap.pdf"
  ),

  plot = p_loss,

  width = 9.5,

  height = 7,

  units = "in"

)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_LOSO_pathway_gene_loss_heatmap.tiff"
  ),

  plot = p_loss,

  width = 9.5,

  height = 7,

  units = "in",

  dpi = 600,

  compression = "lzw"

)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_LOSO_pathway_gene_loss_heatmap.png"
  ),

  plot = p_loss,

  width = 9.5,

  height = 7,

  units = "in",

  dpi = 300

)


############################################################
# 28. Figure 2
# Heatmap: unique gene contribution
############################################################

p_unique <- ggplot(

  plot_df,

  aes(
    x = Dataset_removed,
    y = Description,
    fill = Unique_contribution
  )

) +

  geom_tile(
    colour = "white",
    linewidth = 0.6
  ) +

  geom_text(

    aes(
      label = Unique_contribution
    ),

    size = 3.3

  ) +

  scale_fill_viridis_c(

    option = "D",

    name = "Unique genes"

  ) +

  scale_y_discrete(

    labels = function(x) {

      stringr::str_wrap(
        x,
        width = 38
      )

    }

  ) +

  labs(
    x = "Dataset excluded",
    y = NULL
  ) +

  theme_classic(
    base_size = 11
  ) +

  theme(

    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      size = 10,
      colour = "black"
    ),

    axis.text.y = element_text(
      size = 9.5,
      colour = "black"
    ),

    axis.ticks = element_blank(),

    axis.line = element_blank(),

    legend.position = "right"

  )


############################################################
# Save unique contribution heatmap
############################################################

ggsave(

  filename = file.path(
    output_dir,
    "Figure_LOSO_unique_gene_contribution_heatmap.pdf"
  ),

  plot = p_unique,

  width = 9.5,

  height = 7,

  units = "in"

)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_LOSO_unique_gene_contribution_heatmap.tiff"
  ),

  plot = p_unique,

  width = 9.5,

  height = 7,

  units = "in",

  dpi = 600,

  compression = "lzw"

)


############################################################
# 29. Figure 3
# Overall dataset contribution
############################################################

dataset_plot_df <- dataset_summary %>%

  dplyr::mutate(

    Dataset_removed = factor(

      Dataset_removed,

      levels = Dataset_removed[
        order(
          Mean_loss_percent
        )
      ]

    )

  )


p_dataset <- ggplot(

  dataset_plot_df,

  aes(
    x = Dataset_removed,
    y = Mean_loss_percent
  )

) +

  geom_col(
    width = 0.65
  ) +

  geom_text(

    aes(
      label = sprintf(
        "%.1f%%",
        Mean_loss_percent
      )
    ),

    vjust = -0.4,

    size = 3.7

  ) +

  labs(

    x = NULL,

    y = "Mean pathway gene loss (%)"

  ) +

  theme_classic(
    base_size = 11
  ) +

  theme(

    axis.text.x = element_text(
      angle = 25,
      hjust = 1,
      colour = "black"
    ),

    axis.text.y = element_text(
      colour = "black"
    ),

    axis.title.y = element_text(
      size = 11
    ),

    plot.margin = margin(
      10,
      10,
      10,
      10
    )

  ) +

  expand_limits(
    y = max(
      dataset_plot_df$Mean_loss_percent
    ) * 1.12
  )


############################################################
# Save overall barplot
############################################################

ggsave(

  filename = file.path(
    output_dir,
    "Figure_LOSO_overall_dataset_contribution.pdf"
  ),

  plot = p_dataset,

  width = 7,

  height = 5.5,

  units = "in"

)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_LOSO_overall_dataset_contribution.tiff"
  ),

  plot = p_dataset,

  width = 7,

  height = 5.5,

  units = "in",

  dpi = 600,

  compression = "lzw"

)


############################################################
# 30. Final output
############################################################

cat(
  "\n\n========================================\n",
  "LEAVE-ONE-STUDY-OUT ANALYSIS COMPLETED\n",
  "========================================\n"
)


cat(
  "\nOutput directory:\n",
  output_dir,
  "\n"
)


cat(
  "\nMain Excel:\n",
  excel_file,
  "\n"
)


cat(
  "\nDataset contribution summary:\n"
)

print(
  dataset_summary
)


cat(
  "\nMost influential dataset per pathway:\n"
)

print(
  most_influential
)


cat(
  "\nPublication-quality figure:\n",
  "Figure_LOSO_pathway_gene_loss_heatmap.tiff\n"
)
