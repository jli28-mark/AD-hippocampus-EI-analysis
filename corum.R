############################################################
# ==========================================================
# CORUM protein-complex enrichment analysis
# Whole PPI network
# ==========================================================
#
# Input:
#   CORUM.xlsx
#   PPI.txt
#
# Query:
#   Unique genes in the PPI network
#
# Reference universe:
#   All human genes represented in CORUM
#
# Statistics:
#   One-sided Fisher's exact test
#   Benjamini-Hochberg FDR correction
#
# Significance:
#   FDR < 0.05
#
# Main reporting:
#   FDR < 0.05 and overlap >= 2 genes
#
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
  "openxlsx",
  "scales"
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
  library(scales)

})


############################################################
# 2. File paths
############################################################

base_dir <- "D:/My Desktop/09.09返修/CORUM"


corum_file <- file.path(
  base_dir,
  "CORUM.xlsx"
)


ppi_file <- file.path(
  base_dir,
  "PPI.txt"
)


output_dir <- file.path(
  base_dir,
  "CORUM_enrichment_PPI"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 3. Analysis settings
############################################################

fdr_cutoff <- 0.05

minimum_complex_size <- 2

minimum_overlap_for_reporting <- 2

top_n_plot <- 20


############################################################
# 4. Check input files
############################################################

if (!file.exists(corum_file)) {

  stop(
    paste0(
      "CORUM file not found:\n",
      corum_file
    )
  )

}


if (!file.exists(ppi_file)) {

  stop(
    paste0(
      "PPI gene file not found:\n",
      ppi_file
    )
  )

}


############################################################
# 5. Read gene list
#
# PPI.txt:
# one gene symbol per line
# no header required
############################################################

read_gene_list <- function(file) {

  x <- readLines(
    file,
    warn = FALSE
  )

  x <- trimws(
    as.character(x)
  )

  x <- x[
    !is.na(x) &
      x != ""
  ]

  x <- toupper(x)

  unique(x)

}


############################################################
# 6. Read PPI genes
############################################################

ppi_raw_lines <- readLines(
  ppi_file,
  warn = FALSE
)


ppi_raw_trimmed <- trimws(
  ppi_raw_lines
)


ppi_raw_nonempty <- ppi_raw_trimmed[
  ppi_raw_trimmed != ""
]


ppi_genes <- read_gene_list(
  ppi_file
)


cat(
  "\n========================================\n"
)

cat(
  "PPI INPUT SUMMARY\n"
)

cat(
  "========================================\n"
)

cat(
  "Raw lines in PPI.txt:",
  length(ppi_raw_lines),
  "\n"
)

cat(
  "Non-empty lines:",
  length(ppi_raw_nonempty),
  "\n"
)

cat(
  "Unique PPI genes:",
  length(ppi_genes),
  "\n"
)


############################################################
# 7. Detect duplicated genes in PPI.txt
############################################################

duplicated_genes <- unique(
  toupper(
    ppi_raw_nonempty[
      duplicated(
        toupper(ppi_raw_nonempty)
      )
    ]
  )
)


if (length(duplicated_genes) > 0) {

  cat(
    "\nDuplicated genes detected:\n"
  )

  print(
    duplicated_genes
  )

} else {

  cat(
    "\nNo duplicated genes detected.\n"
  )

}


############################################################
# 8. Read CORUM database
############################################################

corum_raw <- readxl::read_excel(
  corum_file,
  sheet = 1
)


corum_raw <- as.data.frame(
  corum_raw
)


cat(
  "\n========================================\n"
)

cat(
  "CORUM DATABASE\n"
)

cat(
  "========================================\n"
)

cat(
  "Total rows:",
  nrow(corum_raw),
  "\n"
)

cat(
  "Total columns:",
  ncol(corum_raw),
  "\n"
)


############################################################
# 9. Required columns
############################################################

required_cols <- c(
  "complex_id",
  "complex_name",
  "organism",
  "subunits_gene_name"
)


missing_cols <- setdiff(
  required_cols,
  colnames(corum_raw)
)


if (length(missing_cols) > 0) {

  stop(
    paste0(
      "Missing required columns:\n",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  )

}


############################################################
# 10. Keep human CORUM complexes
############################################################

corum_human <- corum_raw %>%

  dplyr::filter(

    grepl(
      "human|homo sapiens",
      organism,
      ignore.case = TRUE
    )

  )


cat(
  "Human CORUM records:",
  nrow(corum_human),
  "\n"
)


############################################################
# 11. Convert CORUM to complex-gene long format
############################################################

corum_long <- corum_human %>%

  dplyr::transmute(

    complex_id =
      as.character(complex_id),

    complex_name =
      as.character(complex_name),

    subunits_gene_name =
      as.character(subunits_gene_name)

  ) %>%

  tidyr::separate_rows(
    subunits_gene_name,
    sep = ";"
  ) %>%

  dplyr::mutate(

    Gene = toupper(
      trimws(
        subunits_gene_name
      )
    )

  ) %>%

  dplyr::filter(

    !is.na(Gene),

    Gene != "",

    Gene != "NA"

  ) %>%

  dplyr::select(

    complex_id,
    complex_name,
    Gene

  ) %>%

  dplyr::distinct()


############################################################
# 12. Human CORUM gene universe
############################################################

corum_human_genes <- sort(
  unique(
    corum_long$Gene
  )
)


cat(
  "Unique genes represented in human CORUM:",
  length(corum_human_genes),
  "\n"
)


############################################################
# 13. Define analysis universe
#
# Current analysis:
# all human genes represented in CORUM
############################################################

universe <- corum_human_genes


############################################################
# 14. Map PPI genes to CORUM
############################################################

query_genes <- intersect(
  ppi_genes,
  universe
)


query_unmapped <- setdiff(
  ppi_genes,
  universe
)


cat(
  "\n========================================\n"
)

cat(
  "GENE MAPPING\n"
)

cat(
  "========================================\n"
)

cat(
  "Unique PPI genes:",
  length(ppi_genes),
  "\n"
)

cat(
  "PPI genes represented in CORUM:",
  length(query_genes),
  "\n"
)

cat(
  "PPI genes not represented in CORUM:",
  length(query_unmapped),
  "\n"
)

cat(
  "CORUM background genes:",
  length(universe),
  "\n"
)


############################################################
# 15. Build CORUM complex table
############################################################

complex_table <- corum_long %>%

  dplyr::filter(
    Gene %in% universe
  ) %>%

  dplyr::group_by(
    complex_id,
    complex_name
  ) %>%

  dplyr::summarise(

    Complex_size =
      dplyr::n_distinct(
        Gene
      ),

    Complex_genes =
      paste(
        sort(
          unique(Gene)
        ),
        collapse = ";"
      ),

    .groups = "drop"

  ) %>%

  dplyr::filter(
    Complex_size >= minimum_complex_size
  )


cat(
  "CORUM complexes tested:",
  nrow(complex_table),
  "\n"
)


############################################################
# 16. Fisher enrichment function
############################################################

run_corum_enrichment <- function(
    complex_id_value,
    complex_name_value,
    complex_gene_string
) {

  complex_genes <- unlist(
    strsplit(
      complex_gene_string,
      ";",
      fixed = TRUE
    )
  )


  complex_genes <- unique(
    complex_genes
  )


  overlap_genes <- intersect(
    query_genes,
    complex_genes
  )


  ##########################################################
  # 2 x 2 contingency table
  ##########################################################

  a <- length(
    overlap_genes
  )


  b <- length(
    query_genes
  ) - a


  c <- length(
    complex_genes
  ) - a


  d <- length(
    universe
  ) - a - b - c


  if (d < 0) {

    stop(
      paste0(
        "Invalid contingency table for complex ID: ",
        complex_id_value
      )
    )

  }


  contingency_table <- matrix(

    c(
      a,
      b,
      c,
      d
    ),

    nrow = 2,

    byrow = TRUE

  )


  ##########################################################
  # One-sided Fisher's exact test
  ##########################################################

  fisher_result <- fisher.test(

    contingency_table,

    alternative = "greater"

  )


  ##########################################################
  # Enrichment statistics
  ##########################################################

  gene_ratio <- a /
    length(query_genes)


  background_ratio <- length(complex_genes) /
    length(universe)


  fold_enrichment <- gene_ratio /
    background_ratio


  data.frame(

    complex_id =
      complex_id_value,

    complex_name =
      complex_name_value,

    Query_size =
      length(query_genes),

    Universe_size =
      length(universe),

    Complex_size =
      length(complex_genes),

    Overlap =
      a,

    GeneRatio =
      gene_ratio,

    BackgroundRatio =
      background_ratio,

    Fold_enrichment =
      fold_enrichment,

    P_value =
      fisher_result$p.value,

    Overlap_genes =
      ifelse(
        a == 0,
        "",
        paste(
          sort(overlap_genes),
          collapse = ";"
        )
      ),

    stringsAsFactors = FALSE

  )

}


############################################################
# 17. Run enrichment for all CORUM complexes
############################################################

result_list <- vector(
  mode = "list",
  length = nrow(complex_table)
)


for (i in seq_len(
  nrow(complex_table)
)) {

  result_list[[i]] <- run_corum_enrichment(

    complex_id_value =
      complex_table$complex_id[i],

    complex_name_value =
      complex_table$complex_name[i],

    complex_gene_string =
      complex_table$Complex_genes[i]

  )

}


corum_result <- dplyr::bind_rows(
  result_list
)


############################################################
# 18. Multiple-testing correction
############################################################

corum_result$FDR <- p.adjust(

  corum_result$P_value,

  method = "BH"

)


############################################################
# 19. Add derived variables
############################################################

corum_result <- corum_result %>%

  dplyr::mutate(

    Significant_FDR_0.05 =
      ifelse(
        FDR < fdr_cutoff,
        "Yes",
        "No"
      ),

    Main_reporting =
      ifelse(
        FDR < fdr_cutoff &
          Overlap >= minimum_overlap_for_reporting,
        "Yes",
        "No"
      ),

    minus_log10_FDR =
      -log10(
        pmax(
          FDR,
          1e-300
        )
      )

  ) %>%

  dplyr::arrange(

    FDR,

    dplyr::desc(
      Fold_enrichment
    ),

    dplyr::desc(
      Overlap
    )

  )


############################################################
# 20. Significant results
############################################################

corum_significant <- corum_result %>%

  dplyr::filter(
    FDR < fdr_cutoff
  )


############################################################
# 21. Main reporting results
#
# FDR < 0.05
# overlap >= 2 genes
############################################################

corum_main <- corum_result %>%

  dplyr::filter(

    FDR < fdr_cutoff,

    Overlap >=
      minimum_overlap_for_reporting

  )


############################################################
# 22. Analysis summary
############################################################

analysis_summary <- data.frame(

  Metric = c(

    "Raw_lines_in_PPI_file",

    "Non_empty_lines_in_PPI_file",

    "Input_unique_PPI_genes",

    "Duplicated_gene_symbols",

    "CORUM_mappable_PPI_genes",

    "PPI_genes_not_in_CORUM",

    "Human_CORUM_background_genes",

    "Human_CORUM_complexes_tested",

    "Significant_complexes_FDR_lt_0.05",

    "Main_reported_complexes_FDR_lt_0.05_overlap_ge_2"

  ),

  Value = c(

    length(
      ppi_raw_lines
    ),

    length(
      ppi_raw_nonempty
    ),

    length(
      ppi_genes
    ),

    length(
      duplicated_genes
    ),

    length(
      query_genes
    ),

    length(
      query_unmapped
    ),

    length(
      universe
    ),

    nrow(
      corum_result
    ),

    nrow(
      corum_significant
    ),

    nrow(
      corum_main
    )

  ),

  stringsAsFactors = FALSE

)


############################################################
# 23. Save CSV files
############################################################

write.csv(

  corum_result,

  file.path(
    output_dir,
    "CORUM_enrichment_ALL_complexes.csv"
  ),

  row.names = FALSE

)


write.csv(

  corum_significant,

  file.path(
    output_dir,
    "CORUM_enrichment_significant_FDR005.csv"
  ),

  row.names = FALSE

)


write.csv(

  corum_main,

  file.path(
    output_dir,
    "CORUM_enrichment_main_results_FDR005_overlap2.csv"
  ),

  row.names = FALSE

)


write.csv(

  analysis_summary,

  file.path(
    output_dir,
    "CORUM_analysis_summary.csv"
  ),

  row.names = FALSE

)


write.csv(

  data.frame(
    Gene = query_genes
  ),

  file.path(
    output_dir,
    "PPI_genes_mapped_to_CORUM.csv"
  ),

  row.names = FALSE

)


write.csv(

  data.frame(
    Gene = query_unmapped
  ),

  file.path(
    output_dir,
    "PPI_genes_not_in_CORUM.csv"
  ),

  row.names = FALSE

)


write.csv(

  data.frame(
    Gene = duplicated_genes
  ),

  file.path(
    output_dir,
    "Duplicated_genes_in_PPI_input.csv"
  ),

  row.names = FALSE

)


############################################################
# 24. Create publication-ready Excel workbook
############################################################

excel_file <- file.path(

  output_dir,

  "CORUM_PPI_enrichment_complete_results.xlsx"

)


wb <- openxlsx::createWorkbook()


############################################################
# Excel styles
############################################################

header_style <- openxlsx::createStyle(

  textDecoration = "bold",

  halign = "center",

  valign = "center",

  border = "Bottom"

)


############################################################
# Summary sheet
############################################################

openxlsx::addWorksheet(
  wb,
  "Summary"
)


openxlsx::writeData(
  wb,
  "Summary",
  analysis_summary
)


############################################################
# Main results
############################################################

openxlsx::addWorksheet(
  wb,
  "Main_results"
)


openxlsx::writeData(
  wb,
  "Main_results",
  corum_main
)


############################################################
# Significant results
############################################################

openxlsx::addWorksheet(
  wb,
  "Significant_FDR005"
)


openxlsx::writeData(
  wb,
  "Significant_FDR005",
  corum_significant
)


############################################################
# Complete enrichment
############################################################

openxlsx::addWorksheet(
  wb,
  "ALL_complexes"
)


openxlsx::writeData(
  wb,
  "ALL_complexes",
  corum_result
)


############################################################
# Mapped PPI genes
############################################################

openxlsx::addWorksheet(
  wb,
  "Mapped_PPI_genes"
)


openxlsx::writeData(

  wb,

  "Mapped_PPI_genes",

  data.frame(
    Gene = query_genes
  )

)


############################################################
# Unmapped PPI genes
############################################################

openxlsx::addWorksheet(
  wb,
  "Unmapped_PPI_genes"
)


openxlsx::writeData(

  wb,

  "Unmapped_PPI_genes",

  data.frame(
    Gene = query_unmapped
  )

)


############################################################
# Duplicate input genes
############################################################

openxlsx::addWorksheet(
  wb,
  "Duplicated_input_genes"
)


openxlsx::writeData(

  wb,

  "Duplicated_input_genes",

  data.frame(
    Gene = duplicated_genes
  )

)


############################################################
# CORUM gene membership
############################################################

openxlsx::addWorksheet(
  wb,
  "CORUM_complex_genes"
)


openxlsx::writeData(
  wb,
  "CORUM_complex_genes",
  corum_long
)


############################################################
# Format worksheets
############################################################

for (sheet_name in names(wb)) {

  openxlsx::freezePane(

    wb,

    sheet = sheet_name,

    firstRow = TRUE

  )


  sheet_data <- openxlsx::readWorkbook(

    wb,

    sheet = sheet_name

  )


  if (ncol(sheet_data) > 0) {

    openxlsx::addStyle(

      wb,

      sheet = sheet_name,

      style = header_style,

      rows = 1,

      cols = seq_len(
        ncol(sheet_data)
      ),

      gridExpand = TRUE

    )


    openxlsx::setColWidths(

      wb,

      sheet = sheet_name,

      cols = seq_len(
        ncol(sheet_data)
      ),

      widths = "auto"

    )

  }

}


openxlsx::saveWorkbook(

  wb,

  excel_file,

  overwrite = TRUE

)


############################################################
# ==========================================================
# 25. Publication-quality CORUM dot plot
# ==========================================================
############################################################

plot_data <- corum_main %>%

  dplyr::slice_head(
    n = top_n_plot
  )


if (nrow(plot_data) > 0) {


  ##########################################################
  # IMPORTANT:
  # Keep complex names on ONE LINE
  ##########################################################

  plot_data <- plot_data %>%

    dplyr::mutate(

      complex_name_plot =
        as.character(
          complex_name
        )

    )


  ##########################################################
  # Preserve ranking
  ##########################################################

  plot_data$complex_name_plot <- factor(

    plot_data$complex_name_plot,

    levels = rev(
      plot_data$complex_name_plot
    )

  )


  ##########################################################
  # Dot plot
  ##########################################################

  p_dot <- ggplot(

    plot_data,

    aes(
      x = Fold_enrichment,
      y = complex_name_plot
    )

  ) +

    geom_point(

      aes(

        size = Overlap,

        colour = minus_log10_FDR

      ),

      alpha = 0.90

    ) +

    scale_size_continuous(

      name = "Overlap genes",

      range = c(
        3,
        10
      )

    ) +

    scale_colour_viridis_c(

      option = "D",

      name = expression(
        -log[10](FDR)
      )

    ) +

    labs(

      x = "Fold enrichment",

      y = NULL

    ) +

    theme_classic(
      base_size = 11
    ) +

    theme(

      axis.text.x = element_text(

        colour = "black",

        size = 9.5

      ),

      axis.text.y = element_text(

        colour = "black",

        size = 8.5,

        lineheight = 1

      ),

      axis.title.x = element_text(

        size = 10.5,

        margin = margin(
          t = 8
        )

      ),

      axis.ticks.y = element_line(),

      legend.title = element_text(
        size = 9.5
      ),

      legend.text = element_text(
        size = 9
      ),

      legend.key.height = unit(
        0.45,
        "cm"
      ),

      legend.spacing.y = unit(
        0.15,
        "cm"
      ),

      plot.margin = margin(

        t = 10,

        r = 20,

        b = 10,

        l = 10

      )

    )


  ##########################################################
  # Save PDF
  ##########################################################

  ggsave(

    filename = file.path(
      output_dir,
      "Figure_CORUM_PPI_top_complexes_dotplot.pdf"
    ),

    plot = p_dot,

    width = 13,

    height = 8,

    units = "in",

    limitsize = FALSE

  )


  ##########################################################
  # Save high-resolution TIFF
  ##########################################################

  ggsave(

    filename = file.path(
      output_dir,
      "Figure_CORUM_PPI_top_complexes_dotplot.tiff"
    ),

    plot = p_dot,

    width = 13,

    height = 8,

    units = "in",

    dpi = 600,

    compression = "lzw",

    limitsize = FALSE

  )


  ##########################################################
  # Save PNG
  ##########################################################

  ggsave(

    filename = file.path(
      output_dir,
      "Figure_CORUM_PPI_top_complexes_dotplot.png"
    ),

    plot = p_dot,

    width = 13,

    height = 8,

    units = "in",

    dpi = 300,

    limitsize = FALSE

  )

}


############################################################
# ==========================================================
# 26. Supplementary overlap bar plot
# ==========================================================
############################################################

if (nrow(plot_data) > 0) {


  p_bar <- ggplot(

    plot_data,

    aes(
      x = Overlap,
      y = complex_name_plot
    )

  ) +

    geom_col(
      width = 0.70
    ) +

    geom_text(

      aes(
        label = Overlap
      ),

      hjust = -0.25,

      size = 3.1

    ) +

    scale_x_continuous(

      expand = expansion(

        mult = c(
          0,
          0.12
        )

      )

    ) +

    labs(

      x = "Number of overlapping PPI genes",

      y = NULL

    ) +

    theme_classic(
      base_size = 11
    ) +

    theme(

      axis.text.x = element_text(

        colour = "black",

        size = 9.5

      ),

      axis.text.y = element_text(

        colour = "black",

        size = 8.5

      ),

      axis.title.x = element_text(

        size = 10.5,

        margin = margin(
          t = 8
        )

      ),

      plot.margin = margin(

        t = 10,

        r = 25,

        b = 10,

        l = 10

      )

    )


  ##########################################################
  # Save bar plot PDF
  ##########################################################

  ggsave(

    filename = file.path(
      output_dir,
      "Figure_CORUM_PPI_overlap_barplot.pdf"
    ),

    plot = p_bar,

    width = 13,

    height = 8,

    units = "in",

    limitsize = FALSE

  )


  ##########################################################
  # Save bar plot TIFF
  ##########################################################

  ggsave(

    filename = file.path(
      output_dir,
      "Figure_CORUM_PPI_overlap_barplot.tiff"
    ),

    plot = p_bar,

    width = 13,

    height = 8,

    units = "in",

    dpi = 600,

    compression = "lzw",

    limitsize = FALSE

  )

}


############################################################
# 27. Save plotted data
############################################################

if (nrow(plot_data) > 0) {

  write.csv(

    plot_data,

    file.path(
      output_dir,
      "CORUM_top20_complexes_used_for_plot.csv"
    ),

    row.names = FALSE

  )

}


############################################################
# 28. Save session information
############################################################

writeLines(

  capture.output(
    sessionInfo()
  ),

  file.path(
    output_dir,
    "sessionInfo.txt"
  )

)


############################################################
# 29. Final report
############################################################

cat(
  "\n\n========================================\n"
)

cat(
  "CORUM ENRICHMENT COMPLETED\n"
)

cat(
  "========================================\n\n"
)


print(
  analysis_summary
)


cat(
  "\nOutput directory:\n",
  output_dir,
  "\n"
)


cat(
  "\nSupplementary Excel file:\n",
  excel_file,
  "\n"
)


cat(
  "\nMain enrichment table:\n",
  file.path(
    output_dir,
    "CORUM_enrichment_main_results_FDR005_overlap2.csv"
  ),
  "\n"
)


cat(
  "\nPublication-quality dot plot:\n",
  file.path(
    output_dir,
    "Figure_CORUM_PPI_top_complexes_dotplot.tiff"
  ),
  "\n"
)


cat(
  "\nVector PDF:\n",
  file.path(
    output_dir,
    "Figure_CORUM_PPI_top_complexes_dotplot.pdf"
  ),
  "\n"
)
