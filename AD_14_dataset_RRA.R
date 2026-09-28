# AD hippocampus: audit 16 contrasts from 14 GEO series and perform RRA.
# Save this file on the Windows computer and run it in RStudio with Source.
# Original DEG, KEGG, STRING and manuscript files are read-only in this script.
#
# One-time installation, if needed:
# install.packages(c("readxl", "RobustRankAggreg"))
# If 00_file_audit.csv says BLOCK, correct that input before interpreting RRA.
# REVIEW rows remain in RRA; verify their full-table provenance before publication.
# A file with only significant DEGs cannot be expanded to all genes here;
# export the COMPLETE results from the original DESeq2/limma analysis.

options(stringsAsFactors = FALSE, scipen = 999)

# --------------------------- USER SETTINGS -------------------------------
input_dir <- "D:/My Desktop/09.09返修/14_data"
fdr_cutoff <- 0.05
log2fc_cutoff <- 0.58
low_coverage_flag_below <- 8000L # REVIEW flag, never a hard platform-size cutoff
min_rankable_fraction <- 0.80 # Block if most genes have no raw P/log2FC
min_overlap_genes <- 100L # Technical guard for study and cross-study ranking
min_series_coverage <- 7L  # Measured in at least half of the 14 GSE series

# Edit a column override or sheet below only when the audit reports that
# automatic column/sheet detection is incorrect. Column names are case-sensitive.
file_names <- c(
  "GSE173954_03_all_genes_limma.csv",
  "DEG_GSE278723_AD_vs_NC_blockBySubject_regionAdj.csv",
  "GSE280268_all_genes_limma.csv",
  "GSE173955.top.table.xlsx",
  "GSE184942.top.table.xlsx",
  "GSE84422.top.table.xlsx",
  "GSE67333.top.table.xlsx",
  "GSE48350.top.table.xlsx",
  "GSE29378.CA1.top.table.xlsx",
  "GSE29378.CA3.top.table.xlsx",
  "GSE36980.top.table.xlsx",
  "GSE13214.rep1.top.table.xlsx",
  "GSE13214.rep2.top.table.xlsx",
  "GSE28146.top.table.xlsx",
  "13.GSE5281.xlsx",
  "14.GSE1297.xlsx"
)

contrast_names <- c(
  "overall", "region_adjusted", "overall", "overall",
  "overall", "overall", "overall", "overall",
  "CA1", "CA3", "overall", "rep1", "rep2",
  "overall", "overall", "overall"
)

get_gse <- function(s) {
  found <- regmatches(s, regexpr("GSE[0-9]+", s))
  if (length(found) != 1L || !nzchar(found)) stop("No GSE ID: ", s)
  found
}

config <- data.frame(
  file_name = file_names,
  gse = vapply(file_names, get_gse, character(1)),
  contrast = contrast_names,
  sheet = rep("", length(file_names)),
  gene_col = rep("", length(file_names)),
  log2fc_col = rep("", length(file_names)),
  pvalue_col = rep("", length(file_names)),
  padj_col = rep("", length(file_names)),
  stat_col = rep("", length(file_names)),
  mean_expr_col = rep("", length(file_names)),
  probe_id_col = rep("", length(file_names))
)

if (length(unique(config$gse)) != 14L || nrow(config) != 16L) {
  stop("The configuration must contain 16 contrasts from 14 GSE series.")
}
if (!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir)
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Install readxl first: install.packages('readxl')")
}

output_root <- file.path(dirname(input_dir), "14_GSE_rank_integration")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
# Keep every run inside ONE project folder, preserving old runs and audit trails.
run_dir <- file.path(output_root,
                     paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S")))
if (!dir.create(run_dir, showWarnings = FALSE)) {
  stop("Output folder already exists (or could not be created): ", run_dir,
       ". Wait one second before rerunning, or check write permissions.")
}
standard_dir <- file.path(run_dir, "standardized_full_tables")
series_dir <- file.path(run_dir, "series_ranks")
dir.create(standard_dir, showWarnings = FALSE)
dir.create(series_dir, showWarnings = FALSE)

write_table <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

normal_name <- function(x) {
  gsub("[^a-z0-9]", "", tolower(x))
}

aliases <- list(
  gene = c("Gene.symbol", "Gene Symbol", "gene_symbol", "Symbol",
           "hgnc_symbol", "external_gene_name", "gene_name", "gene"),
  log2fc = c("log2FoldChange", "log2FC", "logFC", "log_fold_change", "lfc"),
  pvalue = c("P.Value", "pvalue", "p_value", "pval", "p"),
  padj = c("adj.P.Val", "padj", "p.adjust", "adj_p_val",
           "adjusted_p_value", "FDR", "FDR_BH", "BH_pvalue"),
  stat = c("stat", "t", "t_stat", "t_value", "wald_stat",
           "wald_statistic", "z_score", "F"),
  mean_expr = c("AveExpr", "baseMean", "mean_expression", "average_expression"),
  probe_id = c("probe_id", "ProbeID", "probeset_id", "probe", "ID")
)

find_col <- function(headers, candidates, override = "") {
  if (!is.na(override) && nzchar(override)) {
    if (!(override %in% headers)) stop("Column override not found: ", override)
    return(override)
  }
  h <- normal_name(headers)
  for (candidate in candidates) {
    hit <- which(h == normal_name(candidate))
    if (length(hit) == 1L) return(headers[hit])
    if (length(hit) > 1L) {
      stop("Ambiguous columns named like ", candidate, "; use a column override.")
    }
  }
  NA_character_
}

sheet_score <- function(headers, row_cfg) {
  fields <- c("gene", "log2fc", "pvalue", "padj")
  overrides <- c(row_cfg$gene_col, row_cfg$log2fc_col,
                 row_cfg$pvalue_col, row_cfg$padj_col)
  sum(vapply(seq_along(fields), function(j) {
    !is.na(tryCatch(find_col(headers, aliases[[fields[j]]], overrides[j]),
                    error = function(e) NA_character_))
  }, logical(1)))
}

read_input <- function(path, row_cfg) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    tab <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                           fileEncoding = "UTF-8-BOM")
    names(tab) <- sub("^\ufeff", "", names(tab))
    return(list(data = tab, sheet = "CSV"))
  }
  if (ext != "xlsx") stop("Unsupported extension: ", ext)
  sheets <- readxl::excel_sheets(path)
  if (nzchar(row_cfg$sheet)) {
    if (!(row_cfg$sheet %in% sheets)) stop("Sheet not found: ", row_cfg$sheet)
    chosen <- row_cfg$sheet
  } else {
    scores <- vapply(sheets, function(sh) {
      hdr <- names(readxl::read_excel(path, sheet = sh, n_max = 0,
                                      .name_repair = "unique"))
      sheet_score(hdr, row_cfg)
    }, integer(1))
    chosen <- sheets[which.max(scores)]
  }
  tab <- as.data.frame(readxl::read_excel(path, sheet = chosen,
                                           .name_repair = "unique"),
                       check.names = FALSE)
  list(data = tab, sheet = chosen)
}

as_number <- function(x) suppressWarnings(as.numeric(trimws(as.character(x))))

audit_template <- function(i) {
  data.frame(
    file_name = config$file_name[i], gse = config$gse[i],
    contrast = config$contrast[i], sheet = "",
    gene_col = "", log2fc_col = "", pvalue_col = "", padj_col = "",
    stat_col = "", mean_expr_col = "", probe_id_col = "",
    n_input_rows = NA_integer_,
    n_missing_gene = NA_integer_, n_multigene_rows = NA_integer_,
    n_duplicate_gene_rows = NA_integer_, n_unique_genes = NA_integer_,
    n_with_raw_p = NA_integer_, n_with_adjusted_p = NA_integer_,
    n_with_log2fc = NA_integer_, n_rankable = NA_integer_,
    max_raw_pvalue = NA_real_,
    n_raw_p_ge_0_5 = NA_integer_, n_adjusted_p_ge_0_05 = NA_integer_,
    n_raw_feature_rows_both_DEG_conditions = NA_integer_,
    n_adjusted_p_lt_0_05 = NA_integer_,
    n_abs_log2fc_ge_0_58 = NA_integer_,
    n_both_DEG_conditions = NA_integer_,
    status = "BLOCK", note = "", check.names = FALSE
  )
}

standardize_one <- function(i) {
  a <- audit_template(i)
  path <- file.path(input_dir, config$file_name[i])
  if (!file.exists(path)) {
    a$note <- "File not found. Check file name and input_dir."
    return(list(audit = a, table = NULL))
  }
  attempt <- tryCatch({
    row_cfg <- config[i, , drop = FALSE]
    read_result <- read_input(path, row_cfg)
    raw <- read_result$data
    a$sheet <- read_result$sheet
    a$n_input_rows <- nrow(raw)
    if (nrow(raw) == 0L) stop("The selected sheet contains no rows.")

    selected <- c(
      gene = find_col(names(raw), aliases$gene, row_cfg$gene_col),
      log2fc = find_col(names(raw), aliases$log2fc, row_cfg$log2fc_col),
      pvalue = find_col(names(raw), aliases$pvalue, row_cfg$pvalue_col),
      padj = find_col(names(raw), aliases$padj, row_cfg$padj_col),
      stat = find_col(names(raw), aliases$stat, row_cfg$stat_col),
      mean_expr = find_col(names(raw), aliases$mean_expr, row_cfg$mean_expr_col),
      probe_id = find_col(names(raw), aliases$probe_id, row_cfg$probe_id_col)
    )
    a$gene_col <- ifelse(is.na(selected["gene"]), "", selected["gene"])
    a$log2fc_col <- ifelse(is.na(selected["log2fc"]), "", selected["log2fc"])
    a$pvalue_col <- ifelse(is.na(selected["pvalue"]), "", selected["pvalue"])
    a$padj_col <- ifelse(is.na(selected["padj"]), "", selected["padj"])
    a$stat_col <- ifelse(is.na(selected["stat"]), "", selected["stat"])
    a$mean_expr_col <- ifelse(is.na(selected["mean_expr"]), "", selected["mean_expr"])
    a$probe_id_col <- ifelse(is.na(selected["probe_id"]), "", selected["probe_id"])

    needed <- c("gene", "log2fc", "padj")
    if (anyNA(selected[needed])) {
      missing <- paste(needed[is.na(selected[needed])], collapse = ", ")
      stop("Missing required column(s): ", missing, ". Available: ",
           paste(names(raw), collapse = " | "))
    }
    optional_num <- function(key) {
      if (is.na(selected[key])) {
        rep(NA_real_, nrow(raw))
      } else {
        as_number(raw[[selected[key]]])
      }
    }
    optional_text <- function(key) {
      if (is.na(selected[key])) {
        rep(NA_character_, nrow(raw))
      } else {
        trimws(as.character(raw[[selected[key]]]))
      }
    }
    gene <- toupper(trimws(as.character(raw[[selected["gene"]]])))
    missing_gene <- is.na(gene) | gene %in% c("", "NA", "N/A", "---", "NONE")
    multigene <- !missing_gene & grepl("///|;|,|\\|", gene)
    a$n_missing_gene <- sum(missing_gene)
    a$n_multigene_rows <- sum(multigene)

    tab <- data.frame(
      gse = config$gse[i], contrast = config$contrast[i],
      gene_symbol = gene, log2FC = optional_num("log2fc"),
      pvalue = optional_num("pvalue"), padj = optional_num("padj"),
      stat = optional_num("stat"), mean_expression = optional_num("mean_expr"),
      probe_id = optional_text("probe_id"),
      source_file = config$file_name[i], source_row = seq_len(nrow(raw)),
      stringsAsFactors = FALSE
    )
    tab <- tab[!(missing_gene | multigene), , drop = FALSE]
    if (nrow(tab) == 0L) stop("No unambiguous gene symbols remain.")
    tab$pvalue[!is.na(tab$pvalue) &
                 (tab$pvalue < 0 | tab$pvalue > 1)] <- NA_real_
    tab$padj[!is.na(tab$padj) &
               (tab$padj < 0 | tab$padj > 1)] <- NA_real_
    a$n_raw_feature_rows_both_DEG_conditions <- sum(
      is.finite(tab$padj) & tab$padj < fdr_cutoff &
        is.finite(tab$log2FC) & abs(tab$log2FC) >= log2fc_cutoff
    )

    # Resolve multiple probes without choosing the smallest p-value.
    # Prefer highest mean expression; use a unique lexical probe ID for ties.
    duplicate_genes <- unique(tab$gene_symbol[
      duplicated(tab$gene_symbol) |
        duplicated(tab$gene_symbol, fromLast = TRUE)
    ])
    a$n_duplicate_gene_rows <- nrow(tab) - length(unique(tab$gene_symbol))
    keep_rows <- rep(TRUE, nrow(tab))
    for (g in duplicate_genes) {
      rows <- which(tab$gene_symbol == g)
      expr <- tab$mean_expression[rows]
      candidates <- if (any(is.finite(expr))) {
        rows[is.finite(expr) & expr == max(expr, na.rm = TRUE)]
      } else rows
      if (length(candidates) > 1L) {
        ids <- tab$probe_id[candidates]
        if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) > 0L) {
          stop("Multiple probes for ", g,
               " lack a unique mean-expression winner or unique probe IDs. ",
               "Provide a gene-level table or a probe_id_col override.")
        }
        candidates <- candidates[order(ids)][1]
      }
      keep_rows[setdiff(rows, candidates)] <- FALSE
    }
    tab <- tab[keep_rows, , drop = FALSE]
    rownames(tab) <- NULL

    a$n_unique_genes <- nrow(tab)
    a$n_with_raw_p <- sum(is.finite(tab$pvalue))
    a$n_raw_p_ge_0_5 <- sum(is.finite(tab$pvalue) & tab$pvalue >= 0.5)
    if (a$n_with_raw_p > 0L) {
      a$max_raw_pvalue <- max(tab$pvalue, na.rm = TRUE)
    }
    a$n_with_adjusted_p <- sum(is.finite(tab$padj))
    a$n_adjusted_p_ge_0_05 <- sum(is.finite(tab$padj) &
                                   tab$padj >= fdr_cutoff)
    a$n_with_log2fc <- sum(is.finite(tab$log2FC))
    a$n_rankable <- sum(is.finite(tab$pvalue) & is.finite(tab$log2FC))
    a$n_adjusted_p_lt_0_05 <- sum(is.finite(tab$padj) &
                                   tab$padj < fdr_cutoff)
    a$n_abs_log2fc_ge_0_58 <- sum(is.finite(tab$log2FC) &
                                   abs(tab$log2FC) >= log2fc_cutoff)
    a$n_both_DEG_conditions <- sum(is.finite(tab$padj) &
                                     tab$padj < fdr_cutoff &
                                     is.finite(tab$log2FC) &
                                     abs(tab$log2FC) >= log2fc_cutoff)

    issues <- character()
    if (is.na(selected["pvalue"])) {
      issues <- c(issues, "Raw P-value column missing; full ranking needs raw P-values.")
    }
    if (is.na(selected["stat"])) {
      issues <- c(issues, "Test statistic column missing; export full test statistics.")
    }
    if (a$n_rankable < min_overlap_genes ||
        a$n_rankable / nrow(tab) < min_rankable_fraction) {
      issues <- c(issues, paste0("Only ", a$n_rankable, " of ", nrow(tab),
                   " genes have both raw P-values and log2FC; review missing values."))
    }
    if (a$n_with_raw_p > 0L && a$max_raw_pvalue < 0.5) {
      issues <- c(issues, "Maximum raw P < 0.5; check if the table was prefiltered.")
    }
    if (a$n_with_adjusted_p == 0L) {
      issues <- c(issues, "No usable BH-adjusted P-values.")
    } else if (a$n_adjusted_p_lt_0_05 > 100L &&
               a$n_adjusted_p_lt_0_05 / a$n_with_adjusted_p > 0.98) {
      issues <- c(issues, "Nearly all rows pass FDR; check for prior DEG filtering.")
    }
    id_like <- grepl("^ENSG[0-9]+(\\.[0-9]+)?$|^[0-9]+$", tab$gene_symbol)
    if (mean(id_like) > 0.8) {
      issues <- c(issues, "Gene column appears to contain Ensembl/Entrez IDs, not symbols.")
    }
    reviews <- character()
    if (nrow(tab) < low_coverage_flag_below) {
      reviews <- c(reviews, paste0(
        "Only ", nrow(tab), " unique genes: verify this is the complete tested ",
        "gene universe after platform/QC filtering (e.g., topTable(number=Inf))."
      ))
    }
    a$status <- if (length(issues)) "BLOCK" else if (length(reviews)) "REVIEW" else "OK"
    a$note <- paste(c(issues, reviews), collapse = " | ")
    out_name <- paste0(config$gse[i], "_", config$contrast[i], "_all_genes.csv")
    write_table(tab, file.path(standard_dir, out_name))
    list(audit = a, table = tab)
  }, error = function(e) {
    a$status <- "BLOCK"
    a$note <- conditionMessage(e)
    list(audit = a, table = NULL)
  })
  attempt
}

# ------------------------ STAGES 1 AND 2 -------------------------------
processed <- lapply(seq_len(nrow(config)), standardize_one)
audit <- do.call(rbind, lapply(processed, `[[`, "audit"))
rownames(audit) <- NULL
write_table(audit, file.path(run_dir, "00_file_audit.csv"))
write_table(
  audit[, c("gse", "contrast", "file_name", "n_unique_genes",
            "n_adjusted_p_lt_0_05", "n_abs_log2fc_ge_0_58",
            "n_both_DEG_conditions", "status", "note")],
  file.path(run_dir, "01_three_counts_by_contrast.csv")
)

write_table(audit[audit$status == "REVIEW", , drop = FALSE],
            file.path(run_dir, "00_review_flags.csv"))
blocked_rows <- which(audit$status == "BLOCK")
if (length(blocked_rows) > 0L) {
  issue_lines <- vapply(blocked_rows, function(i) {
    paste0(
      "- ", audit$file_name[i], " [sheet: ", audit$sheet[i], "] ",
      audit$note[i], " | rows=", audit$n_input_rows[i],
      ", unique_genes=", audit$n_unique_genes[i],
      ", raw_P_col=", audit$pvalue_col[i],
      ", adjusted_P_col=", audit$padj_col[i],
      ", statistic_col=", audit$stat_col[i]
    )
  }, character(1))
  message_lines <- c(
    paste0("RRA blocked: ", length(blocked_rows), " of ", nrow(config),
           " input tables require review."),
    issue_lines,
    paste0("Complete audit: ", file.path(run_dir, "00_file_audit.csv")),
    "Correct the source table, selected sheet or column override; then rerun."
  )
  writeLines(message_lines, file.path(run_dir, "00_blocking_issues.txt"))
  writeLines(message_lines, file.path(run_dir, "RUN_STATUS.txt"))
  cat(paste(message_lines, collapse = "\n"), "\n")
  stop("Audit completed; RRA blocked. Results folder: ", run_dir)
}
if (any(audit$status == "REVIEW")) {
  message("Continuing RRA with ", sum(audit$status == "REVIEW"),
          " lower-coverage input table(s). See 00_review_flags.csv in: ",
          run_dir,
          ". Verify each file contains all tested genes before publication.")
}

if (!requireNamespace("RobustRankAggreg", quietly = TRUE)) {
  writeLines("Install RobustRankAggreg, then rerun this script: install.packages('RobustRankAggreg')",
             file.path(run_dir, "RUN_STATUS.txt"))
  stop("Install RobustRankAggreg first. Audit files are in: ", run_dir)
}

# ----------------------- STAGE 3: 14 GSE RRA ----------------------------
tabs <- lapply(processed, `[[`, "table")
names(tabs) <- paste(config$gse, config$contrast, sep = "_")
all_gene_statistics <- do.call(rbind, tabs)
rownames(all_gene_statistics) <- NULL
write_table(all_gene_statistics,
            file.path(run_dir, "02_all_gene_statistics_long.csv"))
groups <- split(seq_along(tabs), config$gse)

combine_within_gse <- function(gse, indices) {
  sub_tabs <- lapply(tabs[indices], function(x) {
    x <- x[is.finite(x$pvalue) & is.finite(x$log2FC),
           c("gene_symbol", "log2FC", "pvalue"), drop = FALSE]
    x$signed_score <- sign(x$log2FC) *
      (-log10(pmax(x$pvalue, .Machine$double.xmin)))
    x$up_pct <- rank(-x$signed_score, ties.method = "average") / nrow(x)
    x$down_pct <- rank(x$signed_score, ties.method = "average") / nrow(x)
    x
  })

  # Multiple regions/replicates in the same GSE become ONE series ranking.
  # Use only genes measured in every contrast of that series.
  common <- sort(Reduce(intersect, lapply(sub_tabs, `[[`, "gene_symbol")))
  if (length(common) < min_overlap_genes) {
    stop(gse, " has only ", length(common),
         " rankable genes shared across its contrasts; review input tables.")
  }
  up <- down <- fc <- matrix(NA_real_, nrow = length(common),
                            ncol = length(sub_tabs))
  for (j in seq_along(sub_tabs)) {
    pos <- match(common, sub_tabs[[j]]$gene_symbol)
    up[, j] <- sub_tabs[[j]]$up_pct[pos]
    down[, j] <- sub_tabs[[j]]$down_pct[pos]
    fc[, j] <- sub_tabs[[j]]$log2FC[pos]
  }
  result <- data.frame(
    gse = gse, gene_symbol = common,
    mean_up_rank_fraction = rowMeans(up),
    mean_down_rank_fraction = rowMeans(down),
    mean_log2FC_descriptive = rowMeans(fc),
    n_contrasts = length(indices), stringsAsFactors = FALSE
  )
  result <- result[order(result$gene_symbol), , drop = FALSE]
  write_table(result, file.path(series_dir, paste0(gse, "_series_rank.csv")))
  result
}

series <- lapply(names(groups), function(gse) combine_within_gse(gse, groups[[gse]]))
names(series) <- names(groups)
if (length(series) != 14L) stop("Expected one ranking per 14 GSE series.")

coverage <- table(unlist(lapply(series, `[[`, "gene_symbol"), use.names = FALSE))
eligible_genes <- names(coverage)[coverage >= min_series_coverage]
if (length(eligible_genes) < min_overlap_genes) {
  stop("Too few genes are measured in at least ", min_series_coverage,
       " GSE series: ", length(eligible_genes))
}

rank_lists <- function(direction) {
  lapply(series, function(x) {
    x <- x[x$gene_symbol %in% eligible_genes, , drop = FALSE]
    score_col <- if (direction == "UP") {
      "mean_up_rank_fraction"
    } else {
      "mean_down_rank_fraction"
    }
    x$gene_symbol[order(x[[score_col]], x$gene_symbol)]
  })
}
up_lists <- rank_lists("UP")
down_lists <- rank_lists("DOWN")

# These are COMPLETE rankings of tested genes, but platforms have different
# gene coverage. full=TRUE treats structurally unmeasured genes as NA.
up_matrix <- RobustRankAggreg::rankMatrix(up_lists, full = TRUE)
down_matrix <- RobustRankAggreg::rankMatrix(down_lists, full = TRUE)
up_rra <- RobustRankAggreg::aggregateRanks(rmat = up_matrix,
                                            method = "RRA", exact = FALSE)
down_rra <- RobustRankAggreg::aggregateRanks(rmat = down_matrix,
                                              method = "RRA", exact = FALSE)
up_result <- data.frame(gene_symbol = as.character(up_rra$Name),
                        rra_p_up = as.numeric(up_rra$Score))
down_result <- data.frame(gene_symbol = as.character(down_rra$Name),
                          rra_p_down = as.numeric(down_rra$Score))

result <- merge(up_result, down_result, by = "gene_symbol", all = TRUE,
                sort = FALSE)
result <- result[order(result$gene_symbol), , drop = FALSE]

# Correct across BOTH directional families (2 tests per gene).
all_p <- c(result$rra_p_up, result$rra_p_down)
all_fdr <- p.adjust(all_p, method = "BH")
nn <- nrow(result)
result$rra_fdr_up <- all_fdr[seq_len(nn)]
result$rra_fdr_down <- all_fdr[nn + seq_len(nn)]
result$n_series_measured <- as.integer(coverage[result$gene_symbol])
result$n_series_positive <- 0L
result$n_series_negative <- 0L
for (x in series) {
  pos <- match(x$gene_symbol, result$gene_symbol)
  keep <- !is.na(pos)
  result$n_series_positive[pos[keep]] <- result$n_series_positive[pos[keep]] +
    as.integer(x$mean_log2FC_descriptive[keep] > 0)
  result$n_series_negative[pos[keep]] <- result$n_series_negative[pos[keep]] +
    as.integer(x$mean_log2FC_descriptive[keep] < 0)
}
significant_up <- is.finite(result$rra_fdr_up) & result$rra_fdr_up < 0.05
significant_down <- is.finite(result$rra_fdr_down) & result$rra_fdr_down < 0.05
result$direction_at_fdr_0_05 <- ifelse(
  significant_up & significant_down, "BOTH_DIRECTIONS",
  ifelse(significant_up, "UP",
         ifelse(significant_down, "DOWN", "NOT_SIGNIFICANT"))
)
result$directionally_supported <- (
  significant_up & result$n_series_positive >= 2L &
    result$n_series_positive > result$n_series_negative & !significant_down
) | (
  significant_down & result$n_series_negative >= 2L &
    result$n_series_negative > result$n_series_positive & !significant_up
)
result <- result[order(pmin(result$rra_fdr_up, result$rra_fdr_down),
                       result$gene_symbol), , drop = FALSE]
rownames(result) <- NULL

write_table(result, file.path(run_dir, "03_RRA_all_genes_up_down.csv"))
write_table(result[result$directionally_supported, , drop = FALSE],
            file.path(run_dir, "04_RRA_FDR_0.05_direction_supported.csv"))
write_table(data.frame(
  gse = names(series),
  n_contrasts = vapply(series, function(x) x$n_contrasts[1], integer(1)),
  n_shared_genes_within_series = vapply(series, nrow, integer(1)),
  n_genes_in_RRA = vapply(up_lists, length, integer(1))
), file.path(run_dir, "05_series_coverage.csv"))

writeLines(c(
  "Completed: input audits, standardized input tables and directional RRA for 14 GEO series.",
  "RRA findings remain provisional until REVIEW-flagged tables are confirmed complete.",
  paste("Output directory:", run_dir),
  paste("Input tables flagged for full-table provenance review:",
        sum(audit$status == "REVIEW")),
  "Low gene count is a review flag, not proof of significant-only filtering.",
  "Check 00_review_flags.csv and confirm that each source exported all tested genes.",
  paste("FDR cutoff for DEG audit:", fdr_cutoff),
  paste("Absolute log2FC cutoff for DEG audit:", log2fc_cutoff),
  paste("Genes measured in at least", min_series_coverage, "series:",
        length(eligible_genes)),
  paste("Directional RRA results at BH FDR < 0.05:",
        sum(result$direction_at_fdr_0_05 != "NOT_SIGNIFICANT")),
  paste("Also supported by at least two series in the same direction:",
        sum(result$directionally_supported)),
  "Rank metric within each contrast: sign(log2FC) * [-log10(raw P)].",
  "Multiple contrasts of one GSE: mean fractional rank over genes shared by them.",
  "RRA input: one full UP list and one full DOWN list per GSE (14 each).",
  "RRA full=TRUE: structural absence from another platform is not a worst rank.",
  "BH adjustment: both directions together across all included genes.",
  "RRA does not weight cohorts by sample size and does not validate pathways or modules.",
  "Confirm that distinct GSE accessions do not reuse the same participants/samples.",
  "The original four-dataset DEG, KEGG and STRING files were not modified."
), file.path(run_dir, "RUN_STATUS.txt"))
capture.output(sessionInfo(), file = file.path(run_dir, "R_sessionInfo.txt"))
message("Completed. Results folder: ", run_dir)
