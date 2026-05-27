#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(readxl)
})

# --------------------------------------------------
# Settings
# --------------------------------------------------
runs_dir <- Sys.getenv(
  "RUNS_DIR",
  "/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260428_GE_GEUVADIS_v2/GeneExpression/runs"
)
run_suffix <- Sys.getenv("RUN_SUFFIX", "_peerauto_pmg0_npc5")

shet_xlsx <- Sys.getenv("SHET_XLSX", "/gpfs/data/mostafavilab/shared_data/gene_information/s_het_info.xlsx")
shet_sheet <- Sys.getenv("SHET_SHEET", "Supplementary Table 1")

method_raw <- Sys.getenv("METHOD_RAW", "all_snps_tmm_raw")
method_irnt <- Sys.getenv("METHOD_IRNT", "all_snps_tmm_irnt")
force_tmm_h2 <- tolower(Sys.getenv("FORCE_TMM_H2", "true")) %in% c("1", "true", "yes", "y")
method_h2_tpm_raw <- Sys.getenv("METHOD_H2_TPM_RAW", "")
method_h2_tpm_irnt <- Sys.getenv("METHOD_H2_TPM_IRNT", "")

method_raw_parts <- str_match(method_raw, "^(all_snps|hm3_no_mhc)_(tmm|tpm)_(raw|irnt|inverse_normal)$")
method_irnt_parts <- str_match(method_irnt, "^(all_snps|hm3_no_mhc)_(tmm|tpm)_(raw|irnt|inverse_normal)$")

method_h2_raw <- method_raw
method_h2_irnt <- method_irnt

if (force_tmm_h2) {
  if (!is.na(method_raw_parts[1, 1])) {
    method_h2_raw <- paste(method_raw_parts[1, 2], "tmm", "raw", sep = "_")
  }
  if (!is.na(method_irnt_parts[1, 1])) {
    method_h2_irnt <- paste(method_irnt_parts[1, 2], "tmm", "irnt", sep = "_")
  }
}

if (!nzchar(method_h2_tpm_raw)) {
  if (!is.na(method_raw_parts[1, 1])) {
    method_h2_tpm_raw <- paste(method_raw_parts[1, 2], "tpm", "raw", sep = "_")
  } else {
    method_h2_tpm_raw <- "all_snps_tpm_raw"
  }
}

if (!nzchar(method_h2_tpm_irnt)) {
  if (!is.na(method_irnt_parts[1, 1])) {
    method_h2_tpm_irnt <- paste(method_irnt_parts[1, 2], "tpm", "irnt", sep = "_")
  } else {
    method_h2_tpm_irnt <- "all_snps_tpm_irnt"
  }
}

tpm_mean_method <- Sys.getenv("TPM_MEAN_METHOD", "")
if (!nzchar(tpm_mean_method)) {
  m <- str_match(method_raw, "^(all_snps|hm3_no_mhc)_(tmm|tpm)_(raw|irnt|inverse_normal)$")
  if (!is.na(m[1, 1])) {
    tpm_mean_method <- paste(m[1, 2], "tpm", "raw", sep = "_")
  } else {
    tpm_mean_method <- "all_snps_tpm_raw"
  }
}

shet_x_mode <- "decile"

n_deciles <- suppressWarnings(as.integer(Sys.getenv("N_DECILES", "10")))
if (!is.finite(n_deciles) || n_deciles != 10L) {
  stop("N_DECILES must be 10 for decile-based plots.")
}

cor_method <- tolower(Sys.getenv("COR_METHOD", "spearman"))
if (!cor_method %in% c("pearson", "spearman", "kendall")) {
  stop("COR_METHOD must be one of: pearson, spearman, kendall")
}

h2_frac_cutoff <- as.numeric(Sys.getenv("H2_FRAC_CUTOFF", "0.0005"))
if (!is.finite(h2_frac_cutoff) || h2_frac_cutoff < 0) {
  stop("H2_FRAC_CUTOFF must be a non-negative numeric value.")
}

analysis_root <- file.path(runs_dir, "_analysis", "shet_tpm_h2_plots")
dir.create(analysis_root, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------
# Helpers
# --------------------------------------------------
normalize_gene_id <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace("\\.[0-9]+$", "")
}

method_candidates <- function(method) {
  out <- method
  if (str_detect(method, "_irnt$")) out <- c(out, str_replace(method, "_irnt$", "_inverse_normal"))
  if (str_detect(method, "_inverse_normal$")) out <- c(out, str_replace(method, "_inverse_normal$", "_irnt"))
  unique(out)
}

find_run_dir <- function(method, require_summary = TRUE) {
  all_dirs <- list.dirs(runs_dir, recursive = FALSE, full.names = TRUE)
  if (length(all_dirs) == 0L) stop("No run directories found under: ", runs_dir)

  hits <- map_dfr(method_candidates(method), function(m) {
    patt <- paste0("^", m, run_suffix, "([_].+)?$")
    tibble(run_dir = all_dirs[str_detect(basename(all_dirs), patt)])
  })

  if (nrow(hits) == 0L) {
    stop("No run directory found for method: ", method, " (suffix=", run_suffix, ")")
  }

  valid <- hits %>%
    distinct(run_dir) %>%
    mutate(
      data_dir = file.path(run_dir, "results", "data"),
      summary_dir = file.path(run_dir, "results", "summary"),
      has_pheno = map_lgl(data_dir, ~ length(list.files(.x, pattern = "^pheno_.*\\.phenotypes\\.tsv$", full.names = TRUE)) > 0),
      has_map = map_lgl(data_dir, ~ length(list.files(.x, pattern = "^pheno_.*\\.gene_index_map\\.txt$", full.names = TRUE)) > 0),
      has_summary = map_lgl(summary_dir, ~ length(list.files(.x, pattern = "^final_heritability_summary.*\\.tsv$", full.names = TRUE)) > 0)
    ) %>%
    filter(has_pheno, has_map, if (require_summary) has_summary else TRUE)

  if (nrow(valid) == 0L) {
    stop("Run directory exists but required files are missing for method: ", method)
  }

  exact <- valid %>% filter(basename(run_dir) == paste0(method, run_suffix))
  if (nrow(exact) > 0L) return(exact$run_dir[[1]])

  valid %>%
    mutate(mtime = file.info(run_dir)$mtime) %>%
    arrange(desc(mtime)) %>%
    pull(run_dir) %>%
    .[[1]]
}

resolve_gene_ids_for_summary <- function(gene_raw, run_dir) {
  gene_chr <- as.character(gene_raw)
  idx_mask <- str_detect(str_trim(gene_chr), "^[0-9]+$")
  if (mean(idx_mask, na.rm = TRUE) < 0.8) return(normalize_gene_id(gene_chr))

  map_file <- list.files(file.path(run_dir, "results", "data"), pattern = "^pheno_.*\\.gene_index_map\\.txt$", full.names = TRUE)
  if (length(map_file) == 0) return(normalize_gene_id(gene_chr))

  map_dt <- fread(map_file[1], sep = "\t", data.table = TRUE)
  if (!all(c("gene_name", "mpheno_index") %in% names(map_dt))) return(normalize_gene_id(gene_chr))

  map_tbl <- map_dt %>%
    transmute(mpheno_index = as.integer(mpheno_index), gene_name = as.character(gene_name))

  idx_tbl <- tibble(
    row_id = seq_along(gene_chr),
    mpheno_index = suppressWarnings(as.integer(gene_chr))
  ) %>%
    left_join(map_tbl, by = "mpheno_index")

  resolved <- gene_chr
  replace_mask <- idx_mask & !is.na(idx_tbl$gene_name)
  resolved[replace_mask] <- idx_tbl$gene_name[replace_mask]
  normalize_gene_id(resolved)
}

read_h2_by_gene <- function(run_dir, h2_col_name) {
  sum_file <- list.files(file.path(run_dir, "results", "summary"), pattern = "^final_heritability_summary.*\\.tsv$", full.names = TRUE)
  if (length(sum_file) == 0L) stop("Missing summary file in: ", run_dir)

  dt <- fread(sum_file[1], sep = "\t", data.table = TRUE)
  req <- c("Gene", "Status", "h2_GREML")
  if (!all(req %in% names(dt))) {
    stop("Summary file missing required columns (Gene, Status, h2_GREML): ", sum_file[1])
  }

  tibble(
    Gene = resolve_gene_ids_for_summary(dt$Gene, run_dir),
    Status = toupper(str_trim(as.character(dt$Status))),
    h2 = suppressWarnings(as.numeric(dt$h2_GREML))
  ) %>%
    filter(Status == "PASS", is.finite(h2)) %>%
    group_by(Gene) %>%
    summarise(!!h2_col_name := mean(h2, na.rm = TRUE), .groups = "drop")
}

read_tpm_metrics_by_gene <- function(run_dir) {
  data_dir <- file.path(run_dir, "results", "data")
  pheno_file <- list.files(data_dir, pattern = "^pheno_.*\\.phenotypes\\.tsv$", full.names = TRUE) %>% .[[1]]
  map_file <- list.files(data_dir, pattern = "^pheno_.*\\.gene_index_map\\.txt$", full.names = TRUE) %>% .[[1]]

  pheno <- fread(pheno_file, header = FALSE, sep = "\t", data.table = TRUE)
  map_dt <- fread(map_file, sep = "\t", data.table = TRUE)

  if (!all(c("gene_name", "mpheno_index") %in% names(map_dt))) {
    stop("Unexpected gene_index_map format in: ", map_file)
  }

  map_dt <- map_dt[order(as.integer(mpheno_index))]
  genes <- normalize_gene_id(map_dt$gene_name)
  expr <- as.matrix(pheno[, -(1:2), with = FALSE])
  storage.mode(expr) <- "numeric"

  if ((ncol(pheno) - 2L) != length(genes)) {
    stop("Phenotype/map dimension mismatch in: ", run_dir)
  }

  if (anyDuplicated(genes) > 0L) {
    idx <- split(seq_along(genes), genes)
    expr <- do.call(
      cbind,
      lapply(idx, function(ii) {
        if (length(ii) == 1L) expr[, ii] else rowMeans(expr[, ii, drop = FALSE], na.rm = TRUE)
      })
    )
    genes <- names(idx)
  } else {
    colnames(expr) <- genes
  }

  tibble(
    Gene = normalize_gene_id(colnames(expr)),
    mean_tpm = colMeans(expr, na.rm = TRUE),
    median_tpm = apply(expr, 2, median, na.rm = TRUE)
  ) %>%
    filter(is.finite(mean_tpm), is.finite(median_tpm)) %>%
    group_by(Gene) %>%
    summarise(
      mean_tpm = mean(mean_tpm, na.rm = TRUE),
      median_tpm = mean(median_tpm, na.rm = TRUE),
      .groups = "drop"
    )
}

plot_h2_shet <- function(df, out_file) {
  h2_long <- df %>%
    select(post_mean, post_mean_bin, h2_tmm_raw, h2_tmm_irnt) %>%
    pivot_longer(
      cols = c(h2_tmm_raw, h2_tmm_irnt),
      names_to = "h2_type",
      values_to = "h2"
    ) %>%
    mutate(h2_type = recode(h2_type, h2_tmm_raw = "RAW", h2_tmm_irnt = "IRNT"))

  decile_tbl <- h2_long %>%
    group_by(post_mean_bin, h2_type) %>%
    summarise(
      mean_h2 = mean(h2, na.rm = TRUE),
      n_genes = n(),
      .groups = "drop"
    )

  p <- ggplot(
    decile_tbl,
    aes(
      x = factor(post_mean_bin),
      y = mean_h2,
      color = h2_type,
      group = h2_type
    )
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.8) +
    labs(
      title = "TMM h2 (RAW vs IRNT): mean by s_het decile",
      subtitle = "Deciles recalculated after merge (1 = lowest, 10 = highest s_het)",
      x = "s_het post_mean decile",
      y = "h2_GREML",
      color = "Normalization"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(out_file, p, width = 8, height = 5, dpi = 300)
}

run_plot_suite_pair <- function(df, suite_name, suite_subtitle) {
  suite_root <- file.path(analysis_root, suite_name)
  plots_dir <- file.path(suite_root, "plots")
  tables_dir <- file.path(suite_root, "tables")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

  suite_tbl <- df %>%
    mutate(
      post_mean_bin = ntile(post_mean, n_deciles),
      mean_tpm_decile = ntile(mean_tpm, n_deciles),
      median_tpm_decile = ntile(median_tpm, n_deciles)
    )

  # Plot 1: Shet decile vs mean/median TPM
  p1_decile_tbl <- suite_tbl %>%
    group_by(post_mean_bin) %>%
    summarise(
      n_genes = n(),
      mean_tpm = mean(mean_tpm, na.rm = TRUE),
      median_tpm = median(median_tpm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(post_mean_bin)

  p1_decile_long <- p1_decile_tbl %>%
    pivot_longer(
      cols = c(mean_tpm, median_tpm),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = recode(metric, mean_tpm = "Mean TPM", median_tpm = "Median TPM")
    )

  p1_decile <- ggplot(
    p1_decile_long,
    aes(x = factor(post_mean_bin), y = value, color = metric, group = metric)
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    labs(
      title = "TPM by s_het decile",
      subtitle = suite_subtitle,
      x = "s_het post_mean decile (1-10)",
      y = "TPM",
      color = "TPM summary"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(plots_dir, "plot1_shet_decile_vs_tpm_mean_median.png"), p1_decile, width = 8, height = 5, dpi = 300)
  fwrite(as.data.table(p1_decile_tbl), file.path(tables_dir, "plot1_shet_decile_vs_tpm_mean_median.tsv"), sep = "\t")

  # Plot 2: TPM decile vs correlation(h2 RAW, h2 IRNT) for mean+median TPM deciles
  p2_mean_tbl <- suite_tbl %>%
    group_by(tpm_decile = mean_tpm_decile) %>%
    summarise(
      metric = "Mean TPM decile",
      n_genes = n(),
      cor_h2_raw_irnt = if_else(
        n() >= 3,
        cor(h2_tmm_raw, h2_tmm_irnt, method = cor_method, use = "pairwise.complete.obs"),
        as.numeric(NA)
      ),
      .groups = "drop"
    )

  p2_median_tbl <- suite_tbl %>%
    group_by(tpm_decile = median_tpm_decile) %>%
    summarise(
      metric = "Median TPM decile",
      n_genes = n(),
      cor_h2_raw_irnt = if_else(
        n() >= 3,
        cor(h2_tmm_raw, h2_tmm_irnt, method = cor_method, use = "pairwise.complete.obs"),
        as.numeric(NA)
      ),
      .groups = "drop"
    )

  p2_tbl <- bind_rows(p2_mean_tbl, p2_median_tbl) %>%
    arrange(metric, tpm_decile)

  p2 <- ggplot(
    p2_tbl,
    aes(x = factor(tpm_decile), y = cor_h2_raw_irnt, color = metric, group = metric)
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    labs(
      title = "Correlation(TMM h2 RAW, TMM h2 IRNT) by TPM decile",
      subtitle = paste0(suite_subtitle, " | correlation method: ", toupper(cor_method), " | h2 source: TMM"),
      x = "TPM decile (1-10)",
      y = "Correlation(TMM h2 RAW, TMM h2 IRNT)",
      color = "Decile definition"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(plots_dir, "plot2_tpm_decile_vs_h2_raw_irnt_correlation_mean_median.png"), p2, width = 8, height = 5, dpi = 300)
  fwrite(as.data.table(p2_tbl), file.path(tables_dir, "plot2_tpm_decile_vs_h2_raw_irnt_correlation_mean_median.tsv"), sep = "\t")

  # Plot 3: Shet decile vs RAW+IRNT h2 (decile only)
  plot_h2_shet(
    df = suite_tbl,
    out_file = file.path(plots_dir, "plot3_shet_decile_vs_h2_raw_irnt_mean_median_lines.png")
  )

  fwrite(as.data.table(suite_tbl), file.path(tables_dir, "shet_tpm_h2_merged_gene_level.tsv"), sep = "\t")

  message("Saved outputs in: ", suite_root)
  message("Plots:")
  message(" - ", file.path(plots_dir, "plot1_shet_decile_vs_tpm_mean_median.png"))
  message(" - ", file.path(plots_dir, "plot2_tpm_decile_vs_h2_raw_irnt_correlation_mean_median.png"))
  message(" - ", file.path(plots_dir, "plot3_shet_decile_vs_h2_raw_irnt_mean_median_lines.png"))
  message("Tables:")
  message(" - ", file.path(tables_dir, "shet_tpm_h2_merged_gene_level.tsv"))
  message(" - ", file.path(tables_dir, "plot1_shet_decile_vs_tpm_mean_median.tsv"))
  message(" - ", file.path(tables_dir, "plot2_tpm_decile_vs_h2_raw_irnt_correlation_mean_median.tsv"))
}

make_decile_map <- function(reference_tbl, source_label, direction = c("ascending", "descending")) {
  direction <- match.arg(direction)
  out <- reference_tbl %>%
    select(Gene, post_mean) %>%
    filter(!is.na(Gene), Gene != "", is.finite(post_mean)) %>%
    group_by(Gene) %>%
    summarise(post_mean = mean(post_mean, na.rm = TRUE), .groups = "drop") %>%
    mutate(decile_raw = ntile(post_mean, n_deciles))

  if (direction == "descending") {
    out <- out %>% mutate(decile = n_deciles + 1L - decile_raw)
  } else {
    out <- out %>% mutate(decile = decile_raw)
  }

  out %>%
    transmute(
      Gene = Gene,
      post_mean_bin = as.integer(decile),
      decile_source = source_label,
      decile_direction = direction
    )
}

run_decile_sensitivity_pair <- function(
  df,
  suite_name,
  suite_subtitle,
  shet_reference_tbl,
  shet_tpm_overlap_tbl
) {
  suite_root <- file.path(analysis_root, suite_name)
  plots_dir <- file.path(suite_root, "plots")
  tables_dir <- file.path(suite_root, "tables")

  decile_maps <- bind_rows(
    make_decile_map(df, "merged_h2_genes", "ascending"),
    make_decile_map(shet_reference_tbl, "global_shet_genes", "ascending"),
    make_decile_map(shet_tpm_overlap_tbl, "shet_tpm_overlap_genes", "ascending")
  )

  sens_tbl <- decile_maps %>%
    inner_join(df, by = "Gene") %>%
    filter(is.finite(post_mean_bin), post_mean_bin >= 1L, post_mean_bin <= n_deciles)

  if (nrow(sens_tbl) == 0L) {
    message("Skipping decile sensitivity: no overlap after joining decile maps.")
    return(invisible(NULL))
  }

  sens_decile_summary <- sens_tbl %>%
    group_by(decile_source, decile_direction, post_mean_bin) %>%
    summarise(
      n_genes = n(),
      mean_post_mean = mean(post_mean, na.rm = TRUE),
      mean_tpm = mean(mean_tpm, na.rm = TRUE),
      median_tpm = median(median_tpm, na.rm = TRUE),
      mean_h2_raw = mean(h2_tmm_raw, na.rm = TRUE),
      median_h2_raw = median(h2_tmm_raw, na.rm = TRUE),
      mean_h2_irnt = mean(h2_tmm_irnt, na.rm = TRUE),
      median_h2_irnt = median(h2_tmm_irnt, na.rm = TRUE),
      prop_h2_raw_gt_frac_cutoff = mean(h2_tmm_raw > h2_frac_cutoff, na.rm = TRUE),
      prop_h2_irnt_gt_frac_cutoff = mean(h2_tmm_irnt > h2_frac_cutoff, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(decile_source, decile_direction, post_mean_bin)

  sens_stats <- sens_decile_summary %>%
    group_by(decile_source, decile_direction) %>%
    summarise(
      rho_mean_h2_raw = cor(post_mean_bin, mean_h2_raw, method = "spearman", use = "pairwise.complete.obs"),
      rho_median_h2_raw = cor(post_mean_bin, median_h2_raw, method = "spearman", use = "pairwise.complete.obs"),
      rho_mean_h2_irnt = cor(post_mean_bin, mean_h2_irnt, method = "spearman", use = "pairwise.complete.obs"),
      rho_median_h2_irnt = cor(post_mean_bin, median_h2_irnt, method = "spearman", use = "pairwise.complete.obs"),
      rho_mean_tpm = cor(post_mean_bin, mean_tpm, method = "spearman", use = "pairwise.complete.obs"),
      rho_median_tpm = cor(post_mean_bin, median_tpm, method = "spearman", use = "pairwise.complete.obs"),
      decile10_minus_decile1_mean_h2_raw = mean_h2_raw[post_mean_bin == 10] - mean_h2_raw[post_mean_bin == 1],
      decile10_minus_decile1_mean_h2_irnt = mean_h2_irnt[post_mean_bin == 10] - mean_h2_irnt[post_mean_bin == 1],
      .groups = "drop"
    ) %>%
    arrange(decile_source, decile_direction)

  fwrite(as.data.table(sens_decile_summary), file.path(tables_dir, "decile_sensitivity_summary.tsv"), sep = "\t")
  fwrite(as.data.table(sens_stats), file.path(tables_dir, "decile_sensitivity_stats.tsv"), sep = "\t")

  # Recreation-style plot: left axis mean h2 raw, right axis fraction h2 raw > cutoff
  dual_tbl <- sens_decile_summary %>%
    group_by(decile_source, decile_direction) %>%
    mutate(
      h2_min = min(mean_h2_raw, na.rm = TRUE),
      h2_max = max(mean_h2_raw, na.rm = TRUE),
      frac_min = min(prop_h2_raw_gt_frac_cutoff, na.rm = TRUE),
      frac_max = max(prop_h2_raw_gt_frac_cutoff, na.rm = TRUE),
      frac_scaled = if_else(
        is.finite(frac_max - frac_min) & (frac_max - frac_min) > 0,
        h2_min + (prop_h2_raw_gt_frac_cutoff - frac_min) * (h2_max - h2_min) / (frac_max - frac_min),
        h2_min
      )
    ) %>%
    ungroup()

  p_dual <- ggplot(dual_tbl, aes(x = post_mean_bin)) +
    geom_line(aes(y = mean_h2_raw), color = "#2A6F9E", linewidth = 1) +
    geom_point(aes(y = mean_h2_raw), color = "#2A6F9E", size = 1.8, alpha = 0.8) +
    geom_line(aes(y = frac_scaled), color = "#B22222", linewidth = 1) +
    geom_point(aes(y = frac_scaled), color = "#B22222", size = 1.8, alpha = 0.8) +
    scale_x_continuous(breaks = 1:10) +
    facet_wrap(~decile_source, scales = "free_y", ncol = 3) +
    labs(
      title = "Decile sensitivity recreation view",
      subtitle = paste0(suite_subtitle, " | blue = mean h2 RAW, red = fraction h2 TMM RAW > ", h2_frac_cutoff),
      x = "Decile of selective constraint",
      y = "Expression h2 (RAW); decile 10 = highest constraint"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(
    filename = file.path(plots_dir, "decile_sensitivity_recreation_dual_axis_style.png"),
    plot = p_dual,
    width = 14,
    height = 7,
    dpi = 300
  )

  # TPM recreation-style plot: blue = mean TPM, red = fraction high-h2 (rescaled)
  dual_tpm_tbl <- sens_decile_summary %>%
    group_by(decile_source) %>%
    mutate(
      tpm_min = min(mean_tpm, na.rm = TRUE),
      tpm_max = max(mean_tpm, na.rm = TRUE),
      frac_min = min(prop_h2_raw_gt_frac_cutoff, na.rm = TRUE),
      frac_max = max(prop_h2_raw_gt_frac_cutoff, na.rm = TRUE),
      frac_scaled_to_tpm = if_else(
        is.finite(frac_max - frac_min) & (frac_max - frac_min) > 0,
        tpm_min + (prop_h2_raw_gt_frac_cutoff - frac_min) * (tpm_max - tpm_min) / (frac_max - frac_min),
        tpm_min
      )
    ) %>%
    ungroup()

  p_dual_tpm <- ggplot(dual_tpm_tbl, aes(x = post_mean_bin)) +
    geom_line(aes(y = mean_tpm), color = "#2A6F9E", linewidth = 1) +
    geom_point(aes(y = mean_tpm), color = "#2A6F9E", size = 1.8, alpha = 0.8) +
    geom_line(aes(y = frac_scaled_to_tpm), color = "#B22222", linewidth = 1) +
    geom_point(aes(y = frac_scaled_to_tpm), color = "#B22222", size = 1.8, alpha = 0.8) +
    scale_x_continuous(breaks = 1:10) +
    facet_wrap(~decile_source, scales = "free_y", ncol = 3) +
    labs(
      title = "Decile sensitivity recreation view (TPM)",
      subtitle = paste0(suite_subtitle, " | blue = mean TPM, red = fraction h2 TMM RAW > ", h2_frac_cutoff),
      x = "Decile of selective constraint",
      y = "Mean TPM; decile 10 = highest constraint"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(
    filename = file.path(plots_dir, "decile_sensitivity_recreation_dual_axis_style_tpm.png"),
    plot = p_dual_tpm,
    width = 14,
    height = 6.5,
    dpi = 300
  )

  h2_lines <- sens_decile_summary %>%
    select(decile_source, decile_direction, post_mean_bin, mean_h2_raw, mean_h2_irnt) %>%
    pivot_longer(
      cols = c(mean_h2_raw, mean_h2_irnt),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = recode(
        metric,
        mean_h2_raw = "Mean h2 RAW",
        mean_h2_irnt = "Mean h2 IRNT"
      )
    )

  p_h2_sens <- ggplot(h2_lines, aes(x = factor(post_mean_bin), y = value, color = metric, group = metric)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.7) +
    facet_wrap(~decile_source, scales = "free_y", ncol = 3) +
    labs(
      title = "Mean h2 trend by decile definition",
      subtitle = suite_subtitle,
      x = "s_het decile (1-10)",
      y = "h2",
      color = "Metric"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(
    filename = file.path(plots_dir, "decile_sensitivity_mean_h2_lines.png"),
    plot = p_h2_sens,
    width = 14,
    height = 7,
    dpi = 300
  )
}

plot_h2_shet_single <- function(df, h2_col, h2_label, out_file) {
  plot_df <- df %>%
    transmute(
      post_mean = post_mean,
      post_mean_bin = post_mean_bin,
      h2 = .data[[h2_col]]
    ) %>%
    filter(is.finite(post_mean), is.finite(h2))

  decile_tbl <- plot_df %>%
    group_by(post_mean_bin) %>%
    summarise(
      mean_h2 = mean(h2, na.rm = TRUE),
      .groups = "drop"
    )

  p <- ggplot(decile_tbl, aes(x = factor(post_mean_bin), y = mean_h2, group = 1)) +
    geom_line(linewidth = 1, color = "#1D3557") +
    geom_point(size = 1.8, color = "#1D3557") +
    labs(
      title = paste0("TMM h2 ", h2_label, ": mean by s_het decile"),
      subtitle = "Deciles recalculated after merge (1 = lowest, 10 = highest s_het)",
      x = "s_het post_mean decile",
      y = "h2_GREML"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(out_file, p, width = 8, height = 5, dpi = 300)
}

run_plot_suite_single <- function(df, h2_col, h2_label, suite_name, suite_subtitle) {
  suite_root <- file.path(analysis_root, suite_name)
  plots_dir <- file.path(suite_root, "plots")
  tables_dir <- file.path(suite_root, "tables")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

  suite_tbl <- df %>%
    mutate(
      post_mean_bin = ntile(post_mean, n_deciles),
      mean_tpm_decile = ntile(mean_tpm, n_deciles),
      median_tpm_decile = ntile(median_tpm, n_deciles)
    )

  p1_decile_tbl <- suite_tbl %>%
    group_by(post_mean_bin) %>%
    summarise(
      n_genes = n(),
      mean_tpm = mean(mean_tpm, na.rm = TRUE),
      median_tpm = median(median_tpm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(post_mean_bin)

  p1_decile_long <- p1_decile_tbl %>%
    pivot_longer(cols = c(mean_tpm, median_tpm), names_to = "metric", values_to = "value") %>%
    mutate(metric = recode(metric, mean_tpm = "Mean TPM", median_tpm = "Median TPM"))

  p1_decile <- ggplot(p1_decile_long, aes(x = factor(post_mean_bin), y = value, color = metric, group = metric)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    labs(
      title = "TPM by s_het decile",
      subtitle = paste0(suite_subtitle, " | h2 branch: ", h2_label),
      x = "s_het post_mean decile (1-10)",
      y = "TPM",
      color = "TPM summary"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(plots_dir, "plot1_shet_decile_vs_tpm_mean_median.png"), p1_decile, width = 8, height = 5, dpi = 300)
  fwrite(as.data.table(p1_decile_tbl), file.path(tables_dir, "plot1_shet_decile_vs_tpm_mean_median.tsv"), sep = "\t")

  plot_h2_shet_single(
    df = suite_tbl,
    h2_col = h2_col,
    h2_label = h2_label,
    out_file = file.path(plots_dir, "plot3_shet_decile_vs_h2_single_norm_mean_median_lines.png")
  )

  fwrite(as.data.table(suite_tbl), file.path(tables_dir, "shet_tpm_h2_merged_gene_level.tsv"), sep = "\t")
}

run_decile_sensitivity_single <- function(df, h2_col, h2_label, suite_name, suite_subtitle, shet_reference_tbl, shet_tpm_overlap_tbl) {
  suite_root <- file.path(analysis_root, suite_name)
  plots_dir <- file.path(suite_root, "plots")
  tables_dir <- file.path(suite_root, "tables")

  decile_maps <- bind_rows(
    make_decile_map(df, "merged_h2_genes", "ascending"),
    make_decile_map(shet_reference_tbl, "global_shet_genes", "ascending"),
    make_decile_map(shet_tpm_overlap_tbl, "shet_tpm_overlap_genes", "ascending")
  )

  sens_tbl <- decile_maps %>%
    inner_join(df, by = "Gene") %>%
    filter(is.finite(post_mean_bin), post_mean_bin >= 1L, post_mean_bin <= n_deciles)

  if (nrow(sens_tbl) == 0L) return(invisible(NULL))

  sens_decile_summary <- sens_tbl %>%
    group_by(decile_source, decile_direction, post_mean_bin) %>%
    summarise(
      n_genes = n(),
      mean_tpm = mean(mean_tpm, na.rm = TRUE),
      mean_h2 = mean(.data[[h2_col]], na.rm = TRUE),
      median_h2 = median(.data[[h2_col]], na.rm = TRUE),
      prop_h2_gt_frac_cutoff = mean(.data[[h2_col]] > h2_frac_cutoff, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(decile_source, decile_direction, post_mean_bin)

  fwrite(as.data.table(sens_decile_summary), file.path(tables_dir, "decile_sensitivity_summary_single_norm.tsv"), sep = "\t")

  dual_tbl <- sens_decile_summary %>%
    group_by(decile_source) %>%
    mutate(
      h2_min = min(mean_h2, na.rm = TRUE),
      h2_max = max(mean_h2, na.rm = TRUE),
      frac_min = min(prop_h2_gt_frac_cutoff, na.rm = TRUE),
      frac_max = max(prop_h2_gt_frac_cutoff, na.rm = TRUE),
      frac_scaled = if_else(
        is.finite(frac_max - frac_min) & (frac_max - frac_min) > 0,
        h2_min + (prop_h2_gt_frac_cutoff - frac_min) * (h2_max - h2_min) / (frac_max - frac_min),
        h2_min
      )
    ) %>%
    ungroup()

  p_dual <- ggplot(dual_tbl, aes(x = post_mean_bin)) +
    geom_line(aes(y = mean_h2), color = "#2A6F9E", linewidth = 1) +
    geom_point(aes(y = mean_h2), color = "#2A6F9E", size = 1.8, alpha = 0.8) +
    geom_line(aes(y = frac_scaled), color = "#B22222", linewidth = 1) +
    geom_point(aes(y = frac_scaled), color = "#B22222", size = 1.8, alpha = 0.8) +
    scale_x_continuous(breaks = 1:10) +
    facet_wrap(~decile_source, scales = "free_y", ncol = 3) +
    labs(
      title = paste0("Decile sensitivity recreation view (", h2_label, ")"),
      subtitle = paste0(suite_subtitle, " | blue = mean h2, red = fraction h2 > ", h2_frac_cutoff),
      x = "Decile of selective constraint",
      y = "Expression h2; decile 10 = highest constraint"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(plots_dir, "decile_sensitivity_recreation_single_norm.png"), p_dual, width = 14, height = 6.5, dpi = 300)
}

run_h2_source_comparison <- function(df, h2_tmm_col, h2_tpm_col, norm_label, suite_name, suite_subtitle) {
  suite_root <- file.path(analysis_root, suite_name)
  plots_dir <- file.path(suite_root, "plots")
  tables_dir <- file.path(suite_root, "tables")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

  tbl <- df %>%
    mutate(post_mean_bin = ntile(post_mean, n_deciles)) %>%
    filter(is.finite(.data[[h2_tmm_col]]), is.finite(.data[[h2_tpm_col]]))

  if (nrow(tbl) == 0L) {
    message("Skipping h2 source comparison for ", norm_label, ": no overlap between TPM and TMM h2.")
    return(invisible(NULL))
  }

  long_decile <- tbl %>%
    select(Gene, post_mean_bin, all_of(c(h2_tmm_col, h2_tpm_col))) %>%
    pivot_longer(
      cols = all_of(c(h2_tmm_col, h2_tpm_col)),
      names_to = "h2_source",
      values_to = "h2"
    ) %>%
    mutate(
      h2_source = recode(h2_source, !!h2_tmm_col := "TMM", !!h2_tpm_col := "TPM")
    )

  decile_tbl <- long_decile %>%
    group_by(post_mean_bin, h2_source) %>%
    summarise(
      mean_h2 = mean(h2, na.rm = TRUE),
      n_genes = n(),
      .groups = "drop"
    )

  p_decile <- ggplot(
    decile_tbl,
    aes(
      x = factor(post_mean_bin),
      y = mean_h2,
      color = h2_source,
      group = h2_source
    )
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.8) +
    labs(
      title = paste0("h2 source comparison by s_het decile (", norm_label, ")"),
      subtitle = suite_subtitle,
      x = "s_het post_mean decile (1-10)",
      y = "h2_GREML",
      color = "Expression source"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggsave(file.path(plots_dir, paste0("plot_h2_source_decile_", tolower(norm_label), ".png")), p_decile, width = 8, height = 5, dpi = 300)

  fwrite(as.data.table(tbl), file.path(tables_dir, paste0("h2_source_comparison_", tolower(norm_label), "_gene_level.tsv")), sep = "\t")
}

run_shet_h2_triplet <- function(
  tmm_tbl,
  tpm_tbl,
  overlap_tbl,
  h2_tmm_col,
  h2_tpm_col,
  norm_label,
  suite_name
) {
  suite_root <- file.path(analysis_root, suite_name)
  plots_dir <- file.path(suite_root, "plots")
  tables_dir <- file.path(suite_root, "tables")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  norm_key <- tolower(norm_label)

  tmm_plot_tbl <- tmm_tbl %>%
    mutate(post_mean_bin = ntile(post_mean, n_deciles)) %>%
    transmute(Gene, post_mean, post_mean_bin, h2_tmm = .data[[h2_tmm_col]]) %>%
    filter(is.finite(post_mean), is.finite(h2_tmm))

  tpm_plot_tbl <- tpm_tbl %>%
    mutate(post_mean_bin = ntile(post_mean, n_deciles)) %>%
    transmute(Gene, post_mean, post_mean_bin, h2_tpm = .data[[h2_tpm_col]]) %>%
    filter(is.finite(post_mean), is.finite(h2_tpm))

  overlap_long <- overlap_tbl %>%
    mutate(post_mean_bin = ntile(post_mean, n_deciles)) %>%
    transmute(Gene, post_mean, post_mean_bin, h2_tmm = .data[[h2_tmm_col]], h2_tpm = .data[[h2_tpm_col]]) %>%
    filter(is.finite(post_mean), is.finite(h2_tmm), is.finite(h2_tpm)) %>%
    pivot_longer(
      cols = c(h2_tmm, h2_tpm),
      names_to = "source",
      values_to = "h2"
    ) %>%
    mutate(source = recode(source, h2_tmm = "TMM", h2_tpm = "TPM"))

  tmm_decile <- tmm_plot_tbl %>%
    group_by(post_mean_bin) %>%
    summarise(
      mean_h2 = mean(h2_tmm, na.rm = TRUE),
      .groups = "drop"
    )

  p_tmm <- ggplot(tmm_decile, aes(x = factor(post_mean_bin), y = mean_h2, group = 1)) +
    geom_line(linewidth = 1, color = "#1D3557") +
    geom_point(size = 1.8, color = "#1D3557") +
    labs(
      title = paste0("s_het decile vs h2 (TMM ", toupper(norm_label), ")"),
      subtitle = "Deciles recalculated after TMM merge",
      x = "s_het post_mean decile (1-10)",
      y = "h2_GREML"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggsave(file.path(plots_dir, paste0("plot1_shet_vs_h2_tmm_", norm_key, ".png")), p_tmm, width = 8, height = 5, dpi = 300)

  tpm_decile <- tpm_plot_tbl %>%
    group_by(post_mean_bin) %>%
    summarise(
      mean_h2 = mean(h2_tpm, na.rm = TRUE),
      .groups = "drop"
    )

  p_tpm <- ggplot(tpm_decile, aes(x = factor(post_mean_bin), y = mean_h2, group = 1)) +
    geom_line(linewidth = 1, color = "#2A6F9E") +
    geom_point(size = 1.8, color = "#2A6F9E") +
    labs(
      title = paste0("s_het decile vs h2 (TPM ", toupper(norm_label), ")"),
      subtitle = "Deciles recalculated after TPM merge",
      x = "s_het post_mean decile (1-10)",
      y = "h2_GREML"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggsave(file.path(plots_dir, paste0("plot2_shet_vs_h2_tpm_", norm_key, ".png")), p_tpm, width = 8, height = 5, dpi = 300)

  overlap_decile <- overlap_long %>%
    group_by(post_mean_bin, source) %>%
    summarise(
      mean_h2 = mean(h2, na.rm = TRUE),
      .groups = "drop"
    )

  p_overlap <- ggplot(
    overlap_decile,
    aes(x = factor(post_mean_bin), y = mean_h2, color = source, group = source)
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.8) +
    labs(
      title = paste0("s_het decile vs h2 (overlap TPM+TMM ", toupper(norm_label), ")"),
      subtitle = "Deciles recalculated on overlapping genes only",
      x = "s_het post_mean decile (1-10)",
      y = "h2_GREML",
      color = "h2 source"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggsave(file.path(plots_dir, paste0("plot3_shet_vs_h2_overlap_tpm_tmm_", norm_key, ".png")), p_overlap, width = 8, height = 5, dpi = 300)

  fwrite(as.data.table(tmm_plot_tbl), file.path(tables_dir, paste0("shet_vs_h2_tmm_", norm_key, "_gene_level.tsv")), sep = "\t")
  fwrite(as.data.table(tpm_plot_tbl), file.path(tables_dir, paste0("shet_vs_h2_tpm_", norm_key, "_gene_level.tsv")), sep = "\t")
  fwrite(as.data.table(overlap_long), file.path(tables_dir, paste0("shet_vs_h2_overlap_tpm_tmm_", norm_key, "_gene_level.tsv")), sep = "\t")
}

# --------------------------------------------------
# Main
# --------------------------------------------------
raw_run_dir <- find_run_dir(method_h2_raw, require_summary = TRUE)
irnt_run_dir <- find_run_dir(method_h2_irnt, require_summary = TRUE)
tpm_h2_raw_run_dir <- find_run_dir(method_h2_tpm_raw, require_summary = TRUE)
tpm_h2_irnt_run_dir <- find_run_dir(method_h2_tpm_irnt, require_summary = TRUE)
tpm_run_dir <- find_run_dir(tpm_mean_method, require_summary = FALSE)

message("Using TMM-RAW h2 run:  ", raw_run_dir, " (method=", method_h2_raw, ")")
message("Using TMM-IRNT h2 run: ", irnt_run_dir, " (method=", method_h2_irnt, ")")
message("Using TPM-RAW h2 run:  ", tpm_h2_raw_run_dir, " (method=", method_h2_tpm_raw, ")")
message("Using TPM-IRNT h2 run: ", tpm_h2_irnt_run_dir, " (method=", method_h2_tpm_irnt, ")")
message("Using TPM run:  ", tpm_run_dir)

h2_method_info <- tibble(
  setting = c(
    "force_tmm_h2",
    "h2_frac_cutoff",
    "input_method_raw",
    "input_method_irnt",
    "h2_method_raw_used",
    "h2_method_irnt_used",
    "h2_method_tpm_raw_used",
    "h2_method_tpm_irnt_used",
    "h2_raw_run_dir",
    "h2_irnt_run_dir",
    "h2_tpm_raw_run_dir",
    "h2_tpm_irnt_run_dir"
  ),
  value = c(
    as.character(force_tmm_h2),
    as.character(h2_frac_cutoff),
    method_raw,
    method_irnt,
    method_h2_raw,
    method_h2_irnt,
    method_h2_tpm_raw,
    method_h2_tpm_irnt,
    raw_run_dir,
    irnt_run_dir,
    tpm_h2_raw_run_dir,
    tpm_h2_irnt_run_dir
  )
)
fwrite(as.data.table(h2_method_info), file.path(analysis_root, "h2_method_info.tsv"), sep = "\t")

if (!file.exists(shet_xlsx)) {
  stop("SHET_XLSX not found: ", shet_xlsx)
}

shet_tbl <- read_xlsx(shet_xlsx, sheet = shet_sheet) %>%
  as_tibble()
if (!all(c("ensg", "post_mean") %in% names(shet_tbl))) {
  stop("S_het sheet must include columns named 'ensg' and 'post_mean'.")
}

shet_tbl <- shet_tbl %>%
  transmute(
    Gene = normalize_gene_id(ensg),
    post_mean = suppressWarnings(as.numeric(post_mean))
  ) %>%
  filter(!is.na(Gene), Gene != "", is.finite(post_mean)) %>%
  group_by(Gene) %>%
  summarise(post_mean = mean(post_mean, na.rm = TRUE), .groups = "drop")

tpm_metrics_tbl <- read_tpm_metrics_by_gene(tpm_run_dir)
h2_tmm_raw_tbl <- read_h2_by_gene(raw_run_dir, h2_col_name = "h2_tmm_raw")
h2_tmm_irnt_tbl <- read_h2_by_gene(irnt_run_dir, h2_col_name = "h2_tmm_irnt")
h2_tpm_raw_tbl <- read_h2_by_gene(tpm_h2_raw_run_dir, h2_col_name = "h2_tpm_raw")
h2_tpm_irnt_tbl <- read_h2_by_gene(tpm_h2_irnt_run_dir, h2_col_name = "h2_tpm_irnt")
shet_tpm_overlap_tbl <- shet_tbl %>%
  inner_join(tpm_metrics_tbl %>% select(Gene), by = "Gene")

merged_raw_tbl <- shet_tbl %>%
  inner_join(tpm_metrics_tbl, by = "Gene") %>%
  inner_join(h2_tmm_raw_tbl, by = "Gene") %>%
  mutate(post_mean_bin_merged = ntile(post_mean, n_deciles))

merged_irnt_tbl <- shet_tbl %>%
  inner_join(tpm_metrics_tbl, by = "Gene") %>%
  inner_join(h2_tmm_irnt_tbl, by = "Gene") %>%
  mutate(post_mean_bin_merged = ntile(post_mean, n_deciles))

merged_source_raw_tbl <- merged_raw_tbl %>%
  inner_join(h2_tpm_raw_tbl, by = "Gene") %>%
  mutate(post_mean_bin_merged = ntile(post_mean, n_deciles))

merged_source_irnt_tbl <- merged_irnt_tbl %>%
  inner_join(h2_tpm_irnt_tbl, by = "Gene") %>%
  mutate(post_mean_bin_merged = ntile(post_mean, n_deciles))

merged_tpm_raw_tbl <- shet_tbl %>%
  inner_join(tpm_metrics_tbl, by = "Gene") %>%
  inner_join(h2_tpm_raw_tbl, by = "Gene") %>%
  mutate(post_mean_bin_merged = ntile(post_mean, n_deciles))

merged_tpm_irnt_tbl <- shet_tbl %>%
  inner_join(tpm_metrics_tbl, by = "Gene") %>%
  inner_join(h2_tpm_irnt_tbl, by = "Gene") %>%
  mutate(post_mean_bin_merged = ntile(post_mean, n_deciles))

merged_pair_tbl <- merged_raw_tbl %>%
  inner_join(
    merged_irnt_tbl %>% select(Gene, h2_tmm_irnt),
    by = "Gene"
  ) %>%
  mutate(post_mean_bin_merged = ntile(post_mean, n_deciles))

gene_count_audit <- tibble(
  set = c(
    "shet_tbl",
    "tpm_metrics_tbl",
    "h2_tmm_raw_tbl",
    "h2_tmm_irnt_tbl",
    "h2_tpm_raw_tbl",
    "h2_tpm_irnt_tbl",
    "merged_raw_tbl",
    "merged_irnt_tbl",
    "merged_tpm_raw_tbl",
    "merged_tpm_irnt_tbl",
    "merged_source_raw_tbl",
    "merged_source_irnt_tbl",
    "merged_pair_tbl"
  ),
  n_genes = c(
    nrow(shet_tbl),
    nrow(tpm_metrics_tbl),
    nrow(h2_tmm_raw_tbl),
    nrow(h2_tmm_irnt_tbl),
    nrow(h2_tpm_raw_tbl),
    nrow(h2_tpm_irnt_tbl),
    nrow(merged_raw_tbl),
    nrow(merged_irnt_tbl),
    nrow(merged_tpm_raw_tbl),
    nrow(merged_tpm_irnt_tbl),
    nrow(merged_source_raw_tbl),
    nrow(merged_source_irnt_tbl),
    nrow(merged_pair_tbl)
  )
)
fwrite(as.data.table(gene_count_audit), file.path(analysis_root, "gene_count_audit.tsv"), sep = "\t")

if (nrow(merged_raw_tbl) == 0L) {
  stop("No overlapping genes after merging Shet + TPM + TMM RAW h2.")
}

if (nrow(merged_irnt_tbl) == 0L) {
  stop("No overlapping genes after merging Shet + TPM + TMM IRNT h2.")
}

if (nrow(merged_tpm_raw_tbl) == 0L) {
  stop("No overlapping genes after merging Shet + TPM + TPM RAW h2.")
}

if (nrow(merged_tpm_irnt_tbl) == 0L) {
  stop("No overlapping genes after merging Shet + TPM + TPM IRNT h2.")
}

if (nrow(merged_source_raw_tbl) > 0L) {
  run_shet_h2_triplet(
    tmm_tbl = merged_raw_tbl,
    tpm_tbl = merged_tpm_raw_tbl,
    overlap_tbl = merged_source_raw_tbl,
    h2_tmm_col = "h2_tmm_raw",
    h2_tpm_col = "h2_tpm_raw",
    norm_label = "RAW",
    suite_name = "shet_vs_h2_triplet_raw"
  )
} else {
  message("Skipping Shet-vs-h2 RAW triplet: no overlap between TMM and TPM RAW h2 tables.")
}

if (nrow(merged_source_irnt_tbl) > 0L) {
  run_shet_h2_triplet(
    tmm_tbl = merged_irnt_tbl,
    tpm_tbl = merged_tpm_irnt_tbl,
    overlap_tbl = merged_source_irnt_tbl,
    h2_tmm_col = "h2_tmm_irnt",
    h2_tpm_col = "h2_tpm_irnt",
    norm_label = "IRNT",
    suite_name = "shet_vs_h2_triplet_irnt"
  )
} else {
  message("Skipping Shet-vs-h2 IRNT triplet: no overlap between TMM and TPM IRNT h2 tables.")
}

run_plot_suite_single(
  df = merged_raw_tbl,
  h2_col = "h2_tmm_raw",
  h2_label = "RAW",
  suite_name = "raw_only_all_genes",
  suite_subtitle = "All PASS genes with Shet+TPM+TMM RAW h2 overlap"
)

run_decile_sensitivity_single(
  df = merged_raw_tbl,
  h2_col = "h2_tmm_raw",
  h2_label = "RAW",
  suite_name = "raw_only_all_genes",
  suite_subtitle = "All PASS genes with Shet+TPM+TMM RAW h2 overlap",
  shet_reference_tbl = shet_tbl,
  shet_tpm_overlap_tbl = shet_tpm_overlap_tbl
)

run_plot_suite_single(
  df = merged_irnt_tbl,
  h2_col = "h2_tmm_irnt",
  h2_label = "IRNT",
  suite_name = "irnt_only_all_genes",
  suite_subtitle = "All PASS genes with Shet+TPM+TMM IRNT h2 overlap"
)

run_decile_sensitivity_single(
  df = merged_irnt_tbl,
  h2_col = "h2_tmm_irnt",
  h2_label = "IRNT",
  suite_name = "irnt_only_all_genes",
  suite_subtitle = "All PASS genes with Shet+TPM+TMM IRNT h2 overlap",
  shet_reference_tbl = shet_tbl,
  shet_tpm_overlap_tbl = shet_tpm_overlap_tbl
)

if (nrow(merged_source_raw_tbl) > 0L) {
  run_h2_source_comparison(
    df = merged_source_raw_tbl,
    h2_tmm_col = "h2_tmm_raw",
    h2_tpm_col = "h2_tpm_raw",
    norm_label = "RAW",
    suite_name = "h2_source_compare_raw_all_genes",
    suite_subtitle = "Genes with Shet+TPM+both TMM/TPM RAW h2 overlap"
  )
} else {
  message("Skipping RAW h2 source comparison: no overlap between TMM and TPM RAW h2.")
}

if (nrow(merged_source_irnt_tbl) > 0L) {
  run_h2_source_comparison(
    df = merged_source_irnt_tbl,
    h2_tmm_col = "h2_tmm_irnt",
    h2_tpm_col = "h2_tpm_irnt",
    norm_label = "IRNT",
    suite_name = "h2_source_compare_irnt_all_genes",
    suite_subtitle = "Genes with Shet+TPM+both TMM/TPM IRNT h2 overlap"
  )
} else {
  message("Skipping IRNT h2 source comparison: no overlap between TMM and TPM IRNT h2.")
}

if (nrow(merged_pair_tbl) > 0L) {
  run_plot_suite_pair(
    df = merged_pair_tbl,
    suite_name = "raw_irnt_overlap_all_genes",
    suite_subtitle = "Genes with Shet+TPM+both TMM RAW/IRNT h2 overlap"
  )

  run_decile_sensitivity_pair(
    df = merged_pair_tbl,
    suite_name = "raw_irnt_overlap_all_genes",
    suite_subtitle = "Genes with Shet+TPM+both TMM RAW/IRNT h2 overlap",
    shet_reference_tbl = shet_tbl,
    shet_tpm_overlap_tbl = shet_tpm_overlap_tbl
  )
} else {
  message("Skipping overlap RAW-vs-IRNT suite: no genes shared between RAW and IRNT merges.")
}
