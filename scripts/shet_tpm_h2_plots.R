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

tpm_mean_method <- Sys.getenv("TPM_MEAN_METHOD", "")
if (!nzchar(tpm_mean_method)) {
  m <- str_match(method_raw, "^(all_snps|hm3_no_mhc)_(tmm|tpm)_(raw|irnt|inverse_normal)$")
  if (!is.na(m[1, 1])) {
    tpm_mean_method <- paste(m[1, 2], "tpm", "raw", sep = "_")
  } else {
    tpm_mean_method <- "all_snps_tpm_raw"
  }
}

shet_x_mode <- tolower(Sys.getenv("SHET_X_MODE", "both")) # one of: actual, decile, both
if (!shet_x_mode %in% c("actual", "decile", "both")) {
  stop("SHET_X_MODE must be one of: actual, decile, both")
}

n_deciles <- suppressWarnings(as.integer(Sys.getenv("N_DECILES", "10")))
if (!is.finite(n_deciles) || n_deciles != 10L) {
  stop("N_DECILES must be 10 for decile-based plots.")
}

cor_method <- tolower(Sys.getenv("COR_METHOD", "spearman"))
if (!cor_method %in% c("pearson", "spearman", "kendall")) {
  stop("COR_METHOD must be one of: pearson, spearman, kendall")
}

h2_low_cutoff <- as.numeric(Sys.getenv("H2_LO_CUTOFF", "0.005"))
if (!is.finite(h2_low_cutoff) || h2_low_cutoff < 0) {
  stop("H2_LO_CUTOFF must be a non-negative numeric value.")
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

plot_h2_shet <- function(df, x_mode, out_file) {
  h2_long <- df %>%
    select(post_mean, post_mean_bin, h2_raw, h2_irnt) %>%
    pivot_longer(
      cols = c(h2_raw, h2_irnt),
      names_to = "h2_type",
      values_to = "h2"
    ) %>%
    mutate(h2_type = recode(h2_type, h2_raw = "RAW", h2_irnt = "IRNT"))

  if (x_mode == "actual") {
    p <- ggplot(h2_long, aes(x = post_mean, y = h2)) +
      geom_point(aes(color = h2_type), alpha = 0.35, size = 0.7) +
      geom_smooth(aes(color = h2_type), method = "loess", se = FALSE, linewidth = 1) +
      facet_wrap(~h2_type, ncol = 2, scales = "free_y") +
      labs(
        title = "h2 (RAW vs IRNT) across s_het post_mean",
        subtitle = "Gene-level scatter with LOESS trend",
        x = "s_het post_mean",
        y = "h2_GREML",
        color = "Normalization"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold")
      ) +
      scale_color_manual(values = c("RAW" = "#1D3557", "IRNT" = "#E76F51"))
  } else {
    decile_tbl <- h2_long %>%
      group_by(post_mean_bin, h2_type) %>%
      summarise(
        mean_h2 = mean(h2, na.rm = TRUE),
        median_h2 = median(h2, na.rm = TRUE),
        n_genes = n(),
        .groups = "drop"
      )

    decile_long <- decile_tbl %>%
      pivot_longer(
        cols = c(mean_h2, median_h2),
        names_to = "stat_type",
        values_to = "h2_value"
      ) %>%
      mutate(
        stat_type = recode(stat_type, mean_h2 = "Mean h2", median_h2 = "Median h2")
      )

    p <- ggplot(
      decile_long,
      aes(
        x = factor(post_mean_bin),
        y = h2_value,
        color = h2_type,
        linetype = stat_type,
        group = interaction(h2_type, stat_type)
      )
    ) +
      geom_line(linewidth = 1) +
      geom_point(size = 1.8) +
      labs(
        title = "h2 (RAW vs IRNT): mean and median by s_het decile",
        subtitle = "Each point is a decile-level summary",
        x = "s_het post_mean decile",
        y = "h2_GREML",
        color = "Normalization",
        linetype = "Summary"
      ) +
      theme_minimal(base_size = 11) +
      theme(panel.grid.minor = element_blank())
  }

  ggsave(out_file, p, width = 8, height = 5, dpi = 300)
}

run_plot_suite <- function(df, suite_name, suite_subtitle) {
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
        cor(h2_raw, h2_irnt, method = cor_method, use = "pairwise.complete.obs"),
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
        cor(h2_raw, h2_irnt, method = cor_method, use = "pairwise.complete.obs"),
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
      title = "Correlation(h2 RAW, h2 IRNT) by TPM decile",
      subtitle = paste0(suite_subtitle, " | correlation method: ", toupper(cor_method)),
      x = "TPM decile (1-10)",
      y = "Correlation(h2 RAW, h2 IRNT)",
      color = "Decile definition"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(plots_dir, "plot2_tpm_decile_vs_h2_raw_irnt_correlation_mean_median.png"), p2, width = 8, height = 5, dpi = 300)
  fwrite(as.data.table(p2_tbl), file.path(tables_dir, "plot2_tpm_decile_vs_h2_raw_irnt_correlation_mean_median.tsv"), sep = "\t")

  # Plot 3: Shet vs RAW+IRNT h2
  if (shet_x_mode %in% c("actual", "both")) {
    plot_h2_shet(
      df = suite_tbl,
      x_mode = "actual",
      out_file = file.path(plots_dir, "plot3_shet_actual_vs_h2_raw_irnt_lines.png")
    )
  }

  if (shet_x_mode %in% c("decile", "both")) {
    plot_h2_shet(
      df = suite_tbl,
      x_mode = "decile",
      out_file = file.path(plots_dir, "plot3_shet_decile_vs_h2_raw_irnt_mean_median_lines.png")
    )
  }

  fwrite(as.data.table(suite_tbl), file.path(tables_dir, "shet_tpm_h2_merged_gene_level.tsv"), sep = "\t")

  message("Saved outputs in: ", suite_root)
  message("Plots:")
  message(" - ", file.path(plots_dir, "plot1_shet_decile_vs_tpm_mean_median.png"))
  message(" - ", file.path(plots_dir, "plot2_tpm_decile_vs_h2_raw_irnt_correlation_mean_median.png"))
  message(" - ", file.path(plots_dir, "plot3_shet_actual_vs_h2_raw_irnt_lines.png"))
  message(" - ", file.path(plots_dir, "plot3_shet_decile_vs_h2_raw_irnt_mean_median_lines.png"))
  message("Tables:")
  message(" - ", file.path(tables_dir, "shet_tpm_h2_merged_gene_level.tsv"))
  message(" - ", file.path(tables_dir, "plot1_shet_decile_vs_tpm_mean_median.tsv"))
  message(" - ", file.path(tables_dir, "plot2_tpm_decile_vs_h2_raw_irnt_correlation_mean_median.tsv"))
}

# --------------------------------------------------
# Main
# --------------------------------------------------
raw_run_dir <- find_run_dir(method_raw, require_summary = TRUE)
irnt_run_dir <- find_run_dir(method_irnt, require_summary = TRUE)
tpm_run_dir <- find_run_dir(tpm_mean_method, require_summary = FALSE)

message("Using RAW run:  ", raw_run_dir)
message("Using IRNT run: ", irnt_run_dir)
message("Using TPM run:  ", tpm_run_dir)

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

h2_filter_note <- paste0("Low-h2 filter reference cutoff from prior scripts: h2 > ", format(h2_low_cutoff, scientific = FALSE))
message(h2_filter_note)

tpm_metrics_tbl <- read_tpm_metrics_by_gene(tpm_run_dir)
h2_raw_tbl <- read_h2_by_gene(raw_run_dir, h2_col_name = "h2_raw")
h2_irnt_tbl <- read_h2_by_gene(irnt_run_dir, h2_col_name = "h2_irnt")

merged_tbl <- shet_tbl %>%
  inner_join(tpm_metrics_tbl, by = "Gene") %>%
  inner_join(h2_raw_tbl, by = "Gene") %>%
  inner_join(h2_irnt_tbl, by = "Gene")

if (nrow(merged_tbl) == 0L) {
  stop("No overlapping genes after merging Shet, mean/median TPM, h2 RAW, and h2 IRNT.")
}

run_plot_suite(
  df = merged_tbl,
  suite_name = "all_genes",
  suite_subtitle = "All PASS genes with Shet+TPM+h2 overlap"
)

filtered_tbl <- merged_tbl %>%
  filter(h2_raw > h2_low_cutoff, h2_irnt > h2_low_cutoff)

if (nrow(filtered_tbl) < 100L) {
  message(
    "Low-h2 filtered set is small (n=", nrow(filtered_tbl), "). ",
    "Still writing outputs, but interpret trends with caution."
  )
}

if (nrow(filtered_tbl) >= 10L) {
  run_plot_suite(
    df = filtered_tbl,
    suite_name = paste0("h2_filtered_gt_", gsub("\\.", "p", format(h2_low_cutoff, scientific = FALSE))),
    suite_subtitle = paste0("Genes with h2 RAW > ", h2_low_cutoff, " and h2 IRNT > ", h2_low_cutoff)
  )
} else {
  message(
    "Skipping low-h2 filtered suite: too few genes after filtering (n=",
    nrow(filtered_tbl),
    ")."
  )
}
