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

tpm_cor_bins <- suppressWarnings(as.integer(Sys.getenv("TPM_COR_BINS", "10")))
if (!is.finite(tpm_cor_bins) || tpm_cor_bins < 3L) {
  stop("TPM_COR_BINS must be an integer >= 3")
}

cor_method <- tolower(Sys.getenv("COR_METHOD", "spearman"))
if (!cor_method %in% c("pearson", "spearman", "kendall")) {
  stop("COR_METHOD must be one of: pearson, spearman, kendall")
}

analysis_root <- file.path(runs_dir, "_analysis", "shet_tpm_h2_plots")
plots_dir <- file.path(analysis_root, "plots")
tables_dir <- file.path(analysis_root, "tables")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

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

read_mean_expression_by_gene <- function(run_dir, out_col_name = "mean_tpm") {
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
    value = colMeans(expr, na.rm = TRUE)
  ) %>%
    filter(is.finite(value)) %>%
    group_by(Gene) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    rename(!!out_col_name := value)
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
      stat_density_2d_filled(
        aes(fill = after_stat(level)),
        contour_var = "ndensity",
        alpha = 0.4,
        show.legend = FALSE
      ) +
      geom_point(aes(color = h2_type), alpha = 0.35, size = 0.7) +
      facet_wrap(~h2_type, ncol = 2, scales = "free_y") +
      labs(
        title = "h2 (RAW vs IRNT) across s_het post_mean",
        subtitle = "Dots + density clouds on gene-level values",
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

    p <- ggplot(decile_tbl, aes(x = factor(post_mean_bin), y = mean_h2, color = h2_type, group = h2_type)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      labs(
        title = "Mean h2 (RAW vs IRNT) by s_het decile",
        subtitle = "Each point is a decile-level mean",
        x = "s_het post_mean decile",
        y = "Mean h2_GREML",
        color = "Normalization"
      ) +
      theme_minimal(base_size = 11) +
      theme(panel.grid.minor = element_blank())
  }

  ggsave(out_file, p, width = 8, height = 5, dpi = 300)
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

mean_tpm_tbl <- read_mean_expression_by_gene(tpm_run_dir, out_col_name = "mean_tpm")
h2_raw_tbl <- read_h2_by_gene(raw_run_dir, h2_col_name = "h2_raw")
h2_irnt_tbl <- read_h2_by_gene(irnt_run_dir, h2_col_name = "h2_irnt")

merged_tbl <- shet_tbl %>%
  inner_join(mean_tpm_tbl, by = "Gene") %>%
  inner_join(h2_raw_tbl, by = "Gene") %>%
  inner_join(h2_irnt_tbl, by = "Gene") %>%
  mutate(
    post_mean_bin = ntile(post_mean, 10),
    tpm_bin = ntile(mean_tpm, tpm_cor_bins)
  )

if (nrow(merged_tbl) == 0L) {
  stop("No overlapping genes after merging Shet, mean TPM, h2 RAW, and h2 IRNT.")
}

if (n_distinct(merged_tbl$post_mean_bin) < 3L) {
  stop("Too few Shet bins after merge; check overlap and inputs.")
}

# Plot 1: Shet vs mean TPM
if (shet_x_mode %in% c("actual", "both")) {
  p1_actual <- ggplot(merged_tbl, aes(x = post_mean, y = mean_tpm)) +
    stat_density_2d_filled(
      aes(fill = after_stat(level)),
      contour_var = "ndensity",
      alpha = 0.45,
      show.legend = FALSE
    ) +
    geom_point(alpha = 0.3, size = 0.7, color = "#1D3557") +
    labs(
      title = "Mean TPM across s_het post_mean",
      subtitle = "Dots + density cloud",
      x = "s_het post_mean",
      y = "Mean TPM"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(plots_dir, "plot1_shet_actual_vs_mean_tpm.png"), p1_actual, width = 8, height = 5, dpi = 300)
}

if (shet_x_mode %in% c("decile", "both")) {
  p1_decile_tbl <- merged_tbl %>%
    group_by(post_mean_bin) %>%
    summarise(
      n_genes = n(),
      mean_tpm = mean(mean_tpm, na.rm = TRUE),
      median_tpm = median(mean_tpm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(post_mean_bin)

  p1_decile <- ggplot(p1_decile_tbl, aes(x = factor(post_mean_bin), y = mean_tpm, group = 1)) +
    geom_line(linewidth = 1, color = "#1D3557") +
    geom_point(size = 2, color = "#1D3557") +
    labs(
      title = "Mean TPM by s_het decile",
      subtitle = "Each point is a decile-level mean",
      x = "s_het post_mean decile",
      y = "Mean TPM"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(plots_dir, "plot1_shet_decile_vs_mean_tpm.png"), p1_decile, width = 8, height = 5, dpi = 300)
  fwrite(as.data.table(p1_decile_tbl), file.path(tables_dir, "plot1_shet_decile_vs_mean_tpm.tsv"), sep = "\t")
}

# Plot 2: mean TPM vs correlation(h2 RAW, h2 IRNT)
p2_tbl <- merged_tbl %>%
  group_by(tpm_bin) %>%
  summarise(
    n_genes = n(),
    mean_tpm = mean(mean_tpm, na.rm = TRUE),
    median_tpm = median(mean_tpm, na.rm = TRUE),
    cor_h2_raw_irnt = if_else(
      n() >= 3,
      cor(h2_raw, h2_irnt, method = cor_method, use = "pairwise.complete.obs"),
      as.numeric(NA)
    ),
    .groups = "drop"
  ) %>%
  arrange(mean_tpm)

p2 <- ggplot(p2_tbl, aes(x = mean_tpm, y = cor_h2_raw_irnt)) +
  geom_line(linewidth = 1, color = "#E76F51") +
  geom_point(size = 2, color = "#E76F51") +
  labs(
    title = "Correlation between h2 RAW and h2 IRNT across TPM bins",
    subtitle = paste0("Correlation method: ", toupper(cor_method), " | each point is a TPM bin"),
    x = "Mean TPM (bin mean)",
    y = "Correlation(h2 RAW, h2 IRNT)"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(plots_dir, "plot2_mean_tpm_vs_h2_raw_irnt_correlation.png"), p2, width = 8, height = 5, dpi = 300)
fwrite(as.data.table(p2_tbl), file.path(tables_dir, "plot2_mean_tpm_vs_h2_raw_irnt_correlation.tsv"), sep = "\t")

# Plot 3: Shet vs RAW+IRNT h2 lines
if (shet_x_mode %in% c("actual", "both")) {
  plot_h2_shet(
    df = merged_tbl,
    x_mode = "actual",
    out_file = file.path(plots_dir, "plot3_shet_actual_vs_h2_raw_irnt_lines.png")
  )
}

if (shet_x_mode %in% c("decile", "both")) {
  plot_h2_shet(
    df = merged_tbl,
    x_mode = "decile",
    out_file = file.path(plots_dir, "plot3_shet_decile_vs_h2_raw_irnt_lines.png")
  )
}

# Save merged gene-level table
fwrite(as.data.table(merged_tbl), file.path(tables_dir, "shet_tpm_h2_merged_gene_level.tsv"), sep = "\t")

message("Saved outputs in: ", analysis_root)
message("Plots:")
message(" - ", file.path(plots_dir, "plot1_shet_actual_vs_mean_tpm.png"))
message(" - ", file.path(plots_dir, "plot1_shet_decile_vs_mean_tpm.png"))
message(" - ", file.path(plots_dir, "plot2_mean_tpm_vs_h2_raw_irnt_correlation.png"))
message(" - ", file.path(plots_dir, "plot3_shet_actual_vs_h2_raw_irnt_lines.png"))
message(" - ", file.path(plots_dir, "plot3_shet_decile_vs_h2_raw_irnt_lines.png"))
message("Tables:")
message(" - ", file.path(tables_dir, "shet_tpm_h2_merged_gene_level.tsv"))
message(" - ", file.path(tables_dir, "plot1_shet_decile_vs_mean_tpm.tsv"))
message(" - ", file.path(tables_dir, "plot2_mean_tpm_vs_h2_raw_irnt_correlation.tsv"))
