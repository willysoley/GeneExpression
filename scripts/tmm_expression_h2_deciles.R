suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(magrittr)
})

# --------------------------------------------------
# Settings
# --------------------------------------------------
runs_dir <- Sys.getenv(
  "RUNS_DIR",
  "/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260428_GE_GEUVADIS_v2/GeneExpression/runs"
)
run_suffix <- Sys.getenv("RUN_SUFFIX", "_peerauto_pmg0_npc5")
min_log2_mean_tmm_plus1 <- as.numeric(Sys.getenv("MIN_LOG2_MEAN_TMM_PLUS1", "1.5"))
if (!is.finite(min_log2_mean_tmm_plus1)) {
  stop("MIN_LOG2_MEAN_TMM_PLUS1 must be numeric.")
}
min_log2_median_tpm_plus1 <- as.numeric(Sys.getenv("MIN_LOG2_MEDIAN_TPM_PLUS1", "0"))
if (!is.finite(min_log2_median_tpm_plus1)) {
  stop("MIN_LOG2_MEDIAN_TPM_PLUS1 must be numeric.")
}

# Comma-separated method IDs, e.g.:
# METHOD_IDS="all_snps_tmm_raw,hm3_no_mhc_tmm_raw"
method_ids_env <- Sys.getenv("METHOD_IDS", "all_snps_tmm_raw,all_snps_tmm_irnt")
method_ids <- method_ids_env %>%
  str_split(",") %>%
  unlist() %>%
  str_trim() %>%
  discard(~ .x == "") %>%
  unique()

if (length(method_ids) == 0L) stop("No METHOD_IDS provided.")

analysis_root <- file.path(runs_dir, "_analysis", "tmm_expression_h2_deciles")
plots_dir <- file.path(analysis_root, "plots")
tables_dir <- file.path(analysis_root, "tables")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

out_gene_tsv <- file.path(tables_dir, "tmm_expression_h2_deciles_gene_level.tsv")
out_summary_tsv <- file.path(tables_dir, "tmm_expression_h2_deciles_summary.tsv")
out_plot <- file.path(plots_dir, "tmm_expression_h2_deciles_boxplot.png")
out_median_tpm_gene_tsv <- file.path(tables_dir, "median_tpm_expression_h2_deciles_gene_level.tsv")
out_median_tpm_summary_tsv <- file.path(tables_dir, "median_tpm_expression_h2_deciles_summary.tsv")
out_median_tpm_plot <- file.path(plots_dir, "median_tpm_expression_h2_deciles_boxplot.png")
out_skewness_by_bin_tsv <- file.path(tables_dir, "expression_skewness_by_decile_gene_level.tsv")
out_skewness_by_bin_summary_tsv <- file.path(tables_dir, "expression_skewness_by_decile_summary.tsv")
out_skewness_by_bin_plot <- file.path(plots_dir, "expression_skewness_by_decile_boxplot.png")

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
  method_opts <- method_candidates(method)
  all_dirs <- list.dirs(runs_dir, recursive = FALSE, full.names = TRUE)

  hits <- map_dfr(method_opts, function(m) {
    patt <- paste0("^", m, run_suffix, "([_].+)?$")
    tibble(run_dir = all_dirs[str_detect(basename(all_dirs), patt)])
  })
  if (nrow(hits) == 0) stop("No run directory found for method: ", method)

  valid <- hits %>%
    mutate(
      data_dir = file.path(run_dir, "results", "data"),
      summary_dir = file.path(run_dir, "results", "summary"),
      has_pheno = map_lgl(data_dir, ~ length(list.files(.x, pattern = "^pheno_.*\\.phenotypes\\.tsv$", full.names = TRUE)) > 0),
      has_map = map_lgl(data_dir, ~ length(list.files(.x, pattern = "^pheno_.*\\.gene_index_map\\.txt$", full.names = TRUE)) > 0),
      has_summary = map_lgl(summary_dir, ~ length(list.files(.x, pattern = "^final_heritability_summary.*\\.tsv$", full.names = TRUE)) > 0)
    ) %>%
    filter(has_pheno, has_map, if (require_summary) has_summary else TRUE)

  if (nrow(valid) == 0) stop("Run directory exists but required files are missing for method: ", method)

  exact <- valid %>% filter(basename(run_dir) == paste0(method, run_suffix))
  if (nrow(exact) > 0) return(exact$run_dir[[1]])

  valid %>%
    mutate(mtime = file.info(run_dir)$mtime) %>%
    arrange(desc(mtime)) %>%
    pull(run_dir) %>%
    .[[1]]
}

read_expression_by_gene <- function(run_dir) {
  data_dir <- file.path(run_dir, "results", "data")
  pheno_file <- list.files(data_dir, pattern = "^pheno_.*\\.phenotypes\\.tsv$", full.names = TRUE) %>% .[[1]]
  map_file <- list.files(data_dir, pattern = "^pheno_.*\\.gene_index_map\\.txt$", full.names = TRUE) %>% .[[1]]

  pheno <- fread(pheno_file, header = FALSE, sep = "\t", data.table = TRUE)
  map_dt <- fread(map_file, sep = "\t", data.table = TRUE)
  if (!all(c("gene_name", "mpheno_index") %in% names(map_dt))) stop("Unexpected map format: ", map_file)

  map_dt <- map_dt[order(as.integer(mpheno_index))]
  genes <- normalize_gene_id(map_dt$gene_name)
  expr <- as.matrix(pheno[, -(1:2), with = FALSE])
  storage.mode(expr) <- "numeric"
  if ((ncol(pheno) - 2L) != length(genes)) stop("Phenotype/map dimension mismatch in: ", run_dir)

  # Collapse duplicate gene IDs by mean if needed.
  if (anyDuplicated(genes) > 0L) {
    idx <- split(seq_along(genes), genes)
    expr_collapsed <- do.call(
      cbind,
      lapply(idx, function(ii) {
        if (length(ii) == 1L) expr[, ii] else rowMeans(expr[, ii, drop = FALSE], na.rm = TRUE)
      })
    )
    gene_names <- names(idx)
    mean_expr <- colMeans(expr_collapsed, na.rm = TRUE)
  } else {
    mean_expr <- colMeans(expr, na.rm = TRUE)
    names(mean_expr) <- genes
    gene_names <- genes
  }

  tibble(
    Gene = normalize_gene_id(gene_names),
    mean_tmm_expression = as.numeric(mean_expr)
  ) %>%
    group_by(Gene) %>%
    summarise(mean_tmm_expression = mean(mean_tmm_expression, na.rm = TRUE), .groups = "drop")
}

read_expression_metrics_by_gene <- function(run_dir, value_label) {
  data_dir <- file.path(run_dir, "results", "data")
  pheno_file <- list.files(data_dir, pattern = "^pheno_.*\\.phenotypes\\.tsv$", full.names = TRUE) %>% .[[1]]
  map_file <- list.files(data_dir, pattern = "^pheno_.*\\.gene_index_map\\.txt$", full.names = TRUE) %>% .[[1]]

  pheno <- fread(pheno_file, header = FALSE, sep = "\t", data.table = TRUE)
  map_dt <- fread(map_file, sep = "\t", data.table = TRUE)
  if (!all(c("gene_name", "mpheno_index") %in% names(map_dt))) stop("Unexpected map format: ", map_file)

  map_dt <- map_dt[order(as.integer(mpheno_index))]
  genes <- normalize_gene_id(map_dt$gene_name)
  expr <- as.matrix(pheno[, -(1:2), with = FALSE])
  storage.mode(expr) <- "numeric"
  if ((ncol(pheno) - 2L) != length(genes)) stop("Phenotype/map dimension mismatch in: ", run_dir)

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

  skewness_third_moment <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 3L) return(NA_real_)
    s <- sd(x)
    if (!is.finite(s) || s == 0) return(NA_real_)
    mean(((x - mean(x)) / s)^3)
  }

  tibble(
    Gene = normalize_gene_id(colnames(expr)),
    value_type = value_label,
    mean_expr = colMeans(expr, na.rm = TRUE),
    median_expr = apply(expr, 2, median, na.rm = TRUE),
    log2_mean_expr_plus1 = log2(mean_expr + 1),
    log2_median_expr_plus1 = log2(median_expr + 1),
    skewness = apply(expr, 2, skewness_third_moment),
    zero_fraction = colMeans(expr == 0, na.rm = TRUE)
  ) %>%
    group_by(Gene, value_type) %>%
    summarise(
      mean_expr = mean(mean_expr, na.rm = TRUE),
      median_expr = mean(median_expr, na.rm = TRUE),
      log2_mean_expr_plus1 = mean(log2_mean_expr_plus1, na.rm = TRUE),
      log2_median_expr_plus1 = mean(log2_median_expr_plus1, na.rm = TRUE),
      skewness = mean(skewness, na.rm = TRUE),
      zero_fraction = mean(zero_fraction, na.rm = TRUE),
      .groups = "drop"
    )
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

read_h2_summary <- function(run_dir) {
  sum_file <- list.files(file.path(run_dir, "results", "summary"), pattern = "^final_heritability_summary.*\\.tsv$", full.names = TRUE)
  if (length(sum_file) == 0) stop("Missing summary in run dir: ", run_dir)

  dt <- fread(sum_file[1], sep = "\t", data.table = TRUE)
  tibble(
    Gene = resolve_gene_ids_for_summary(dt$Gene, run_dir),
    Status = toupper(str_trim(as.character(dt$Status))),
    h2_GREML = suppressWarnings(as.numeric(dt$h2_GREML))
  ) %>%
    filter(Status == "PASS", is.finite(h2_GREML))
}

make_method_label <- function(method_id) {
  m <- str_match(method_id, "^(all_snps|hm3_no_mhc)_tmm_(raw|irnt|inverse_normal)$")
  if (is.na(m[1, 1])) return(method_id)
  snp <- ifelse(m[1, 2] == "all_snps", "ALL", "HM3")
  norm <- ifelse(m[1, 3] == "inverse_normal", "IRNT", toupper(m[1, 3]))
  paste(snp, "TMM", norm, sep = " | ")
}

method_to_tmm_raw_baseline <- function(method_id) {
  m <- str_match(method_id, "^(all_snps|hm3_no_mhc)_(tpm|tmm)_(raw|irnt|inverse_normal)$")
  if (is.na(m[1, 1])) return(method_id)
  paste(m[1, 2], "tmm", "raw", sep = "_")
}

method_to_tpm_raw_baseline <- function(method_id) {
  m <- str_match(method_id, "^(all_snps|hm3_no_mhc)_(tpm|tmm)_(raw|irnt|inverse_normal)$")
  if (is.na(m[1, 1])) return(method_id)
  paste(m[1, 2], "tpm", "raw", sep = "_")
}

plot_h2_by_decile <- function(plot_tbl, summary_tbl, x_var, title_text, subtitle_text, x_text, out_file) {
  y_min <- min(plot_tbl$h2_GREML, na.rm = TRUE)
  y_max <- max(plot_tbl$h2_GREML, na.rm = TRUE)
  y_pad <- max((y_max - y_min) * 0.05, 0.02)

  p <- ggplot(plot_tbl, aes(x = .data[[x_var]], y = h2_GREML, fill = .data[[x_var]])) +
    geom_boxplot(width = 0.68, outlier.alpha = 0.15, size = 0.28) +
    geom_point(
      data = summary_tbl,
      aes(x = .data[[x_var]], y = mean_h2),
      inherit.aes = FALSE,
      shape = 23,
      fill = "white",
      color = "black",
      size = 1.9,
      stroke = 0.35
    ) +
    geom_point(
      data = summary_tbl,
      aes(x = .data[[x_var]], y = median_h2),
      inherit.aes = FALSE,
      shape = 21,
      fill = "#F4A261",
      color = "black",
      size = 1.8,
      stroke = 0.3
    ) +
    geom_text(
      data = summary_tbl,
      aes(
        x = .data[[x_var]],
        y = pmax(mean_h2, median_h2) + y_pad,
        label = sprintf("mean=%.3f\nmedian=%.3f", mean_h2, median_h2)
      ),
      inherit.aes = FALSE,
      size = 2.5,
      vjust = 0
    ) +
    scale_fill_manual(
      values = setNames(
        colorRampPalette(RColorBrewer::brewer.pal(9, "GnBu"))(10),
        paste0("D", 1:10)
      )
    ) +
    facet_wrap(~ method_label, ncol = 2, scales = "fixed") +
    coord_cartesian(ylim = c(y_min, y_max + y_pad * 2.4), clip = "off") +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = x_text,
      y = "GREML h2 estimate"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.margin = margin(10, 16, 6, 6)
    )

  ggsave(out_file, p, width = 13, height = 7.5, dpi = 300)
}

# --------------------------------------------------
# Build analysis table
# --------------------------------------------------
analysis_tbl <- map_dfr(method_ids, function(mid) {
  run_dir <- find_run_dir(mid)
  expr_method <- method_to_tmm_raw_baseline(mid)
  expr_run_dir <- find_run_dir(expr_method, require_summary = FALSE)
  expr_tbl <- read_expression_by_gene(expr_run_dir)
  h2_tbl <- read_h2_summary(run_dir)

  merged <- h2_tbl %>%
    inner_join(expr_tbl, by = "Gene") %>%
    mutate(log2_mean_tmm_plus1 = log2(mean_tmm_expression + 1)) %>%
    filter(
      is.finite(mean_tmm_expression),
      is.finite(log2_mean_tmm_plus1),
      log2_mean_tmm_plus1 > min_log2_mean_tmm_plus1,
      is.finite(h2_GREML)
    )

  if (nrow(merged) < 50) {
    stop("Too few genes after merge for method ", mid, " (n=", nrow(merged), ").")
  }

  merged %>%
    mutate(
      method_id = mid,
      expression_method_id = expr_method,
      expression_filter_min_log2_mean_tmm_plus1 = min_log2_mean_tmm_plus1,
      run_dir = run_dir,
      expression_run_dir = expr_run_dir,
      method_label = make_method_label(mid)
    )
})

median_tpm_analysis_tbl <- map_dfr(method_ids, function(mid) {
  run_dir <- find_run_dir(mid)
  tpm_method <- method_to_tpm_raw_baseline(mid)
  tpm_run_dir <- find_run_dir(tpm_method, require_summary = FALSE)
  tpm_metrics <- read_expression_metrics_by_gene(tpm_run_dir, "TPM RAW") %>%
    select(Gene, median_tpm_expression = median_expr, log2_median_tpm_plus1 = log2_median_expr_plus1)
  h2_tbl <- read_h2_summary(run_dir)

  merged <- h2_tbl %>%
    inner_join(tpm_metrics, by = "Gene") %>%
    filter(
      is.finite(median_tpm_expression),
      is.finite(log2_median_tpm_plus1),
      log2_median_tpm_plus1 > min_log2_median_tpm_plus1,
      is.finite(h2_GREML)
    )

  if (nrow(merged) < 50) {
    stop("Too few genes after TPM median merge for method ", mid, " (n=", nrow(merged), ").")
  }

  merged %>%
    mutate(
      method_id = mid,
      expression_method_id = tpm_method,
      expression_filter_min_log2_median_tpm_plus1 = min_log2_median_tpm_plus1,
      run_dir = run_dir,
      expression_run_dir = tpm_run_dir,
      method_label = make_method_label(mid)
    )
})

analysis_tbl <- analysis_tbl %>%
  group_by(method_id) %>%
  mutate(
    expr_decile = ntile(mean_tmm_expression, 10),
    expr_decile = factor(expr_decile, levels = 1:10, labels = paste0("D", 1:10))
  ) %>%
  ungroup()

median_tpm_analysis_tbl <- median_tpm_analysis_tbl %>%
  group_by(method_id) %>%
  mutate(
    median_tpm_decile = ntile(median_tpm_expression, 10),
    median_tpm_decile = factor(median_tpm_decile, levels = 1:10, labels = paste0("D", 1:10))
  ) %>%
  ungroup()

summary_tbl <- analysis_tbl %>%
  group_by(method_id, method_label, expr_decile) %>%
  summarise(
    n_genes = n(),
    mean_h2 = mean(h2_GREML, na.rm = TRUE),
    median_h2 = median(h2_GREML, na.rm = TRUE),
    q25_h2 = quantile(h2_GREML, 0.25, na.rm = TRUE),
    q75_h2 = quantile(h2_GREML, 0.75, na.rm = TRUE),
    mean_tmm_expression = mean(mean_tmm_expression, na.rm = TRUE),
    mean_log2_tmm_plus1 = mean(log2_mean_tmm_plus1, na.rm = TRUE),
    expression_filter_min_log2_mean_tmm_plus1 = first(expression_filter_min_log2_mean_tmm_plus1),
    .groups = "drop"
  )

median_tpm_summary_tbl <- median_tpm_analysis_tbl %>%
  group_by(method_id, method_label, median_tpm_decile) %>%
  summarise(
    n_genes = n(),
    mean_h2 = mean(h2_GREML, na.rm = TRUE),
    median_h2 = median(h2_GREML, na.rm = TRUE),
    q25_h2 = quantile(h2_GREML, 0.25, na.rm = TRUE),
    q75_h2 = quantile(h2_GREML, 0.75, na.rm = TRUE),
    median_tpm_expression = median(median_tpm_expression, na.rm = TRUE),
    mean_log2_median_tpm_plus1 = mean(log2_median_tpm_plus1, na.rm = TRUE),
    expression_filter_min_log2_median_tpm_plus1 = first(expression_filter_min_log2_median_tpm_plus1),
    .groups = "drop"
  )

skewness_tbl <- map_dfr(method_ids, function(mid) {
  tmm_method <- method_to_tmm_raw_baseline(mid)
  tpm_method <- method_to_tpm_raw_baseline(mid)

  tmm_run_dir <- find_run_dir(tmm_method, require_summary = FALSE)
  tpm_run_dir <- find_run_dir(tpm_method, require_summary = FALSE)
  h2_run_dir <- find_run_dir(mid)

  h2_genes <- read_h2_summary(h2_run_dir) %>%
    select(Gene, h2_GREML)

  bind_rows(
    read_expression_metrics_by_gene(tmm_run_dir, "TMM RAW"),
    read_expression_metrics_by_gene(tpm_run_dir, "TPM RAW")
  ) %>%
    inner_join(h2_genes, by = "Gene") %>%
    filter(is.finite(skewness), is.finite(h2_GREML)) %>%
    group_by(value_type) %>%
    mutate(
      expr_decile = ntile(median_expr, 10),
      expr_decile = factor(expr_decile, levels = 1:10, labels = paste0("D", 1:10))
    ) %>%
    ungroup() %>%
    mutate(
      method_id = mid,
      method_label = make_method_label(mid),
      h2_run_dir = h2_run_dir,
      expression_run_dir = if_else(value_type == "TMM RAW", tmm_run_dir, tpm_run_dir)
    )
})

skewness_summary_tbl <- skewness_tbl %>%
  group_by(method_id, method_label, value_type, expr_decile) %>%
  summarise(
    n_genes = n(),
    mean_skewness = mean(skewness, na.rm = TRUE),
    median_skewness = median(skewness, na.rm = TRUE),
    q25_skewness = quantile(skewness, 0.25, na.rm = TRUE),
    q75_skewness = quantile(skewness, 0.75, na.rm = TRUE),
    mean_h2 = mean(h2_GREML, na.rm = TRUE),
    median_h2 = median(h2_GREML, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(analysis_tbl, out_gene_tsv)
write_tsv(summary_tbl, out_summary_tsv)
write_tsv(median_tpm_analysis_tbl, out_median_tpm_gene_tsv)
write_tsv(median_tpm_summary_tbl, out_median_tpm_summary_tsv)
write_tsv(skewness_tbl, out_skewness_by_bin_tsv)
write_tsv(skewness_summary_tbl, out_skewness_by_bin_summary_tsv)

# --------------------------------------------------
# Plots
# --------------------------------------------------
plot_h2_by_decile(
  analysis_tbl,
  summary_tbl,
  "expr_decile",
  "GREML h2 by mean TMM expression decile",
  sprintf(
    "Genes are binned by matching TMM-RAW mean expression; PASS genes only; filtered to log2(mean TMM + 1) > %.2f\nWhite diamond = mean, orange circle = median",
    min_log2_mean_tmm_plus1
  ),
  "TMM-RAW mean expression decile after low-expression filter",
  out_plot
)

plot_h2_by_decile(
  median_tpm_analysis_tbl,
  median_tpm_summary_tbl,
  "median_tpm_decile",
  "TMM GREML h2 by median TPM expression decile",
  sprintf(
    "Genes are binned by median TPM-RAW expression; h2 estimates come from the selected TMM method summaries; filtered to log2(median TPM + 1) > %.2f\nWhite diamond = mean, orange circle = median",
    min_log2_median_tpm_plus1
  ),
  "TPM-RAW median expression decile",
  out_median_tpm_plot
)

p_skewness <- ggplot(skewness_tbl, aes(x = expr_decile, y = skewness, fill = expr_decile)) +
  geom_boxplot(width = 0.68, outlier.alpha = 0.15, size = 0.28) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#E76F51", linewidth = 0.35) +
  scale_fill_manual(
    values = setNames(
      colorRampPalette(RColorBrewer::brewer.pal(9, "GnBu"))(10),
      paste0("D", 1:10)
    )
  ) +
  facet_grid(value_type ~ method_label, scales = "free_y") +
  labs(
    title = "Expression skewness distribution by expression decile",
    subtitle = "Deciles are based on median expression separately within TMM RAW and TPM RAW; each point/gene contributes one skewness value",
    x = "Median expression decile",
    y = "Skewness (third standardized moment)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(out_skewness_by_bin_plot, p_skewness, width = 14, height = 8, dpi = 300)

cat("Saved:\n")
cat("- ", out_gene_tsv, "\n", sep = "")
cat("- ", out_summary_tsv, "\n", sep = "")
cat("- ", out_plot, "\n", sep = "")
cat("- ", out_median_tpm_gene_tsv, "\n", sep = "")
cat("- ", out_median_tpm_summary_tsv, "\n", sep = "")
cat("- ", out_median_tpm_plot, "\n", sep = "")
cat("- ", out_skewness_by_bin_tsv, "\n", sep = "")
cat("- ", out_skewness_by_bin_summary_tsv, "\n", sep = "")
cat("- ", out_skewness_by_bin_plot, "\n", sep = "")
