suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# -----------------------------
# Settings
# -----------------------------
runs_dir <- "/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260428_GE_GEUVADIS_v2/GeneExpression/runs"

# 1-vs-all mode: set method id here, e.g. "all_snps_tmm_raw_peerauto_pmg0_npc5"
# all-pairs mode: keep NULL
focus_method <- NULL
include_mixed <- FALSE

# PI-driven additions
skew_reference_method <- "all_snps_tmm_raw"  # TMM (RAW)
n_random_genes <- 200L
n_skew_genes_per_tail <- 4L

# -----------------------------
# Helpers
# -----------------------------
normalize_norm <- function(x) if_else(x == "inverse_normal", "irnt", x)

normalize_focus <- function(x) {
  x %>%
    str_replace("_inverse_normal_", "_irnt_") %>%
    str_remove("_peerauto_pmg0_npc5$")
}

method_to_label <- function(method_id) {
  m <- str_match(method_id, "^(all_snps|hm3_no_mhc)_(tpm|tmm)_(irnt|inverse_normal|raw)$")
  snp_lbl <- if_else(m[, 2] == "all_snps", "ALL", "HM3")
  expr_lbl <- toupper(m[, 3])
  norm_lbl <- toupper(normalize_norm(m[, 4]))
  paste(snp_lbl, expr_lbl, norm_lbl, sep = " | ")
}

parse_method <- function(run_name) {
  m <- str_match(run_name, "^(all_snps|hm3_no_mhc)_(tpm|tmm)_(irnt|inverse_normal|raw)_peerauto_pmg0_npc5$")
  tibble(
    run_name = run_name,
    snp_set = m[, 2],
    expr = m[, 3],
    norm = normalize_norm(m[, 4])
  ) %>%
    mutate(method = paste(snp_set, expr, norm, sep = "_"))
}

method_to_run_name <- function(method_id) {
  paste0(method_id, "_peerauto_pmg0_npc5")
}

method_to_raw_baseline <- function(method_id) {
  str_replace(method_id, "_(irnt|raw)$", "_raw")
}

# read one run's phenotype matrix and return per-gene log expression
read_log_expr_for_run <- function(run_name) {
  data_dir <- file.path(runs_dir, run_name, "results", "data")
  pheno_file <- list.files(data_dir, pattern = "^pheno_.*\\.phenotypes\\.tsv$", full.names = TRUE)
  map_file <- list.files(data_dir, pattern = "^pheno_.*\\.gene_index_map\\.txt$", full.names = TRUE)

  if (length(pheno_file) == 0 || length(map_file) == 0) {
    warning("Missing phenotype/map file for run: ", run_name)
    return(tibble(Gene = character(), log_expr = numeric(), mean_expr = numeric()))
  }

  pheno <- fread(pheno_file[1], header = FALSE, sep = "\t", data.table = TRUE)
  map_dt <- fread(map_file[1], sep = "\t", data.table = TRUE)

  if (!all(c("gene_name", "mpheno_index") %in% names(map_dt))) {
    warning("Unexpected gene index map columns in run: ", run_name)
    return(tibble(Gene = character(), log_expr = numeric(), mean_expr = numeric()))
  }

  map_dt <- map_dt[order(as.integer(mpheno_index))]
  gene_names <- as.character(map_dt$gene_name)

  if (ncol(pheno) < 3 || (ncol(pheno) - 2L) != length(gene_names)) {
    warning("Phenotype/map dimension mismatch in run: ", run_name)
    return(tibble(Gene = character(), log_expr = numeric(), mean_expr = numeric()))
  }

  setnames(pheno, c("FID", "IID", gene_names))
  expr_df <- pheno[, ..gene_names]
  expr_mat <- as.matrix(expr_df)
  storage.mode(expr_mat) <- "numeric"

  mean_expr <- colMeans(expr_mat, na.rm = TRUE)
  tibble(
    Gene = gene_names,
    mean_expr = mean_expr,
    log_expr = log2(pmax(mean_expr, 0) + 1)
  )
}

calc_skewness <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 3L) return(NA_real_)
  m <- mean(x)
  m2 <- mean((x - m)^2)
  if (!is.finite(m2) || m2 <= 0) return(0)
  m3 <- mean((x - m)^3)
  g1 <- m3 / (m2^(3 / 2))
  # Bias-corrected Fisher-Pearson skewness (still third-moment based)
  sqrt(n * (n - 1)) / (n - 2) * g1
}

build_qq_table <- function(mat, genes) {
  map_dfr(genes, function(g) {
    vals <- mat[, g]
    vals <- vals[is.finite(vals)]
    n <- length(vals)
    if (n < 3L) return(tibble())
    mu <- mean(vals)
    sig <- sd(vals)
    if (!is.finite(sig) || sig == 0) return(tibble())

    tibble(
      Gene = g,
      theo_q = qnorm(ppoints(n)),
      obs_q = sort((vals - mu) / sig)
    )
  })
}

# -----------------------------
# 1) Read all summary files
# -----------------------------
files <- list.files(
  runs_dir,
  pattern = "final_heritability_summary.*\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(files) == 0) stop("No summary files found under: ", runs_dir)

h2_long <- map_dfr(files, function(f) {
  run_name <- str_match(f, "runs/([^/]+)/results/summary/")[, 2]
  meta <- parse_method(run_name)

  if (any(is.na(meta$snp_set))) return(tibble())

  read_tsv(f, show_col_types = FALSE) %>%
    transmute(
      Gene,
      Status,
      h2 = suppressWarnings(as.numeric(h2_GREML)),
      se = suppressWarnings(as.numeric(SE_GREML)),
      pval = suppressWarnings(as.numeric(Pval_GREML)),
      method = meta$method
    )
}) %>%
  filter(Status == "PASS", is.finite(h2)) %>%
  mutate(
    z = if_else(is.finite(se) & se > 0, h2 / se, NA_real_)
  )

if (nrow(h2_long) == 0) stop("No PASS h2 values found.")

# Collapse duplicates created by inverse_normal -> irnt harmonization
h2_long <- h2_long %>%
  group_by(Gene, method) %>%
  summarise(
    h2 = mean(h2, na.rm = TRUE),
    se = mean(se, na.rm = TRUE),
    pval = mean(pval, na.rm = TRUE),
    z = mean(z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(method) %>%
  mutate(qval = p.adjust(pval, method = "BH")) %>%
  ungroup()

# -----------------------------
# 2) Gene x method matrix
# -----------------------------
h2_wide <- h2_long %>%
  pivot_wider(
    names_from = method,
    values_from = c(h2, z, qval),
    names_sep = "__"
  )

methods <- h2_long %>% distinct(method) %>% pull(method) %>% sort()
if (length(methods) < 2) stop("Need at least 2 methods for pairwise comparison.")

if (!is.null(focus_method)) {
  focus_method <- normalize_focus(focus_method)
  if (!focus_method %in% methods) stop("focus_method not found after normalization: ", focus_method)
}

# -----------------------------
# 3) Log-expression map from RAW baselines
# -----------------------------
raw_methods <- methods %>% map_chr(method_to_raw_baseline) %>% unique()
log_expr_map <- map_dfr(raw_methods, function(m_raw) {
  run_name <- method_to_run_name(m_raw)
  read_log_expr_for_run(run_name) %>%
    mutate(raw_method = m_raw)
})

# -----------------------------
# 4) Build pairwise table
# -----------------------------
pair_df <- combn(methods, 2, simplify = FALSE) %>%
  map_dfr(function(p) {
    h2_wide %>%
      transmute(
        Gene,
        m1 = p[1],
        m2 = p[2],
        x = .data[[paste0("h2__", p[1])]],
        y = .data[[paste0("h2__", p[2])]],
        x_z = .data[[paste0("z__", p[1])]],
        y_z = .data[[paste0("z__", p[2])]],
        x_q = .data[[paste0("qval__", p[1])]],
        y_q = .data[[paste0("qval__", p[2])]]
      ) %>%
      filter(!is.na(x), !is.na(y))
  }) %>%
  mutate(
    m1_snp = str_match(m1, "^(all_snps|hm3_no_mhc)_(tpm|tmm)_(irnt|inverse_normal|raw)$")[, 2],
    m1_expr = str_match(m1, "^(all_snps|hm3_no_mhc)_(tpm|tmm)_(irnt|inverse_normal|raw)$")[, 3],
    m1_norm = normalize_norm(str_match(m1, "^(all_snps|hm3_no_mhc)_(tpm|tmm)_(irnt|inverse_normal|raw)$")[, 4]),
    m2_snp = str_match(m2, "^(all_snps|hm3_no_mhc)_(tpm|tmm)_(irnt|inverse_normal|raw)$")[, 2],
    m2_expr = str_match(m2, "^(all_snps|hm3_no_mhc)_(tpm|tmm)_(irnt|inverse_normal|raw)$")[, 3],
    m2_norm = normalize_norm(str_match(m2, "^(all_snps|hm3_no_mhc)_(tpm|tmm)_(irnt|inverse_normal|raw)$")[, 4]),
    section = case_when(
      m1_snp != m2_snp & m1_expr == m2_expr & m1_norm == m2_norm ~ "SNP set: ALL vs HM3",
      m1_snp == m2_snp & m1_expr != m2_expr & m1_norm == m2_norm ~ "Expression: TPM vs TMM",
      m1_snp == m2_snp & m1_expr == m2_expr & m1_norm != m2_norm ~ "Normalization: RAW vs IRNT",
      TRUE ~ "Mixed changes"
    ),
    m1_label = method_to_label(m1),
    m2_label = method_to_label(m2)
  )

if (!is.null(focus_method)) {
  pair_df <- pair_df %>% filter(m1 == focus_method | m2 == focus_method)
}

cat("Methods detected:\n")
print(methods)
cat("Section counts before mixed-filter:\n")
print(table(pair_df$section, useNA = "ifany"))

if (!include_mixed) {
  pair_df <- pair_df %>% filter(section != "Mixed changes")
}

if (nrow(pair_df) == 0) stop("No rows left after filtering.")

# -----------------------------
# 5) Orient x/y for grouped sections
# -----------------------------
snp_df <- pair_df %>%
  filter(section == "SNP set: ALL vs HM3") %>%
  mutate(
    x_plot = if_else(m1_snp == "all_snps", x, y),
    y_plot = if_else(m1_snp == "all_snps", y, x),
    x_z_plot = if_else(m1_snp == "all_snps", x_z, y_z),
    y_z_plot = if_else(m1_snp == "all_snps", y_z, x_z),
    x_q_plot = if_else(m1_snp == "all_snps", x_q, y_q),
    y_q_plot = if_else(m1_snp == "all_snps", y_q, x_q),
    x_method = if_else(m1_snp == "all_snps", m1, m2),
    y_method = if_else(m1_snp == "all_snps", m2, m1),
    expr_fix = toupper(if_else(m1_snp == "all_snps", m1_expr, m2_expr)),
    norm_fix = toupper(if_else(m1_snp == "all_snps", m1_norm, m2_norm)),
    facet_label = paste(expr_fix, norm_fix, sep = " | ")
  )

expr_df <- pair_df %>%
  filter(section == "Expression: TPM vs TMM") %>%
  mutate(
    x_plot = if_else(m1_expr == "tmm", x, y),
    y_plot = if_else(m1_expr == "tmm", y, x),
    x_z_plot = if_else(m1_expr == "tmm", x_z, y_z),
    y_z_plot = if_else(m1_expr == "tmm", y_z, x_z),
    x_q_plot = if_else(m1_expr == "tmm", x_q, y_q),
    y_q_plot = if_else(m1_expr == "tmm", y_q, x_q),
    x_method = if_else(m1_expr == "tmm", m1, m2),
    y_method = if_else(m1_expr == "tmm", m2, m1),
    snp_fix = if_else(if_else(m1_expr == "tmm", m1_snp, m2_snp) == "all_snps", "ALL", "HM3"),
    norm_fix = toupper(if_else(m1_expr == "tmm", m1_norm, m2_norm)),
    facet_label = paste(snp_fix, norm_fix, sep = " | ")
  )

norm_df <- pair_df %>%
  filter(section == "Normalization: RAW vs IRNT") %>%
  mutate(
    x_plot = if_else(m1_norm == "raw", x, y),
    y_plot = if_else(m1_norm == "raw", y, x),
    x_z_plot = if_else(m1_norm == "raw", x_z, y_z),
    y_z_plot = if_else(m1_norm == "raw", y_z, x_z),
    x_q_plot = if_else(m1_norm == "raw", x_q, y_q),
    y_q_plot = if_else(m1_norm == "raw", y_q, x_q),
    x_method = if_else(m1_norm == "raw", m1, m2),
    y_method = if_else(m1_norm == "raw", m2, m1),
    snp_fix = if_else(if_else(m1_norm == "raw", m1_snp, m2_snp) == "all_snps", "ALL", "HM3"),
    expr_fix = toupper(if_else(m1_norm == "raw", m1_expr, m2_expr)),
    facet_label = paste(snp_fix, expr_fix, sep = " | ")
  )

add_plot_covariates <- function(df) {
  if (nrow(df) == 0) return(df)

  expr_lookup <- log_expr_map %>%
    rename(raw_method_x = raw_method, log_expr_x = log_expr, mean_expr_x = mean_expr)

  out <- df %>%
    mutate(
      raw_method_x = method_to_raw_baseline(x_method),
      raw_method_y = method_to_raw_baseline(y_method),
      z_plot = (x_z_plot + y_z_plot) / 2
    ) %>%
    left_join(expr_lookup %>% select(Gene, raw_method_x, log_expr_x, mean_expr_x),
      by = c("Gene", "raw_method_x")
    )

  expr_lookup_y <- log_expr_map %>%
    rename(raw_method_y = raw_method, log_expr_y = log_expr, mean_expr_y = mean_expr)

  out <- out %>%
    left_join(expr_lookup_y %>% select(Gene, raw_method_y, log_expr_y, mean_expr_y),
      by = c("Gene", "raw_method_y")
    ) %>%
    mutate(
      log_expr_plot = rowMeans(cbind(log_expr_x, log_expr_y), na.rm = TRUE),
      mean_expr_plot = rowMeans(cbind(mean_expr_x, mean_expr_y), na.rm = TRUE)
    )

  out
}

snp_df <- add_plot_covariates(snp_df)
expr_df <- add_plot_covariates(expr_df)
norm_df <- add_plot_covariates(norm_df)

# -----------------------------
# 6) Plot helpers + stats
# -----------------------------
calc_stats <- function(df) {
  df %>%
    group_by(facet_label) %>%
    summarise(
      n = n(),
      r = cor(x_plot, y_plot, use = "complete.obs"),
      n_sig_x = sum(!is.na(x_q_plot) & x_q_plot < 0.05),
      n_sig_y = sum(!is.na(y_q_plot) & y_q_plot < 0.05),
      n_abs_x_hi_y_lo = sum(x_plot > 0.05 & y_plot < 0.005, na.rm = TRUE),
      n_abs_y_hi_x_lo = sum(y_plot > 0.05 & x_plot < 0.005, na.rm = TRUE),
      n_fold_x_gt10y = sum(x_plot > 0 & y_plot > 0 & x_plot / y_plot >= 10, na.rm = TRUE),
      n_fold_y_gt10x = sum(x_plot > 0 & y_plot > 0 & y_plot / x_plot >= 10, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      lbl_top = sprintf(
        "n=%d, r=%.2f\nq<0.05: X=%d, Y=%d\nabs: X>0.05,Y<0.005=%d\nabs: Y>0.05,X<0.005=%d\nfold(>=10x): X>>Y=%d\nfold(>=10x): Y>>X=%d",
        n, r,
        n_sig_x, n_sig_y,
        n_abs_x_hi_y_lo, n_abs_y_hi_x_lo,
        n_fold_x_gt10y, n_fold_y_gt10x
      )
    )
}

plot_section <- function(df, title, x_lab, y_lab, out_file, color_var, color_label, color_mode = c("seq", "div")) {
  if (nrow(df) == 0) return(invisible(NULL))
  color_mode <- match.arg(color_mode)
  stats <- calc_stats(df)

  p <- df %>%
    ggplot(aes(x = x_plot, y = y_plot, color = .data[[color_var]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E76F51", linewidth = 0.35) +
    geom_point(alpha = 0.4, size = 0.7) +
    geom_text(
      data = stats,
      aes(x = -Inf, y = Inf, label = lbl_top),
      inherit.aes = FALSE,
      hjust = -0.1,
      vjust = 1.1,
      size = 2.8,
      color = "#264653"
    ) +
    facet_wrap(~facet_label, scales = "free", ncol = 4) +
    labs(
      title = title,
      subtitle = "PASS genes only",
      x = x_lab,
      y = y_lab,
      color = color_label
    ) +
    theme_minimal(base_size = 10.5) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 8.5)
    )

  if (color_mode == "seq") {
    p <- p + scale_color_viridis_c(option = "magma", na.value = "grey80")
  } else {
    p <- p + scale_color_gradient2(
      low = "#3B4CC0",
      mid = "#F7F7F7",
      high = "#B40426",
      midpoint = 0,
      na.value = "grey80"
    )
  }

  ggsave(out_file, p, width = 14, height = 8, dpi = 300)
  invisible(stats)
}

prefix <- ifelse(is.null(focus_method), "all_pairs", paste0("focus_", focus_method))
out_stats <- file.path(runs_dir, paste0("pairwise_h2_stats_", prefix, ".tsv"))
out_discrepancy_types <- file.path(runs_dir, paste0("pairwise_h2_discrepancy_types_", prefix, ".tsv"))

stats_all <- bind_rows(
  calc_stats(snp_df) %>% mutate(section = "SNP set: ALL vs HM3"),
  calc_stats(expr_df) %>% mutate(section = "Expression: TPM vs TMM"),
  calc_stats(norm_df) %>% mutate(section = "Normalization: RAW vs IRNT")
)
write_tsv(stats_all, out_stats)

discrepancy_type_counts <- stats_all %>%
  select(section, facet_label, n_abs_x_hi_y_lo, n_abs_y_hi_x_lo, n_fold_x_gt10y, n_fold_y_gt10x) %>%
  pivot_longer(
    cols = c(n_abs_x_hi_y_lo, n_abs_y_hi_x_lo, n_fold_x_gt10y, n_fold_y_gt10x),
    names_to = "type",
    values_to = "n_genes"
  ) %>%
  mutate(
    type = recode(
      type,
      n_abs_x_hi_y_lo = "ABS_X_gt_0.05_and_Y_lt_0.005",
      n_abs_y_hi_x_lo = "ABS_Y_gt_0.05_and_X_lt_0.005",
      n_fold_x_gt10y = "FOLD_X_over_Y_ge_10",
      n_fold_y_gt10x = "FOLD_Y_over_X_ge_10"
    )
  )
write_tsv(discrepancy_type_counts, out_discrepancy_types)

# 6A) Color by log expression
plot_section(
  snp_df,
  title = "Pairwise GREML h2: SNP effect (ALL vs HM3)",
  x_lab = "X-axis h2 (ALL)",
  y_lab = "Y-axis h2 (HM3)",
  out_file = file.path(runs_dir, paste0("pairwise_h2_scatter_", prefix, "_snp_set_all_vs_hm3_color_logexpr.png")),
  color_var = "log_expr_plot",
  color_label = "log2(mean expr + 1)",
  color_mode = "seq"
)

plot_section(
  expr_df,
  title = "Pairwise GREML h2: Expression effect (TPM vs TMM)",
  x_lab = "X-axis h2 (TMM)",
  y_lab = "Y-axis h2 (TPM)",
  out_file = file.path(runs_dir, paste0("pairwise_h2_scatter_", prefix, "_expression_tpm_vs_tmm_color_logexpr.png")),
  color_var = "log_expr_plot",
  color_label = "log2(mean expr + 1)",
  color_mode = "seq"
)

plot_section(
  norm_df,
  title = "Pairwise GREML h2: Normalization effect (RAW vs IRNT)",
  x_lab = "X-axis h2 (RAW)",
  y_lab = "Y-axis h2 (IRNT)",
  out_file = file.path(runs_dir, paste0("pairwise_h2_scatter_", prefix, "_normalization_raw_vs_irnt_color_logexpr.png")),
  color_var = "log_expr_plot",
  color_label = "log2(mean expr + 1)",
  color_mode = "seq"
)

# 6B) Color by heritability significance proxy (Wald Z)
plot_section(
  snp_df,
  title = "Pairwise GREML h2: SNP effect (ALL vs HM3)",
  x_lab = "X-axis h2 (ALL)",
  y_lab = "Y-axis h2 (HM3)",
  out_file = file.path(runs_dir, paste0("pairwise_h2_scatter_", prefix, "_snp_set_all_vs_hm3_color_zscore.png")),
  color_var = "z_plot",
  color_label = "mean Wald Z (X,Y)",
  color_mode = "div"
)

plot_section(
  expr_df,
  title = "Pairwise GREML h2: Expression effect (TPM vs TMM)",
  x_lab = "X-axis h2 (TMM)",
  y_lab = "Y-axis h2 (TPM)",
  out_file = file.path(runs_dir, paste0("pairwise_h2_scatter_", prefix, "_expression_tpm_vs_tmm_color_zscore.png")),
  color_var = "z_plot",
  color_label = "mean Wald Z (X,Y)",
  color_mode = "div"
)

plot_section(
  norm_df,
  title = "Pairwise GREML h2: Normalization effect (RAW vs IRNT)",
  x_lab = "X-axis h2 (RAW)",
  y_lab = "Y-axis h2 (IRNT)",
  out_file = file.path(runs_dir, paste0("pairwise_h2_scatter_", prefix, "_normalization_raw_vs_irnt_color_zscore.png")),
  color_var = "z_plot",
  color_label = "mean Wald Z (X,Y)",
  color_mode = "div"
)

# -----------------------------
# 7) Skewness + individual-level raw TMM distributions
# -----------------------------
skew_run_name <- method_to_run_name(skew_reference_method)
skew_data_dir <- file.path(runs_dir, skew_run_name, "results", "data")
skew_pheno_file <- list.files(skew_data_dir, pattern = "^pheno_.*\\.phenotypes\\.tsv$", full.names = TRUE)
skew_map_file <- list.files(skew_data_dir, pattern = "^pheno_.*\\.gene_index_map\\.txt$", full.names = TRUE)

if (length(skew_pheno_file) > 0 && length(skew_map_file) > 0) {
  skew_pheno <- fread(skew_pheno_file[1], header = FALSE, sep = "\t", data.table = TRUE)
  skew_map <- fread(skew_map_file[1], sep = "\t", data.table = TRUE)

  if (all(c("gene_name", "mpheno_index") %in% names(skew_map))) {
    skew_map <- skew_map[order(as.integer(mpheno_index))]
    skew_genes <- as.character(skew_map$gene_name)

    if (ncol(skew_pheno) >= 3 && (ncol(skew_pheno) - 2L) == length(skew_genes)) {
      setnames(skew_pheno, c("FID", "IID", skew_genes))

      skew_mat <- as.matrix(skew_pheno[, ..skew_genes])
      storage.mode(skew_mat) <- "numeric"

      skew_tbl <- tibble(
        Gene = skew_genes,
        n_individuals = colSums(is.finite(skew_mat)),
        mean_expr = colMeans(skew_mat, na.rm = TRUE),
        log_expr = log2(pmax(mean_expr, 0) + 1),
        skewness = apply(skew_mat, 2, calc_skewness)
      ) %>%
        arrange(desc(abs(skewness)))

      write_tsv(
        skew_tbl,
        file.path(runs_dir, paste0("tmm_raw_gene_skewness_", prefix, ".tsv"))
      )

      set.seed(1)
      random_genes <- sample(skew_genes, size = min(n_random_genes, length(skew_genes)), replace = FALSE)
      random_long <- as_tibble(skew_mat[, random_genes, drop = FALSE]) %>%
        pivot_longer(cols = everything(), names_to = "Gene", values_to = "expr") %>%
        filter(is.finite(expr))

      p_random <- random_long %>%
        ggplot(aes(x = expr)) +
        geom_histogram(bins = 80, fill = "#2A9D8F", color = "white", linewidth = 0.1) +
        labs(
          title = "TMM (RAW): distribution of expression values",
          subtitle = paste0("Randomly sampled ", length(random_genes), " genes across individuals"),
          x = "Raw TMM expression",
          y = "Count"
        ) +
        theme_minimal(base_size = 11)

      ggsave(
        file.path(runs_dir, paste0("tmm_raw_expression_random_distribution_", prefix, ".png")),
        p_random,
        width = 10,
        height = 5.5,
        dpi = 300
      )

      selected_skew <- bind_rows(
        skew_tbl %>% filter(is.finite(skewness)) %>% slice_max(skewness, n = n_skew_genes_per_tail, with_ties = FALSE),
        skew_tbl %>% filter(is.finite(skewness)) %>% slice_min(skewness, n = n_skew_genes_per_tail, with_ties = FALSE)
      ) %>%
        distinct(Gene, .keep_all = TRUE)

      qq_tbl <- build_qq_table(skew_mat, selected_skew$Gene) %>%
        left_join(selected_skew %>% select(Gene, skewness), by = "Gene") %>%
        mutate(
          facet_label = paste0(Gene, "\n", "Skew=", sprintf("%.2f", skewness))
        )

      if (nrow(qq_tbl) > 0) {
        p_qq <- qq_tbl %>%
          ggplot(aes(x = theo_q, y = obs_q)) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E76F51", linewidth = 0.35) +
          geom_point(color = "#2A9D8F", alpha = 0.5, size = 0.8) +
          facet_wrap(~facet_label, scales = "free", ncol = 4) +
          labs(
            title = "TMM (RAW): QQ plots for skewness-selected genes",
            subtitle = "X = theoretical normal quantiles, Y = standardized observed quantiles",
            x = "Theoretical quantiles",
            y = "Observed quantiles (standardized)"
          ) +
          theme_minimal(base_size = 11) +
          theme(
            panel.grid.minor = element_blank(),
            strip.text = element_text(face = "bold", size = 8.5)
          )

        ggsave(
          file.path(runs_dir, paste0("tmm_raw_selected_gene_qq_", prefix, ".png")),
          p_qq,
          width = 12,
          height = 8,
          dpi = 300
        )
      }
    }
  }
}

cat("Saved grouped plots and stats:\n")
cat("- ", file.path(runs_dir, paste0("pairwise_h2_scatter_", prefix, "_snp_set_all_vs_hm3_color_logexpr.png")), "\n", sep = "")
cat("- ", file.path(runs_dir, paste0("pairwise_h2_scatter_", prefix, "_expression_tpm_vs_tmm_color_logexpr.png")), "\n", sep = "")
cat("- ", file.path(runs_dir, paste0("pairwise_h2_scatter_", prefix, "_normalization_raw_vs_irnt_color_logexpr.png")), "\n", sep = "")
cat("- ", file.path(runs_dir, paste0("pairwise_h2_scatter_", prefix, "_snp_set_all_vs_hm3_color_zscore.png")), "\n", sep = "")
cat("- ", file.path(runs_dir, paste0("pairwise_h2_scatter_", prefix, "_expression_tpm_vs_tmm_color_zscore.png")), "\n", sep = "")
cat("- ", file.path(runs_dir, paste0("pairwise_h2_scatter_", prefix, "_normalization_raw_vs_irnt_color_zscore.png")), "\n", sep = "")
cat("- ", out_stats, "\n", sep = "")
cat("- ", out_discrepancy_types, "\n", sep = "")
cat("- ", file.path(runs_dir, paste0("tmm_raw_gene_skewness_", prefix, ".tsv")), "\n", sep = "")
cat("- ", file.path(runs_dir, paste0("tmm_raw_expression_random_distribution_", prefix, ".png")), "\n", sep = "")
cat("- ", file.path(runs_dir, paste0("tmm_raw_selected_gene_qq_", prefix, ".png")), "\n", sep = "")
