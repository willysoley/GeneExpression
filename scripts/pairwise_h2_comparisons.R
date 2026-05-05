suppressPackageStartupMessages({
  library(tidyverse)
})

# -----------------------------
# Settings
# -----------------------------
runs_dir <- "/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260428_GE_GEUVADIS_v2/GeneExpression/runs"

# 1-vs-all mode: set method id here, e.g. "all_snps_tmm_raw_peerauto_pmg0_npc5"
# all-pairs mode: keep NULL
focus_method <- NULL
include_mixed <- FALSE

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
  m <- str_match(method_id, "^(all_snps|hm3_no_mhc)_(tpm|tmm)_(irnt|raw)$")
  snp_lbl <- if_else(m[, 2] == "all_snps", "ALL", "HM3")
  expr_lbl <- toupper(m[, 3])
  norm_lbl <- toupper(m[, 4])
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

h2_long <- map_dfr(files, \(f) {
  run_name <- str_match(f, "runs/([^/]+)/results/summary/")[, 2]
  meta <- parse_method(run_name)

  if (any(is.na(meta$snp_set))) return(tibble())

  read_tsv(f, show_col_types = FALSE) %>%
    transmute(
      Gene,
      Status,
      h2 = as.numeric(h2_GREML),
      method = meta$method
    )
}) %>%
  filter(Status == "PASS", is.finite(h2))

if (nrow(h2_long) == 0) stop("No PASS h2 values found.")

# Collapse duplicates created by inverse_normal -> irnt harmonization
h2_long <- h2_long %>%
  group_by(Gene, method) %>%
  summarise(h2 = mean(h2, na.rm = TRUE), .groups = "drop")

# -----------------------------
# 2) Gene x method matrix
# -----------------------------
h2_wide <- h2_long %>%
  pivot_wider(names_from = method, values_from = h2)

methods <- setdiff(names(h2_wide), "Gene")
if (length(methods) < 2) stop("Need at least 2 methods for pairwise comparison.")

if (!is.null(focus_method)) {
  focus_method <- normalize_focus(focus_method)
  if (!focus_method %in% methods) stop("focus_method not found after normalization: ", focus_method)
}

# -----------------------------
# 3) Build pairwise table
# -----------------------------
pair_df <- combn(methods, 2, simplify = FALSE) %>%
  map_dfr(\(p) {
    h2_wide %>%
      transmute(
        Gene,
        m1 = p[1],
        m2 = p[2],
        x = .data[[p[1]]],
        y = .data[[p[2]]]
      ) %>%
      filter(!is.na(x), !is.na(y))
  }) %>%
  mutate(
    m1_snp = str_match(m1, "^(all_snps|hm3_no_mhc)_(tpm|tmm)_(irnt|raw)$")[, 2],
    m1_expr = str_match(m1, "^(all_snps|hm3_no_mhc)_(tpm|tmm)_(irnt|raw)$")[, 3],
    m1_norm = str_match(m1, "^(all_snps|hm3_no_mhc)_(tpm|tmm)_(irnt|raw)$")[, 4],
    m2_snp = str_match(m2, "^(all_snps|hm3_no_mhc)_(tpm|tmm)_(irnt|raw)$")[, 2],
    m2_expr = str_match(m2, "^(all_snps|hm3_no_mhc)_(tpm|tmm)_(irnt|raw)$")[, 3],
    m2_norm = str_match(m2, "^(all_snps|hm3_no_mhc)_(tpm|tmm)_(irnt|raw)$")[, 4],
    section = case_when(
      m1_snp != m2_snp & m1_expr == m2_expr & m1_norm == m2_norm ~ "SNP set: ALL vs HM3",
      m1_snp == m2_snp & m1_expr != m2_expr & m1_norm == m2_norm ~ "Expression: TPM vs TMM",
      m1_snp == m2_snp & m1_expr == m2_expr & m1_norm != m2_norm ~ "Normalization: RAW vs IRNT",
      TRUE ~ "Mixed changes"
    ),
    m1_label = method_to_label(m1),
    m2_label = method_to_label(m2),
    pair = paste(m1_label, "vs", m2_label),
    facet_label = paste0("X: ", m1_label, "\nY: ", m2_label)
  )

if (!is.null(focus_method)) {
  pair_df <- pair_df %>% filter(m1 == focus_method | m2 == focus_method)
}

if (!include_mixed) {
  pair_df <- pair_df %>% filter(section != "Mixed changes")
}

if (nrow(pair_df) == 0) stop("No rows left after filtering.")

# -----------------------------
# 4) Grouped section plots
# -----------------------------
calc_stats <- function(df) {
  df %>%
    group_by(facet_label) %>%
    summarise(
      n = n(),
      r = cor(x_plot, y_plot, use = "complete.obs"),
      n_abs_x_hi_y_lo = sum(x_plot > 0.05 & y_plot < 0.005, na.rm = TRUE),
      n_abs_y_hi_x_lo = sum(y_plot > 0.05 & x_plot < 0.005, na.rm = TRUE),
      n_fold_x_gt10y = sum(x_plot > 0 & y_plot > 0 & x_plot / y_plot >= 10, na.rm = TRUE),
      n_fold_y_gt10x = sum(x_plot > 0 & y_plot > 0 & y_plot / x_plot >= 10, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      lbl_top = sprintf(
        "n=%d, r=%.2f\nabs: X>0.05,Y<0.005=%d\nabs: Y>0.05,X<0.005=%d\nfold(>=10x): X>>Y=%d\nfold(>=10x): Y>>X=%d",
        n, r,
        n_abs_x_hi_y_lo, n_abs_y_hi_x_lo,
        n_fold_x_gt10y, n_fold_y_gt10x
      )
    )
}

plot_section <- function(df, title, x_lab, y_lab, out_file) {
  if (nrow(df) == 0) return(invisible(NULL))
  stats <- calc_stats(df)

  p <- df %>%
    ggplot(aes(x = x_plot, y = y_plot)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E76F51", linewidth = 0.35) +
    geom_point(color = "#2A9D8F", alpha = 0.35, size = 0.75) +
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
      y = y_lab
    ) +
    theme_minimal(base_size = 10.5) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 8.5)
    )

  ggsave(out_file, p, width = 14, height = 8, dpi = 300)
  print(p)
  invisible(stats)
}

snp_df <- pair_df %>%
  filter(section == "SNP set: ALL vs HM3") %>%
  mutate(
    x_plot = if_else(m1_snp == "all_snps", x, y),
    y_plot = if_else(m1_snp == "all_snps", y, x),
    expr_fix = toupper(if_else(m1_snp == "all_snps", m1_expr, m2_expr)),
    norm_fix = toupper(if_else(m1_snp == "all_snps", m1_norm, m2_norm)),
    facet_label = paste(expr_fix, norm_fix, sep = " | ")
  )

expr_df <- pair_df %>%
  filter(section == "Expression: TPM vs TMM") %>%
  mutate(
    x_plot = if_else(m1_expr == "tmm", x, y),
    y_plot = if_else(m1_expr == "tmm", y, x),
    snp_fix = if_else(if_else(m1_expr == "tmm", m1_snp, m2_snp) == "all_snps", "ALL", "HM3"),
    norm_fix = toupper(if_else(m1_expr == "tmm", m1_norm, m2_norm)),
    facet_label = paste(snp_fix, norm_fix, sep = " | ")
  )

norm_df <- pair_df %>%
  filter(section == "Normalization: RAW vs IRNT") %>%
  mutate(
    x_plot = if_else(m1_norm == "raw", x, y),
    y_plot = if_else(m1_norm == "raw", y, x),
    snp_fix = if_else(if_else(m1_norm == "raw", m1_snp, m2_snp) == "all_snps", "ALL", "HM3"),
    expr_fix = toupper(if_else(m1_norm == "raw", m1_expr, m2_expr)),
    facet_label = paste(snp_fix, expr_fix, sep = " | ")
  )

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

plot_section(
  snp_df,
  title = "Pairwise GREML h2: SNP effect (ALL vs HM3)",
  x_lab = "X-axis h2 (ALL)",
  y_lab = "Y-axis h2 (HM3)",
  out_file = file.path(runs_dir, paste0("pairwise_h2_scatter_", prefix, "_snp_set_all_vs_hm3.png"))
)

plot_section(
  expr_df,
  title = "Pairwise GREML h2: Expression effect (TPM vs TMM)",
  x_lab = "X-axis h2 (TMM)",
  y_lab = "Y-axis h2 (TPM)",
  out_file = file.path(runs_dir, paste0("pairwise_h2_scatter_", prefix, "_expression_tpm_vs_tmm.png"))
)

plot_section(
  norm_df,
  title = "Pairwise GREML h2: Normalization effect (RAW vs IRNT)",
  x_lab = "X-axis h2 (RAW)",
  y_lab = "Y-axis h2 (IRNT)",
  out_file = file.path(runs_dir, paste0("pairwise_h2_scatter_", prefix, "_normalization_raw_vs_irnt.png"))
)

cat("Saved grouped plots and stats:\n")
cat("- ", file.path(runs_dir, paste0("pairwise_h2_scatter_", prefix, "_snp_set_all_vs_hm3.png")), "\n", sep = "")
cat("- ", file.path(runs_dir, paste0("pairwise_h2_scatter_", prefix, "_expression_tpm_vs_tmm.png")), "\n", sep = "")
cat("- ", file.path(runs_dir, paste0("pairwise_h2_scatter_", prefix, "_normalization_raw_vs_irnt.png")), "\n", sep = "")
cat("- ", out_stats, "\n", sep = "")
cat("- ", out_discrepancy_types, "\n", sep = "")
