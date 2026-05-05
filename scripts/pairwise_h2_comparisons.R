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
      filter(!is.na(x), !is.na(y), !(x == 0 & y == 0))
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

pair_stats <- pair_df %>%
  group_by(section, pair, facet_label, m1, m2, m1_label, m2_label) %>%
  summarise(
    n = n(),
    r = cor(x, y, use = "complete.obs"),
    slope = ifelse(n() >= 2 && sd(x, na.rm = TRUE) > 0, coef(lm(y ~ x))[2], NA_real_),
    intercept = ifelse(n() >= 2 && sd(x, na.rm = TRUE) > 0, coef(lm(y ~ x))[1], NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    slope_flag = case_when(
      is.na(slope) ~ "m=NA",
      abs(slope - 1) < 0.05 ~ "m≈1",
      slope > 1 ~ "m>1",
      TRUE ~ "m<1"
    ),
    eq_label = if_else(
      is.na(slope) | is.na(intercept),
      "y = NA",
      sprintf("y = %.3fx %+.3f", slope, intercept)
    ),
    lbl_top = sprintf("n=%d, r=%.2f\n%s (%s)", n, r, eq_label, slope_flag)
  )

# -----------------------------
# 4) Plot by section (tidyverse + non-default colors)
# -----------------------------
prefix <- ifelse(is.null(focus_method), "all_pairs", paste0("focus_", focus_method))
out_stats <- file.path(runs_dir, paste0("pairwise_h2_stats_", prefix, ".tsv"))
write_tsv(pair_stats, out_stats)

section_levels <- c("SNP set: ALL vs HM3", "Expression: TPM vs TMM", "Normalization: RAW vs IRNT", "Mixed changes")
plot_sections <- intersect(section_levels, unique(pair_df$section))
out_root <- file.path(runs_dir, paste0("pairwise_h2_plots_", prefix))
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

for (sec in plot_sections) {
  sec_df <- pair_df %>% filter(section == sec)
  sec_stats <- pair_stats %>% filter(section == sec)
  if (nrow(sec_df) == 0) next

  sec_slug <- sec %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
  sec_dir <- file.path(out_root, sec_slug)
  dir.create(sec_dir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_len(nrow(sec_stats))) {
    row_i <- sec_stats[i, ]
    pair_i <- sec_df %>% filter(pair == row_i$pair)
    if (nrow(pair_i) == 0) next

    p <- pair_i %>%
      ggplot(aes(x = x, y = y)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E76F51", linewidth = 0.4) +
      geom_point(color = "#2A9D8F", alpha = 0.4, size = 0.9) +
      geom_text(
        data = row_i,
        aes(x = -Inf, y = Inf, label = lbl_top),
        inherit.aes = FALSE,
        hjust = -0.1,
        vjust = 1.1,
        size = 3.0,
        color = "#264653"
      ) +
      labs(
        title = sec,
        subtitle = "PASS genes only; dropped rows where both estimates are zero",
        x = paste0("X: ", row_i$m1_label),
        y = paste0("Y: ", row_i$m2_label)
      ) +
      theme_minimal(base_size = 11) +
      theme(panel.grid.minor = element_blank())

    file_i <- file.path(sec_dir, paste0("scatter_", row_i$m1, "_vs_", row_i$m2, ".png"))
    ggsave(file_i, p, width = 6.8, height = 5.2, dpi = 300)
    cat("Saved plot:\n", file_i, "\n")
  }
}

cat("Saved stats:\n", out_stats, "\n")
