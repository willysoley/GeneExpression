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
    m1_label = method_to_label(m1),
    m2_label = method_to_label(m2),
    pair = paste(m1_label, "vs", m2_label)
  )

if (!is.null(focus_method)) {
  pair_df <- pair_df %>% filter(m1 == focus_method | m2 == focus_method)
}

if (nrow(pair_df) == 0) stop("No rows left after filtering.")

pair_stats <- pair_df %>%
  group_by(pair, m1, m2) %>%
  summarise(
    n = n(),
    r = cor(x, y, use = "complete.obs"),
    .groups = "drop"
  )

# -----------------------------
# 4) Plot (tidyverse + non-default colors)
# -----------------------------
p <- pair_df %>%
  ggplot(aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E76F51", linewidth = 0.35) +
  geom_point(color = "#2A9D8F", alpha = 0.35, size = 0.6) +
  geom_text(
    data = pair_stats,
    aes(x = -Inf, y = Inf, label = sprintf("n=%d, r=%.2f", n, r)),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.1,
    size = 2.3,
    color = "#264653"
  ) +
  facet_wrap(~pair, scales = "free", ncol = ifelse(is.null(focus_method), 5, 4)) +
  labs(
    title = ifelse(
      is.null(focus_method),
      "All pairwise GREML h2 comparisons",
      paste0("GREML h2 comparisons vs ", method_to_label(focus_method))
    ),
    subtitle = "PASS genes only; dropped rows where both estimates are zero",
    x = "h2 estimate",
    y = "h2 estimate"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 7)
  )

print(p)

prefix <- ifelse(is.null(focus_method), "all_pairs", paste0("focus_", focus_method))
out_plot <- file.path(runs_dir, paste0("pairwise_h2_scatter_", prefix, ".png"))
out_stats <- file.path(runs_dir, paste0("pairwise_h2_stats_", prefix, ".tsv"))

ggsave(out_plot, p, width = 16, height = ifelse(is.null(focus_method), 12, 8), dpi = 300)
write_tsv(pair_stats, out_stats)

cat("Saved:\n", out_plot, "\n", out_stats, "\n")
