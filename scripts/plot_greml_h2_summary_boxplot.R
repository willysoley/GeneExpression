suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(magrittr)
})

# -----------------------------
# 1) Inputs
# -----------------------------
base_dir <- "/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260428_GE_GEUVADIS_v2/GeneExpression/runs"

# If you want to manually specify folders, replace this block with your vector.
run_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE) %>%
  keep(~ str_detect(basename(.x), "^(all_snps|hm3_no_mhc)_(tpm|tmm)_(irnt|inverse_normal|raw)_peerauto_pmg0_npc5$")) %>%
  sort()

if (length(run_dirs) == 0) {
  stop("No run directories found under: ", base_dir)
}

# -----------------------------
# 2) Find summary files
# -----------------------------
summary_index <- tibble(run_dir = run_dirs) %>%
  mutate(
    run_name = basename(run_dir),
    summary_file = map_chr(run_dir, ~ {
      f <- list.files(
        path = file.path(.x, "results", "summary"),
        pattern = "^final_heritability_summary.*\\.tsv$",
        full.names = TRUE
      )
      if (length(f) == 0) NA_character_ else f[1]
    })
  )

if (any(is.na(summary_index$summary_file))) {
  missing_runs <- summary_index %>% filter(is.na(summary_file)) %>% pull(run_name)
  warning("Missing summary file for: ", paste(missing_runs, collapse = ", "))
}

summary_index <- summary_index %>% filter(!is.na(summary_file))

if (nrow(summary_index) == 0) {
  stop("No summary files found.")
}

# -----------------------------
# 3) Parse run type from folder name
# -----------------------------
parse_run_name <- function(x) {
  m <- str_match(x, "^(all_snps|hm3_no_mhc)_(tpm|tmm)_(irnt|inverse_normal|raw)_peerauto_pmg0_npc5$")

  out <- tibble(
    run_name = x,
    snp_set = m[, 2],
    expr_method = m[, 3],
    norm_method = m[, 4]
  )

  bad <- which(!complete.cases(out))
  if (length(bad) > 0) {
    stop("Run-name parse failed for: ", paste(out$run_name[bad], collapse = ", "))
  }
  out
}

run_meta <- parse_run_name(summary_index$run_name)

# -----------------------------
# 4) Read and combine summary tables
# -----------------------------
all_results <- summary_index %>%
  left_join(run_meta, by = "run_name") %>%
  mutate(dt = map(summary_file, fread)) %>%
  transmute(run_name, snp_set, expr_method, norm_method, dt) %>%
  mutate(
    tbl = pmap(
      list(dt, run_name, snp_set, expr_method, norm_method),
      ~ as_tibble(..1) %>%
        mutate(
          run_name = ..2,
          snp_set = ..3,
          expr_method = ..4,
          norm_method = ..5
        )
    )
  ) %>%
  select(tbl) %>%
  unnest(tbl) %>%
  mutate(
    h2_GREML = suppressWarnings(as.numeric(h2_GREML)),
    snp_set = recode(snp_set,
      "hm3_no_mhc" = "HM3 no MHC",
      "all_snps" = "All SNPs"
    ),
    expr_method = toupper(expr_method),
    norm_method = recode(norm_method,
      "raw" = "Raw",
      "irnt" = "IRNT",
      "inverse_normal" = "Inverse normal"
    ),
    analysis_type = paste(snp_set, expr_method, norm_method, sep = " | ")
  )

# -----------------------------
# 5) Summary-of-summaries table
# -----------------------------
summary_of_summaries <- all_results %>%
  group_by(analysis_type, snp_set, expr_method, norm_method) %>%
  summarise(
    n_total = n(),
    n_pass = sum(Status == "PASS", na.rm = TRUE),
    n_fail = sum(Status == "FAIL", na.rm = TRUE),
    median_h2 = median(h2_GREML[Status == "PASS"], na.rm = TRUE),
    mean_h2 = mean(h2_GREML[Status == "PASS"], na.rm = TRUE),
    q1_h2 = quantile(h2_GREML[Status == "PASS"], 0.25, na.rm = TRUE),
    q3_h2 = quantile(h2_GREML[Status == "PASS"], 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(snp_set, expr_method, norm_method)

print(summary_of_summaries, n = nrow(summary_of_summaries))

fwrite(
  summary_of_summaries,
  file.path(base_dir, "summary_of_summaries.tsv"),
  sep = "\t"
)

# -----------------------------
# 6) Horizontal boxplot
# -----------------------------
# Explicit hierarchical order for y-axis:
# SNP set -> expression method -> normalization
snp_levels <- c("All SNPs", "HM3 no MHC")
expr_levels <- c("TPM", "TMM")
norm_levels <- c("Raw", "IRNT", "Inverse normal")

analysis_levels <- crossing(
  snp_set = snp_levels,
  expr_method = expr_levels,
  norm_method = norm_levels
) %>%
  transmute(analysis_type = paste(snp_set, expr_method, norm_method, sep = " | ")) %>%
  pull(analysis_type)

plot_df <- all_results %>%
  filter(Status == "PASS", is.finite(h2_GREML)) %>%
  mutate(
    snp_set = factor(snp_set, levels = snp_levels),
    expr_method = factor(expr_method, levels = expr_levels),
    norm_method = factor(norm_method, levels = norm_levels),
    analysis_type = factor(analysis_type, levels = rev(analysis_levels))
  )

plot_stats <- plot_df %>%
  group_by(analysis_type, snp_set) %>%
  summarise(
    q1 = quantile(h2_GREML, 0.25, na.rm = TRUE),
    mean_h2 = mean(h2_GREML, na.rm = TRUE),
    median_h2 = median(h2_GREML, na.rm = TRUE),
    q3 = quantile(h2_GREML, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

x_min <- min(plot_df$h2_GREML, na.rm = TRUE)
x_max <- max(plot_df$h2_GREML, na.rm = TRUE)
x_span <- max(x_max - x_min, 1e-6)

plot_stats <- plot_stats %>%
  mutate(
    label_x = x_max + (0.25 * x_span),
    label = sprintf(
      "Q1=%.3f | Mean=%.3f | Median=%.3f | Q3=%.3f",
      q1, mean_h2, median_h2, q3
    )
  )

p <- plot_df %>%
  ggplot(aes(x = h2_GREML, y = analysis_type, fill = snp_set)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_boxplot(width = 0.65, outlier.alpha = 0.2, size = 0.3) +
  geom_point(
    data = plot_stats,
    aes(x = mean_h2, y = analysis_type),
    inherit.aes = FALSE,
    shape = 23,
    size = 2.2,
    fill = "white",
    color = "black",
    stroke = 0.4
  ) +
  geom_text(
    data = plot_stats,
    aes(x = label_x, y = analysis_type, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3
  ) +
  coord_cartesian(
    xlim = c(x_min, x_max + (0.9 * x_span)),
    clip = "off"
  ) +
  labs(
    title = "GREML h2 by analysis setting",
    subtitle = "PASS genes only (white diamond = mean; right text = Q1/Mean/Median/Q3)",
    x = "GREML h2 estimate",
    y = "SNP set | expression method | normalization",
    fill = "SNP set"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 300, 5.5, 5.5)
  )

print(p)

ggsave(
  filename = file.path(base_dir, "greml_h2_boxplot_by_setting.png"),
  plot = p,
  width = 18,
  height = 7,
  dpi = 300
)

cat("\nSaved:\n")
cat("- ", file.path(base_dir, "summary_of_summaries.tsv"), "\n", sep = "")
cat("- ", file.path(base_dir, "greml_h2_boxplot_by_setting.png"), "\n", sep = "")
