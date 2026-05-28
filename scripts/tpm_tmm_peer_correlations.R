# Purpose:
# - Helper or downstream analysis script for the GEUVADIS v2 GeneExpression workflow.
# - Used for phenotype preparation, result summarization, or plotting; see the filename and `nf/main.nf` for context.

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(magrittr)
})

# -----------------------------
# Settings
# -----------------------------
runs_dir <- "/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260428_GE_GEUVADIS_v2/GeneExpression/runs"
run_suffix <- "_peerauto_pmg0_npc5"
analysis_root <- file.path(runs_dir, "_analysis", "tpm_tmm_peer_correlations")
plots_dir <- file.path(analysis_root, "plots")
tables_dir <- file.path(analysis_root, "tables")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

out_gene <- file.path(tables_dir, "tpm_tmm_peer_gene_correlations.tsv")
out_summary <- file.path(tables_dir, "tpm_tmm_peer_correlation_summary.tsv")
out_plot <- file.path(plots_dir, "tpm_tmm_peer_correlation_boxplot.png")
out_scatter <- file.path(plots_dir, "tpm_tmm_peer_mean_scatter.png")

# -----------------------------
# Helpers
# -----------------------------
normalize_gene_id <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace("\\.[0-9]+$", "") %>%
    str_extract("ENSG[0-9]+|.+")
}

method_candidates <- function(method_id) {
  out <- method_id
  if (str_detect(method_id, "_irnt$")) {
    out <- c(out, str_replace(method_id, "_irnt$", "_inverse_normal"))
  }
  if (str_detect(method_id, "_inverse_normal$")) {
    out <- c(out, str_replace(method_id, "_inverse_normal$", "_irnt"))
  }
  unique(out)
}

find_run_dir <- function(method_id) {
  method_opts <- method_candidates(method_id)
  dir_tbl <- map_dfr(method_opts, function(mid) {
    patt <- paste0("^", mid, run_suffix, "([_].+)?$")
    dirs <- list.dirs(runs_dir, recursive = FALSE, full.names = TRUE) %>%
      keep(~ str_detect(basename(.x), patt))
    tibble(method_option = mid, run_dir = dirs)
  })

  if (nrow(dir_tbl) == 0) {
    stop("No run dir found for method: ", method_id)
  }

  valid_tbl <- dir_tbl %>%
    mutate(
      data_dir = file.path(run_dir, "results", "data"),
      has_pheno = map_lgl(data_dir, ~ length(list.files(.x, pattern = "^pheno_.*\\.phenotypes\\.tsv$", full.names = TRUE)) > 0),
      has_map = map_lgl(data_dir, ~ length(list.files(.x, pattern = "^pheno_.*\\.gene_index_map\\.txt$", full.names = TRUE)) > 0)
    ) %>%
    filter(has_pheno, has_map)

  if (nrow(valid_tbl) == 0) {
    stop("Run dir exists but phenotype/map files are missing for method: ", method_id)
  }

  exact <- valid_tbl %>%
    filter(basename(run_dir) == paste0(method_id, run_suffix))
  if (nrow(exact) > 0) {
    return(exact$run_dir[[1]])
  }

  valid_tbl %>%
    mutate(mtime = file.info(run_dir)$mtime) %>%
    arrange(desc(mtime)) %>%
    pull(run_dir) %>%
    .[[1]]
}

read_method_expression <- function(method_id) {
  run_dir <- find_run_dir(method_id)
  data_dir <- file.path(run_dir, "results", "data")

  pheno_file <- list.files(data_dir, pattern = "^pheno_.*\\.phenotypes\\.tsv$", full.names = TRUE) %>% .[[1]]
  map_file <- list.files(data_dir, pattern = "^pheno_.*\\.gene_index_map\\.txt$", full.names = TRUE) %>% .[[1]]

  pheno <- fread(pheno_file, header = FALSE, sep = "\t", data.table = TRUE)
  map_dt <- fread(map_file, sep = "\t", data.table = TRUE)

  if (!all(c("gene_name", "mpheno_index") %in% names(map_dt))) {
    stop("Unexpected gene map format: ", map_file)
  }

  map_dt <- map_dt[order(as.integer(mpheno_index))]
  genes <- normalize_gene_id(map_dt$gene_name)

  if ((ncol(pheno) - 2L) != length(genes)) {
    stop("Phenotype/map mismatch in run: ", run_dir)
  }

  expr_mat <- as.matrix(pheno[, -(1:2), with = FALSE])
  storage.mode(expr_mat) <- "numeric"
  rownames(expr_mat) <- paste0(as.character(pheno[[1]]), "::", as.character(pheno[[2]]))

  if (anyDuplicated(genes) > 0L) {
    idx <- split(seq_along(genes), genes)
    expr_mat <- do.call(
      cbind,
      lapply(idx, function(ii) {
        if (length(ii) == 1L) expr_mat[, ii] else rowMeans(expr_mat[, ii, drop = FALSE], na.rm = TRUE)
      })
    )
    colnames(expr_mat) <- names(idx)
  } else {
    colnames(expr_mat) <- genes
  }

  list(
    method = method_id,
    run_dir = run_dir,
    expr = expr_mat
  )
}

compute_pair <- function(snp_set, norm) {
  tmm_method <- paste(snp_set, "tmm", norm, sep = "_")
  tpm_method <- paste(snp_set, "tpm", norm, sep = "_")

  tmm <- read_method_expression(tmm_method)
  tpm <- read_method_expression(tpm_method)

  common_samples <- intersect(rownames(tmm$expr), rownames(tpm$expr))
  common_genes <- intersect(colnames(tmm$expr), colnames(tpm$expr))

  if (length(common_samples) < 3L || length(common_genes) == 0L) {
    return(tibble())
  }

  tmm_mat <- tmm$expr[common_samples, common_genes, drop = FALSE]
  tpm_mat <- tpm$expr[common_samples, common_genes, drop = FALSE]

  gene_cor <- map_dbl(seq_along(common_genes), ~ cor(tmm_mat[, .x], tpm_mat[, .x], use = "pairwise.complete.obs"))

  tibble(
    Gene = common_genes,
    snp_set = snp_set,
    norm = norm,
    panel = paste(if_else(snp_set == "all_snps", "ALL", "HM3"), toupper(norm), sep = " | "),
    cor_tpm_vs_tmm = gene_cor,
    mean_tmm = colMeans(tmm_mat, na.rm = TRUE),
    mean_tpm = colMeans(tpm_mat, na.rm = TRUE),
    n_samples = length(common_samples),
    n_genes_overlap = length(common_genes),
    run_tmm = basename(tmm$run_dir),
    run_tpm = basename(tpm$run_dir)
  )
}

# -----------------------------
# Run
# -----------------------------
combos <- crossing(
  snp_set = c("all_snps", "hm3_no_mhc"),
  norm = c("raw", "irnt")
)

cor_df <- pmap_dfr(combos, compute_pair) %>%
  filter(is.finite(cor_tpm_vs_tmm))

if (nrow(cor_df) == 0) {
  stop("No TPM/TMM gene-wise correlations were computed.")
}

summary_df <- cor_df %>%
  group_by(panel) %>%
  summarise(
    n_genes = n(),
    n_samples = first(n_samples),
    mean_cor = mean(cor_tpm_vs_tmm, na.rm = TRUE),
    median_cor = median(cor_tpm_vs_tmm, na.rm = TRUE),
    q25_cor = quantile(cor_tpm_vs_tmm, 0.25, na.rm = TRUE),
    q75_cor = quantile(cor_tpm_vs_tmm, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(match(panel, c("ALL | RAW", "ALL | IRNT", "HM3 | RAW", "HM3 | IRNT")))

write_tsv(cor_df, out_gene)
write_tsv(summary_df, out_summary)

panel_order <- c("ALL | RAW", "ALL | IRNT", "HM3 | RAW", "HM3 | IRNT")
plot_df <- cor_df %>%
  mutate(panel = factor(panel, levels = panel_order))

p <- ggplot(plot_df, aes(x = cor_tpm_vs_tmm, y = panel, fill = panel)) +
  geom_boxplot(width = 0.62, outlier.alpha = 0.15) +
  scale_fill_manual(values = c(
    "ALL | RAW" = "#1b9e77",
    "ALL | IRNT" = "#66a61e",
    "HM3 | RAW" = "#7570b3",
    "HM3 | IRNT" = "#e7298a"
  )) +
  labs(
    title = "TPM vs TMM Correlation Per Gene (post-PEER phenotype values)",
    x = "Per-gene Pearson correlation across individuals",
    y = "SNP set | normalization",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

ggsave(out_plot, p, width = 10, height = 4.8, dpi = 300)

# Per-gene mean expression scatter: TMM on X, TPM on Y
scatter_df <- cor_df %>%
  mutate(panel = factor(panel, levels = panel_order))

p_scatter <- ggplot(scatter_df, aes(x = mean_tmm, y = mean_tpm)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E76F51", linewidth = 0.4) +
  geom_point(alpha = 0.18, size = 0.6, color = "#2A9D8F") +
  geom_smooth(method = "lm", se = FALSE, color = "#264653", linewidth = 0.55) +
  facet_wrap(~panel, scales = "free", ncol = 2) +
  labs(
    title = "Per-gene Mean Expression: TPM vs TMM (post-PEER)",
    subtitle = "X-axis: mean TMM across individuals | Y-axis: mean TPM across individuals",
    x = "Mean TMM (per gene)",
    y = "Mean TPM (per gene)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(out_scatter, p_scatter, width = 10, height = 8, dpi = 300)

cat("Saved:\n")
cat("- ", out_gene, "\n", sep = "")
cat("- ", out_summary, "\n", sep = "")
cat("- ", out_plot, "\n", sep = "")
cat("- ", out_scatter, "\n", sep = "")
