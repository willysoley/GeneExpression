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

# Example:
# METHOD_ID=all_snps_tmm_raw N_SAMPLE_GENES=2000 RNG_SEED=1 N_EXAMPLES_EACH=5 Rscript script.R
method_id <- Sys.getenv("METHOD_ID", "all_snps_tmm_raw")
n_sample_genes <- as.integer(Sys.getenv("N_SAMPLE_GENES", "2000"))
rng_seed <- as.integer(Sys.getenv("RNG_SEED", "1"))
n_examples_each <- as.integer(Sys.getenv("N_EXAMPLES_EACH", "5"))

h2_hi_cutoff <- as.numeric(Sys.getenv("H2_HI_CUTOFF", "0.05"))
h2_lo_cutoff <- as.numeric(Sys.getenv("H2_LO_CUTOFF", "0.005"))

skew_hist_bins <- as.integer(Sys.getenv("SKEW_HIST_BINS", "320"))
expr_hist_bins <- as.integer(Sys.getenv("EXPR_HIST_BINS", "160"))

# Stricter low-zero filter for high-skew examples
max_zero_frac_for_nonzero_high_skew <- as.numeric(Sys.getenv("MAX_ZERO_FRAC_NONZERO_HIGH_SKEW", "0.05"))
# Majority of counts >10
min_frac_gt10_for_high_skew <- as.numeric(Sys.getenv("MIN_FRAC_GT10_FOR_HIGH_SKEW", "0.50"))

if (!is.finite(n_sample_genes) || n_sample_genes < 10) stop("N_SAMPLE_GENES must be >= 10.")
if (!is.finite(n_examples_each) || n_examples_each < 1) stop("N_EXAMPLES_EACH must be >= 1.")
if (!is.finite(skew_hist_bins) || skew_hist_bins < 20) stop("SKEW_HIST_BINS must be >= 20.")
if (!is.finite(expr_hist_bins) || expr_hist_bins < 20) stop("EXPR_HIST_BINS must be >= 20.")
if (!is.finite(max_zero_frac_for_nonzero_high_skew) || max_zero_frac_for_nonzero_high_skew < 0 || max_zero_frac_for_nonzero_high_skew > 1) {
  stop("MAX_ZERO_FRAC_NONZERO_HIGH_SKEW must be between 0 and 1.")
}
if (!is.finite(min_frac_gt10_for_high_skew) || min_frac_gt10_for_high_skew < 0 || min_frac_gt10_for_high_skew > 1) {
  stop("MIN_FRAC_GT10_FOR_HIGH_SKEW must be between 0 and 1.")
}

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

method_candidates <- function(method) {
  out <- method
  if (str_detect(method, "_irnt$")) out <- c(out, str_replace(method, "_irnt$", "_inverse_normal"))
  if (str_detect(method, "_inverse_normal$")) out <- c(out, str_replace(method, "_inverse_normal$", "_irnt"))
  unique(out)
}

find_run_dir <- function(method) {
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
      has_pheno = map_lgl(data_dir, ~ length(list.files(.x, pattern = "^pheno_.*\\.phenotypes\\.tsv$", full.names = TRUE)) > 0),
      has_map = map_lgl(data_dir, ~ length(list.files(.x, pattern = "^pheno_.*\\.gene_index_map\\.txt$", full.names = TRUE)) > 0)
    ) %>%
    filter(has_pheno, has_map)
  if (nrow(valid) == 0) stop("Run directory exists but phenotype/map files are missing for method: ", method)

  exact <- valid %>% filter(basename(run_dir) == paste0(method, run_suffix))
  if (nrow(exact) > 0) return(exact$run_dir[[1]])

  valid %>%
    mutate(mtime = file.info(run_dir)$mtime) %>%
    arrange(desc(mtime)) %>%
    pull(run_dir) %>%
    .[[1]]
}

read_expression_matrix <- function(method) {
  run_dir <- find_run_dir(method)
  data_dir <- file.path(run_dir, "results", "data")

  pheno_file <- list.files(data_dir, pattern = "^pheno_.*\\.phenotypes\\.tsv$", full.names = TRUE) %>% .[[1]]
  map_file <- list.files(data_dir, pattern = "^pheno_.*\\.gene_index_map\\.txt$", full.names = TRUE) %>% .[[1]]

  pheno <- fread(pheno_file, header = FALSE, sep = "\t", data.table = TRUE)
  map_dt <- fread(map_file, sep = "\t", data.table = TRUE)
  if (!all(c("gene_name", "mpheno_index") %in% names(map_dt))) stop("Unexpected map format: ", map_file)

  map_dt <- map_dt[order(as.integer(mpheno_index))]
  genes <- normalize_gene_id(map_dt$gene_name)
  if ((ncol(pheno) - 2L) != length(genes)) stop("Phenotype/map dimensions mismatch for run: ", run_dir)

  expr <- as.matrix(pheno[, -(1:2), with = FALSE])
  storage.mode(expr) <- "numeric"
  rownames(expr) <- paste0(as.character(pheno[[1]]), "::", as.character(pheno[[2]]))

  if (anyDuplicated(genes) > 0L) {
    idx <- split(seq_along(genes), genes)
    expr <- do.call(cbind, lapply(idx, function(ii) {
      if (length(ii) == 1L) expr[, ii] else rowMeans(expr[, ii, drop = FALSE], na.rm = TRUE)
    }))
    colnames(expr) <- names(idx)
  } else {
    colnames(expr) <- genes
  }

  list(run_dir = run_dir, expr = expr)
}

skewness_third_moment <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  mu <- mean(x)
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  mean(((x - mu) / s)^3)
}

resolve_gene_ids_for_run <- function(gene_raw, run_dir) {
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

read_h2_summary <- function(method) {
  run_dir <- find_run_dir(method)
  sum_file <- list.files(file.path(run_dir, "results", "summary"), pattern = "^final_heritability_summary.*\\.tsv$", full.names = TRUE)
  if (length(sum_file) == 0) stop("Missing summary for method: ", method)

  dt <- fread(sum_file[1], sep = "\t", data.table = TRUE)
  tibble(
    Gene = resolve_gene_ids_for_run(dt$Gene, run_dir),
    Status = toupper(str_trim(as.character(dt$Status))),
    h2 = suppressWarnings(as.numeric(dt$h2_GREML))
  ) %>%
    filter(is.finite(h2))
}

# -----------------------------
# Compute skewness on random sample
# -----------------------------
tmm_obj <- read_expression_matrix(method_id)
tmm_expr <- tmm_obj$expr
all_genes <- colnames(tmm_expr)
if (length(all_genes) < 50) stop("Too few genes found in matrix.")

set.seed(rng_seed)
n_take <- min(n_sample_genes, length(all_genes))
sampled_genes <- sample(all_genes, size = n_take, replace = FALSE)
tmm_sub <- tmm_expr[, sampled_genes, drop = FALSE]

skew_tbl <- tibble(
  Gene = sampled_genes,
  mean_expr = colMeans(tmm_sub, na.rm = TRUE),
  median_expr = apply(tmm_sub, 2, median, na.rm = TRUE),
  sd_expr = apply(tmm_sub, 2, sd, na.rm = TRUE),
  zero_fraction = colMeans(tmm_sub == 0, na.rm = TRUE),
  frac_gt10 = colMeans(tmm_sub > 10, na.rm = TRUE),
  skewness = apply(tmm_sub, 2, skewness_third_moment)
) %>%
  filter(is.finite(skewness))
if (nrow(skew_tbl) == 0) stop("No finite skewness values were computed.")

summary_tbl <- skew_tbl %>%
  summarise(
    method = method_id,
    run_dir = tmm_obj$run_dir,
    n_individuals = nrow(tmm_sub),
    n_genes_requested = n_sample_genes,
    n_genes_used = nrow(skew_tbl),
    mean_skewness = mean(skewness, na.rm = TRUE),
    median_skewness = median(skewness, na.rm = TRUE),
    q25_skewness = quantile(skewness, 0.25, na.rm = TRUE),
    q75_skewness = quantile(skewness, 0.75, na.rm = TRUE),
    frac_right_skew = mean(skewness > 0, na.rm = TRUE)
  )

# -----------------------------
# Output paths
# -----------------------------
slug <- str_replace_all(method_id, "[^A-Za-z0-9_]+", "_")
out_gene <- file.path(runs_dir, paste0("tmm_raw_skewness_gene_table_", slug, ".tsv"))
out_summary <- file.path(runs_dir, paste0("tmm_raw_skewness_summary_", slug, ".tsv"))
out_hist <- file.path(runs_dir, paste0("tmm_raw_skewness_hist_density_", slug, ".png"))
out_scatter <- file.path(runs_dir, paste0("tmm_raw_skewness_mean_vs_skew_", slug, ".png"))
out_skew_examples_tsv <- file.path(runs_dir, paste0("tmm_raw_skewness_selected_genes_", slug, ".tsv"))
out_skew_examples_plot <- file.path(runs_dir, paste0("tmm_raw_skewness_selected_gene_distributions_", slug, ".png"))
out_h2_examples_tsv <- file.path(runs_dir, paste0("tmm_vs_tpm_h2_discrepancy_selected_genes_", slug, ".tsv"))
out_h2_examples_plot <- file.path(runs_dir, paste0("tmm_vs_tpm_h2_discrepancy_expr_distributions_", slug, ".png"))

write_tsv(skew_tbl, out_gene)
write_tsv(summary_tbl, out_summary)

# 1) Skewness distribution
p_hist <- ggplot(skew_tbl, aes(x = skewness)) +
  geom_histogram(aes(y = after_stat(density)), bins = skew_hist_bins, fill = "#2A9D8F", color = "white", alpha = 0.8) +
  geom_density(color = "#264653", linewidth = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#E76F51", linewidth = 0.4) +
  labs(
    title = paste0("Gene-wise skewness (third moment) - ", method_id),
    subtitle = paste0("Random sample of ", nrow(skew_tbl), " genes"),
    x = "Skewness",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
ggsave(out_hist, p_hist, width = 10, height = 5, dpi = 300)

# 2) Mean expression vs skewness
p_scatter <- ggplot(skew_tbl, aes(x = log2(mean_expr + 1), y = skewness)) +
  geom_point(alpha = 0.35, size = 1.0, color = "#1F78B4") +
  geom_smooth(method = "loess", se = TRUE, color = "#D95F02", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  labs(
    title = paste0("Mean expression vs skewness - ", method_id),
    subtitle = "Each point is one gene",
    x = "log2(mean expression + 1)",
    y = "Skewness (third moment)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
ggsave(out_scatter, p_scatter, width = 9, height = 5.5, dpi = 300)

# 3) Combined requested skewness examples
sel_high <- skew_tbl %>%
  arrange(desc(skewness)) %>%
  slice_head(n = n_examples_each) %>%
  mutate(group = "High positive skew")

sel_low <- skew_tbl %>%
  arrange(skewness) %>%
  slice_head(n = n_examples_each) %>%
  mutate(group = "Low negative skew")

sel_zero <- skew_tbl %>%
  mutate(abs_skew = abs(skewness)) %>%
  arrange(abs_skew) %>%
  slice_head(n = n_examples_each) %>%
  select(-abs_skew) %>%
  mutate(group = "Near zero skew")

sel_high_low_zero <- skew_tbl %>%
  filter(zero_fraction <= max_zero_frac_for_nonzero_high_skew) %>%
  arrange(desc(skewness)) %>%
  slice_head(n = n_examples_each) %>%
  mutate(group = paste0("High positive skew | zero<=", sprintf("%.2f", max_zero_frac_for_nonzero_high_skew)))

if (nrow(sel_high_low_zero) == 0) {
  warning("No genes met strict zero_fraction cutoff. Using fallback selection.")
  sel_high_low_zero <- skew_tbl %>%
    arrange(desc(skewness), zero_fraction) %>%
    slice_head(n = n_examples_each) %>%
    mutate(group = paste0("High positive skew | zero<=", sprintf("%.2f", max_zero_frac_for_nonzero_high_skew), " (fallback)"))
}

sel_high_major_gt10 <- skew_tbl %>%
  filter(frac_gt10 >= min_frac_gt10_for_high_skew) %>%
  arrange(desc(skewness)) %>%
  slice_head(n = n_examples_each) %>%
  mutate(group = paste0("High positive skew | frac(>10)>=", sprintf("%.2f", min_frac_gt10_for_high_skew)))

if (nrow(sel_high_major_gt10) == 0) {
  warning("No genes met frac_gt10 cutoff. Using fallback selection.")
  sel_high_major_gt10 <- skew_tbl %>%
    arrange(desc(skewness), desc(frac_gt10)) %>%
    slice_head(n = n_examples_each) %>%
    mutate(group = paste0("High positive skew | frac(>10)>=", sprintf("%.2f", min_frac_gt10_for_high_skew), " (fallback)"))
}

skew_examples <- bind_rows(
  sel_high, sel_low, sel_zero, sel_high_low_zero, sel_high_major_gt10
) %>%
  group_by(group) %>%
  mutate(rank_within_group = row_number()) %>%
  ungroup() %>%
  mutate(rank_id = row_number())

write_tsv(skew_examples, out_skew_examples_tsv)

skew_plot_df <- as_tibble(tmm_sub[, unique(skew_examples$Gene), drop = FALSE]) %>%
  mutate(sample_id = rownames(tmm_sub)) %>%
  pivot_longer(cols = -sample_id, names_to = "Gene", values_to = "expr") %>%
  left_join(skew_examples %>% select(Gene, group, rank_id, skewness, zero_fraction, frac_gt10), by = "Gene") %>%
  mutate(
    facet_label = paste0(
      group, "\n", Gene,
      "\nskew=", sprintf("%.2f", skewness),
      ", zero=", sprintf("%.2f", zero_fraction),
      ", frac>10=", sprintf("%.2f", frac_gt10)
    ),
    facet_label = factor(facet_label, levels = unique(facet_label[order(rank_id)]))
  )

p_skew_examples <- ggplot(skew_plot_df, aes(x = expr)) +
  geom_histogram(bins = expr_hist_bins, fill = "#4DB6AC", color = "white", alpha = 0.9) +
  facet_wrap(~ facet_label, scales = "free_y", ncol = 4) +
  labs(
    title = paste0("Selected skewness example genes - ", method_id),
    subtitle = "Includes high+/low-/near-zero skew, strict low-zero high-skew, and high-skew majority >10",
    x = "Expression value (post-PEER phenotype values)",
    y = "Number of individuals"
  ) +
  theme_minimal(base_size = 10.5) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 8))
ggsave(out_skew_examples_plot, p_skew_examples, width = 18, height = 10, dpi = 300)

# 4) TPM vs TMM h2 discrepancy examples + skewness
method_match <- str_match(method_id, "^(all_snps|hm3_no_mhc)_tmm_raw$")
if (!is.na(method_match[1, 2])) {
  snp_set <- method_match[1, 2]
  tpm_method <- paste0(snp_set, "_tpm_raw")

  h2_tmm <- read_h2_summary(method_id) %>% filter(Status == "PASS") %>% select(Gene, h2_tmm = h2)
  h2_tpm <- read_h2_summary(tpm_method) %>% filter(Status == "PASS") %>% select(Gene, h2_tpm = h2)
  h2_join <- h2_tmm %>% inner_join(h2_tpm, by = "Gene")

  tpm_high <- h2_join %>%
    filter(h2_tpm > h2_hi_cutoff, h2_tmm < h2_lo_cutoff) %>%
    arrange(desc(h2_tpm - h2_tmm)) %>%
    slice_head(n = n_examples_each) %>%
    mutate(group = "TPM high, TMM near-zero")

  tmm_high <- h2_join %>%
    filter(h2_tmm > h2_hi_cutoff, h2_tpm < h2_lo_cutoff) %>%
    arrange(desc(h2_tmm - h2_tpm)) %>%
    slice_head(n = n_examples_each) %>%
    mutate(group = "TMM high, TPM near-zero")

  h2_examples <- bind_rows(tpm_high, tmm_high) %>%
    distinct(Gene, .keep_all = TRUE) %>%
    left_join(skew_tbl %>% select(Gene, skewness, zero_fraction, frac_gt10, mean_expr), by = "Gene")

  write_tsv(h2_examples, out_h2_examples_tsv)

  if (nrow(h2_examples) > 0) {
    tpm_obj <- read_expression_matrix(tpm_method)
    common_samples <- intersect(rownames(tmm_obj$expr), rownames(tpm_obj$expr))
    common_genes <- intersect(h2_examples$Gene, intersect(colnames(tmm_obj$expr), colnames(tpm_obj$expr)))

    if (length(common_samples) >= 3 && length(common_genes) > 0) {
      skew_tmm_vals <- map_dbl(common_genes, ~ skewness_third_moment(tmm_obj$expr[common_samples, .x]))
      skew_tpm_vals <- map_dbl(common_genes, ~ skewness_third_moment(tpm_obj$expr[common_samples, .x]))

      h2_examples <- h2_examples %>%
        left_join(
          tibble(
            Gene = common_genes,
            skewness_tmm_raw = skew_tmm_vals,
            skewness_tpm_raw = skew_tpm_vals,
            skewness_delta_tpm_minus_tmm = skew_tpm_vals - skew_tmm_vals
          ),
          by = "Gene"
        )
      write_tsv(h2_examples, out_h2_examples_tsv)

      tmm_long <- as_tibble(tmm_obj$expr[common_samples, common_genes, drop = FALSE]) %>%
        mutate(sample_id = common_samples) %>%
        pivot_longer(cols = -sample_id, names_to = "Gene", values_to = "expr") %>%
        mutate(method = "TMM RAW")

      tpm_long <- as_tibble(tpm_obj$expr[common_samples, common_genes, drop = FALSE]) %>%
        mutate(sample_id = common_samples) %>%
        pivot_longer(cols = -sample_id, names_to = "Gene", values_to = "expr") %>%
        mutate(method = "TPM RAW")

      h2_plot_df <- bind_rows(tmm_long, tpm_long) %>%
        left_join(h2_examples, by = "Gene") %>%
        mutate(
          facet_label = paste0(
            group, "\n", Gene,
            "\nTMM h2=", sprintf("%.3f", h2_tmm), ", TPM h2=", sprintf("%.3f", h2_tpm),
            "\nskew(TMM)=", sprintf("%.2f", skewness_tmm_raw), ", skew(TPM)=", sprintf("%.2f", skewness_tpm_raw)
          )
        )

      p_h2_examples <- ggplot(h2_plot_df, aes(x = expr, color = method, fill = method)) +
        geom_density(alpha = 0.15, linewidth = 0.7) +
        facet_wrap(~ facet_label, scales = "free_y", ncol = 2) +
        scale_color_manual(values = c("TMM RAW" = "#1F78B4", "TPM RAW" = "#D95F02")) +
        scale_fill_manual(values = c("TMM RAW" = "#1F78B4", "TPM RAW" = "#D95F02")) +
        labs(
          title = paste0("Expression distributions for h2-discrepant genes - ", snp_set, " RAW"),
          subtitle = "Near-zero vs high h2 between TPM and TMM",
          x = "Expression value (post-PEER phenotype values)",
          y = "Density",
          color = NULL,
          fill = NULL
        ) +
        theme_minimal(base_size = 10.5) +
        theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 8), legend.position = "top")

      ggsave(out_h2_examples_plot, p_h2_examples, width = 13, height = 8.5, dpi = 300)
    }
  }
} else {
  write_tsv(tibble(), out_h2_examples_tsv)
  message("Skipping h2-discrepancy example selection because METHOD_ID is not *_tmm_raw.")
}

cat("Saved:\n")
cat("- ", out_gene, "\n", sep = "")
cat("- ", out_summary, "\n", sep = "")
cat("- ", out_hist, "\n", sep = "")
cat("- ", out_scatter, "\n", sep = "")
cat("- ", out_skew_examples_tsv, "\n", sep = "")
cat("- ", out_skew_examples_plot, "\n", sep = "")
cat("- ", out_h2_examples_tsv, "\n", sep = "")
cat("- ", out_h2_examples_plot, "\n", sep = "")

