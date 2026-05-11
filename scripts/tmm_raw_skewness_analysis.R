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

# Override via env vars if needed:
# METHOD_ID=hm3_no_mhc_tmm_raw N_SAMPLE_GENES=2000 RNG_SEED=1 Rscript script.R
method_id <- Sys.getenv("METHOD_ID", "all_snps_tmm_raw")
n_sample_genes <- as.integer(Sys.getenv("N_SAMPLE_GENES", "2000"))
rng_seed <- as.integer(Sys.getenv("RNG_SEED", "1"))
n_example_genes <- 12L

if (!is.finite(n_sample_genes) || n_sample_genes < 10) {
  stop("N_SAMPLE_GENES must be >= 10.")
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
  if (str_detect(method, "_irnt$")) {
    out <- c(out, str_replace(method, "_irnt$", "_inverse_normal"))
  }
  if (str_detect(method, "_inverse_normal$")) {
    out <- c(out, str_replace(method, "_inverse_normal$", "_irnt"))
  }
  unique(out)
}

find_run_dir <- function(method) {
  method_opts <- method_candidates(method)
  all_dirs <- list.dirs(runs_dir, recursive = FALSE, full.names = TRUE)

  hits <- map_dfr(method_opts, function(m) {
    patt <- paste0("^", m, run_suffix, "([_].+)?$")
    tibble(run_dir = all_dirs[str_detect(basename(all_dirs), patt)])
  })

  if (nrow(hits) == 0) {
    stop("No run directory found for method: ", method)
  }

  valid <- hits %>%
    mutate(
      data_dir = file.path(run_dir, "results", "data"),
      has_pheno = map_lgl(data_dir, ~ length(list.files(.x, pattern = "^pheno_.*\\.phenotypes\\.tsv$", full.names = TRUE)) > 0),
      has_map = map_lgl(data_dir, ~ length(list.files(.x, pattern = "^pheno_.*\\.gene_index_map\\.txt$", full.names = TRUE)) > 0)
    ) %>%
    filter(has_pheno, has_map)

  if (nrow(valid) == 0) {
    stop("Run directory exists but phenotype/map files are missing for method: ", method)
  }

  exact <- valid %>%
    filter(basename(run_dir) == paste0(method, run_suffix))
  if (nrow(exact) > 0) {
    return(exact$run_dir[[1]])
  }

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

  if (!all(c("gene_name", "mpheno_index") %in% names(map_dt))) {
    stop("Unexpected map format: ", map_file)
  }

  map_dt <- map_dt[order(as.integer(mpheno_index))]
  genes <- normalize_gene_id(map_dt$gene_name)

  if ((ncol(pheno) - 2L) != length(genes)) {
    stop("Phenotype/map dimensions mismatch for run: ", run_dir)
  }

  expr <- as.matrix(pheno[, -(1:2), with = FALSE])
  storage.mode(expr) <- "numeric"
  rownames(expr) <- paste0(as.character(pheno[[1]]), "::", as.character(pheno[[2]]))

  # Collapse duplicate gene IDs by averaging columns.
  if (anyDuplicated(genes) > 0L) {
    idx <- split(seq_along(genes), genes)
    expr <- do.call(
      cbind,
      lapply(idx, function(ii) {
        if (length(ii) == 1L) expr[, ii] else rowMeans(expr[, ii, drop = FALSE], na.rm = TRUE)
      })
    )
    colnames(expr) <- names(idx)
  } else {
    colnames(expr) <- genes
  }

  list(run_dir = run_dir, expr = expr)
}

skewness_third_moment <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 3) return(NA_real_)
  mu <- mean(x)
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  mean(((x - mu) / s)^3)
}

pick_quantile_genes <- function(skew_tbl, n_pick = 12L) {
  q_targets <- quantile(
    skew_tbl$skewness,
    probs = seq(0.05, 0.95, length.out = n_pick),
    na.rm = TRUE,
    names = FALSE
  )

  chosen <- character(0)
  out <- vector("list", length(q_targets))

  for (i in seq_along(q_targets)) {
    qv <- q_targets[i]
    cand <- skew_tbl %>%
      filter(!Gene %in% chosen) %>%
      mutate(dist_q = abs(skewness - qv)) %>%
      arrange(dist_q) %>%
      slice(1)
    chosen <- c(chosen, cand$Gene)
    out[[i]] <- cand
  }

  bind_rows(out) %>%
    mutate(rank_id = row_number())
}

# -----------------------------
# Compute skewness on random sample of genes
# -----------------------------
obj <- read_expression_matrix(method_id)
expr <- obj$expr
all_genes <- colnames(expr)

if (length(all_genes) < 50) {
  stop("Too few genes found in matrix.")
}

set.seed(rng_seed)
n_take <- min(n_sample_genes, length(all_genes))
sampled_genes <- sample(all_genes, size = n_take, replace = FALSE)
expr_sub <- expr[, sampled_genes, drop = FALSE]

skew_tbl <- tibble(
  Gene = sampled_genes,
  mean_expr = colMeans(expr_sub, na.rm = TRUE),
  median_expr = apply(expr_sub, 2, median, na.rm = TRUE),
  sd_expr = apply(expr_sub, 2, sd, na.rm = TRUE),
  skewness = apply(expr_sub, 2, skewness_third_moment)
) %>%
  filter(is.finite(skewness))

if (nrow(skew_tbl) == 0) {
  stop("No finite skewness values were computed.")
}

summary_tbl <- skew_tbl %>%
  summarise(
    method = method_id,
    run_dir = obj$run_dir,
    n_individuals = nrow(expr_sub),
    n_genes_requested = n_sample_genes,
    n_genes_used = nrow(skew_tbl),
    mean_skewness = mean(skewness, na.rm = TRUE),
    median_skewness = median(skewness, na.rm = TRUE),
    q25_skewness = quantile(skewness, 0.25, na.rm = TRUE),
    q75_skewness = quantile(skewness, 0.75, na.rm = TRUE),
    frac_right_skew = mean(skewness > 0, na.rm = TRUE)
  )

# -----------------------------
# Outputs
# -----------------------------
slug <- str_replace_all(method_id, "[^A-Za-z0-9_]+", "_")
out_gene <- file.path(runs_dir, paste0("tmm_raw_skewness_gene_table_", slug, ".tsv"))
out_summary <- file.path(runs_dir, paste0("tmm_raw_skewness_summary_", slug, ".tsv"))
out_hist <- file.path(runs_dir, paste0("tmm_raw_skewness_hist_density_", slug, ".png"))
out_scatter <- file.path(runs_dir, paste0("tmm_raw_skewness_mean_vs_skew_", slug, ".png"))
out_examples <- file.path(runs_dir, paste0("tmm_raw_skewness_example_gene_distributions_", slug, ".png"))

write_tsv(skew_tbl, out_gene)
write_tsv(summary_tbl, out_summary)

# 1) Skewness distribution
p_hist <- ggplot(skew_tbl, aes(x = skewness)) +
  geom_histogram(aes(y = after_stat(density)), bins = 70, fill = "#2A9D8F", color = "white", alpha = 0.8) +
  geom_density(color = "#264653", linewidth = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#E76F51", linewidth = 0.4) +
  labs(
    title = paste0("Gene-wise skewness (third moment) - ", method_id),
    subtitle = paste0("Random sample of ", nrow(skew_tbl), " genes; values across individuals"),
    x = "Skewness",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

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
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(out_scatter, p_scatter, width = 9, height = 5.5, dpi = 300)

# 3) Distribution across individuals for representative genes
example_genes <- pick_quantile_genes(skew_tbl, n_pick = n_example_genes)
expr_example <- as_tibble(expr_sub[, example_genes$Gene, drop = FALSE]) %>%
  mutate(sample_id = rownames(expr_sub)) %>%
  pivot_longer(cols = -sample_id, names_to = "Gene", values_to = "expr") %>%
  left_join(example_genes %>% select(Gene, skewness, mean_expr, rank_id), by = "Gene") %>%
  mutate(
    facet_label = paste0(
      "Q", rank_id,
      "\n", Gene,
      "\nskew=", sprintf("%.2f", skewness)
    ),
    facet_label = factor(facet_label, levels = unique(facet_label[order(rank_id)]))
  )

p_examples <- ggplot(expr_example, aes(x = expr)) +
  geom_histogram(bins = 40, fill = "#4DB6AC", color = "white", alpha = 0.9) +
  facet_wrap(~ facet_label, scales = "free_y", ncol = 4) +
  labs(
    title = paste0("Expression distribution across individuals - ", method_id),
    subtitle = "Representative genes spanning skewness quantiles",
    x = "Expression value (post-PEER phenotype values)",
    y = "Number of individuals"
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 8)
  )

ggsave(out_examples, p_examples, width = 13, height = 8, dpi = 300)

cat("Saved:\n")
cat("- ", out_gene, "\n", sep = "")
cat("- ", out_summary, "\n", sep = "")
cat("- ", out_hist, "\n", sep = "")
cat("- ", out_scatter, "\n", sep = "")
cat("- ", out_examples, "\n", sep = "")

