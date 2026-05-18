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

# Comma-separated method IDs, e.g.:
# METHOD_IDS="all_snps_tmm_raw,hm3_no_mhc_tmm_raw"
method_ids_env <- Sys.getenv("METHOD_IDS", "all_snps_tmm_raw,hm3_no_mhc_tmm_raw")
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
      summary_dir = file.path(run_dir, "results", "summary"),
      has_pheno = map_lgl(data_dir, ~ length(list.files(.x, pattern = "^pheno_.*\\.phenotypes\\.tsv$", full.names = TRUE)) > 0),
      has_map = map_lgl(data_dir, ~ length(list.files(.x, pattern = "^pheno_.*\\.gene_index_map\\.txt$", full.names = TRUE)) > 0),
      has_summary = map_lgl(summary_dir, ~ length(list.files(.x, pattern = "^final_heritability_summary.*\\.tsv$", full.names = TRUE)) > 0)
    ) %>%
    filter(has_pheno, has_map, has_summary)

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

# --------------------------------------------------
# Build analysis table
# --------------------------------------------------
analysis_tbl <- map_dfr(method_ids, function(mid) {
  run_dir <- find_run_dir(mid)
  expr_tbl <- read_expression_by_gene(run_dir)
  h2_tbl <- read_h2_summary(run_dir)

  merged <- h2_tbl %>%
    inner_join(expr_tbl, by = "Gene") %>%
    filter(is.finite(mean_tmm_expression), is.finite(h2_GREML))

  if (nrow(merged) < 50) {
    stop("Too few genes after merge for method ", mid, " (n=", nrow(merged), ").")
  }

  merged %>%
    mutate(
      method_id = mid,
      run_dir = run_dir,
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

summary_tbl <- analysis_tbl %>%
  group_by(method_id, method_label, expr_decile) %>%
  summarise(
    n_genes = n(),
    mean_h2 = mean(h2_GREML, na.rm = TRUE),
    median_h2 = median(h2_GREML, na.rm = TRUE),
    q25_h2 = quantile(h2_GREML, 0.25, na.rm = TRUE),
    q75_h2 = quantile(h2_GREML, 0.75, na.rm = TRUE),
    mean_tmm_expression = mean(mean_tmm_expression, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(analysis_tbl, out_gene_tsv)
write_tsv(summary_tbl, out_summary_tsv)

# --------------------------------------------------
# Plot: h2 distribution by mean-expression decile
# --------------------------------------------------
y_min <- min(analysis_tbl$h2_GREML, na.rm = TRUE)
y_max <- max(analysis_tbl$h2_GREML, na.rm = TRUE)
y_pad <- max((y_max - y_min) * 0.05, 0.02)

p <- ggplot(analysis_tbl, aes(x = expr_decile, y = h2_GREML, fill = expr_decile)) +
  geom_boxplot(width = 0.68, outlier.alpha = 0.15, size = 0.28) +
  geom_point(
    data = summary_tbl,
    aes(x = expr_decile, y = mean_h2),
    inherit.aes = FALSE,
    shape = 23,
    fill = "white",
    color = "black",
    size = 1.9,
    stroke = 0.35
  ) +
  geom_text(
    data = summary_tbl,
    aes(x = expr_decile, y = mean_h2 + y_pad, label = sprintf("mean=%.3f", mean_h2)),
    inherit.aes = FALSE,
    size = 2.7,
    vjust = 0
  ) +
  scale_fill_brewer(palette = "GnBu") +
  facet_wrap(~ method_label, ncol = 2, scales = "fixed") +
  coord_cartesian(ylim = c(y_min, y_max + y_pad * 2.4), clip = "off") +
  labs(
    title = "GREML h2 by mean TMM expression decile",
    subtitle = "Genes are binned into 10 deciles by per-gene mean TMM expression (PASS genes only)",
    x = "Mean TMM expression decile (D1 = lowest, D10 = highest)",
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

ggsave(out_plot, p, width = 13, height = 7.5, dpi = 300)

cat("Saved:\n")
cat("- ", out_gene_tsv, "\n", sep = "")
cat("- ", out_summary_tsv, "\n", sep = "")
cat("- ", out_plot, "\n", sep = "")
