suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# -----------------------------
# Settings
# -----------------------------
runs_dir <- "/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260428_GE_GEUVADIS_v2/GeneExpression/runs"
analysis_root <- file.path(runs_dir, "_analysis", "pairwise_h2_comparisons")
plots_dir <- file.path(analysis_root, "plots")
tables_dir <- file.path(analysis_root, "tables")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# 1-vs-all mode: set method id here, e.g. "all_snps_tmm_raw_peerauto_pmg0_npc5"
# all-pairs mode: keep NULL
focus_method <- NULL
include_mixed <- FALSE
# fallback strategy when PASS-only intersections are empty:
# "auto" -> try PASS finite first, then fall back to any finite h2 (including non-PASS)
# "strict" -> require PASS finite intersections
# "any_finite" -> always use any finite h2
fallback_strategy <- "auto"
log_expr_caps <- c(1.5, 5, 8)

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

normalize_gene_id <- function(x) {
  x_norm <- x %>%
    as.character() %>%
    str_trim() %>%
    str_replace("\\.[0-9]+$", "")
  ensg <- str_extract(x_norm, "ENSG[0-9]+")
  if_else(!is.na(ensg), ensg, x_norm)
}

resolve_gene_ids_for_run <- function(gene_raw, run_name) {
  gene_chr <- as.character(gene_raw)
  idx_mask <- str_detect(str_trim(gene_chr), "^[0-9]+$")
  if (mean(idx_mask, na.rm = TRUE) < 0.8) {
    return(gene_chr)
  }

  data_dir <- file.path(runs_dir, run_name, "results", "data")
  map_file <- list.files(data_dir, pattern = "^pheno_.*\\.gene_index_map\\.txt$", full.names = TRUE)
  if (length(map_file) == 0) {
    warning("Gene IDs look numeric but no gene_index_map found for run: ", run_name)
    return(gene_chr)
  }

  map_dt <- fread(map_file[1], sep = "\t", data.table = TRUE)
  if (!all(c("gene_name", "mpheno_index") %in% names(map_dt))) {
    warning("Unexpected gene_index_map columns for run: ", run_name)
    return(gene_chr)
  }

  map_tbl <- map_dt %>%
    transmute(
      mpheno_index = as.integer(mpheno_index),
      gene_name = as.character(gene_name)
    )

  idx_tbl <- tibble(
    row_id = seq_along(gene_chr),
    mpheno_index = suppressWarnings(as.integer(gene_chr))
  ) %>%
    left_join(map_tbl, by = "mpheno_index")

  resolved <- gene_chr
  replace_mask <- idx_mask & !is.na(idx_tbl$gene_name)
  resolved[replace_mask] <- idx_tbl$gene_name[replace_mask]
  resolved
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

read_expression_phenotype_for_method <- function(method_id) {
  run_name <- method_to_run_name(method_id)
  data_dir <- file.path(runs_dir, run_name, "results", "data")
  pheno_file <- list.files(data_dir, pattern = "^pheno_.*\\.phenotypes\\.tsv$", full.names = TRUE)
  map_file <- list.files(data_dir, pattern = "^pheno_.*\\.gene_index_map\\.txt$", full.names = TRUE)

  if (length(pheno_file) == 0 || length(map_file) == 0) {
    warning("Missing phenotype/map file for method: ", method_id, " (run=", run_name, ")")
    return(NULL)
  }

  pheno <- fread(pheno_file[1], header = FALSE, sep = "\t", data.table = TRUE)
  map_dt <- fread(map_file[1], sep = "\t", data.table = TRUE)
  if (!all(c("gene_name", "mpheno_index") %in% names(map_dt))) {
    warning("Unexpected gene index map columns for method: ", method_id)
    return(NULL)
  }

  map_dt <- map_dt[order(as.integer(mpheno_index))]
  genes <- normalize_gene_id(as.character(map_dt$gene_name))
  if (ncol(pheno) < 3 || (ncol(pheno) - 2L) != length(genes)) {
    warning("Phenotype/map dimension mismatch for method: ", method_id)
    return(NULL)
  }

  expr_mat <- as.matrix(pheno[, -(1:2), with = FALSE])
  storage.mode(expr_mat) <- "numeric"
  sample_id <- paste0(as.character(pheno[[1]]), "::", as.character(pheno[[2]]))
  rownames(expr_mat) <- sample_id

  # Collapse duplicate gene IDs (rare) by averaging duplicate columns.
  if (anyDuplicated(genes) > 0L) {
    idx <- split(seq_along(genes), genes)
    expr_mat <- do.call(cbind, lapply(idx, function(ii) {
      if (length(ii) == 1L) {
        expr_mat[, ii]
      } else {
        rowMeans(expr_mat[, ii, drop = FALSE], na.rm = TRUE)
      }
    }))
    colnames(expr_mat) <- names(idx)
  } else {
    colnames(expr_mat) <- genes
  }

  expr_mat <- as.matrix(expr_mat)
  rownames(expr_mat) <- sample_id

  list(
    method = method_id,
    run_name = run_name,
    expr = expr_mat
  )
}

compute_tpm_tmm_expression_correlations <- function() {
  combos <- tidyr::crossing(
    snp_set = c("all_snps", "hm3_no_mhc"),
    norm = c("raw", "irnt")
  )

  map_dfr(seq_len(nrow(combos)), function(i) {
    snp <- combos$snp_set[i]
    norm <- combos$norm[i]
    tmm_method <- paste(snp, "tmm", norm, sep = "_")
    tpm_method <- paste(snp, "tpm", norm, sep = "_")

    tmm_obj <- read_expression_phenotype_for_method(tmm_method)
    tpm_obj <- read_expression_phenotype_for_method(tpm_method)
    if (is.null(tmm_obj) || is.null(tpm_obj)) return(tibble())

    common_samples <- intersect(rownames(tmm_obj$expr), rownames(tpm_obj$expr))
    common_genes <- intersect(colnames(tmm_obj$expr), colnames(tpm_obj$expr))
    if (length(common_samples) < 3L || length(common_genes) == 0L) return(tibble())

    tmm_mat <- tmm_obj$expr[common_samples, common_genes, drop = FALSE]
    tpm_mat <- tpm_obj$expr[common_samples, common_genes, drop = FALSE]

    gene_cor <- vapply(seq_along(common_genes), function(j) {
      cor(tmm_mat[, j], tpm_mat[, j], use = "pairwise.complete.obs")
    }, numeric(1))

    tibble(
      Gene = common_genes,
      snp_set = snp,
      norm = norm,
      cor_expr_tpm_vs_tmm = gene_cor,
      mean_tmm = colMeans(tmm_mat, na.rm = TRUE),
      mean_tpm = colMeans(tpm_mat, na.rm = TRUE),
      n_samples = length(common_samples),
      n_genes_overlap = length(common_genes),
      panel = paste(if_else(snp == "all_snps", "ALL", "HM3"), toupper(norm), sep = " | ")
    )
  })
}

summarise_h2_long <- function(df) {
  df %>%
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
}

build_overlap_matrix <- function(source_tbl, method_vec, value_col = "h2") {
  out <- matrix(0L, nrow = length(method_vec), ncol = length(method_vec), dimnames = list(method_vec, method_vec))
  if (length(method_vec) == 0) return(out)

  gene_sets <- setNames(
    lapply(method_vec, function(m) {
      source_tbl %>%
        filter(method == m, is.finite(.data[[value_col]])) %>%
        pull(Gene) %>%
        unique()
    }),
    method_vec
  )

  for (i in seq_along(method_vec)) {
    for (j in seq_along(method_vec)) {
      out[i, j] <- length(intersect(gene_sets[[i]], gene_sets[[j]]))
    }
  }
  out
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

h2_raw <- map_dfr(files, function(f) {
  run_name <- str_match(f, "runs/([^/]+)/results/summary/")[, 2]
  meta <- parse_method(run_name)

  if (any(is.na(meta$snp_set))) return(tibble())

  sum_tbl <- read_tsv(f, show_col_types = FALSE)
  gene_resolved <- resolve_gene_ids_for_run(sum_tbl$Gene, run_name)

  sum_tbl %>%
    transmute(
      Gene_raw = Gene,
      Gene = normalize_gene_id(gene_resolved),
      Status = toupper(str_trim(Status)),
      h2 = suppressWarnings(as.numeric(h2_GREML)),
      se = suppressWarnings(as.numeric(SE_GREML)),
      pval = suppressWarnings(as.numeric(Pval_GREML)),
      run_name = run_name,
      method = meta$method
    )
})

if (nrow(h2_raw) == 0) {
  stop("No parseable summary rows were found. Check run directory naming and parse_method regex.")
}

method_run_counts <- h2_raw %>%
  distinct(method, run_name) %>%
  count(method, name = "n_runs")

if (any(method_run_counts$n_runs > 1L)) {
  warning(
    "Multiple runs detected for at least one method. ",
    "This can create disjoint gene universes if runs came from different cohorts/sweeps.\n",
    paste0(
      apply(
        as.data.frame(method_run_counts %>% filter(n_runs > 1L)),
        1,
        function(r) paste0("  - ", r[["method"]], ": ", r[["n_runs"]], " runs")
      ),
      collapse = "\n"
    )
  )
}

method_status_summary <- h2_raw %>%
  group_by(method) %>%
  summarise(
    n_rows = n(),
    n_genes_total = n_distinct(Gene),
    n_pass_rows = sum(Status == "PASS"),
    n_pass_finite_h2 = sum(Status == "PASS" & is.finite(h2)),
    n_nonpass_finite_h2 = sum(Status != "PASS" & is.finite(h2)),
    n_any_finite_h2 = sum(is.finite(h2)),
    n_fail_rows = sum(Status != "PASS"),
    pct_ensg_like = mean(str_detect(Gene, "^ENSG[0-9]+$"), na.rm = TRUE),
    .groups = "drop"
  )

cat("Method-level status summary:\n")
print(method_status_summary)

h2_long_pass <- h2_raw %>%
  filter(Status == "PASS", is.finite(h2)) %>%
  mutate(
    z = if_else(is.finite(se) & se > 0, h2 / se, NA_real_)
  ) %>%
  summarise_h2_long()

h2_long_any <- h2_raw %>%
  filter(is.finite(h2)) %>%
  mutate(
    z = if_else(is.finite(se) & se > 0, h2 / se, NA_real_)
  ) %>%
  summarise_h2_long()

methods <- h2_raw %>% distinct(method) %>% pull(method) %>% sort()
if (length(methods) < 2) stop("Need at least 2 methods for pairwise comparison.")

if (!is.null(focus_method)) {
  focus_method <- normalize_focus(focus_method)
  if (!focus_method %in% methods) stop("focus_method not found after normalization: ", focus_method)
}
fallback_strategy <- match.arg(fallback_strategy, c("auto", "strict", "any_finite"))

build_overlap_table <- function(source_tbl, method_vec, value_col = "h2") {
  combn(method_vec, 2, simplify = FALSE) %>%
    map_dfr(function(p) {
      m1 <- p[1]
      m2 <- p[2]
      g1 <- source_tbl %>%
        filter(method == m1, is.finite(.data[[value_col]])) %>%
        pull(Gene) %>%
        unique()
      g2 <- source_tbl %>%
        filter(method == m2, is.finite(.data[[value_col]])) %>%
        pull(Gene) %>%
        unique()
      tibble(
        m1 = m1,
        m2 = m2,
        n_genes_m1 = length(g1),
        n_genes_m2 = length(g2),
        n_overlap = length(intersect(g1, g2))
      )
    })
}

build_pair_df_from_long <- function(h2_long_tbl, method_vec) {
  h2_wide_local <- h2_long_tbl %>%
    select(Gene, method, h2, z, qval) %>%
    pivot_wider(
      id_cols = Gene,
      names_from = method,
      values_from = c(h2, z, qval),
      names_sep = "__"
    )

  combn(method_vec, 2, simplify = FALSE) %>%
    map_dfr(function(p) {
      h2_wide_local %>%
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
}

prefix <- ifelse(is.null(focus_method), "all_pairs", paste0("focus_", focus_method))
method_status_file <- file.path(tables_dir, paste0("pairwise_method_status_", prefix, ".tsv"))
write_tsv(method_status_summary, method_status_file)

overlap_all_tbl <- build_overlap_table(h2_raw %>% filter(is.finite(h2)), methods)
overlap_pass_tbl <- build_overlap_table(h2_long_pass, methods)

overlap_file <- file.path(
  tables_dir,
  paste0("pairwise_gene_overlap_diagnostics_", prefix, ".tsv")
)
overlap_diag <- overlap_pass_tbl %>%
  left_join(
    overlap_all_tbl %>% select(m1, m2, n_overlap_all = n_overlap),
    by = c("m1", "m2")
  ) %>%
  mutate(
    n_overlap_fail_only = pmax(n_overlap_all - n_overlap, 0L)
  )
write_tsv(overlap_diag, overlap_file)

overlap_mat_pass <- build_overlap_matrix(h2_long_pass, methods, value_col = "h2")
overlap_mat_any <- build_overlap_matrix(h2_raw %>% filter(is.finite(h2)), methods, value_col = "h2")
overlap_matrix_file <- file.path(tables_dir, paste0("pairwise_overlap_matrix_", prefix, ".tsv"))
bind_rows(
  as_tibble(as.data.frame(as.table(overlap_mat_pass))) %>%
    transmute(mode = "pass_finite", method_x = Var1, method_y = Var2, n_overlap = as.integer(Freq)),
  as_tibble(as.data.frame(as.table(overlap_mat_any))) %>%
    transmute(mode = "any_finite", method_x = Var1, method_y = Var2, n_overlap = as.integer(Freq))
) %>%
  write_tsv(overlap_matrix_file)
cat("Pairwise overlap matrix (PASS + finite h2):\n")
print(overlap_mat_pass, quote = FALSE)
cat("Pairwise overlap matrix (any finite h2; PASS + non-PASS):\n")
print(overlap_mat_any, quote = FALSE)

data_mode <- if (identical(fallback_strategy, "any_finite")) "any_finite" else "pass_finite"
h2_long_active <- if (identical(data_mode, "pass_finite")) h2_long_pass else h2_long_any

if (nrow(h2_long_active) == 0 &&
    identical(data_mode, "pass_finite") &&
    identical(fallback_strategy, "auto")) {
  message("No PASS-finite rows available; falling back to any finite h2 across statuses.")
  data_mode <- "any_finite"
  h2_long_active <- h2_long_any
}

if (nrow(h2_long_active) == 0) {
  stop(
    "No finite h2 rows available for the selected fallback strategy. ",
    "See method-specific report: ", method_status_file
  )
}

methods_active <- h2_long_active %>% distinct(method) %>% pull(method) %>% sort()
if (length(methods_active) < 2) {
  stop("Need at least 2 methods with finite h2 in active mode: ", data_mode)
}

# -----------------------------
# 3) Log-expression map from RAW baselines
# -----------------------------
raw_methods <- methods_active %>% map_chr(method_to_raw_baseline) %>% unique()
log_expr_map <- map_dfr(raw_methods, function(m_raw) {
  run_name <- method_to_run_name(m_raw)
  read_log_expr_for_run(run_name) %>%
    mutate(raw_method = m_raw)
})

# -----------------------------
# 4) Build pairwise table
# -----------------------------
pair_df <- build_pair_df_from_long(h2_long_active, methods_active)

if (!is.null(focus_method)) {
  pair_df <- pair_df %>% filter(m1 == focus_method | m2 == focus_method)
  overlap_diag <- overlap_diag %>% filter(m1 == focus_method | m2 == focus_method)
}

cat("Methods detected in active mode:\n")
print(methods_active)
cat("Pairwise overlap diagnostics (PASS vs all finite h2):\n")
print(overlap_diag)
cat("Section counts before mixed-filter:\n")
print(table(pair_df$section, useNA = "ifany"))

if (!include_mixed) {
  pair_df <- pair_df %>% filter(section != "Mixed changes")
}

if (nrow(pair_df) == 0 &&
    identical(data_mode, "pass_finite") &&
    identical(fallback_strategy, "auto")) {
  message("No PASS pair rows after filters; retrying with any finite h2 across statuses.")
  data_mode <- "any_finite"
  h2_long_active <- h2_long_any
  methods_active <- h2_long_active %>% distinct(method) %>% pull(method) %>% sort()
  raw_methods <- methods_active %>% map_chr(method_to_raw_baseline) %>% unique()
  log_expr_map <- map_dfr(raw_methods, function(m_raw) {
    run_name <- method_to_run_name(m_raw)
    read_log_expr_for_run(run_name) %>%
      mutate(raw_method = m_raw)
  })
  pair_df <- build_pair_df_from_long(h2_long_active, methods_active)
  if (!is.null(focus_method)) {
    pair_df <- pair_df %>% filter(m1 == focus_method | m2 == focus_method)
  }
  if (!include_mixed) {
    pair_df <- pair_df %>% filter(section != "Mixed changes")
  }
  cat("Section counts after fallback to any_finite:\n")
  print(table(pair_df$section, useNA = "ifany"))
}

if (nrow(pair_df) == 0) {
  stop(
    "No rows left after filtering (mode=", data_mode, "). ",
    "See method-specific report: ", method_status_file, " and overlap diagnostics: ", overlap_file,
    ". Set fallback_strategy='any_finite' to force non-PASS finite fallback."
  )
}

plot_subtitle <- if (identical(data_mode, "pass_finite")) {
  "PASS genes only"
} else {
  "Finite h2 genes (PASS + non-PASS fallback)"
}

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

plot_section <- function(df, title, subtitle_text, x_lab, y_lab, out_file, color_var, color_label, color_mode = c("seq", "div"), color_cap = NA_real_) {
  if (nrow(df) == 0) return(invisible(NULL))
  color_mode <- match.arg(color_mode)
  stats <- calc_stats(df)
  plot_df <- df
  color_col <- color_var
  color_label_use <- color_label

  if (color_mode == "seq" && is.finite(color_cap)) {
    plot_df <- plot_df %>%
      mutate(.color_capped = pmin(pmax(.data[[color_var]], 0), color_cap))
    color_col <- ".color_capped"
    color_label_use <- paste0(color_label, "\n(capped at ", color_cap, ")")
  }

  p <- plot_df %>%
    ggplot(aes(x = x_plot, y = y_plot, color = .data[[color_col]])) +
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
      subtitle = subtitle_text,
      x = x_lab,
      y = y_lab,
      color = color_label_use
    ) +
    theme_minimal(base_size = 10.5) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 8.5)
    )

  if (color_mode == "seq") {
    if (is.finite(color_cap)) {
      p <- p + scale_color_viridis_c(
        option = "magma",
        limits = c(0, color_cap),
        oob = scales::squish,
        na.value = "grey80"
      )
    } else {
      p <- p + scale_color_viridis_c(option = "magma", na.value = "grey80")
    }
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

plot_tpm_tmm_correlation <- function(df, subtitle_text, out_file) {
  if (nrow(df) == 0) return(invisible(NULL))

  corr_stats <- df %>%
    group_by(facet_label) %>%
    summarise(
      n = n(),
      r = cor(x_plot, y_plot, use = "complete.obs"),
      .groups = "drop"
    ) %>%
    mutate(lbl = sprintf("n=%d, r=%.2f", n, r))

  p <- df %>%
    ggplot(aes(x = x_plot, y = y_plot)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E76F51", linewidth = 0.35) +
    geom_point(alpha = 0.35, color = "#2A9D8F", size = 0.8) +
    geom_smooth(method = "lm", se = FALSE, color = "#264653", linewidth = 0.55) +
    geom_text(
      data = corr_stats,
      aes(x = -Inf, y = Inf, label = lbl),
      inherit.aes = FALSE,
      hjust = -0.1,
      vjust = 1.1,
      size = 3.0,
      color = "#264653"
    ) +
    facet_wrap(~facet_label, scales = "free", ncol = 2) +
    labs(
      title = "TPM vs TMM Correlation (by normalization)",
      subtitle = subtitle_text,
      x = "X-axis h2 (TMM)",
      y = "Y-axis h2 (TPM)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 9)
    )

  ggsave(out_file, p, width = 10, height = 7, dpi = 300)
  invisible(corr_stats)
}

out_stats <- file.path(tables_dir, paste0("pairwise_h2_stats_", prefix, ".tsv"))
out_discrepancy_types <- file.path(tables_dir, paste0("pairwise_h2_discrepancy_types_", prefix, ".tsv"))

stats_all <- bind_rows(
  calc_stats(snp_df) %>% mutate(section = "SNP set: ALL vs HM3"),
  calc_stats(expr_df) %>% mutate(section = "Expression: TPM vs TMM"),
  calc_stats(norm_df) %>% mutate(section = "Normalization: RAW vs IRNT")
) %>%
  mutate(data_mode = data_mode)
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
  subtitle_text = plot_subtitle,
  x_lab = "X-axis h2 (ALL)",
  y_lab = "Y-axis h2 (HM3)",
  out_file = file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_snp_set_all_vs_hm3_color_logexpr.png")),
  color_var = "log_expr_plot",
  color_label = "log2(mean expr + 1)",
  color_mode = "seq"
)

for (cap in log_expr_caps) {
  plot_section(
    expr_df,
    title = "Pairwise GREML h2: Expression effect (TPM vs TMM)",
    subtitle_text = paste0(plot_subtitle, " | color cap=", cap),
    x_lab = "X-axis h2 (TMM)",
    y_lab = "Y-axis h2 (TPM)",
    out_file = file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_expression_tpm_vs_tmm_color_logexpr_cap", cap, ".png")),
    color_var = "log_expr_plot",
    color_label = "log2(mean expr + 1)",
    color_mode = "seq",
    color_cap = cap
  )

  plot_section(
    norm_df,
    title = "Pairwise GREML h2: Normalization effect (RAW vs IRNT)",
    subtitle_text = paste0(plot_subtitle, " | color cap=", cap),
    x_lab = "X-axis h2 (RAW)",
    y_lab = "Y-axis h2 (IRNT)",
    out_file = file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_normalization_raw_vs_irnt_color_logexpr_cap", cap, ".png")),
    color_var = "log_expr_plot",
    color_label = "log2(mean expr + 1)",
    color_mode = "seq",
    color_cap = cap
  )

  # Extra expression-only views with explicit color source
  # X-axis is TMM h2 and Y-axis is TPM h2, so:
  # - log_expr_x corresponds to TMM RAW expression
  # - log_expr_y corresponds to TPM RAW expression
  plot_section(
    expr_df,
    title = "Pairwise GREML h2: Expression effect (TPM vs TMM)",
    subtitle_text = paste0(plot_subtitle, " | color = TMM RAW log2(mean expr + 1), cap=", cap),
    x_lab = "X-axis h2 (TMM)",
    y_lab = "Y-axis h2 (TPM)",
    out_file = file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_expression_tpm_vs_tmm_color_logexpr_tmm_cap", cap, ".png")),
    color_var = "log_expr_x",
    color_label = "log2(mean expr + 1) | TMM RAW",
    color_mode = "seq",
    color_cap = cap
  )

  plot_section(
    expr_df,
    title = "Pairwise GREML h2: Expression effect (TPM vs TMM)",
    subtitle_text = paste0(plot_subtitle, " | color = TPM RAW log2(mean expr + 1), cap=", cap),
    x_lab = "X-axis h2 (TMM)",
    y_lab = "Y-axis h2 (TPM)",
    out_file = file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_expression_tpm_vs_tmm_color_logexpr_tpm_cap", cap, ".png")),
    color_var = "log_expr_y",
    color_label = "log2(mean expr + 1) | TPM RAW",
    color_mode = "seq",
    color_cap = cap
  )
}

# 6B) Color by heritability significance proxy (Wald Z)
plot_section(
  snp_df,
  title = "Pairwise GREML h2: SNP effect (ALL vs HM3)",
  subtitle_text = plot_subtitle,
  x_lab = "X-axis h2 (ALL)",
  y_lab = "Y-axis h2 (HM3)",
  out_file = file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_snp_set_all_vs_hm3_color_zscore.png")),
  color_var = "z_plot",
  color_label = "mean Wald Z (X,Y)",
  color_mode = "div"
)

plot_section(
  expr_df,
  title = "Pairwise GREML h2: Expression effect (TPM vs TMM)",
  subtitle_text = plot_subtitle,
  x_lab = "X-axis h2 (TMM)",
  y_lab = "Y-axis h2 (TPM)",
  out_file = file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_expression_tpm_vs_tmm_color_zscore.png")),
  color_var = "z_plot",
  color_label = "mean Wald Z (X,Y)",
  color_mode = "div"
)

plot_section(
  norm_df,
  title = "Pairwise GREML h2: Normalization effect (RAW vs IRNT)",
  subtitle_text = plot_subtitle,
  x_lab = "X-axis h2 (RAW)",
  y_lab = "Y-axis h2 (IRNT)",
  out_file = file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_normalization_raw_vs_irnt_color_zscore.png")),
  color_var = "z_plot",
  color_label = "mean Wald Z (X,Y)",
  color_mode = "div"
)

# Extra: dedicated TPM vs TMM correlation figure, split by RAW vs IRNT and ALL vs HM3
tpm_tmm_corr_file <- file.path(plots_dir, paste0("pairwise_h2_tpm_vs_tmm_correlation_", prefix, ".png"))
plot_tpm_tmm_correlation(
  expr_df,
  subtitle_text = plot_subtitle,
  out_file = tpm_tmm_corr_file
)

cat("Saved grouped plots and stats:\n")
cat("- ", file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_snp_set_all_vs_hm3_color_logexpr.png")), "\n", sep = "")
cat("- ", file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_expression_tpm_vs_tmm_color_logexpr_cap1.5.png")), "\n", sep = "")
cat("- ", file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_expression_tpm_vs_tmm_color_logexpr_cap5.png")), "\n", sep = "")
cat("- ", file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_expression_tpm_vs_tmm_color_logexpr_cap8.png")), "\n", sep = "")
cat("- ", file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_expression_tpm_vs_tmm_color_logexpr_tmm_cap1.5.png")), "\n", sep = "")
cat("- ", file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_expression_tpm_vs_tmm_color_logexpr_tmm_cap5.png")), "\n", sep = "")
cat("- ", file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_expression_tpm_vs_tmm_color_logexpr_tmm_cap8.png")), "\n", sep = "")
cat("- ", file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_expression_tpm_vs_tmm_color_logexpr_tpm_cap1.5.png")), "\n", sep = "")
cat("- ", file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_expression_tpm_vs_tmm_color_logexpr_tpm_cap5.png")), "\n", sep = "")
cat("- ", file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_expression_tpm_vs_tmm_color_logexpr_tpm_cap8.png")), "\n", sep = "")
cat("- ", file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_normalization_raw_vs_irnt_color_logexpr_cap1.5.png")), "\n", sep = "")
cat("- ", file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_normalization_raw_vs_irnt_color_logexpr_cap5.png")), "\n", sep = "")
cat("- ", file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_normalization_raw_vs_irnt_color_logexpr_cap8.png")), "\n", sep = "")
cat("- ", file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_snp_set_all_vs_hm3_color_zscore.png")), "\n", sep = "")
cat("- ", file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_expression_tpm_vs_tmm_color_zscore.png")), "\n", sep = "")
cat("- ", file.path(plots_dir, paste0("pairwise_h2_scatter_", prefix, "_normalization_raw_vs_irnt_color_zscore.png")), "\n", sep = "")
cat("- ", out_stats, "\n", sep = "")
cat("- ", out_discrepancy_types, "\n", sep = "")
cat("- ", method_status_file, "\n", sep = "")
cat("- ", overlap_file, "\n", sep = "")
cat("- ", overlap_matrix_file, "\n", sep = "")
cat("- ", tpm_tmm_corr_file, "\n", sep = "")
cat("Data mode used: ", data_mode, "\n", sep = "")
