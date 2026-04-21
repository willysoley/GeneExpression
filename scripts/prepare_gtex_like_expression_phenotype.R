#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

if (!requireNamespace("edgeR", quietly = TRUE)) {
  stop(
    "Package 'edgeR' is required for TMM normalization.\n",
    "Install with: BiocManager::install('edgeR')"
  )
}

# ------------------------------ USER CONFIG -----------------------------------
# Update these paths before running.
cfg <- list(
  tpm_tsv = "/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260125_salmon_nextflow/results/matrices/gene_tpm.tsv",
  counts_tsv = "/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260125_salmon_nextflow/results/matrices/gene_counts.tsv",
  output_dir = "results/gtex_like_expression_phenotype",
  sdrf_url = "https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt",
  eur_pops = c("British", "Finnish", "Tuscan", "Utah"),
  expression_filter_mode = "gtex_tpm_only", # options: gtex_tpm_only, gtex_tpm_and_counts
  gtex_tpm_threshold = 0.1,
  gtex_count_threshold = 6L,
  gtex_sample_frac_threshold = 0.2,
  sample_id_map_tsv = NULL # optional 2-column TSV: sample_id, participant_id
)


# ------------------------------- HELPERS --------------------------------------
strip_gene_version <- function(x) {
  x %>%
    as.character() %>%
    str_remove("\\.[0-9]+$")
}

safe_message <- function(...) {
  message(sprintf(...))
}

stop_if_missing <- function(path, label) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop("Missing ", label, ": ", path)
  }
}

sanity_check <- function(condition, pass_msg, fail_msg) {
  if (isTRUE(condition)) {
    message("[Sanity] ", pass_msg)
  } else {
    stop("[Sanity failed] ", fail_msg)
  }
}

detect_gene_col <- function(dt) {
  candidates <- c("gene_id", "Gene", "feature_id", "gene")
  col <- candidates[candidates %in% names(dt)][1]
  if (is.na(col) || !nzchar(col)) {
    stop("Could not identify gene ID column. Tried: ", paste(candidates, collapse = ", "))
  }
  col
}

select_and_rename_samples <- function(dt, gene_col, keep_run_ids, label) {
  sample_cols <- setdiff(names(dt), gene_col)
  sample_map <- tibble(
    raw_col = sample_cols,
    run_id = str_remove(sample_cols, "_\\d+$")
  ) %>%
    filter(run_id %in% keep_run_ids)

  if (nrow(sample_map) == 0L) {
    stop("No matching sample columns found in ", label, " after ENA_RUN cleanup.")
  }

  dup_tbl <- sample_map %>%
    count(run_id, name = "n") %>%
    filter(n > 1L)

  if (nrow(dup_tbl) > 0L) {
    safe_message(
      "%s: dropping duplicate columns after suffix-cleanup for %d run IDs (keeping first).",
      label,
      nrow(dup_tbl)
    )
  }

  sample_map <- sample_map %>%
    group_by(run_id) %>%
    slice_head(n = 1L) %>%
    ungroup()

  out <- dt %>%
    as_tibble() %>%
    select(all_of(c(gene_col, sample_map$raw_col))) %>%
    rename_with(
      .fn = ~ sample_map$run_id[match(.x, sample_map$raw_col)],
      .cols = all_of(sample_map$raw_col)
    ) %>%
    mutate(gene_id_clean = strip_gene_version(.data[[gene_col]])) %>%
    select(gene_id_clean, everything(), -all_of(gene_col))

  out
}

collapse_duplicate_genes <- function(tbl, sample_cols, label, collapse_fn) {
  dup_n <- tbl %>%
    count(gene_id_clean, name = "n") %>%
    filter(n > 1L) %>%
    nrow()

  if (dup_n > 0L) {
    safe_message(
      "%s: collapsing %d duplicate gene IDs after version stripping.",
      label,
      dup_n
    )
  }

  tbl %>%
    group_by(gene_id_clean) %>%
    summarise(
      across(all_of(sample_cols), collapse_fn),
      .groups = "drop"
    )
}

inverse_normal_transform_vec <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)

  if (sum(ok) <= 1L) {
    out <- rep(NA_real_, length(x))
    if (sum(ok) == 1L) {
      out[ok] <- 0
    }
    return(out)
  }

  ranks <- rank(x[ok], ties.method = "average")
  n_ok <- sum(ok)
  transformed <- qnorm((ranks - 0.5) / n_ok)

  out <- rep(NA_real_, length(x))
  out[ok] <- transformed
  out
}

write_matrix_tsv <- function(mat, out_tsv_gz) {
  out_dt <- as.data.table(mat, keep.rownames = "gene_id_clean")
  fwrite(out_dt, out_tsv_gz, sep = "\t")
}


# ------------------------------- WORKFLOW -------------------------------------
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing(cfg$tpm_tsv, "TPM matrix")
stop_if_missing(cfg$counts_tsv, "Counts matrix")

safe_message("Step 1: loading SDRF and selecting European GEUVADIS runs")
sdrf <- fread(cfg$sdrf_url) %>%
  {
    names(.) <- make.unique(names(.))
    .
  } %>%
  as_tibble()

ancestry_col <- "Characteristics[ancestry category]"
run_col <- "Comment[ENA_RUN]"

if (!ancestry_col %in% names(sdrf) || !run_col %in% names(sdrf)) {
  stop("SDRF is missing required columns: ", ancestry_col, " and/or ", run_col)
}

eur_runs <- sdrf %>%
  filter(.data[[ancestry_col]] %in% cfg$eur_pops) %>%
  pull(.data[[run_col]]) %>%
  unique()

sanity_check(
  condition = length(eur_runs) > 0L,
  pass_msg = paste0("European run IDs selected: n = ", length(eur_runs)),
  fail_msg = "No EUR run IDs selected from SDRF."
)

safe_message("Step 2: loading TPM and counts matrices")
tpm_raw <- fread(cfg$tpm_tsv)
counts_raw <- fread(cfg$counts_tsv)

tpm_gene_col <- detect_gene_col(tpm_raw)
counts_gene_col <- detect_gene_col(counts_raw)

tpm_tbl <- select_and_rename_samples(tpm_raw, tpm_gene_col, eur_runs, label = "TPM matrix")
counts_tbl <- select_and_rename_samples(counts_raw, counts_gene_col, eur_runs, label = "Counts matrix")

shared_samples <- intersect(
  setdiff(names(tpm_tbl), "gene_id_clean"),
  setdiff(names(counts_tbl), "gene_id_clean")
)

sanity_check(
  condition = length(shared_samples) >= 10L,
  pass_msg = paste0("Shared EUR samples between TPM and counts: n = ", length(shared_samples)),
  fail_msg = "Too few shared samples between TPM and counts."
)

tpm_tbl <- tpm_tbl %>%
  select(gene_id_clean, all_of(shared_samples)) %>%
  mutate(across(all_of(shared_samples), as.numeric))

counts_tbl <- counts_tbl %>%
  select(gene_id_clean, all_of(shared_samples)) %>%
  mutate(across(all_of(shared_samples), as.numeric))

tpm_tbl <- collapse_duplicate_genes(
  tbl = tpm_tbl,
  sample_cols = shared_samples,
  label = "TPM matrix",
  collapse_fn = ~ mean(.x, na.rm = TRUE)
)
counts_tbl <- collapse_duplicate_genes(
  tbl = counts_tbl,
  sample_cols = shared_samples,
  label = "Counts matrix",
  collapse_fn = ~ sum(.x, na.rm = TRUE)
)

shared_genes <- intersect(tpm_tbl$gene_id_clean, counts_tbl$gene_id_clean)
sanity_check(
  condition = length(shared_genes) >= 1000L,
  pass_msg = paste0("Shared genes between TPM and counts: n = ", length(shared_genes)),
  fail_msg = "Too few shared genes between TPM and counts."
)

tpm_tbl <- tpm_tbl %>%
  filter(gene_id_clean %in% shared_genes) %>%
  arrange(gene_id_clean)
counts_tbl <- counts_tbl %>%
  filter(gene_id_clean %in% shared_genes) %>%
  arrange(gene_id_clean)

stopifnot(identical(tpm_tbl$gene_id_clean, counts_tbl$gene_id_clean))

safe_message("Step 3: applying GTEx-style expression filter")
n_samples <- length(shared_samples)
min_samples <- ceiling(cfg$gtex_sample_frac_threshold * n_samples)

tpm_mat <- tpm_tbl %>%
  select(all_of(shared_samples)) %>%
  as.matrix()

counts_mat <- counts_tbl %>%
  select(all_of(shared_samples)) %>%
  as.matrix()

tpm_mat[!is.finite(tpm_mat)] <- 0
counts_mat[!is.finite(counts_mat)] <- 0
counts_mat[counts_mat < 0] <- 0

tpm_keep <- rowSums(tpm_mat >= cfg$gtex_tpm_threshold, na.rm = TRUE) >= min_samples

if (identical(cfg$expression_filter_mode, "gtex_tpm_and_counts")) {
  count_keep <- rowSums(counts_mat >= cfg$gtex_count_threshold, na.rm = TRUE) >= min_samples
  keep_mask <- tpm_keep & count_keep
} else if (identical(cfg$expression_filter_mode, "gtex_tpm_only")) {
  keep_mask <- tpm_keep
} else {
  stop("Unknown cfg$expression_filter_mode: ", cfg$expression_filter_mode)
}

filtered_genes <- tpm_tbl$gene_id_clean[keep_mask]
sanity_check(
  condition = length(filtered_genes) > 0L,
  pass_msg = paste0("Genes retained after expression filter: n = ", length(filtered_genes)),
  fail_msg = "No genes passed the expression filter."
)

counts_filt <- counts_mat[keep_mask, , drop = FALSE]
rownames(counts_filt) <- filtered_genes
colnames(counts_filt) <- shared_samples

safe_message("Step 4: TMM-normalizing counts")
dge <- edgeR::DGEList(counts = counts_filt)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
tmm_cpm <- edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = FALSE, prior.count = 0)

sanity_check(
  condition = all(dim(tmm_cpm) == dim(counts_filt)),
  pass_msg = paste0(
    "TMM CPM matrix dimensions: ",
    nrow(tmm_cpm),
    " genes x ",
    ncol(tmm_cpm),
    " samples"
  ),
  fail_msg = "TMM matrix dimensions are inconsistent."
)

safe_message("Step 5: inverse-normal transform per gene (across samples)")
phenotype_mat <- t(apply(tmm_cpm, 1L, inverse_normal_transform_vec))
rownames(phenotype_mat) <- rownames(tmm_cpm)
colnames(phenotype_mat) <- colnames(tmm_cpm)

sanity_check(
  condition = all(dim(phenotype_mat) == dim(tmm_cpm)),
  pass_msg = "INT phenotype matrix dimensions match TMM matrix.",
  fail_msg = "INT phenotype matrix dimensions are inconsistent."
)

sample_ids_out <- colnames(phenotype_mat)

if (!is.null(cfg$sample_id_map_tsv) && nzchar(cfg$sample_id_map_tsv)) {
  stop_if_missing(cfg$sample_id_map_tsv, "sample ID map TSV")

  sample_map <- fread(cfg$sample_id_map_tsv) %>%
    as_tibble()

  if (ncol(sample_map) < 2L) {
    stop("sample_id_map_tsv must have at least two columns: sample_id and participant_id.")
  }

  names(sample_map)[1:2] <- c("sample_id", "participant_id")
  sample_map <- sample_map %>%
    mutate(
      sample_id = as.character(sample_id),
      participant_id = as.character(participant_id)
    )

  missing_map <- setdiff(sample_ids_out, sample_map$sample_id)
  if (length(missing_map) > 0L) {
    stop(
      "sample_id_map_tsv is missing mappings for ",
      length(missing_map),
      " samples."
    )
  }

  sample_ids_out <- sample_map$participant_id[match(sample_ids_out, sample_map$sample_id)]
  colnames(tmm_cpm) <- sample_ids_out
  colnames(phenotype_mat) <- sample_ids_out
}

safe_message("Step 6: writing outputs")
writeLines(filtered_genes, file.path(cfg$output_dir, "filtered_gene_ids.txt"))
writeLines(colnames(phenotype_mat), file.path(cfg$output_dir, "phenotype_sample_ids.txt"))

write_matrix_tsv(
  mat = tmm_cpm,
  out_tsv_gz = file.path(cfg$output_dir, "tmm_cpm_gene_by_sample.tsv.gz")
)
write_matrix_tsv(
  mat = phenotype_mat,
  out_tsv_gz = file.path(cfg$output_dir, "phenotype_int_tmm_gene_by_sample.tsv.gz")
)

run_summary <- tibble(
  metric = c(
    "n_samples_input",
    "n_genes_shared_tpm_counts",
    "n_genes_retained",
    "gtex_tpm_threshold",
    "gtex_count_threshold",
    "gtex_sample_frac_threshold",
    "min_samples_required",
    "expression_filter_mode"
  ),
  value = c(
    n_samples,
    length(shared_genes),
    length(filtered_genes),
    cfg$gtex_tpm_threshold,
    cfg$gtex_count_threshold,
    cfg$gtex_sample_frac_threshold,
    min_samples,
    cfg$expression_filter_mode
  )
)

fwrite(
  as.data.table(run_summary),
  file.path(cfg$output_dir, "phenotype_preparation_summary.tsv"),
  sep = "\t"
)

message("Done. Outputs written to: ", normalizePath(cfg$output_dir))
