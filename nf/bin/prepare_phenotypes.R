#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(peer)
})

SCRIPT_VERSION <- "prepare_phenotypes.R 2026-04-29a"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 11 || length(args) > 14) {
  stop(
    "Usage: prepare_phenotypes.R ",
    "tpm_file counts_file map_file pca_file fam_file out_prefix peer_input ",
    "tpm_threshold count_threshold sample_frac_threshold n_pcs [normalization_type] [peer_max_genes] ",
    "OR tpm_file ... n_pcs expression_source normalization_type [peer_max_genes]"
  )
}

tpm_file <- args[1]
counts_file <- args[2]
map_file <- args[3]
pca_file <- args[4]
fam_file <- args[5]
out_prefix <- args[6]
peer_input <- args[7]
tpm_threshold <- as.numeric(args[8])
count_threshold <- as.numeric(args[9])
sample_frac_threshold <- as.numeric(args[10])
n_pcs <- as.integer(args[11])
allowed_expression_sources <- c("tmm", "tpm")
expression_source <- "tmm"
normalization_type <- "irnt"
peer_max_genes <- 0L

if (length(args) == 12L) {
  token12 <- tolower(trimws(args[12]))
  if (token12 %in% allowed_expression_sources) {
    expression_source <- token12
  } else {
    normalization_type <- token12
  }
} else if (length(args) == 13L) {
  token12 <- tolower(trimws(args[12]))
  token13 <- tolower(trimws(args[13]))
  if (token12 %in% allowed_expression_sources) {
    expression_source <- token12
    normalization_type <- token13
  } else {
    normalization_type <- token12
    peer_max_genes <- suppressWarnings(as.integer(args[13]))
  }
} else if (length(args) >= 14L) {
  expression_source <- tolower(trimws(args[12]))
  normalization_type <- tolower(trimws(args[13]))
  peer_max_genes <- suppressWarnings(as.integer(args[14]))
}
if (identical(normalization_type, "ukb_irnt")) {
  normalization_type <- "irnt"
}
if (identical(normalization_type, "tmm_only")) {
  normalization_type <- "raw"
}
if (!expression_source %in% allowed_expression_sources) {
  stop(
    "Invalid expression_source: '",
    expression_source,
    "'. Allowed: ",
    paste(allowed_expression_sources, collapse = ", ")
  )
}
allowed_normalization_types <- c("irnt", "inverse_normal", "raw")
if (!normalization_type %in% allowed_normalization_types) {
  stop(
    "Invalid normalization_type: '",
    normalization_type,
    "'. Allowed: ",
    paste(allowed_normalization_types, collapse = ", ")
  )
}
if (!is.finite(peer_max_genes) || is.na(peer_max_genes) || peer_max_genes < 0L) {
  peer_max_genes <- 0L
}

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

strip_gene_version <- function(x) {
  sub("\\.[0-9]+$", "", as.character(x))
}

detect_gene_col <- function(dt) {
  cands <- c("gene_id", "Gene", "feature_id", "gene")
  gene_col <- cands[cands %in% names(dt)][1]
  if (is.na(gene_col) || !nzchar(gene_col)) {
    stop("Could not identify gene ID column in input table.")
  }
  gene_col
}

require_cols <- function(dt, required, label) {
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0L) {
    stop(
      label,
      " is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }
}

dedup_cols_by_run <- function(cols) {
  map_dt <- data.table(raw_col = as.character(cols))
  map_dt[, raw_col_base := sub("\\.[0-9]+$", "", raw_col)]
  map_dt[, run_id := sub("_[0-9]+$", "", raw_col_base)]
  map_dt <- map_dt[!duplicated(run_id)]
  map_dt[, raw_col_base := NULL]
  map_dt
}

numeric_sample_dt <- function(dt, cols, label) {
  if (length(cols) == 0L) {
    stop("No sample columns selected for ", label, ".")
  }

  out <- vector("list", length(cols))
  bad_total <- 0L
  bad_cols <- character()
  bad_examples <- character()

  for (i in seq_along(cols)) {
    col_name <- cols[[i]]
    raw_vals <- dt[[col_name]]
    if (is.numeric(raw_vals)) {
      num_vals <- as.numeric(raw_vals)
      bad_mask <- rep(FALSE, length(num_vals))
      chr_vals <- rep("", length(num_vals))
    } else {
      chr_vals <- trimws(as.character(raw_vals))
      num_vals <- suppressWarnings(as.numeric(chr_vals))
      bad_mask <- is.na(num_vals) & !is.na(raw_vals) & nzchar(chr_vals)
    }

    if (any(bad_mask)) {
      bad_total <- bad_total + sum(bad_mask)
      bad_cols <- c(bad_cols, col_name)
      bad_examples <- unique(c(bad_examples, unique(chr_vals[bad_mask])))
      if (length(bad_examples) > 5L) {
        bad_examples <- bad_examples[seq_len(5L)]
      }
    }

    out[[i]] <- num_vals
  }

  out_dt <- as.data.table(out)
  setnames(out_dt, cols)

  if (bad_total > 0L) {
    log_msg(sprintf(
      "%s: coerced %d non-numeric values to NA across %d sample columns (examples: %s)",
      label,
      bad_total,
      length(unique(bad_cols)),
      paste(bad_examples, collapse = ", ")
    ))
  }

  out_dt
}

collapse_gene_matrix <- function(mat, gene_ids, mode = c("mean", "sum")) {
  mode <- match.arg(mode)
  mat <- as.matrix(mat)
  if (is.null(dim(mat))) {
    mat <- matrix(mat, ncol = 1L)
  }
  if (length(gene_ids) != nrow(mat)) {
    stop(
      "Gene ID length/matrix row mismatch in collapse_gene_matrix: ",
      length(gene_ids),
      " gene IDs for ",
      nrow(mat),
      " matrix rows."
    )
  }
  if (anyDuplicated(gene_ids) == 0L) {
    rownames(mat) <- gene_ids
    return(mat)
  }

  g <- as.character(gene_ids)
  sum_mat <- rowsum(mat, group = g, reorder = FALSE, na.rm = TRUE)
  if (mode == "sum") {
    return(sum_mat)
  }

  non_missing <- rowsum((!is.na(mat)) * 1.0, group = g, reorder = FALSE, na.rm = TRUE)
  mean_mat <- sum_mat / non_missing
  mean_mat[!is.finite(mean_mat)] <- NA_real_
  mean_mat
}

rank_inverse_normal <- function(x, offset = 0.5) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  if (sum(ok) <= 1L) {
    out <- rep(NA_real_, length(x))
    if (sum(ok) == 1L) out[ok] <- 0
    return(out)
  }
  r <- rank(x[ok], ties.method = "average")
  n <- sum(ok)
  out <- rep(NA_real_, length(x))
  out[ok] <- qnorm((r - offset) / (n - 2 * offset + 1))
  out
}

if (!file.exists(tpm_file)) stop("Missing tpm_file: ", tpm_file)
if (!file.exists(counts_file)) stop("Missing counts_file: ", counts_file)
if (!file.exists(map_file)) stop("Missing map_file: ", map_file)
if (!file.exists(pca_file)) stop("Missing pca_file: ", pca_file)
if (!file.exists(fam_file)) stop("Missing fam_file: ", fam_file)

log_msg("Loading inputs")
log_msg(sprintf("Script version: %s", SCRIPT_VERSION))
log_msg(sprintf("Expression source: %s", expression_source))
log_msg(sprintf("Normalization mode: %s", normalization_type))
log_msg(sprintf("PEER max genes (0 = all): %d", peer_max_genes))
eur_map <- fread(map_file) # RNA_ID, IID, FID
pca <- fread(pca_file, header = FALSE)
tpm_dt <- fread(tpm_file)
counts_dt <- fread(counts_file)
require_cols(eur_map, c("RNA_ID", "IID", "FID"), "map_file")
eur_map[, RNA_ID := as.character(RNA_ID)]
eur_map[, IID := as.character(IID)]
eur_map[, FID := as.character(FID)]
names(tpm_dt) <- make.unique(names(tpm_dt))
names(counts_dt) <- make.unique(names(counts_dt))

tpm_gene_col <- detect_gene_col(tpm_dt)
counts_gene_col <- detect_gene_col(counts_dt)

tpm_sample_cols <- setdiff(names(tpm_dt), tpm_gene_col)
counts_sample_cols <- setdiff(names(counts_dt), counts_gene_col)

tpm_lookup <- dedup_cols_by_run(tpm_sample_cols)
counts_lookup <- dedup_cols_by_run(counts_sample_cols)

common_runs <- intersect(tpm_lookup$run_id, counts_lookup$run_id)
common_runs <- intersect(common_runs, eur_map$RNA_ID)

if (length(common_runs) == 0L) {
  stop("No overlapping RNA IDs across TPM, counts, and EUR map.")
}

sample_map <- data.table(run_id = common_runs)
sample_map <- merge(sample_map, eur_map[, .(RNA_ID, IID, FID)], by.x = "run_id", by.y = "RNA_ID")
sample_map <- sample_map[!duplicated(IID)]
sample_map <- merge(sample_map, tpm_lookup[, .(run_id, tpm_col = raw_col)], by = "run_id")
sample_map <- merge(sample_map, counts_lookup[, .(run_id, counts_col = raw_col)], by = "run_id")
sample_map <- sample_map[!is.na(IID) & !is.na(FID) & !is.na(tpm_col) & !is.na(counts_col)]
sample_map[, IID := as.character(IID)]
sample_map[, FID := as.character(FID)]
sample_map[, tpm_col := as.character(tpm_col)]
sample_map[, counts_col := as.character(counts_col)]
setorder(sample_map, IID)

if (nrow(sample_map) == 0L) {
  stop("No mapped EUR samples after IID de-duplication.")
}
log_msg(sprintf("Mapped EUR samples: n = %d", nrow(sample_map)))

tpm_cols <- as.character(sample_map[["tpm_col"]])
counts_cols <- as.character(sample_map[["counts_col"]])
iid_vec <- as.character(sample_map[["IID"]])

if (length(tpm_cols) != nrow(sample_map) || length(counts_cols) != nrow(sample_map)) {
  stop("Internal sample_map mismatch between rows and mapped column vectors.")
}
if (anyDuplicated(iid_vec) > 0L) {
  dup_iid <- unique(iid_vec[duplicated(iid_vec)])
  stop("Duplicated IID values after mapping: ", paste(head(dup_iid, 10), collapse = ", "))
}
if (anyDuplicated(tpm_cols) > 0L) {
  dup_cols <- unique(tpm_cols[duplicated(tpm_cols)])
  stop("Duplicated TPM column mappings detected: ", paste(head(dup_cols, 10), collapse = ", "))
}
if (anyDuplicated(counts_cols) > 0L) {
  dup_cols <- unique(counts_cols[duplicated(counts_cols)])
  stop("Duplicated count column mappings detected: ", paste(head(dup_cols, 10), collapse = ", "))
}
if (!all(tpm_cols %in% names(tpm_dt))) {
  missing_cols <- setdiff(tpm_cols, names(tpm_dt))
  stop("Missing TPM columns after mapping: ", paste(head(missing_cols, 10), collapse = ", "))
}
if (!all(counts_cols %in% names(counts_dt))) {
  missing_cols <- setdiff(counts_cols, names(counts_dt))
  stop("Missing count columns after mapping: ", paste(head(missing_cols, 10), collapse = ", "))
}

tpm_num_dt <- numeric_sample_dt(tpm_dt, tpm_cols, "TPM matrix")
counts_num_dt <- numeric_sample_dt(counts_dt, counts_cols, "Count matrix")

tpm_mat <- as.matrix(tpm_num_dt)
counts_mat <- as.matrix(counts_num_dt)
storage.mode(tpm_mat) <- "double"
storage.mode(counts_mat) <- "double"
if (!is.matrix(tpm_mat) || !is.matrix(counts_mat)) {
  stop("Failed to build sample matrices from TPM/count tables.")
}
log_msg(sprintf(
  "Matrix dimensions before colnames assignment: TPM=%d x %d, Counts=%d x %d, IID_n=%d",
  nrow(tpm_mat), ncol(tpm_mat), nrow(counts_mat), ncol(counts_mat), length(iid_vec)
))
if (ncol(tpm_mat) != length(iid_vec)) {
  stop(
    "TPM matrix/sample ID mismatch: ncol(tpm_mat)=",
    ncol(tpm_mat),
    " vs length(IID)=",
    length(iid_vec),
    " | first_tpm_cols=",
    paste(head(colnames(tpm_mat), 5), collapse = ", "),
    " | first_iids=",
    paste(head(iid_vec, 5), collapse = ", ")
  )
}
if (ncol(counts_mat) != length(iid_vec)) {
  stop(
    "Count matrix/sample ID mismatch: ncol(counts_mat)=",
    ncol(counts_mat),
    " vs length(IID)=",
    length(iid_vec),
    " | first_count_cols=",
    paste(head(colnames(counts_mat), 5), collapse = ", "),
    " | first_iids=",
    paste(head(iid_vec, 5), collapse = ", ")
  )
}

tpm_mat[!is.finite(tpm_mat)] <- 0
counts_mat[!is.finite(counts_mat)] <- 0
counts_mat[counts_mat < 0] <- 0

colnames(tpm_mat) <- iid_vec
colnames(counts_mat) <- iid_vec

tpm_gene_ids <- strip_gene_version(tpm_dt[[tpm_gene_col]])
counts_gene_ids <- strip_gene_version(counts_dt[[counts_gene_col]])

tpm_mat <- collapse_gene_matrix(tpm_mat, tpm_gene_ids, mode = "mean")
counts_mat <- collapse_gene_matrix(counts_mat, counts_gene_ids, mode = "sum")

shared_genes <- intersect(rownames(tpm_mat), rownames(counts_mat))
if (length(shared_genes) == 0L) {
  stop("No shared genes between TPM and counts after gene ID cleanup.")
}

tpm_mat <- tpm_mat[shared_genes, , drop = FALSE]
counts_mat <- counts_mat[shared_genes, , drop = FALSE]

N <- ncol(tpm_mat)
min_samples <- ceiling(sample_frac_threshold * N)
log_msg(sprintf("Applying mandatory GTEx-style filter: TPM >= %.3f AND counts >= %d in >= %d/%d samples",
                tpm_threshold, as.integer(count_threshold), min_samples, N))

tpm_keep <- rowSums(tpm_mat >= tpm_threshold) >= min_samples
count_keep <- rowSums(counts_mat >= count_threshold) >= min_samples
keep <- tpm_keep & count_keep

tpm_filt <- tpm_mat[keep, , drop = FALSE]
counts_filt <- counts_mat[keep, , drop = FALSE]
log_msg(sprintf("Genes retained after mandatory filter: %d", nrow(counts_filt)))

if (nrow(counts_filt) == 0L) {
  stop("No genes passed mandatory TPM+count filtering.")
}

writeLines(rownames(counts_filt), paste0(out_prefix, ".filtered_gene_ids.txt"))

expr_base_by_gene <- if (identical(expression_source, "tpm")) {
  log_msg("Running TMM normalization on raw counts: skipped (expression_source=tpm)")
  log_msg("Using TPM matrix as phenotype source")
  tpm_filt
} else {
  log_msg("Running TMM normalization on raw counts")
  dge <- DGEList(counts = counts_filt)
  dge <- calcNormFactors(dge, method = "TMM")
  tmm_cpm <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE, prior.count = 0)
  log_msg("Using TMM-normalized CPM matrix as phenotype source")
  tmm_cpm
}

if (identical(normalization_type, "raw")) {
  log_msg("Using raw expression scale (no rank-based normalization)")
  transformed_gene_by_sample <- expr_base_by_gene
} else {
  log_expr <- log2(expr_base_by_gene + 1)
  norm_offset <- if (identical(normalization_type, "irnt")) 3 / 8 else 0.5
  norm_label <- if (identical(normalization_type, "irnt")) {
    "Applying IRNT per gene (Blom offset c=3/8)"
  } else {
    "Applying inverse normal transform per gene (offset c=0.5)"
  }
  log_msg(norm_label)
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    ok <- is.finite(log_expr)
    work <- log_expr
    work[!ok] <- NA_real_
    rr <- matrixStats::rowRanks(work, ties.method = "average", preserveShape = TRUE)
    n_ok <- rowSums(ok)
    denom <- n_ok - 2 * norm_offset + 1
    transformed_gene_by_sample <- matrix(
      NA_real_,
      nrow = nrow(work),
      ncol = ncol(work),
      dimnames = dimnames(work)
    )
    valid_rows <- which(n_ok > 1L)
    for (i in valid_rows) {
      idx <- ok[i, ]
      transformed_gene_by_sample[i, idx] <- qnorm((rr[i, idx] - norm_offset) / denom[i])
    }
    singleton_rows <- which(n_ok == 1L)
    for (i in singleton_rows) {
      idx <- ok[i, ]
      transformed_gene_by_sample[i, idx] <- 0
    }
  } else {
    transformed_gene_by_sample <- t(apply(
      log_expr,
      1,
      function(x) rank_inverse_normal(x, offset = norm_offset)
    ))
    if (!identical(dim(transformed_gene_by_sample), dim(log_expr))) {
      transformed_gene_by_sample <- matrix(
        as.numeric(transformed_gene_by_sample),
        nrow = nrow(log_expr),
        ncol = ncol(log_expr),
        dimnames = list(rownames(log_expr), colnames(log_expr))
      )
    }
  }
}

exp_int <- t(transformed_gene_by_sample) # rows=samples, cols=genes
if (!identical(rownames(exp_int), sample_map$IID)) {
  rownames(exp_int) <- sample_map$IID
}

if (peer_input == "auto") {
  log_msg("PEER auto rule (GTEx): 15 if N<150, 30 if 150<=N<250, 45 if 250<=N<350, 60 if N>=350")
  if (N < 150) {
    requested_nk <- 15
  } else if (N < 250) {
    requested_nk <- 30
  } else if (N < 350) {
    requested_nk <- 45
  } else {
    requested_nk <- 60
  }
  log_msg(sprintf("Auto-selected PEER request K=%d using N=%d samples", requested_nk, N))
} else {
  requested_nk <- as.integer(peer_input)
  log_msg(sprintf("User-specified PEER request K=%d", requested_nk))
}
if (!is.finite(requested_nk) || requested_nk < 0L) {
  stop("Invalid PEER factor request: ", peer_input)
}

max_peer_nk <- min(nrow(exp_int) - 1L, ncol(exp_int) - 1L)
if (!is.finite(max_peer_nk)) {
  max_peer_nk <- 0L
}
max_peer_nk <- as.integer(max(0L, max_peer_nk))

target_nk <- as.integer(min(requested_nk, max_peer_nk))
if (target_nk < requested_nk) {
  log_msg(sprintf(
    "Reducing PEER factors from requested K=%d to feasible K=%d (samples=%d, genes=%d)",
    requested_nk, target_nk, nrow(exp_int), ncol(exp_int)
  ))
}

if (target_nk > 0L) {
  peer_expr <- exp_int
  if (peer_max_genes > 0L && ncol(peer_expr) > peer_max_genes) {
    gene_var <- apply(peer_expr, 2, var, na.rm = TRUE)
    gene_var[!is.finite(gene_var)] <- -Inf
    keep_idx <- order(gene_var, decreasing = TRUE)[seq_len(peer_max_genes)]
    peer_expr <- peer_expr[, keep_idx, drop = FALSE]
    log_msg(sprintf(
      "PEER input genes reduced by variance ranking: %d -> %d genes",
      ncol(exp_int), ncol(peer_expr)
    ))
  }
  subset_max_peer_nk <- min(nrow(peer_expr) - 1L, ncol(peer_expr) - 1L)
  subset_max_peer_nk <- as.integer(max(0L, subset_max_peer_nk))
  if (target_nk > subset_max_peer_nk) {
    log_msg(sprintf(
      "Reducing PEER factors after gene subsetting: K=%d -> K=%d",
      target_nk, subset_max_peer_nk
    ))
    target_nk <- subset_max_peer_nk
  }
  if (target_nk < 1L) {
    log_msg("Skipping PEER factors: not enough dimensions after PEER gene subsetting.")
    peer_df <- data.table()
  } else {

    log_msg(sprintf("Estimating PEER factors: K=%d", target_nk))
    model <- PEER()
    PEER_setPhenoMean(model, peer_expr)
    PEER_setNk(model, target_nk)
    tryCatch(
      PEER_update(model),
      error = function(e) {
        stop("PEER_update failed at K=", target_nk, ": ", conditionMessage(e))
      }
    )

    peer_x <- PEER_getX(model)
    peer_mat <- if (ncol(peer_x) > 1L) peer_x[, -1, drop = FALSE] else matrix(nrow = nrow(peer_x), ncol = 0L)
    peer_df <- as.data.table(peer_mat)
    if (ncol(peer_df) > 0L) {
      setnames(peer_df, paste0("PEER", seq_len(ncol(peer_df))))
    }
  }
} else {
  log_msg("Skipping PEER factors: not enough dimensions after filtering.")
  peer_df <- data.table()
}
peer_df[, IID := rownames(exp_int)]
peer_df <- merge(peer_df, sample_map[, .(IID, FID)], by = "IID", all.x = TRUE)

setnames(pca, c("FID", "IID", paste0("PC", seq_len(ncol(pca) - 2L))))
pca[, FID := as.character(FID)]
pca[, IID := as.character(IID)]
pc_cols <- paste0("PC", seq_len(min(n_pcs, ncol(pca) - 2L)))
pca_subset <- pca[, c("FID", "IID", pc_cols), with = FALSE]
final_covar <- merge(pca_subset, peer_df, by = c("FID", "IID"))
if (nrow(final_covar) == 0L) {
  stop(
    "No samples remained after PCA/PEER merge. ",
    "Check IID/FID consistency across eigenvec/map/FAM."
  )
}

fwrite(final_covar, paste0(out_prefix, ".qcovar"), sep = "\t", col.names = FALSE)

pheno_dt <- as.data.table(exp_int, keep.rownames = "IID")
pheno_dt <- merge(pheno_dt, sample_map[, .(IID, FID)], by = "IID")
setcolorder(pheno_dt, c("FID", "IID", setdiff(names(pheno_dt), c("FID", "IID"))))

ordered_ids <- final_covar[, .(FID, IID)]
pheno_final <- merge(ordered_ids, pheno_dt, by = c("FID", "IID"))
if (nrow(pheno_final) == 0L) {
  stop("No rows in final phenotype table after covariate ordering merge.")
}
fwrite(pheno_final, paste0(out_prefix, ".phenotypes.tsv"), sep = "\t", col.names = FALSE)

genes <- names(pheno_final)[-(1:2)]
map_dt <- data.table(
  gene_name = genes,
  mpheno_index = seq_along(genes)
)
fwrite(map_dt, paste0(out_prefix, ".gene_index_map.txt"), sep = "\t")

log_msg("Phenotype prep completed successfully")
