#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(peer)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 11) {
  stop(
    "Usage: prepare_phenotypes.R ",
    "tpm_file counts_file map_file pca_file fam_file out_prefix peer_input ",
    "tpm_threshold count_threshold sample_frac_threshold n_pcs"
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

dedup_cols_by_run <- function(cols) {
  map_dt <- data.table(raw_col = cols)
  map_dt[, run_id := sub("_[0-9]+$", "", raw_col)]
  map_dt <- map_dt[!duplicated(run_id)]
  map_dt
}

collapse_gene_matrix <- function(mat, gene_ids, mode = c("mean", "sum")) {
  mode <- match.arg(mode)
  dt <- as.data.table(mat)
  dt[, gene_id_clean := gene_ids]

  if (mode == "mean") {
    out <- dt[, lapply(.SD, mean, na.rm = TRUE), by = gene_id_clean]
  } else {
    out <- dt[, lapply(.SD, sum, na.rm = TRUE), by = gene_id_clean]
  }

  gene_out <- out$gene_id_clean
  out[, gene_id_clean := NULL]
  mat_out <- as.matrix(out)
  rownames(mat_out) <- gene_out
  mat_out
}

inv_norm <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  if (sum(ok) <= 1L) {
    out <- rep(NA_real_, length(x))
    if (sum(ok) == 1L) out[ok] <- 0
    return(out)
  }
  r <- rank(x[ok], ties.method = "average")
  out <- rep(NA_real_, length(x))
  out[ok] <- qnorm((r - 0.5) / sum(ok))
  out
}

if (!file.exists(tpm_file)) stop("Missing tpm_file: ", tpm_file)
if (!file.exists(counts_file)) stop("Missing counts_file: ", counts_file)
if (!file.exists(map_file)) stop("Missing map_file: ", map_file)
if (!file.exists(pca_file)) stop("Missing pca_file: ", pca_file)
if (!file.exists(fam_file)) stop("Missing fam_file: ", fam_file)

log_msg("Loading inputs")
eur_map <- fread(map_file) # RNA_ID, IID, FID
pca <- fread(pca_file, header = FALSE)
fam <- fread(fam_file, header = FALSE)
tpm_dt <- fread(tpm_file)
counts_dt <- fread(counts_file)

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
setorder(sample_map, IID)

if (nrow(sample_map) == 0L) {
  stop("No mapped EUR samples after IID de-duplication.")
}
log_msg(sprintf("Mapped EUR samples: n = %d", nrow(sample_map)))

tpm_mat <- as.matrix(tpm_dt[, ..sample_map$tpm_col])
counts_mat <- as.matrix(counts_dt[, ..sample_map$counts_col])
mode(tpm_mat) <- "numeric"
mode(counts_mat) <- "numeric"

tpm_mat[!is.finite(tpm_mat)] <- 0
counts_mat[!is.finite(counts_mat)] <- 0
counts_mat[counts_mat < 0] <- 0

colnames(tpm_mat) <- sample_map$IID
colnames(counts_mat) <- sample_map$IID

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

writeLines(rownames(counts_filt), "filtered_gene_ids.txt")

log_msg("Running TMM normalization on raw counts")
dge <- DGEList(counts = counts_filt)
dge <- calcNormFactors(dge, method = "TMM")
tmm_cpm <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE, prior.count = 0)

log_msg("Applying inverse normal transform per gene")
log_expr <- log2(tmm_cpm + 1)
int_gene_by_sample <- t(apply(log_expr, 1, inv_norm))

exp_int <- t(int_gene_by_sample) # rows=samples, cols=genes
if (!identical(rownames(exp_int), sample_map$IID)) {
  rownames(exp_int) <- sample_map$IID
}

if (peer_input == "auto") {
  if (N < 150) {
    target_nk <- 15
  } else if (N < 250) {
    target_nk <- 30
  } else if (N < 350) {
    target_nk <- 45
  } else {
    target_nk <- 60
  }
} else {
  target_nk <- as.integer(peer_input)
}
log_msg(sprintf("Estimating PEER factors: K=%d", target_nk))

model <- PEER()
PEER_setPhenoMean(model, exp_int)
PEER_setNk(model, target_nk)
PEER_update(model)

peer_mat <- PEER_getX(model)[, -1, drop = FALSE]
peer_df <- as.data.table(peer_mat)
if (ncol(peer_df) > 0L) {
  setnames(peer_df, paste0("PEER", seq_len(ncol(peer_df))))
}
peer_df[, IID := rownames(exp_int)]
peer_df <- merge(peer_df, sample_map[, .(IID, FID)], by = "IID", all.x = TRUE)

setnames(pca, c("FID", "IID", paste0("PC", seq_len(ncol(pca) - 2L))))
pc_cols <- paste0("PC", seq_len(min(n_pcs, ncol(pca) - 2L)))
pca_subset <- pca[, c("FID", "IID", pc_cols), with = FALSE]
final_covar <- merge(pca_subset, peer_df, by = c("FID", "IID"))

fwrite(final_covar, paste0(out_prefix, ".qcovar"), sep = "\t", col.names = FALSE)

pheno_dt <- as.data.table(exp_int, keep.rownames = "IID")
pheno_dt <- merge(pheno_dt, sample_map[, .(IID, FID)], by = "IID")
setcolorder(pheno_dt, c("FID", "IID", setdiff(names(pheno_dt), c("FID", "IID"))))

ordered_ids <- final_covar[, .(FID, IID)]
pheno_final <- merge(ordered_ids, pheno_dt, by = c("FID", "IID"))
fwrite(pheno_final, paste0(out_prefix, ".phenotypes.tsv"), sep = "\t", col.names = FALSE)

genes <- names(pheno_final)[-(1:2)]
map_dt <- data.table(
  gene_name = genes,
  mpheno_index = seq_along(genes)
)
fwrite(map_dt, "gene_index_map.txt", sep = "\t")

log_msg("Phenotype prep completed successfully")
