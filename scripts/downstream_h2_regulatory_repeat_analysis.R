#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readxl)
  library(GenomicRanges)
  library(IRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(ggplot2)
  library(broom)
})

# ------------------------------ USER CONFIG -----------------------------------
# Edit these paths for your environment before running.
cfg <- list(
  summary_tsv = "/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260112_GREML_Nextflow_v2/results/summary/final_heritability_summary.tsv",
  shet_xlsx = "/gpfs/data/mostafavilab/shared_data/gene_information/s_het_info.xlsx",
  shet_sheet = "Supplementary Table 1",
  tpm_tsv = "/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260125_salmon_nextflow/results/matrices/gene_tpm.tsv",
  gene_gtf = "/gpfs/data/mostafavilab/shared_data/GENCODE_data/hg38/gencode.v48.annotation.gtf.gz",
  output_dir = "results/downstream_h2_regulatory_repeat",
  resource_dir = "resources/public_annotations",
  sdrf_url = "https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt",
  eur_pops = c("British", "Finnish", "Tuscan", "Utah"),
  tss_window_bp = 100000L,
  cis_window_bp = 1000000L,
  enhancer_bed = NULL,   # optional: use your own enhancer BED
  open_bed = NULL,       # optional: use your own ATAC/DNase/open BED
  repeat_rmsk = NULL     # optional: use your own UCSC rmsk.txt(.gz)
)

# Defaults are widely used public resources.
defaults <- list(
  enhancer_url = "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.ELS.bed",
  open_url = "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed",
  rmsk_url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
)
standard_chr <- paste0("chr", c(1:22, "X", "Y", "M", "MT"))

# ------------------------------ HELPERS ---------------------------------------
stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("Missing %s at: %s", label, path))
  }
}

resolve_first_existing <- function(paths, label) {
  existing <- paths[file.exists(paths)]

  if (length(existing) == 0L) {
    stop(
      sprintf(
        "Missing %s. Tried:\\n- %s",
        label,
        paste(paths, collapse = "\\n- ")
      )
    )
  }

  existing[[1]]
}

strip_gene_version <- function(x) {
  str_replace(x, "\\..*$", "")
}

safe_download <- function(url, dest_path) {
  if (file.exists(dest_path)) {
    message("Using existing: ", dest_path)
    return(dest_path)
  }

  dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)
  message("Downloading: ", url)
  download.file(url = url, destfile = dest_path, mode = "wb", quiet = FALSE)
  dest_path
}

to_ucsc_style <- function(gr) {
  out <- gr

  if (!all(startsWith(as.character(seqnames(out)), "chr"))) {
    seqlevels(out) <- seqlevels(out) %>%
      ifelse(
        startsWith(., "chr"),
        ., 
        paste0("chr", .)
      )
  }

  out
}

keep_standard_chr <- function(gr) {
  keep <- as.character(seqnames(gr)) %in% standard_chr
  gr[keep]
}

detect_chr_col <- function(dt) {
  scores <- seq_len(ncol(dt)) %>%
    lapply(function(i) {
      x <- as.character(dt[[i]])
      mean(x %in% standard_chr, na.rm = TRUE)
    }) %>%
    unlist()

  if (length(scores) == 0L || max(scores, na.rm = TRUE) < 0.05) {
    return(NA_integer_)
  }

  which.max(scores)
}

detect_start_end_cols <- function(dt, chr_col) {
  num_cols <- seq_len(ncol(dt)) %>%
    lapply(function(i) suppressWarnings(as.numeric(dt[[i]])))

  best <- list(score = -Inf, start_col = NA_integer_, end_col = NA_integer_)

  for (start_col in seq_len(ncol(dt))) {
    for (end_col in seq_len(ncol(dt))) {
      if (start_col == end_col || start_col == chr_col || end_col == chr_col) {
        next
      }

      start_num <- num_cols[[start_col]]
      end_num <- num_cols[[end_col]]
      ok <- is.finite(start_num) & is.finite(end_num) & end_num > start_num
      score <- mean(ok, na.rm = TRUE)

      if (start_col == chr_col + 1L && end_col == chr_col + 2L) {
        score <- score + 0.25
      }

      if (is.finite(score) && score > best$score) {
        best <- list(score = score, start_col = start_col, end_col = end_col)
      }
    }
  }

  if (!is.finite(best$score) || best$score < 0.20) {
    stop("Could not detect start/end coordinate columns in BED-like input.")
  }

  best
}

import_bed_like <- function(path, label) {
  gr_try <- tryCatch(import(path, format = "BED"), error = function(e) NULL)

  if (!is.null(gr_try)) {
    seqn <- as.character(seqnames(gr_try))
    chr_fraction <- mean(seqn %in% standard_chr, na.rm = TRUE)

    if (is.finite(chr_fraction) && chr_fraction > 0.80) {
      return(gr_try)
    }
  }

  message("Re-parsing ", label, " as BED-like table: ", path)
  dt <- fread(path, sep = "\t", header = FALSE, fill = TRUE, comment.char = "#")

  if (nrow(dt) == 0L || ncol(dt) < 3L) {
    stop("BED-like annotation file is empty or has fewer than 3 columns: ", path)
  }

  chr_col <- detect_chr_col(dt)
  if (is.na(chr_col)) {
    stop("Could not detect chromosome column in BED-like file: ", path)
  }

  se_cols <- detect_start_end_cols(dt, chr_col)
  chr <- as.character(dt[[chr_col]])
  start_num <- suppressWarnings(as.numeric(dt[[se_cols$start_col]]))
  end_num <- suppressWarnings(as.numeric(dt[[se_cols$end_col]]))

  keep <- chr %in% standard_chr & is.finite(start_num) & is.finite(end_num) & end_num > start_num
  if (!any(keep)) {
    stop("No valid genomic intervals found after BED-like parsing: ", path)
  }

  GRanges(
    seqnames = chr[keep],
    ranges = IRanges(
      start = as.integer(start_num[keep]) + 1L,
      end = as.integer(end_num[keep])
    ),
    strand = "*"
  )
}

symm_window <- function(gr_tss, flank_bp) {
  starts <- pmax(1L, start(gr_tss) - flank_bp)
  ends <- end(gr_tss) + flank_bp

  GRanges(
    seqnames = seqnames(gr_tss),
    ranges = IRanges(start = starts, end = ends),
    strand = "*",
    gene_id_clean = mcols(gr_tss)$gene_id_clean
  )
}

safe_spearman <- function(x, y, label_x, label_y) {
  ok <- is.finite(x) & is.finite(y)

  if (sum(ok) < 5) {
    return(
      data.table(
        x = label_x,
        y = label_y,
        estimate = NA_real_,
        p_value = NA_real_,
        n = sum(ok)
      )
    )
  }

  test <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))

  data.table(
    x = label_x,
    y = label_y,
    estimate = unname(test$estimate),
    p_value = test$p.value,
    n = sum(ok)
  )
}

count_overlaps_dt <- function(window_gr, feature_gr, col_name) {
  data.table(
    gene_id_clean = mcols(window_gr)$gene_id_clean,
    count = as.integer(countOverlaps(window_gr, feature_gr, ignore.strand = TRUE))
  ) %>%
    setnames("count", col_name)
}

count_by_group <- function(window_gr, feature_gr, feature_group, prefix) {
  group_levels <- feature_group %>%
    unique() %>%
    sort()

  group_levels <- group_levels[!is.na(group_levels) & group_levels != ""]

  out <- data.table(gene_id_clean = mcols(window_gr)$gene_id_clean)
  if (length(group_levels) == 0L) {
    return(out)
  }

  group_map <- data.table(raw = as.character(group_levels)) %>%
    .[
      ,
      clean := gsub("[^A-Za-z0-9]+", "_", raw)
    ] %>%
    .[
      ,
      clean := make.unique(clean)
    ]

  class_tables <- group_map %>%
    split(by = "raw", keep.by = FALSE) %>%
    lapply(function(x) {
      raw_name <- x$raw[[1]]
      clean_name <- x$clean[[1]]
      idx <- which(feature_group == raw_name)

      data.table(
        gene_id_clean = mcols(window_gr)$gene_id_clean,
        count = as.integer(countOverlaps(window_gr, feature_gr[idx], ignore.strand = TRUE))
      ) %>%
        setnames("count", paste0(prefix, "_", clean_name))
    })

  Reduce(
    f = function(left, right) merge(left, right, by = "gene_id_clean", all = TRUE),
    x = class_tables,
    init = out
  )
}

# --------------------------- PREP OUTPUT DIRS ---------------------------------
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$resource_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

# ---------------------- 1) EUROPEAN NON-ZERO FILTER ---------------------------
message("Step 1: load SDRF and filter GEUVADIS European runs")

sdrf <- fread(cfg$sdrf_url) %>%
  {
    names(.) <- make.unique(names(.))
    .
  }

ancestry_col <- "Characteristics[ancestry category]"
run_col <- "Comment[ENA_RUN]"

if (!ancestry_col %in% names(sdrf)) {
  stop("Missing ancestry column in SDRF.")
}
if (!run_col %in% names(sdrf)) {
  stop("Missing ENA run column in SDRF.")
}

eur_runs <- sdrf %>%
  .[get(ancestry_col) %in% cfg$eur_pops, unique(get(run_col))]

if (length(eur_runs) == 0L) {
  stop("No European runs found in SDRF with the configured population labels.")
}

message("Step 2: load TPM and keep genes with TPM > 0 in at least one EU sample")
stop_if_missing(cfg$tpm_tsv, "TPM matrix")

tpm <- fread(cfg$tpm_tsv)

gene_col <- intersect(c("gene_id", "Gene", "feature_id", "gene"), names(tpm))[1]
if (is.na(gene_col)) {
  stop("Could not find gene id column in TPM matrix.")
}

all_cols <- names(tpm)
clean_cols <- str_remove(all_cols, "_\\d+$")
keep_cols <- all_cols[clean_cols %in% eur_runs | all_cols == gene_col]

tpm_eur <- tpm %>%
  .[, ..keep_cols] %>%
  {
    eur_sample_cols <- setdiff(names(.), gene_col)

    if (length(eur_sample_cols) == 0L) {
      stop("No European sample columns matched TPM columns after cleaning suffixes.")
    }

    .[, (eur_sample_cols) := lapply(.SD, as.numeric), .SDcols = eur_sample_cols]

    keep_rows <- rowSums(as.matrix(.[, ..eur_sample_cols]) > 0, na.rm = TRUE) > 0
    . <- .[keep_rows]

    .[, gene_id_clean := strip_gene_version(get(gene_col))]
    .[, mean_tpm := rowMeans(as.matrix(.SD), na.rm = TRUE), .SDcols = eur_sample_cols]

    .
  }

expressing_gene_ids <- tpm_eur %>%
  .[, unique(gene_id_clean)]

writeLines(expressing_gene_ids, con = file.path(cfg$output_dir, "non_zero_eur_gene_ids.txt"))

data.table(
  metric = c("original_genes", "eur_nonzero_genes"),
  value = c(nrow(tpm), nrow(tpm_eur))
) %>%
  fwrite(file.path(cfg$output_dir, "expression_filter_counts.tsv"), sep = "\t")

# ----------------------- 2) HERITABILITY + S_HET TABLE ------------------------
message("Step 3: merge heritability summary and s_het")
stop_if_missing(cfg$summary_tsv, "heritability summary TSV")
stop_if_missing(cfg$shet_xlsx, "s_het xlsx")

sum_dt <- fread(cfg$summary_tsv) %>%
  .[, .(Gene, h2_GREML, SE_GREML, Pval_GREML)] %>%
  .[, gene_id_clean := strip_gene_version(Gene)] %>%
  .[]

shet_dt <- read_xlsx(cfg$shet_xlsx, sheet = cfg$shet_sheet) %>%
  as.data.table() %>%
  {
    if (!all(c("ensg", "post_mean") %in% names(.))) {
      stop("s_het sheet must include columns named 'ensg' and 'post_mean'.")
    }

    .[, gene_id_clean := strip_gene_version(ensg)]
    .
  }

merged <- sum_dt %>%
  merge(
    shet_dt[, .(gene_id_clean, post_mean)],
    by = "gene_id_clean",
    all.x = TRUE
  ) %>%
  .[!is.na(post_mean) & gene_id_clean %in% expressing_gene_ids] %>%
  {
    .[, post_mean_bin := dplyr::ntile(post_mean, 10L)]
    .
  }

if (nrow(merged) == 0L) {
  stop("No genes left after merging h2 + s_het and applying EU non-zero filter.")
}

tpm_means <- tpm_eur %>%
  .[, .(gene_id_clean, mean_tpm)] %>%
  unique(by = "gene_id_clean")

merged <- merged %>%
  merge(tpm_means, by = "gene_id_clean", all.x = TRUE)

# ----------------------- 3) GENE COORDINATES / WINDOWS ------------------------
message("Step 4: load gene annotation and build TSS windows")

gtf_path <- resolve_first_existing(
  unique(c(
    cfg$gene_gtf,
    "/gpfs/data/mostafavilab/shared_data/GENCODE_data/hg38/gencode.v48.annotation.gtf.gz",
    "/gpfs/data/mostafavilab/shared_data/annotations/gencode.v44.annotation.gtf.gz"
  )),
  label = "gene annotation GTF"
)
message("Using gene annotation: ", gtf_path)

gtf <- import(gtf_path)

gene_gr <- gtf[gtf$type == "gene"] %>%
  keepStandardChromosomes(pruning.mode = "coarse") %>%
  to_ucsc_style()

gene_tbl <- data.table(
  gene_id_clean = strip_gene_version(mcols(gene_gr)$gene_id),
  chr = as.character(seqnames(gene_gr)),
  start = start(gene_gr),
  end = end(gene_gr),
  strand = as.character(strand(gene_gr)),
  gene_name = if ("gene_name" %in% colnames(mcols(gene_gr))) {
    as.character(mcols(gene_gr)$gene_name)
  } else {
    strip_gene_version(mcols(gene_gr)$gene_id)
  },
  gene_length = width(gene_gr)
) %>%
  unique(by = "gene_id_clean") %>%
  .[gene_id_clean %in% merged$gene_id_clean]

if (nrow(gene_tbl) == 0L) {
  stop("No gene coordinates matched merged gene IDs. Check gene ID format and GTF build.")
}

tss_pos <- ifelse(gene_tbl$strand == "-", gene_tbl$end, gene_tbl$start)

tss_gr <- GRanges(
  seqnames = gene_tbl$chr,
  ranges = IRanges(start = tss_pos, end = tss_pos),
  strand = "*",
  gene_id_clean = gene_tbl$gene_id_clean
)

tss_win_small <- symm_window(tss_gr, cfg$tss_window_bp)
tss_win_cis <- symm_window(tss_gr, cfg$cis_window_bp)

# ----------------------- 4) REGULATORY FEATURE COUNTS -------------------------
message("Step 5: load enhancer/open annotations and count overlaps")

enhancer_path <- cfg$enhancer_bed
if (is.null(enhancer_path) || enhancer_path == "") {
  enhancer_path <- safe_download(
    defaults$enhancer_url,
    file.path(cfg$resource_dir, basename(defaults$enhancer_url))
  )
}

open_path <- cfg$open_bed
if (is.null(open_path) || open_path == "") {
  open_path <- safe_download(
    defaults$open_url,
    file.path(cfg$resource_dir, basename(defaults$open_url))
  )
}

enh_gr <- import_bed_like(enhancer_path, label = "enhancer annotations") %>%
  to_ucsc_style() %>%
  keep_standard_chr()

open_gr <- import_bed_like(open_path, label = "open-chromatin annotations") %>%
  to_ucsc_style() %>%
  keep_standard_chr()

reg_features <- count_overlaps_dt(tss_win_small, enh_gr, "enh_count_100kb") %>%
  merge(
    count_overlaps_dt(tss_win_cis, enh_gr, "enh_count_1mb"),
    by = "gene_id_clean",
    all = TRUE
  ) %>%
  merge(
    count_overlaps_dt(tss_win_small, open_gr, "open_count_100kb"),
    by = "gene_id_clean",
    all = TRUE
  ) %>%
  merge(
    count_overlaps_dt(tss_win_cis, open_gr, "open_count_1mb"),
    by = "gene_id_clean",
    all = TRUE
  )

# ------------------------ 5) UCSC REPEATMASKER COUNTS -------------------------
message("Step 6: load UCSC RepeatMasker annotations and summarize repeats")

rmsk_path <- cfg$repeat_rmsk
if (is.null(rmsk_path) || rmsk_path == "") {
  rmsk_path <- safe_download(
    defaults$rmsk_url,
    file.path(cfg$resource_dir, basename(defaults$rmsk_url))
  )
}

rmsk <- fread(rmsk_path, sep = "\t", header = FALSE, showProgress = TRUE)
if (ncol(rmsk) < 17L) {
  stop("Unexpected rmsk format. Expected at least 17 columns.")
}

setnames(
  rmsk,
  old = names(rmsk)[1:17],
  new = c(
    "bin", "swScore", "milliDiv", "milliDel", "milliIns", "genoName", "genoStart",
    "genoEnd", "genoLeft", "strand", "repName", "repClass", "repFamily",
    "repStart", "repEnd", "repLeft", "id"
  )
)

rmsk <- rmsk %>%
  .[
    genoName %in% paste0("chr", c(1:22, "X", "Y")) &
      !is.na(repClass) &
      genoEnd > genoStart
  ]

rep_gr <- GRanges(
  seqnames = rmsk$genoName,
  ranges = IRanges(start = rmsk$genoStart + 1L, end = rmsk$genoEnd),
  strand = "*",
  repName = rmsk$repName,
  repClass = rmsk$repClass,
  repFamily = rmsk$repFamily,
  milliDiv = rmsk$milliDiv
)

# Repeat features are window-based burdens near TSS, not gene-body overlap flags.
repeat_features <- count_overlaps_dt(tss_win_small, rep_gr, "repeat_count_100kb") %>%
  merge(
    count_overlaps_dt(tss_win_cis, rep_gr, "repeat_count_1mb"),
    by = "gene_id_clean",
    all = TRUE
  )

repeat_class_dt <- count_by_group(
  window_gr = tss_win_small,
  feature_gr = rep_gr,
  feature_group = mcols(rep_gr)$repClass,
  prefix = "repeat_class"
)

# --------------------------- 6) FINAL ANALYSIS TABLE --------------------------
analysis_dt <- merged %>%
  merge(
    gene_tbl[, .(gene_id_clean, gene_name, chr, start, end, gene_length)],
    by = "gene_id_clean",
    all.x = TRUE
  ) %>%
  merge(reg_features, by = "gene_id_clean", all.x = TRUE) %>%
  merge(repeat_features, by = "gene_id_clean", all.x = TRUE) %>%
  merge(repeat_class_dt, by = "gene_id_clean", all.x = TRUE)

count_cols <- names(analysis_dt) %>%
  str_subset("^(enh_count_|open_count_|repeat_count_|repeat_class_)")

if (length(count_cols) > 0L) {
  analysis_dt[, (count_cols) := lapply(
    .SD,
    function(x) {
      x <- as.numeric(x)
      x[is.na(x)] <- 0
      x
    }
  ), .SDcols = count_cols]
}

fwrite(analysis_dt, file.path(cfg$output_dir, "gene_level_features.tsv"), sep = "\t")

# ----------------------------- 7) SUMMARIES -----------------------------------
decile_summary <- analysis_dt %>%
  .[
    ,
    .(
      n_genes = .N,
      median_h2 = median(h2_GREML, na.rm = TRUE),
      prop_h2_sig = mean(Pval_GREML < 0.05, na.rm = TRUE),
      median_enh_100kb = median(enh_count_100kb, na.rm = TRUE),
      median_open_100kb = median(open_count_100kb, na.rm = TRUE),
      median_repeat_100kb = median(repeat_count_100kb, na.rm = TRUE)
    ),
    by = post_mean_bin
  ] %>%
  .[order(post_mean_bin)]

fwrite(decile_summary, file.path(cfg$output_dir, "decile_feature_summary.tsv"), sep = "\t")

cor_dt <- rbindlist(list(
  safe_spearman(analysis_dt$post_mean, analysis_dt$h2_GREML, "post_mean", "h2_GREML"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$enh_count_100kb, "post_mean", "enh_count_100kb"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$open_count_100kb, "post_mean", "open_count_100kb"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$repeat_count_100kb, "post_mean", "repeat_count_100kb"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$enh_count_100kb, "h2_GREML", "enh_count_100kb"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$open_count_100kb, "h2_GREML", "open_count_100kb"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$repeat_count_100kb, "h2_GREML", "repeat_count_100kb")
))

fwrite(cor_dt, file.path(cfg$output_dir, "spearman_correlations.tsv"), sep = "\t")

# --------------------------- 8) ASSOCIATION MODELS ----------------------------
model_dt <- analysis_dt %>%
  .[
    is.finite(h2_GREML) &
      is.finite(post_mean) &
      is.finite(mean_tpm) &
      is.finite(gene_length)
  ] %>%
  {
    .[, h2_sig := as.integer(Pval_GREML < 0.05)]
    .
  }

post_mean_term <- NULL
post_mean_levels <- model_dt %>%
  .[, uniqueN(post_mean_bin[!is.na(post_mean_bin)])]
if (post_mean_levels >= 2L) {
  post_mean_term <- "factor(post_mean_bin)"
} else {
  message("Skipping factor(post_mean_bin): fewer than 2 levels in model_dt.")
}

lm_terms <- c(
  "scale(log1p(enh_count_100kb))",
  "scale(log1p(open_count_100kb))",
  "scale(log1p(repeat_count_100kb))",
  "scale(log1p(mean_tpm))",
  "scale(log1p(gene_length))",
  post_mean_term
) %>%
  .[!is.na(.)]

lm_formula <- as.formula(
  paste("h2_GREML ~", paste(lm_terms, collapse = " + "))
)
lm_fit <- lm(lm_formula, data = model_dt)

fwrite(
  as.data.table(tidy(lm_fit)),
  file.path(cfg$output_dir, "model_lm_h2.tsv"),
  sep = "\t"
)

glm_terms <- c(
  "scale(log1p(enh_count_100kb))",
  "scale(log1p(open_count_100kb))",
  "scale(log1p(repeat_count_100kb))",
  "scale(log1p(mean_tpm))",
  "scale(log1p(gene_length))",
  post_mean_term
) %>%
  .[!is.na(.)]

glm_formula <- as.formula(
  paste("h2_sig ~", paste(glm_terms, collapse = " + "))
)
glm_fit <- glm(glm_formula, family = binomial(), data = model_dt)

fwrite(
  as.data.table(tidy(glm_fit)),
  file.path(cfg$output_dir, "model_glm_h2sig.tsv"),
  sep = "\t"
)

repeat_cols <- names(model_dt) %>%
  str_subset("^repeat_class_")

if (length(repeat_cols) > 0L) {
  repeat_stats <- data.table(repeat_col = repeat_cols) %>%
    .[
      ,
      `:=`(
        total = vapply(repeat_col, function(col) sum(model_dt[[col]], na.rm = TRUE), numeric(1)),
        n_unique = vapply(
          repeat_col,
          function(col) uniqueN(model_dt[[col]][!is.na(model_dt[[col]])]),
          integer(1)
        )
      )
    ]

  top_repeat_cols <- repeat_stats %>%
    .[total > 0 & n_unique > 1] %>%
    .[order(-total)] %>%
    head(6) %>%
    .[, as.character(repeat_col)]

  top_repeat_cols <- top_repeat_cols %>%
    .[!is.na(.) & nzchar(.)] %>%
    .[. %in% names(model_dt)]

  if (length(top_repeat_cols) == 0L) {
    message("Skipping repeat-class lm: no valid repeat_class predictors with variation.")
  } else {
    repeat_terms <- c(
      sprintf("scale(log1p(`%s`))", top_repeat_cols),
      "scale(log1p(mean_tpm))",
      "scale(log1p(gene_length))",
      post_mean_term
    ) %>%
      .[!is.na(.) & nzchar(.)]

    rhs <- repeat_terms %>%
      paste(collapse = " + ")

    repeat_formula <- as.formula(paste("h2_GREML ~", rhs))
    repeat_fit <- lm(repeat_formula, data = model_dt)

    fwrite(
      as.data.table(tidy(repeat_fit)),
      file.path(cfg$output_dir, "model_lm_repeat_classes.tsv"),
      sep = "\t"
    )
  }
}

# ------------------------------- 9) PLOTS -------------------------------------
message("Step 7: generate summary plots")

p_decile <- ggplot(decile_summary, aes(x = post_mean_bin, y = median_h2)) +
  geom_line(color = "#1f3b4d", linewidth = 1) +
  geom_point(color = "#d1495b", size = 2.5) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    title = "Median h2_GREML across s_het deciles",
    x = "s_het post_mean decile",
    y = "Median h2_GREML"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(cfg$output_dir, "plots", "median_h2_by_decile.pdf"),
  plot = p_decile,
  width = 8,
  height = 5
)

p_reg <- analysis_dt %>%
  as_tibble() %>%
  select(h2_GREML, enh_count_100kb, open_count_100kb) %>%
  pivot_longer(
    cols = c(enh_count_100kb, open_count_100kb),
    names_to = "feature",
    values_to = "count"
  ) %>%
  ggplot(aes(x = log1p(count), y = h2_GREML)) +
  geom_point(alpha = 0.15, size = 0.8, color = "gray35") +
  geom_smooth(
    method = "gam",
    color = "#e07a5f",
    fill = "#e07a5f",
    alpha = 0.2
  ) +
  facet_wrap(~feature, scales = "free_x") +
  labs(
    title = "Regulatory feature burden vs heritability",
    x = "log1p(feature count within +/-100 kb)",
    y = "h2_GREML"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(cfg$output_dir, "plots", "regulatory_burden_vs_h2.pdf"),
  plot = p_reg,
  width = 10,
  height = 5.5
)

repeat_class_cols <- names(analysis_dt) %>%
  str_subset("^repeat_class_")

if (length(repeat_class_cols) > 0L) {
  repeat_heat <- analysis_dt %>%
    as_tibble() %>%
    select(post_mean_bin, all_of(repeat_class_cols)) %>%
    pivot_longer(
      cols = all_of(repeat_class_cols),
      names_to = "rep_class",
      values_to = "count"
    ) %>%
    group_by(post_mean_bin, rep_class) %>%
    summarize(median_count = median(count, na.rm = TRUE), .groups = "drop")

  p_heat <- ggplot(
    repeat_heat,
    aes(x = factor(post_mean_bin), y = rep_class, fill = median_count)
  ) +
    geom_tile() +
    scale_fill_viridis_c(option = "C") +
    labs(
      title = "Repeat class burden by s_het decile",
      x = "s_het decile",
      y = "Repeat class",
      fill = "Median count"
    ) +
    theme_minimal(base_size = 11)

  ggsave(
    filename = file.path(cfg$output_dir, "plots", "repeat_class_heatmap.pdf"),
    plot = p_heat,
    width = 9,
    height = 6
  )
}

capture.output(sessionInfo()) %>%
  writeLines(con = file.path(cfg$output_dir, "sessionInfo.txt"))

message("Done. Outputs written to: ", normalizePath(cfg$output_dir))
