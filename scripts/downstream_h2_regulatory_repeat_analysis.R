#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(dplyr)
  library(purrr)
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

setwd("/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260302_GE_Heritability_and_Repeat/GeneExpression/")

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
  expression_filter_mode = "gtex_tpm",  # options: "gtex_tpm", "any_nonzero"
  gtex_tpm_threshold = 0.1,
  gtex_sample_frac_threshold = 0.2,
  tss_window_bp = 100000L,
  quarter_window_bp = 250000L,
  mid_window_bp = 500000L,
  cis_window_bp = 1000000L,
  enhancer_source = "roadmap_links",  # options: "roadmap_links", "window_count"
  enhancer_fallback_source = "window_count",  # options: "window_count", "none"
  roadmap_links_dir = NULL,           # optional local dir with Roadmap link files
  roadmap_links_url = "https://ernstlab.github.io/roadmaplinking/",
  roadmap_links_urls = c(
    "https://ernstlab.github.io/roadmaplinking/",
    "https://promoter.bx.psu.edu/hi-c/publication/data/Roadmap_links/",
    "http://promoter.bx.psu.edu/hi-c/publication/data/Roadmap_links/"
  ),
  roadmap_links_zip_url = "https://zenodo.org/records/12042965/files/RoadmapLinks.zip?download=1",
  roadmap_links_pattern = "^links_.*_2\\.5\\.txt$",
  download_timeout_sec = 300L,
  enhancer_bed = NULL,   # optional: use your own enhancer BED
  open_bed = NULL,       # optional: use your own ATAC/DNase/open BED
  repeat_rmsk = NULL,    # optional: use your own UCSC rmsk.txt(.gz)
  chrom_info_tsv = NULL, # optional: UCSC chromInfo table with chr sizes
  repeat_bg_method = "explicit_genome", # options: "explicit_genome", "poisson_window"
  repeat_bg_repeat_fraction = 0.5,      # target genome repeat fraction for reference filter set
  repeat_bg_type_col = "repClass_norm", # column used for repeat composition in simulation
  repeat_bg_composition_basis = "bp",   # options: "bp", "count"
  repeat_bg_n_iter = 100L,
  repeat_bg_seed = 20260305L
)

# Defaults are widely used public resources.
defaults <- list(
  enhancer_url = "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.ELS.bed",
  open_url = "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed",
  rmsk_url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz",
  chrom_info_url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/chromInfo.txt.gz"
)
standard_chr <- paste0("chr", c(1:22, "X", "Y"))

# ------------------------------ HELPERS ---------------------------------------
script_file_args <- commandArgs(trailingOnly = FALSE)
script_file_args <- script_file_args[grepl("^--file=", script_file_args)]

if (length(script_file_args) > 0L) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_args[[1]])))
} else {
  script_dir <- getwd()
}

helper_candidates <- c(
  file.path(script_dir, "downstream_h2_regulatory_repeat_helpers.R"),
  file.path(script_dir, "scripts", "downstream_h2_regulatory_repeat_helpers.R"),
  "scripts/downstream_h2_regulatory_repeat_helpers.R",
  "downstream_h2_regulatory_repeat_helpers.R"
)
helper_path <- helper_candidates[file.exists(helper_candidates)][1]

if (is.na(helper_path) || !nzchar(helper_path)) {
  stop("Could not locate downstream_h2_regulatory_repeat_helpers.R from script_dir: ", script_dir)
}

source(helper_path)
message("[Sanity] Loaded helper script: ", normalizePath(helper_path))
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
sanity_rows(sdrf, "SDRF table", min_rows = 10L)

ancestry_col <- "Characteristics[ancestry category]"
run_col <- "Comment[ENA_RUN]"

if (!ancestry_col %in% names(sdrf)) {
  stop("Missing ancestry column in SDRF.")
}
if (!run_col %in% names(sdrf)) {
  stop("Missing ENA run column in SDRF.")
}

eur_runs <- sdrf %>%
  as_tibble() %>%
  filter(.data[[ancestry_col]] %in% cfg$eur_pops) %>%
  pull(.data[[run_col]]) %>%
  unique()

if (length(eur_runs) == 0L) {
  stop("No European runs found in SDRF with the configured population labels.")
}
sanity_check(
  condition = length(eur_runs) >= 10L,
  pass_msg = paste0("European runs detected: n = ", length(eur_runs)),
  fail_msg = "Too few European runs detected; check cfg$eur_pops and SDRF labels."
)

message("Step 2: load TPM and apply expression filter (GTEx-style by default)")
stop_if_missing(cfg$tpm_tsv, "TPM matrix")

tpm <- fread(cfg$tpm_tsv)
sanity_rows(tpm, "TPM matrix", min_rows = 1000L)

gene_col <- intersect(c("gene_id", "Gene", "feature_id", "gene"), names(tpm))[1]
if (is.na(gene_col)) {
  stop("Could not find gene id column in TPM matrix.")
}

all_cols <- names(tpm)
clean_cols <- str_remove(all_cols, "_\\d+$")
keep_cols <- all_cols[clean_cols %in% eur_runs | all_cols == gene_col]

tpm_eur <- tpm %>%
  as_tibble() %>%
  select(all_of(keep_cols))

eur_sample_cols <- setdiff(names(tpm_eur), gene_col)
if (length(eur_sample_cols) == 0L) {
  stop("No European sample columns matched TPM columns after cleaning suffixes.")
}

n_eur_samples <- length(eur_sample_cols)
min_samples_required <- ceiling(cfg$gtex_sample_frac_threshold * n_eur_samples)

tpm_eur <- tpm_eur %>%
  mutate(across(all_of(eur_sample_cols), as.numeric))

if (identical(cfg$expression_filter_mode, "gtex_tpm")) {
  tpm_eur <- tpm_eur %>%
    mutate(
      n_samples_tpm_ge_threshold = rowSums(
        across(
          all_of(eur_sample_cols),
          ~ replace_na(.x, 0) >= cfg$gtex_tpm_threshold
        )
      ),
      pass_expression_filter = n_samples_tpm_ge_threshold >= min_samples_required
    ) %>%
    filter(pass_expression_filter) %>%
    mutate(
      gene_id_clean = strip_gene_version(.data[[gene_col]]),
      mean_tpm = rowMeans(select(., all_of(eur_sample_cols)), na.rm = TRUE)
    )
} else if (identical(cfg$expression_filter_mode, "any_nonzero")) {
  tpm_eur <- tpm_eur %>%
    mutate(
      n_samples_tpm_ge_threshold = rowSums(
        across(
          all_of(eur_sample_cols),
          ~ replace_na(.x, 0) > 0
        )
      ),
      pass_expression_filter = n_samples_tpm_ge_threshold > 0
    ) %>%
    filter(pass_expression_filter) %>%
    mutate(
      gene_id_clean = strip_gene_version(.data[[gene_col]]),
      mean_tpm = rowMeans(select(., all_of(eur_sample_cols)), na.rm = TRUE)
    )
} else {
  stop("Unknown cfg$expression_filter_mode: ", cfg$expression_filter_mode)
}

tpm_eur <- tpm_eur %>%
  as.data.table()
sanity_rows(tpm_eur, "EU expression-filtered gene table", min_rows = 1000L)
sanity_check(
  condition = n_eur_samples >= 10L,
  pass_msg = paste0("EU sample columns in TPM: n = ", n_eur_samples),
  fail_msg = "Too few EU samples matched TPM columns."
)
sanity_check(
  condition = min_samples_required >= 1L && min_samples_required <= n_eur_samples,
  pass_msg = paste0(
    "Expression threshold requires >= ",
    min_samples_required,
    " of ",
    n_eur_samples,
    " EU samples."
  ),
  fail_msg = "Invalid minimum sample count derived for expression filtering."
)
if (identical(cfg$expression_filter_mode, "gtex_tpm")) {
  sanity_check(
    condition = all(tpm_eur$n_samples_tpm_ge_threshold >= min_samples_required),
    pass_msg = paste0(
      "GTEx-style TPM filter applied: TPM >= ",
      cfg$gtex_tpm_threshold,
      " in >= ",
      min_samples_required,
      " samples."
    ),
    fail_msg = "GTEx-style TPM filter sanity check failed."
  )
}

expressing_gene_ids <- tpm_eur$gene_id_clean %>% unique()
sanity_check(
  condition = length(expressing_gene_ids) > 0L,
  pass_msg = paste0("Expressed genes retained: n = ", length(expressing_gene_ids)),
  fail_msg = "No expressed genes were retained after TPM filtering."
)

writeLines(expressing_gene_ids, con = file.path(cfg$output_dir, "non_zero_eur_gene_ids.txt"))

data.table(
  metric = c(
    "original_genes",
    "eur_genes_after_expression_filter",
    "eur_sample_n",
    "expression_min_samples_required",
    "expression_tpm_threshold",
    "expression_sample_fraction_threshold"
  ),
  value = c(
    nrow(tpm),
    nrow(tpm_eur),
    n_eur_samples,
    min_samples_required,
    cfg$gtex_tpm_threshold,
    cfg$gtex_sample_frac_threshold
  )
) %>%
  fwrite(file.path(cfg$output_dir, "expression_filter_counts.tsv"), sep = "\t")

tibble(
  setting = c(
    "expression_filter_mode",
    "gtex_reference_note",
    "gtex_reference_url"
  ),
  value = c(
    as.character(cfg$expression_filter_mode),
    "GTEx eQTL pipeline uses TPM >= 0.1 in >= 20% samples (plus read-count filter when counts are available).",
    "https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/README.md"
  )
) %>%
  as.data.table() %>%
  fwrite(file.path(cfg$output_dir, "expression_filter_method.tsv"), sep = "\t")

# ----------------------- 2) HERITABILITY + S_HET TABLE ------------------------
message("Step 3: merge heritability summary and s_het")
stop_if_missing(cfg$summary_tsv, "heritability summary TSV")
stop_if_missing(cfg$shet_xlsx, "s_het xlsx")

sum_dt <- fread(cfg$summary_tsv) %>%
  as_tibble() %>%
  select(Gene, h2_GREML, SE_GREML, Pval_GREML) %>%
  mutate(gene_id_clean = strip_gene_version(Gene))
sanity_rows(sum_dt, "Heritability summary table", min_rows = 1000L)

shet_dt <- read_xlsx(cfg$shet_xlsx, sheet = cfg$shet_sheet) %>%
  as_tibble()
sanity_rows(shet_dt, "s_het table", min_rows = 1000L)

if (!all(c("ensg", "post_mean") %in% names(shet_dt))) {
  stop("s_het sheet must include columns named 'ensg' and 'post_mean'.")
}

shet_dt <- shet_dt %>%
  mutate(gene_id_clean = strip_gene_version(ensg))

shet_coord_cols <- detect_shet_coord_columns(shet_dt)
shet_coord_dt <- NULL
if (
  !is.na(shet_coord_cols$chr_col) &&
    !is.na(shet_coord_cols$start_col) &&
    !is.na(shet_coord_cols$end_col)
) {
  shet_coord_dt <- shet_dt %>%
    transmute(
      gene_id_clean = gene_id_clean,
      chr_shet = normalize_chr(.data[[shet_coord_cols$chr_col]]),
      start_shet = suppressWarnings(as.numeric(.data[[shet_coord_cols$start_col]])),
      end_shet = suppressWarnings(as.numeric(.data[[shet_coord_cols$end_col]])),
      gene_name_shet = if (!is.na(shet_coord_cols$gene_name_col)) {
        as.character(.data[[shet_coord_cols$gene_name_col]])
      } else {
        NA_character_
      }
    ) %>%
    filter(
      !is.na(gene_id_clean),
      gene_id_clean != "",
      chr_shet %in% standard_chr,
      is.finite(start_shet),
      is.finite(end_shet),
      end_shet > start_shet
    ) %>%
    mutate(
      start_shet = as.integer(start_shet),
      end_shet = as.integer(end_shet)
    ) %>%
    distinct(gene_id_clean, .keep_all = TRUE)

  message("Using S_het coordinates for ", nrow(shet_coord_dt), " genes when available.")
} else {
  message("S_het coordinates not found in sheet; using GTF coordinates.")
}

merged <- sum_dt %>%
  left_join(select(shet_dt, gene_id_clean, post_mean), by = "gene_id_clean") %>%
  filter(!is.na(post_mean), gene_id_clean %in% expressing_gene_ids) %>%
  mutate(post_mean_bin = dplyr::ntile(post_mean, 10L))

if (nrow(merged) == 0L) {
  stop("No genes left after merging h2 + s_het and applying EU non-zero filter.")
}
sanity_rows(merged, "Merged h2 + s_het + expression table", min_rows = 1000L)
sanity_check(
  condition = dplyr::n_distinct(merged$post_mean_bin) >= 5L,
  pass_msg = paste0("post_mean bins present: n = ", dplyr::n_distinct(merged$post_mean_bin)),
  fail_msg = "Too few post_mean bins available after merge; inspect post_mean distribution."
)

tpm_means <- tpm_eur %>%
  as_tibble() %>%
  select(gene_id_clean, mean_tpm) %>%
  distinct(gene_id_clean, .keep_all = TRUE)

merged <- merged %>%
  left_join(tpm_means, by = "gene_id_clean")

# ----------------------- 3) GENE COORDINATES / WINDOWS ------------------------
message("Step 4: load gene annotation and build gene/TSS windows")

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
sanity_check(
  condition = length(gtf) > 0L,
  pass_msg = paste0("Imported GTF entries: n = ", format(length(gtf), big.mark = ",")),
  fail_msg = "Imported GTF is empty."
)

gene_gr <- gtf[gtf$type == "gene"] %>%
  keepStandardChromosomes(pruning.mode = "coarse") %>%
  to_ucsc_style()

gene_tbl <- tibble(
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
  distinct(gene_id_clean, .keep_all = TRUE) %>%
  filter(gene_id_clean %in% merged$gene_id_clean)

if (!is.null(shet_coord_dt) && nrow(shet_coord_dt) > 0L) {
  gene_tbl <- gene_tbl %>%
    full_join(shet_coord_dt, by = "gene_id_clean") %>%
    mutate(
      chr = coalesce(chr_shet, chr),
      start = if_else(is.finite(start_shet), as.integer(start_shet), as.integer(start)),
      end = if_else(is.finite(end_shet), as.integer(end_shet), as.integer(end)),
      gene_name = if_else(
        !is.na(gene_name_shet) & gene_name_shet != "",
        gene_name_shet,
        gene_name
      ),
      strand = if_else(is.na(strand) | strand == "", "+", strand),
      gene_length = as.integer(end - start + 1L)
    ) %>%
    filter(
      !is.na(gene_id_clean),
      gene_id_clean %in% merged$gene_id_clean,
      chr %in% standard_chr,
      is.finite(start),
      is.finite(end),
      end > start
    ) %>%
    select(-any_of(c("chr_shet", "start_shet", "end_shet", "gene_name_shet"))) %>%
    distinct(gene_id_clean, .keep_all = TRUE)
}

if (nrow(gene_tbl) == 0L) {
  stop("No gene coordinates matched merged gene IDs. Check gene ID format and GTF build.")
}
sanity_rows(gene_tbl, "Gene coordinate table", min_rows = 1000L)

tss_pos <- ifelse(gene_tbl$strand == "-", gene_tbl$end, gene_tbl$start)

tss_gr <- GRanges(
  seqnames = gene_tbl$chr,
  ranges = IRanges(start = tss_pos, end = tss_pos),
  strand = "*",
  gene_id_clean = gene_tbl$gene_id_clean
)

tss_win_small <- symm_window(tss_gr, cfg$tss_window_bp)
tss_win_quarter <- symm_window(tss_gr, cfg$quarter_window_bp)
tss_win_mid <- symm_window(tss_gr, cfg$mid_window_bp)
tss_win_cis <- symm_window(tss_gr, cfg$cis_window_bp)

gene_body_gr <- GRanges(
  seqnames = gene_tbl$chr,
  ranges = IRanges(start = gene_tbl$start, end = gene_tbl$end),
  strand = "*",
  gene_id_clean = gene_tbl$gene_id_clean
)

gene_win_small <- expand_interval_window(gene_body_gr, cfg$tss_window_bp)
gene_win_quarter <- expand_interval_window(gene_body_gr, cfg$quarter_window_bp)
gene_win_mid <- expand_interval_window(gene_body_gr, cfg$mid_window_bp)
gene_win_cis <- expand_interval_window(gene_body_gr, cfg$cis_window_bp)

gene_window_widths <- tibble(
  gene_id_clean = mcols(gene_win_small)$gene_id_clean,
  gene_window_bp_100kb = as.numeric(width(gene_win_small)),
  gene_window_bp_250kb = as.numeric(width(gene_win_quarter)),
  gene_window_bp_500kb = as.numeric(width(gene_win_mid)),
  gene_window_bp_1mb = as.numeric(width(gene_win_cis))
)

gene_tbl <- gene_tbl %>%
  left_join(gene_window_widths, by = "gene_id_clean")
sanity_check(
  condition = all(gene_tbl$gene_window_bp_100kb > 0, na.rm = TRUE) &&
    all(gene_tbl$gene_window_bp_250kb > 0, na.rm = TRUE) &&
    all(gene_tbl$gene_window_bp_500kb > 0, na.rm = TRUE) &&
    all(gene_tbl$gene_window_bp_1mb > 0, na.rm = TRUE),
  pass_msg = "Gene window widths are positive for all configured windows.",
  fail_msg = "One or more gene window widths are non-positive."
)

# ----------------------- 4) REGULATORY FEATURE COUNTS -------------------------
message("Step 5: build enhancer and open-chromatin features")

reg_features <- gene_tbl %>%
  select(gene_id_clean) %>%
  distinct(gene_id_clean)

enh_window_obj <- build_window_enhancer_features(
  cfg = cfg,
  gene_win_small = gene_win_small,
  gene_win_quarter = gene_win_quarter,
  gene_win_mid = gene_win_mid,
  gene_win_cis = gene_win_cis,
  tss_win_small = tss_win_small,
  tss_win_quarter = tss_win_quarter,
  tss_win_mid = tss_win_mid,
  tss_win_cis = tss_win_cis,
  gene_body_gr = gene_body_gr
)
enh_count_features <- enh_window_obj$features
enh_gr <- enh_window_obj$enhancer_gr
sanity_check(
  condition = length(enh_gr) > 0L,
  pass_msg = paste0("Enhancer intervals loaded: n = ", format(length(enh_gr), big.mark = ",")),
  fail_msg = "Enhancer interval object is empty."
)
sanity_rows(enh_count_features, "Enhancer feature table", min_rows = 1000L)

if (identical(cfg$enhancer_source, "roadmap_links")) {
  message("Enhancer mode requested: roadmap_links. Active-biosample metrics disabled; using fixed-window enhancer counts.")
  reg_features <- reg_features %>%
    left_join(as_tibble(enh_count_features), by = "gene_id_clean")
} else if (identical(cfg$enhancer_source, "window_count")) {
  message("Enhancer mode: fixed-window counts around gene body and TSS.")

  reg_features <- reg_features %>%
    left_join(as_tibble(enh_count_features), by = "gene_id_clean")
} else {
  stop("Unknown cfg$enhancer_source: ", cfg$enhancer_source)
}

open_path <- cfg$open_bed
if (is.null(open_path) || open_path == "") {
  open_path <- safe_download(
    defaults$open_url,
    file.path(cfg$resource_dir, basename(defaults$open_url))
  )
}

open_gr <- import_bed_like(open_path, label = "open-chromatin annotations") %>%
  to_ucsc_style() %>%
  keep_standard_chr()

open_features <- list(
  count_overlaps_dt(gene_win_small, open_gr, "open_count_100kb") %>% as_tibble(),
  count_overlaps_dt(gene_win_quarter, open_gr, "open_count_250kb") %>% as_tibble(),
  count_overlaps_dt(gene_win_mid, open_gr, "open_count_500kb") %>% as_tibble(),
  count_overlaps_dt(gene_win_cis, open_gr, "open_count_1mb") %>% as_tibble(),
  count_overlaps_dt(tss_win_small, open_gr, "open_count_tss_100kb") %>% as_tibble(),
  count_overlaps_dt(tss_win_quarter, open_gr, "open_count_tss_250kb") %>% as_tibble(),
  count_overlaps_dt(tss_win_mid, open_gr, "open_count_tss_500kb") %>% as_tibble(),
  count_overlaps_dt(tss_win_cis, open_gr, "open_count_tss_1mb") %>% as_tibble(),
  tibble(
    gene_id_clean = mcols(gene_body_gr)$gene_id_clean,
    open_overlap_gene_body_n = as.integer(countOverlaps(gene_body_gr, open_gr, ignore.strand = TRUE))
  )
) %>%
  purrr::reduce(left_join, by = "gene_id_clean")

reg_features <- reg_features %>%
  left_join(open_features, by = "gene_id_clean")
sanity_check(
  condition = length(open_gr) > 0L,
  pass_msg = paste0("Open chromatin intervals loaded: n = ", format(length(open_gr), big.mark = ",")),
  fail_msg = "Open chromatin interval object is empty."
)
sanity_rows(reg_features, "Regulatory feature table", min_rows = 1000L)

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
sanity_rows(rmsk, "Raw rmsk table", min_rows = 10000L)

names(rmsk)[1:17] <- c(
  "bin", "swScore", "milliDiv", "milliDel", "milliIns", "genoName", "genoStart",
  "genoEnd", "genoLeft", "strand", "repName", "repClass", "repFamily",
  "repStart", "repEnd", "repLeft", "id"
)

rmsk <- rmsk %>%
  as_tibble() %>%
  mutate(
    repClass = normalize_repeat_class(repClass),
    repFamily = normalize_repeat_class(repFamily),
    repClass_norm = normalize_repeat_class(repClass),
    repFamily_norm = normalize_repeat_class(repFamily),
    rep_len_bp = as.numeric(genoEnd - genoStart),
    milliDiv = as.numeric(milliDiv)
  ) %>%
  filter(
    genoName %in% paste0("chr", c(1:22, "X", "Y")),
    !is.na(repClass),
    genoEnd > genoStart,
    is.finite(rep_len_bp),
    rep_len_bp > 0
  )
sanity_rows(rmsk, "Filtered rmsk table", min_rows = 10000L)

rep_gr <- GRanges(
  seqnames = rmsk$genoName,
  ranges = IRanges(start = rmsk$genoStart + 1L, end = rmsk$genoEnd),
  strand = "*",
  repName = rmsk$repName,
  repClass = rmsk$repClass,
  repFamily = rmsk$repFamily,
  milliDiv = rmsk$milliDiv
)

# Repeat features include gene-window burdens (primary), TSS-window burdens, and gene-body overlap summaries.
repeat_features <- data.table(
  gene_id_clean = mcols(gene_win_small)$gene_id_clean,
  repeat_count_100kb = as.integer(countOverlaps(gene_win_small, rep_gr, ignore.strand = TRUE)),
  repeat_count_250kb = as.integer(countOverlaps(gene_win_quarter, rep_gr, ignore.strand = TRUE)),
  repeat_count_500kb = as.integer(countOverlaps(gene_win_mid, rep_gr, ignore.strand = TRUE)),
  repeat_count_1mb = as.integer(countOverlaps(gene_win_cis, rep_gr, ignore.strand = TRUE)),
  repeat_count_tss_100kb = as.integer(countOverlaps(tss_win_small, rep_gr, ignore.strand = TRUE)),
  repeat_count_tss_250kb = as.integer(countOverlaps(tss_win_quarter, rep_gr, ignore.strand = TRUE)),
  repeat_count_tss_500kb = as.integer(countOverlaps(tss_win_mid, rep_gr, ignore.strand = TRUE)),
  repeat_count_tss_1mb = as.integer(countOverlaps(tss_win_cis, rep_gr, ignore.strand = TRUE)),
  repeat_overlap_gene_body_n = as.integer(countOverlaps(gene_body_gr, rep_gr, ignore.strand = TRUE))
)

gene_window_list <- list(
  "100kb" = gene_win_small,
  "250kb" = gene_win_quarter,
  "500kb" = gene_win_mid,
  "1mb" = gene_win_cis
)

repeat_filter_profiles <- build_repeat_filter_profiles()
fwrite(
  as.data.table(repeat_filter_profiles),
  file.path(cfg$output_dir, "repeat_filter_criteria.tsv"),
  sep = "\t"
)

repeat_profile_tbls <- repeat_filter_profiles %>%
  split(.$filter_set) %>%
  imap(function(profile_df, filter_set_name) {
    profile_row <- profile_df[1, , drop = FALSE]
    filtered_tbl <- apply_repeat_filter_profile(rmsk, profile_row)
    filtered_tbl %>%
      mutate(filter_set = filter_set_name)
  })

repeat_profile_stats <- repeat_profile_tbls %>%
  imap_dfr(function(tbl, filter_set_name) {
    tibble(
      filter_set = filter_set_name,
      n_elements = nrow(tbl),
      mean_len_bp = mean(tbl$rep_len_bp, na.rm = TRUE),
      median_len_bp = median(tbl$rep_len_bp, na.rm = TRUE),
      mean_milliDiv = mean(tbl$milliDiv, na.rm = TRUE),
      median_milliDiv = median(tbl$milliDiv, na.rm = TRUE)
    )
  }) %>%
  left_join(
    repeat_filter_profiles %>%
      select(filter_set, source_note),
    by = "filter_set"
  ) %>%
  arrange(desc(n_elements))
sanity_check(
  condition = nrow(repeat_profile_stats) > 0L && any(repeat_profile_stats$n_elements > 0, na.rm = TRUE),
  pass_msg = "At least one repeat filter set has non-zero elements.",
  fail_msg = "All repeat filter sets are empty; revise filter criteria."
)

fwrite(
  as.data.table(repeat_profile_stats),
  file.path(cfg$output_dir, "repeat_filter_set_counts.tsv"),
  sep = "\t"
)

repeat_profile_feature_tbls <- repeat_profile_tbls %>%
  imap(function(tbl, filter_set_name) {
    if (nrow(tbl) == 0L) {
      out <- tibble(gene_id_clean = mcols(gene_win_small)$gene_id_clean)
      for (window_label in names(gene_window_list)) {
        out <- out %>%
          mutate(
            !!paste0(
              "repeat_filt_",
              sanitize_name(filter_set_name),
              "_count_",
              window_label
            ) := 0
          )
      }
      return(out)
    }

    rep_gr_subset <- GRanges(
      seqnames = tbl$genoName,
      ranges = IRanges(start = tbl$genoStart + 1L, end = tbl$genoEnd),
      strand = "*"
    )
    build_filtered_repeat_window_counts(
      gene_windows = gene_window_list,
      rep_gr = rep_gr_subset,
      filter_set = filter_set_name
    )
  })

repeat_profile_dt <- repeat_profile_feature_tbls %>%
  purrr::reduce(left_join, by = "gene_id_clean")

repeat_profile_chr_stats <- repeat_profile_tbls %>%
  imap_dfr(function(tbl, filter_set_name) {
    if (nrow(tbl) == 0L) {
      return(tibble(
        filter_set = character(0),
        chr = character(0),
        n_elements = integer(0),
        mean_len_bp = numeric(0)
      ))
    }

    tbl %>%
      group_by(genoName) %>%
      summarise(
        n_elements = as.integer(n()),
        mean_len_bp = as.numeric(mean(rep_len_bp, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      transmute(
        filter_set = filter_set_name,
        chr = normalize_chr(genoName),
        n_elements = as.integer(n_elements),
        mean_len_bp = as.numeric(mean_len_bp)
      )
  }) %>%
  filter(chr %in% paste0("chr", c(1:22, "X", "Y")))

repeat_class_dt <- list(
  count_by_group(
    window_gr = gene_win_small,
    feature_gr = rep_gr,
    feature_group = mcols(rep_gr)$repClass,
    prefix = "repeat_class",
    suffix = "_count_100kb"
  ) %>%
    as_tibble(),
  count_by_group(
    window_gr = gene_win_quarter,
    feature_gr = rep_gr,
    feature_group = mcols(rep_gr)$repClass,
    prefix = "repeat_class",
    suffix = "_count_250kb"
  ) %>%
    as_tibble(),
  count_by_group(
    window_gr = gene_win_mid,
    feature_gr = rep_gr,
    feature_group = mcols(rep_gr)$repClass,
    prefix = "repeat_class",
    suffix = "_count_500kb"
  ) %>%
    as_tibble(),
  count_by_group(
    window_gr = gene_win_cis,
    feature_gr = rep_gr,
    feature_group = mcols(rep_gr)$repClass,
    prefix = "repeat_class",
    suffix = "_count_1mb"
  ) %>%
    as_tibble()
) %>%
  purrr::reduce(left_join, by = "gene_id_clean")

repeat_count_class_cols <- names(repeat_class_dt) %>%
  str_subset("^repeat_class_.*_count_100kb$")

if (length(repeat_count_class_cols) > 0L) {
  count_mat <- repeat_class_dt %>%
    select(all_of(repeat_count_class_cols)) %>%
    as.matrix()
  active_mat <- count_mat > 0
  active_n <- rowSums(active_mat, na.rm = TRUE)

  repeat_class_dt <- repeat_class_dt %>%
    mutate(
      repeat_active_class_n_100kb = active_n,
      repeat_mean_count_active_class_100kb = if_else(
        active_n > 0,
        rowSums(count_mat * active_mat, na.rm = TRUE) / active_n,
        0
      )
    )
} else {
  repeat_class_dt <- repeat_class_dt %>%
    mutate(
      repeat_active_class_n_100kb = 0,
      repeat_mean_count_active_class_100kb = 0
    )
}

for (cls in c("LINE", "SINE", "LTR", "DNA")) {
  count_col <- paste0("repeat_class_", cls, "_count_100kb")
  if (!count_col %in% names(repeat_class_dt)) {
    repeat_class_dt <- repeat_class_dt %>%
      mutate(!!count_col := 0)
  }
}

# Add alias columns requested by user:
# repeat_<repClass>_count from 100kb, plus window-specific aliases
# repeat_<repClass>_count_<window>.
repeat_class_count_cols <- names(repeat_class_dt) %>%
  str_subset("^repeat_class_.*_count_(100kb|250kb|500kb|1mb)$")

if (length(repeat_class_count_cols) > 0L) {
  repeat_alias_map <- tibble(src_col = repeat_class_count_cols) %>%
    mutate(
      rep_class_name = src_col %>%
        str_remove("^repeat_class_") %>%
        str_remove("_count_(100kb|250kb|500kb|1mb)$") %>%
        normalize_repeat_class() %>%
        sanitize_name(),
      window_label = src_col %>%
        str_extract("(100kb|250kb|500kb|1mb)$"),
      alias_col_window = paste0("repeat_", rep_class_name, "_count_", window_label)
    ) %>%
    filter(!is.na(rep_class_name), rep_class_name != "")

  alias_groups <- split(repeat_alias_map$src_col, repeat_alias_map$alias_col_window)
  for (alias_col_window in names(alias_groups)) {
    src_cols <- alias_groups[[alias_col_window]]
    repeat_class_dt <- repeat_class_dt %>%
      mutate(
        !!alias_col_window := rowSums(
          as.matrix(select(., all_of(src_cols))),
          na.rm = TRUE
        )
      )
  }

  base_types <- repeat_alias_map$rep_class_name %>% unique()
  for (rep_type in base_types) {
    for (window_label in c("100kb", "250kb", "500kb", "1mb")) {
      alias_col_window <- paste0("repeat_", rep_type, "_count_", window_label)
      if (!alias_col_window %in% names(repeat_class_dt)) {
        repeat_class_dt <- repeat_class_dt %>%
          mutate(!!alias_col_window := 0)
      }
    }

    alias_col <- paste0("repeat_", rep_type, "_count")
    alias_100kb <- paste0("repeat_", rep_type, "_count_100kb")
    if (!alias_col %in% names(repeat_class_dt) && alias_100kb %in% names(repeat_class_dt)) {
      repeat_class_dt <- repeat_class_dt %>%
        mutate(!!alias_col := .data[[alias_100kb]])
    }
  }
}

repeat_family_dt <- list(
  count_by_group(
    window_gr = gene_win_small,
    feature_gr = rep_gr,
    feature_group = mcols(rep_gr)$repFamily,
    prefix = "repeat_family",
    suffix = "_count_100kb"
  ) %>%
    as_tibble(),
  count_by_group(
    window_gr = gene_win_quarter,
    feature_gr = rep_gr,
    feature_group = mcols(rep_gr)$repFamily,
    prefix = "repeat_family",
    suffix = "_count_250kb"
  ) %>%
    as_tibble(),
  count_by_group(
    window_gr = gene_win_mid,
    feature_gr = rep_gr,
    feature_group = mcols(rep_gr)$repFamily,
    prefix = "repeat_family",
    suffix = "_count_500kb"
  ) %>%
    as_tibble(),
  count_by_group(
    window_gr = gene_win_cis,
    feature_gr = rep_gr,
    feature_group = mcols(rep_gr)$repFamily,
    prefix = "repeat_family",
    suffix = "_count_1mb"
  ) %>%
    as_tibble()
) %>%
  purrr::reduce(left_join, by = "gene_id_clean")

# --------------------------- 6) FINAL ANALYSIS TABLE --------------------------
analysis_dt <- merged %>%
  left_join(
    gene_tbl %>%
      select(
        gene_id_clean,
        gene_name,
        chr,
        start,
        end,
        gene_length,
        gene_window_bp_100kb,
        gene_window_bp_250kb,
        gene_window_bp_500kb,
        gene_window_bp_1mb
      ),
    by = "gene_id_clean"
  ) %>%
  left_join(reg_features, by = "gene_id_clean") %>%
  left_join(as_tibble(repeat_features), by = "gene_id_clean") %>%
  left_join(repeat_profile_dt, by = "gene_id_clean") %>%
  left_join(repeat_class_dt, by = "gene_id_clean") %>%
  left_join(repeat_family_dt, by = "gene_id_clean")

count_cols <- names(analysis_dt) %>%
  str_subset("^(enh_count_|enh_overlap_|open_count_|open_overlap_|repeat_count_|repeat_overlap_|repeat_class_|repeat_family_|repeat_filt_|repeat_[A-Za-z0-9_]+_count(_(100kb|250kb|500kb|1mb))?$)")

if (length(count_cols) > 0L) {
  analysis_dt <- analysis_dt %>%
    mutate(
      across(
        all_of(count_cols),
        ~ replace_na(as.numeric(.x), 0)
      )
    )
}
sanity_check(
  condition = nrow(analysis_dt) == nrow(merged),
  pass_msg = paste0("Final analysis table rows match merged genes: n = ", nrow(analysis_dt)),
  fail_msg = "Final analysis table row count does not match merged gene table."
)

fwrite(as.data.table(analysis_dt), file.path(cfg$output_dir, "gene_level_features.tsv"), sep = "\t")

# ------- 6b) FILTERED-REPEAT DECILE SUMMARIES + GENOMIC BACKGROUND ----------
message("Step 6b: build filtered-repeat summaries and genomic-background simulations")

window_label_map <- c(
  "100kb" = "100 kb",
  "250kb" = "250 kb",
  "500kb" = "500 kb",
  "1mb" = "1 Mb"
)
window_width_col_map <- c(
  "100kb" = "gene_window_bp_100kb",
  "250kb" = "gene_window_bp_250kb",
  "500kb" = "gene_window_bp_500kb",
  "1mb" = "gene_window_bp_1mb"
)

repeat_filter_map <- repeat_filter_profiles %>%
  mutate(filter_set_clean = sanitize_name(filter_set))

repeat_filt_cols <- names(analysis_dt) %>%
  str_subset("^repeat_filt_.*_count_(100kb|250kb|500kb|1mb)$")

repeat_filt_observed_summary <- tibble(
  filter_set = character(0),
  filter_set_clean = character(0),
  window = character(0),
  post_mean_bin = integer(0),
  observed_mean = numeric(0),
  observed_median = numeric(0),
  observed_prop_nonzero = numeric(0),
  n_genes = integer(0)
)

repeat_background_iter <- tibble(
  filter_set = character(0),
  filter_set_clean = character(0),
  window = character(0),
  seed_used = integer(0),
  iter = integer(0),
  post_mean_bin = integer(0),
  bg_mean = numeric(0),
  bg_median = numeric(0),
  bg_prop_nonzero = numeric(0)
)

repeat_bg_seed_map <- tibble(
  filter_set = character(0),
  filter_set_clean = character(0),
  window = character(0),
  seed_used = integer(0),
  n_iter = integer(0),
  seed_base = integer(0),
  seed_strategy = character(0)
)

if (length(repeat_filt_cols) > 0L) {
  repeat_filt_observed_summary <- analysis_dt %>%
    as_tibble() %>%
    select(post_mean_bin, all_of(repeat_filt_cols)) %>%
    pivot_longer(
      cols = -post_mean_bin,
      names_to = "feature_col",
      values_to = "value"
    ) %>%
    mutate(
      filter_set_clean = feature_col %>%
        str_remove("^repeat_filt_") %>%
        str_remove("_count_(100kb|250kb|500kb|1mb)$"),
      window = feature_col %>%
        str_extract("(100kb|250kb|500kb|1mb)$")
    ) %>%
    left_join(
      repeat_filter_map %>% select(filter_set, filter_set_clean),
      by = "filter_set_clean"
    ) %>%
    group_by(filter_set, filter_set_clean, window, post_mean_bin) %>%
    summarise(
      observed_mean = mean(value, na.rm = TRUE),
      observed_median = median(value, na.rm = TRUE),
      observed_prop_nonzero = mean(value > 0, na.rm = TRUE),
      n_genes = n(),
      .groups = "drop"
    ) %>%
    arrange(filter_set, window, post_mean_bin)

  fwrite(
    as.data.table(repeat_filt_observed_summary),
    file.path(cfg$output_dir, "postmean_repeat_filtered_observed_summary.tsv"),
    sep = "\t"
  )

  chrom_sizes_dt <- load_chrom_sizes(cfg, rmsk)
  sanity_rows(chrom_sizes_dt, "Chromosome size table", min_rows = 20L)
  fwrite(
    as.data.table(chrom_sizes_dt),
    file.path(cfg$output_dir, "chrom_sizes_used.tsv"),
    sep = "\t"
  )

  if (isTRUE(as.integer(cfg$repeat_bg_n_iter) > 0L)) {
    bg_method <- as.character(cfg$repeat_bg_method)
    if (!bg_method %in% c("explicit_genome", "poisson_window")) {
      stop("Unknown cfg$repeat_bg_method: ", bg_method)
    }

    message(
      "Running ",
      cfg$repeat_bg_n_iter,
      " genomic-background simulations per filter set and window (method = ",
      bg_method,
      ")."
    )

    repeat_bg_seed_info <- tibble(
      parameter = c(
        "repeat_bg_method",
        "repeat_bg_seed_base",
        "repeat_bg_n_iter",
        "repeat_bg_repeat_fraction",
        "repeat_bg_type_col",
        "repeat_bg_composition_basis",
        "seed_strategy"
      ),
      value = c(
        bg_method,
        as.character(cfg$repeat_bg_seed),
        as.character(cfg$repeat_bg_n_iter),
        as.character(cfg$repeat_bg_repeat_fraction),
        as.character(cfg$repeat_bg_type_col),
        as.character(cfg$repeat_bg_composition_basis),
        "seed_used = repeat_bg_seed_base + simulation_index_in_filter_window_loop"
      )
    )
    fwrite(
      as.data.table(repeat_bg_seed_info),
      file.path(cfg$output_dir, "repeat_background_seed_info.tsv"),
      sep = "\t"
    )

    bg_iter_tbls <- list()
    seed_map_tbls <- list()
    sim_idx <- 0L

    gene_decile_map <- analysis_dt %>%
      as_tibble() %>%
      select(gene_id_clean, post_mean_bin) %>%
      distinct(gene_id_clean, .keep_all = TRUE)

    genome_bp <- sum(chrom_sizes_dt$chrom_size, na.rm = TRUE)
    filter_bp_dt <- repeat_profile_tbls %>%
      imap_dfr(function(tbl, filter_set_name) {
        tbl <- as_tibble(tbl)
        bp <- sum(as.numeric(tbl$rep_len_bp), na.rm = TRUE)
        type_n <- if (cfg$repeat_bg_type_col %in% names(tbl)) {
          n_distinct(tbl[[cfg$repeat_bg_type_col]])
        } else {
          NA_integer_
        }
        tibble(
          filter_set = filter_set_name,
          observed_bp = as.numeric(bp),
          observed_fraction = if_else(genome_bp > 0, as.numeric(bp / genome_bp), 0),
          observed_n = as.integer(nrow(tbl)),
          n_types = as.integer(type_n)
        )
      })

    reference_filter_set <- if ("te_core_interspersed" %in% filter_bp_dt$filter_set) {
      "te_core_interspersed"
    } else {
      filter_bp_dt %>%
        arrange(desc(observed_fraction)) %>%
        dplyr::slice_head(n = 1) %>%
        pull(filter_set)
    }
    reference_fraction <- filter_bp_dt %>%
      filter(filter_set == reference_filter_set) %>%
      pull(observed_fraction)
    if (length(reference_fraction) == 0L || !is.finite(reference_fraction[[1]]) || reference_fraction[[1]] <= 0) {
      reference_fraction <- filter_bp_dt$observed_fraction %>% max(na.rm = TRUE)
    } else {
      reference_fraction <- reference_fraction[[1]]
    }

    sim_plan <- filter_bp_dt %>%
      mutate(
        reference_filter_set = reference_filter_set,
        target_fraction = pmin(
          0.95,
          if_else(
            is.finite(reference_fraction) & reference_fraction > 0,
            as.numeric(cfg$repeat_bg_repeat_fraction) * (observed_fraction / reference_fraction),
            0
          )
        ),
        repeat_bg_method = bg_method,
        type_col = as.character(cfg$repeat_bg_type_col),
        composition_basis = as.character(cfg$repeat_bg_composition_basis)
      )

    fwrite(
      as.data.table(sim_plan),
      file.path(cfg$output_dir, "repeat_background_simulation_plan.tsv"),
      sep = "\t"
    )

    if (identical(bg_method, "explicit_genome")) {
      for (row_idx in seq_len(nrow(repeat_filter_map))) {
        filter_set_name <- repeat_filter_map$filter_set[[row_idx]]
        filter_clean <- repeat_filter_map$filter_set_clean[[row_idx]]
        source_tbl <- repeat_profile_tbls[[filter_set_name]] %>%
          as_tibble()

        plan_row <- sim_plan %>%
          filter(filter_set == filter_set_name) %>%
          dplyr::slice_head(n = 1)
        if (nrow(plan_row) == 0L) {
          next
        }

        target_fraction_set <- as.numeric(plan_row$target_fraction[[1]])
        if (!is.finite(target_fraction_set) || target_fraction_set <= 0 || nrow(source_tbl) == 0L) {
          next
        }

        catalog <- build_repeat_sim_catalog(
          repeat_tbl = source_tbl,
          chrom_sizes_dt = chrom_sizes_dt,
          type_col = as.character(cfg$repeat_bg_type_col),
          composition_basis = as.character(cfg$repeat_bg_composition_basis)
        )
        if (catalog$source_n == 0L) {
          next
        }

        for (iter_idx in seq_len(as.integer(cfg$repeat_bg_n_iter))) {
          sim_idx <- sim_idx + 1L
          seed_used <- as.integer(cfg$repeat_bg_seed) + sim_idx

          sim_dt <- simulate_repeat_genome_dt(
            catalog = catalog,
            chrom_sizes_dt = chrom_sizes_dt,
            repeat_fraction = target_fraction_set,
            seed = as.integer(seed_used)
          )
          if (nrow(sim_dt) == 0L) {
            next
          }

          sim_gr <- GRanges(
            seqnames = sim_dt$chr,
            ranges = IRanges(start = sim_dt$start, end = sim_dt$end),
            strand = "*"
          )

          seed_map_tbls[[length(seed_map_tbls) + 1L]] <- tibble(
            filter_set = filter_set_name,
            filter_set_clean = filter_clean,
            window = "all_windows",
            seed_used = as.integer(seed_used),
            n_iter = as.integer(cfg$repeat_bg_n_iter),
            seed_base = as.integer(cfg$repeat_bg_seed),
            seed_strategy = "seed_used = repeat_bg_seed_base + simulation_index_in_filter_window_loop",
            repeat_bg_method = bg_method,
            target_fraction = as.numeric(target_fraction_set),
            simulated_elements = as.integer(nrow(sim_dt)),
            simulated_bp = as.numeric(sum(sim_dt$rep_len_bp, na.rm = TRUE))
          )

          for (window_label in names(window_width_col_map)) {
            window_gr <- gene_window_list[[window_label]]
            bg_decile <- tibble(
              gene_id_clean = mcols(window_gr)$gene_id_clean,
              bg_count = as.integer(countOverlaps(window_gr, sim_gr, ignore.strand = TRUE))
            ) %>%
              left_join(gene_decile_map, by = "gene_id_clean") %>%
              filter(!is.na(post_mean_bin)) %>%
              group_by(post_mean_bin) %>%
              summarise(
                bg_mean = mean(bg_count, na.rm = TRUE),
                bg_median = median(bg_count, na.rm = TRUE),
                bg_prop_nonzero = mean(bg_count > 0, na.rm = TRUE),
                .groups = "drop"
              ) %>%
              mutate(
                filter_set = filter_set_name,
                filter_set_clean = filter_clean,
                window = window_label,
                seed_used = as.integer(seed_used),
                iter = as.integer(iter_idx)
              )

            bg_iter_tbls[[length(bg_iter_tbls) + 1L]] <- bg_decile
          }
        }
      }
    } else {
      for (row_idx in seq_len(nrow(repeat_filter_map))) {
        filter_set_name <- repeat_filter_map$filter_set[[row_idx]]
        filter_clean <- repeat_filter_map$filter_set_clean[[row_idx]]

        chr_stats <- repeat_profile_chr_stats %>%
          filter(filter_set == filter_set_name) %>%
          select(chr, n_elements, mean_len_bp)

        for (window_label in names(window_width_col_map)) {
          obs_col <- paste0("repeat_filt_", filter_clean, "_count_", window_label)
          width_col <- window_width_col_map[[window_label]]

          if (!obs_col %in% names(analysis_dt) || !width_col %in% names(analysis_dt)) {
            next
          }

          sim_input <- analysis_dt %>%
            as_tibble() %>%
            transmute(
              post_mean_bin = post_mean_bin,
              chr = chr,
              gene_window_bp = as.numeric(.data[[width_col]])
            ) %>%
            left_join(chrom_sizes_dt, by = "chr") %>%
            left_join(chr_stats, by = "chr") %>%
            mutate(
              chrom_size = as.numeric(chrom_size),
              n_elements = replace_na(as.numeric(n_elements), 0),
              mean_len_bp = replace_na(as.numeric(mean_len_bp), 0),
              gene_window_bp = replace_na(as.numeric(gene_window_bp), 0),
              effective_bp = pmin(chrom_size, pmax(gene_window_bp + mean_len_bp, 0)),
              lambda = if_else(
                is.finite(chrom_size) & chrom_size > 0,
                (n_elements / chrom_size) * effective_bp,
                0
              ),
              lambda = pmax(replace_na(lambda, 0), 0)
            ) %>%
            filter(!is.na(post_mean_bin))

          sim_idx <- sim_idx + 1L
          seed_used <- as.integer(cfg$repeat_bg_seed) + sim_idx

          seed_map_tbls[[length(seed_map_tbls) + 1L]] <- tibble(
            filter_set = filter_set_name,
            filter_set_clean = filter_clean,
            window = window_label,
            seed_used = as.integer(seed_used),
            n_iter = as.integer(cfg$repeat_bg_n_iter),
            seed_base = as.integer(cfg$repeat_bg_seed),
            seed_strategy = "seed_used = repeat_bg_seed_base + simulation_index_in_filter_window_loop",
            repeat_bg_method = bg_method,
            target_fraction = NA_real_,
            simulated_elements = NA_integer_,
            simulated_bp = NA_real_
          )

          bg_iter <- simulate_poisson_background(
            sim_input = sim_input,
            n_iter = as.integer(cfg$repeat_bg_n_iter),
            seed = as.integer(seed_used)
          ) %>%
            mutate(
              filter_set = filter_set_name,
              filter_set_clean = filter_clean,
              window = window_label,
              seed_used = as.integer(seed_used)
            )

          bg_iter_tbls[[length(bg_iter_tbls) + 1L]] <- bg_iter
        }
      }
    }

    if (length(seed_map_tbls) > 0L) {
      repeat_bg_seed_map <- bind_rows(seed_map_tbls) %>%
        arrange(filter_set, window, seed_used)

      sanity_check(
        condition = n_distinct(repeat_bg_seed_map$seed_used) == nrow(repeat_bg_seed_map),
        pass_msg = "Background seeds are unique across simulation runs.",
        fail_msg = "Duplicate background seeds detected across simulation runs."
      )

      fwrite(
        as.data.table(repeat_bg_seed_map),
        file.path(cfg$output_dir, "repeat_background_seed_map.tsv"),
        sep = "\t"
      )
    }

    if (length(bg_iter_tbls) > 0L) {
      repeat_background_iter <- bind_rows(bg_iter_tbls) %>%
        arrange(filter_set, window, seed_used, iter, post_mean_bin)

      fwrite(
        as.data.table(repeat_background_iter),
        file.path(cfg$output_dir, "repeat_filtered_background_iter.tsv"),
        sep = "\t"
      )
    }
  }
}

repeat_background_summary <- tibble(
  filter_set = character(0),
  filter_set_clean = character(0),
  window = character(0),
  post_mean_bin = integer(0),
  bg_mean_mean = numeric(0),
  bg_mean_q025 = numeric(0),
  bg_mean_q975 = numeric(0),
  bg_median_mean = numeric(0),
  bg_median_q025 = numeric(0),
  bg_median_q975 = numeric(0),
  bg_prop_nonzero_mean = numeric(0),
  bg_prop_nonzero_q025 = numeric(0),
  bg_prop_nonzero_q975 = numeric(0)
)

repeat_filtered_enrichment <- tibble(
  filter_set = character(0),
  filter_set_clean = character(0),
  window = character(0),
  post_mean_bin = integer(0),
  observed_mean = numeric(0),
  observed_median = numeric(0),
  observed_prop_nonzero = numeric(0),
  bg_mean_mean = numeric(0),
  bg_mean_q025 = numeric(0),
  bg_mean_q975 = numeric(0),
  bg_median_mean = numeric(0),
  bg_median_q025 = numeric(0),
  bg_median_q975 = numeric(0),
  bg_prop_nonzero_mean = numeric(0),
  bg_prop_nonzero_q025 = numeric(0),
  bg_prop_nonzero_q975 = numeric(0),
  enrichment_ratio_median = numeric(0),
  enrichment_ratio_mean = numeric(0)
)

if (nrow(repeat_background_iter) > 0L && nrow(repeat_filt_observed_summary) > 0L) {
  repeat_background_summary <- repeat_background_iter %>%
    group_by(filter_set, filter_set_clean, window, post_mean_bin) %>%
    summarise(
      bg_mean_mean = mean(bg_mean, na.rm = TRUE),
      bg_mean_q025 = quantile(bg_mean, 0.025, na.rm = TRUE, names = FALSE),
      bg_mean_q975 = quantile(bg_mean, 0.975, na.rm = TRUE, names = FALSE),
      bg_median_mean = mean(bg_median, na.rm = TRUE),
      bg_median_q025 = quantile(bg_median, 0.025, na.rm = TRUE, names = FALSE),
      bg_median_q975 = quantile(bg_median, 0.975, na.rm = TRUE, names = FALSE),
      bg_prop_nonzero_mean = mean(bg_prop_nonzero, na.rm = TRUE),
      bg_prop_nonzero_q025 = quantile(bg_prop_nonzero, 0.025, na.rm = TRUE, names = FALSE),
      bg_prop_nonzero_q975 = quantile(bg_prop_nonzero, 0.975, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )

  repeat_filtered_enrichment <- repeat_filt_observed_summary %>%
    left_join(
      repeat_background_summary,
      by = c("filter_set", "filter_set_clean", "window", "post_mean_bin")
    ) %>%
    mutate(
      enrichment_ratio_median = if_else(
        is.finite(bg_median_mean) & bg_median_mean > 0,
        observed_median / bg_median_mean,
        NA_real_
      ),
      enrichment_ratio_mean = if_else(
        is.finite(bg_mean_mean) & bg_mean_mean > 0,
        observed_mean / bg_mean_mean,
        NA_real_
      )
    ) %>%
    arrange(filter_set, window, post_mean_bin)

  fwrite(
    as.data.table(repeat_background_summary),
    file.path(cfg$output_dir, "postmean_repeat_filtered_background_summary.tsv"),
    sep = "\t"
  )
  fwrite(
    as.data.table(repeat_filtered_enrichment),
    file.path(cfg$output_dir, "postmean_repeat_filtered_enrichment.tsv"),
    sep = "\t"
  )
  sanity_rows(repeat_filtered_enrichment, "Filtered repeat enrichment summary", min_rows = 10L)
}

# ----------------------------- 7) SUMMARIES -----------------------------------
expected_feature_cols <- c(
  "enh_count_100kb",
  "enh_count_250kb",
  "enh_count_500kb",
  "enh_count_1mb",
  "enh_count_tss_100kb",
  "enh_count_tss_250kb",
  "enh_count_tss_500kb",
  "enh_count_tss_1mb",
  "enh_overlap_gene_body_n",
  "open_count_100kb",
  "open_count_250kb",
  "open_count_500kb",
  "open_count_1mb",
  "open_count_tss_100kb",
  "open_count_tss_250kb",
  "open_count_tss_500kb",
  "open_count_tss_1mb",
  "open_overlap_gene_body_n",
  "repeat_count_100kb",
  "repeat_count_250kb",
  "repeat_count_500kb",
  "repeat_count_1mb",
  "repeat_count_tss_100kb",
  "repeat_count_tss_250kb",
  "repeat_count_tss_500kb",
  "repeat_count_tss_1mb",
  "repeat_overlap_gene_body_n",
  "repeat_active_class_n_100kb",
  "repeat_mean_count_active_class_100kb",
  "repeat_class_LINE_count_100kb",
  "repeat_class_SINE_count_100kb"
)
analysis_dt <- add_missing_numeric_cols(analysis_dt, expected_feature_cols, default_value = 0)

decile_summary <- analysis_dt %>%
  group_by(post_mean_bin) %>%
  summarise(
    n_genes = n(),
    median_h2 = median(h2_GREML, na.rm = TRUE),
    prop_h2_sig = mean(Pval_GREML < 0.05, na.rm = TRUE),
    median_enh_100kb = median(enh_count_100kb, na.rm = TRUE),
    median_enh_250kb = median(enh_count_250kb, na.rm = TRUE),
    median_enh_500kb = median(enh_count_500kb, na.rm = TRUE),
    median_enh_1mb = median(enh_count_1mb, na.rm = TRUE),
    median_enh_overlap_gene_body_n = median(enh_overlap_gene_body_n, na.rm = TRUE),
    median_open_100kb = median(open_count_100kb, na.rm = TRUE),
    median_open_250kb = median(open_count_250kb, na.rm = TRUE),
    median_open_500kb = median(open_count_500kb, na.rm = TRUE),
    median_open_1mb = median(open_count_1mb, na.rm = TRUE),
    median_open_overlap_gene_body_n = median(open_overlap_gene_body_n, na.rm = TRUE),
    median_repeat_100kb = median(repeat_count_100kb, na.rm = TRUE),
    median_repeat_250kb = median(repeat_count_250kb, na.rm = TRUE),
    median_repeat_500kb = median(repeat_count_500kb, na.rm = TRUE),
    median_repeat_1mb = median(repeat_count_1mb, na.rm = TRUE),
    median_repeat_overlap_gene_body_n = median(repeat_overlap_gene_body_n, na.rm = TRUE),
    median_repeat_active_class_n_100kb = median(repeat_active_class_n_100kb, na.rm = TRUE),
    median_repeat_LINE_count_100kb = median(repeat_class_LINE_count_100kb, na.rm = TRUE),
    median_repeat_SINE_count_100kb = median(repeat_class_SINE_count_100kb, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(post_mean_bin)

fwrite(as.data.table(decile_summary), file.path(cfg$output_dir, "decile_feature_summary.tsv"), sep = "\t")

cor_dt <- rbindlist(list(
  safe_spearman(analysis_dt$post_mean, analysis_dt$h2_GREML, "post_mean", "h2_GREML"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$enh_count_100kb, "post_mean", "enh_count_100kb"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$enh_count_250kb, "post_mean", "enh_count_250kb"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$enh_count_500kb, "post_mean", "enh_count_500kb"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$enh_overlap_gene_body_n, "post_mean", "enh_overlap_gene_body_n"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$open_count_100kb, "post_mean", "open_count_100kb"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$open_count_250kb, "post_mean", "open_count_250kb"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$open_count_500kb, "post_mean", "open_count_500kb"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$open_overlap_gene_body_n, "post_mean", "open_overlap_gene_body_n"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$repeat_count_100kb, "post_mean", "repeat_count_100kb"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$repeat_count_250kb, "post_mean", "repeat_count_250kb"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$repeat_count_500kb, "post_mean", "repeat_count_500kb"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$repeat_overlap_gene_body_n, "post_mean", "repeat_overlap_gene_body_n"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$repeat_class_LINE_count_100kb, "post_mean", "repeat_class_LINE_count_100kb"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$repeat_class_SINE_count_100kb, "post_mean", "repeat_class_SINE_count_100kb"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$enh_count_100kb, "h2_GREML", "enh_count_100kb"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$enh_count_250kb, "h2_GREML", "enh_count_250kb"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$enh_count_500kb, "h2_GREML", "enh_count_500kb"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$enh_overlap_gene_body_n, "h2_GREML", "enh_overlap_gene_body_n"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$open_count_100kb, "h2_GREML", "open_count_100kb"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$open_count_250kb, "h2_GREML", "open_count_250kb"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$open_count_500kb, "h2_GREML", "open_count_500kb"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$open_overlap_gene_body_n, "h2_GREML", "open_overlap_gene_body_n"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$repeat_count_100kb, "h2_GREML", "repeat_count_100kb"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$repeat_count_250kb, "h2_GREML", "repeat_count_250kb"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$repeat_count_500kb, "h2_GREML", "repeat_count_500kb"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$repeat_overlap_gene_body_n, "h2_GREML", "repeat_overlap_gene_body_n"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$repeat_class_LINE_count_100kb, "h2_GREML", "repeat_class_LINE_count_100kb"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$repeat_class_SINE_count_100kb, "h2_GREML", "repeat_class_SINE_count_100kb")
))

fwrite(cor_dt, file.path(cfg$output_dir, "spearman_correlations.tsv"), sep = "\t")

overlap_cols <- c(
  "enh_overlap_gene_body_n",
  "open_overlap_gene_body_n",
  "repeat_overlap_gene_body_n"
)
overlap_cols <- overlap_cols[overlap_cols %in% names(analysis_dt)]

if (length(overlap_cols) > 0L) {
  overlap_summary <- tibble(feature = overlap_cols) %>%
    mutate(
      n_genes_with_overlap = vapply(feature, function(col) sum(analysis_dt[[col]] > 0, na.rm = TRUE), numeric(1)),
      prop_genes_with_overlap = vapply(feature, function(col) mean(analysis_dt[[col]] > 0, na.rm = TRUE), numeric(1)),
      median_overlap_count = vapply(feature, function(col) median(analysis_dt[[col]], na.rm = TRUE), numeric(1))
    )
  fwrite(as.data.table(overlap_summary), file.path(cfg$output_dir, "gene_body_overlap_summary.tsv"), sep = "\t")
}

# --------------------------- 8) ASSOCIATION MODELS ----------------------------
model_dt <- analysis_dt %>%
  filter(
    is.finite(h2_GREML),
    is.finite(post_mean),
    is.finite(mean_tpm),
    is.finite(gene_length)
  ) %>%
  mutate(h2_sig = as.integer(Pval_GREML < 0.05))

post_mean_term <- NULL
post_mean_levels <- model_dt %>%
  pull(post_mean_bin) %>%
  discard(is.na) %>%
  dplyr::n_distinct()
if (post_mean_levels >= 2L) {
  post_mean_term <- "factor(post_mean_bin)"
} else {
  message("Skipping factor(post_mean_bin): fewer than 2 levels in model_dt.")
}

primary_numeric_cols <- c(
  "enh_count_100kb",
  "enh_count_250kb",
  "enh_count_500kb",
  "enh_count_1mb",
  "enh_overlap_gene_body_n",
  "open_count_100kb",
  "open_count_250kb",
  "open_count_500kb",
  "open_count_1mb",
  "open_overlap_gene_body_n",
  "repeat_count_100kb",
  "repeat_count_250kb",
  "repeat_count_500kb",
  "repeat_count_1mb",
  "repeat_overlap_gene_body_n",
  "repeat_class_LINE_count_100kb",
  "repeat_class_SINE_count_100kb",
  "mean_tpm",
  "gene_length"
)

numeric_terms <- build_log_terms(model_dt, primary_numeric_cols)
if (length(numeric_terms) == 0L) {
  stop("No variable numeric predictors available for lm/glm models.")
}

lm_terms <- c(
  numeric_terms,
  post_mean_term
) %>%
  discard(~ is.na(.x) || !nzchar(.x))

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
  numeric_terms,
  post_mean_term
) %>%
  discard(~ is.na(.x) || !nzchar(.x))

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
  str_subset("^repeat_class_.*_count_100kb$")

if (length(repeat_cols) > 0L) {
  repeat_stats <- tibble(repeat_col = repeat_cols) %>%
    mutate(
      total = map_dbl(repeat_col, ~ sum(model_dt[[.x]], na.rm = TRUE)),
      n_unique = map_int(repeat_col, ~ n_distinct(model_dt[[.x]][!is.na(model_dt[[.x]])]))
    )

  top_repeat_cols <- repeat_stats %>%
    filter(total > 0, n_unique > 1) %>%
    arrange(desc(total)) %>%
    dplyr::slice_head(n = 6) %>%
    pull(repeat_col)
  top_repeat_cols <- top_repeat_cols[!is.na(top_repeat_cols) & nzchar(top_repeat_cols)]
  top_repeat_cols <- top_repeat_cols[top_repeat_cols %in% names(model_dt)]

  if (length(top_repeat_cols) == 0L) {
    message("Skipping repeat-class lm: no valid repeat_class predictors with variation.")
  } else {
    repeat_terms <- c(
      sprintf("scale(log1p(`%s`))", top_repeat_cols),
      "scale(log1p(mean_tpm))",
      "scale(log1p(gene_length))",
      post_mean_term
    ) %>%
      discard(~ is.na(.x) || !nzchar(.x))

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

analysis_tbl <- analysis_dt %>%
  as_tibble()

postmean_feature_cols <- names(analysis_tbl) %>%
  str_subset("^(enh_count_|enh_overlap_|open_count_|open_overlap_|repeat_count_|repeat_overlap_|repeat_[A-Za-z0-9_]+_count$)")

postmean_feature_cols <- postmean_feature_cols[postmean_feature_cols %in% names(analysis_tbl)]
postmean_feature_cols <- setdiff(postmean_feature_cols, "post_mean_bin")

if (length(postmean_feature_cols) > 0L) {
  postmean_feature_stats <- analysis_tbl %>%
    select(post_mean_bin, all_of(postmean_feature_cols)) %>%
    pivot_longer(
      cols = -post_mean_bin,
      names_to = "feature",
      values_to = "value"
    ) %>%
    group_by(post_mean_bin, feature) %>%
    summarise(
      n_genes = n(),
      mean_value = mean(value, na.rm = TRUE),
      median_value = median(value, na.rm = TRUE),
      q25_value = quantile(value, 0.25, na.rm = TRUE, names = FALSE),
      q75_value = quantile(value, 0.75, na.rm = TRUE, names = FALSE),
      prop_nonzero = mean(value > 0, na.rm = TRUE),
      .groups = "drop"
    )

  fwrite(
    as.data.table(postmean_feature_stats),
    file.path(cfg$output_dir, "postmean_bin_feature_summary.tsv"),
    sep = "\t"
  )

  postmean_feature_trend_long <- postmean_feature_stats %>%
    select(post_mean_bin, feature, mean_value, median_value) %>%
    pivot_longer(
      cols = c(mean_value, median_value),
      names_to = "stat",
      values_to = "value"
    ) %>%
    mutate(
      feature = factor(feature, levels = ordered_feature_levels(unique(feature))),
      stat = recode(
        stat,
        mean_value = "Mean",
        median_value = "Median"
      )
    )

  p_postmean_all <- ggplot(
    postmean_feature_trend_long,
    aes(x = post_mean_bin, y = value, color = stat, group = stat)
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    scale_x_continuous(breaks = 1:10) +
    scale_color_manual(values = c("Mean" = "#457b9d", "Median" = "#e76f51")) +
    facet_wrap(~feature, scales = "free_y", ncol = 4) +
    labs(
      title = "Feature Trends Across s_het post_mean Deciles",
      subtitle = "Mean and median by decile for all count-based features",
      x = "s_het post_mean decile",
      y = "Feature value",
      color = "Statistic"
    ) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "top")

  ggsave(
    filename = file.path(cfg$output_dir, "plots", "postmean_all_feature_trends_mean_median.pdf"),
    plot = p_postmean_all,
    width = 16,
    height = 12
  )
}

reg_decile_cols <- c(
  "enh_count_100kb",
  "enh_count_250kb",
  "enh_count_500kb",
  "enh_count_1mb",
  "enh_overlap_gene_body_n",
  "open_count_100kb",
  "open_count_250kb",
  "open_count_500kb",
  "open_count_1mb",
  "open_overlap_gene_body_n",
  "repeat_count_100kb",
  "repeat_count_250kb",
  "repeat_count_500kb",
  "repeat_count_1mb",
  "repeat_overlap_gene_body_n"
)
reg_decile_cols <- reg_decile_cols[reg_decile_cols %in% names(analysis_tbl)]

if (length(reg_decile_cols) > 0L) {
  reg_postmean_stats <- analysis_tbl %>%
    select(post_mean_bin, all_of(reg_decile_cols)) %>%
    pivot_longer(
      cols = -post_mean_bin,
      names_to = "feature",
      values_to = "value"
    ) %>%
    group_by(post_mean_bin, feature) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      median_value = median(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(mean_value, median_value),
      names_to = "stat",
      values_to = "value"
    ) %>%
    mutate(
      feature = factor(feature, levels = reg_decile_cols),
      stat = recode(
        stat,
        mean_value = "Mean",
        median_value = "Median"
      )
    )

  p_reg_postmean <- ggplot(
    reg_postmean_stats,
    aes(x = post_mean_bin, y = value, color = stat, group = stat)
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    scale_x_continuous(breaks = 1:10) +
    scale_color_manual(values = c("Mean" = "#2a9d8f", "Median" = "#e76f51")) +
    facet_wrap(~feature, scales = "free_y", ncol = 3) +
    labs(
      title = "Regulatory and Repeat Burden vs s_het Decile",
      subtitle = "Mean and median feature values by post_mean bin",
      x = "s_het post_mean decile",
      y = "Feature value",
      color = "Statistic"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top")

  ggsave(
    filename = file.path(cfg$output_dir, "plots", "postmean_regulatory_repeat_trends_mean_median.pdf"),
    plot = p_reg_postmean,
    width = 13,
    height = 9
  )
}

repeat_type_window_cols <- names(analysis_tbl) %>%
  str_subset("^repeat_[A-Za-z0-9_]+_count_(100kb|250kb|500kb|1mb)$")
repeat_type_window_cols <- setdiff(
  repeat_type_window_cols,
  names(analysis_tbl) %>% str_subset("^repeat_class_|^repeat_filt_")
)

if (length(repeat_type_window_cols) > 0L) {
  repeat_type_postmean_stats <- analysis_tbl %>%
    select(post_mean_bin, all_of(repeat_type_window_cols)) %>%
    pivot_longer(
      cols = -post_mean_bin,
      names_to = "feature_col",
      values_to = "value"
    ) %>%
    mutate(
      repeat_type = feature_col %>%
        str_remove("^repeat_") %>%
        str_remove("_count_(100kb|250kb|500kb|1mb)$"),
      window = feature_col %>%
        str_extract("(100kb|250kb|500kb|1mb)$"),
      window = factor(window, levels = c("100kb", "250kb", "500kb", "1mb"))
    ) %>%
    group_by(post_mean_bin, repeat_type, window) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      median_value = median(value, na.rm = TRUE),
      prop_nonzero = mean(value > 0, na.rm = TRUE),
      .groups = "drop"
    )

  fwrite(
    as.data.table(repeat_type_postmean_stats),
    file.path(cfg$output_dir, "postmean_bin_repeat_type_summary.tsv"),
    sep = "\t"
  )

  repeat_type_trend_long <- repeat_type_postmean_stats %>%
    select(post_mean_bin, repeat_type, window, mean_value, median_value) %>%
    pivot_longer(
      cols = c(mean_value, median_value),
      names_to = "stat",
      values_to = "value"
    ) %>%
    mutate(
      repeat_type = factor(repeat_type, levels = sort(unique(repeat_type))),
      stat = recode(
        stat,
        mean_value = "Mean",
        median_value = "Median"
      )
    )

  window_levels <- c("100kb", "250kb", "500kb", "1mb")
  window_label_map <- c(
    "100kb" = "100 kb",
    "250kb" = "250 kb",
    "500kb" = "500 kb",
    "1mb" = "1 Mb"
  )
  for (window_label in window_levels) {
    plot_dt <- repeat_type_trend_long %>%
      filter(as.character(window) == window_label)

    if (nrow(plot_dt) == 0L) {
      next
    }

    ncol_plot <- 3L
    types_per_page <- 24L
    type_levels <- plot_dt %>%
      pull(repeat_type) %>%
      as.character() %>%
      unique() %>%
      sort()

    page_ids <- ceiling(seq_along(type_levels) / types_per_page)
    type_pages <- split(type_levels, page_ids)

    for (page_idx in seq_along(type_pages)) {
      page_types <- type_pages[[page_idx]]
      page_dt <- plot_dt %>%
        filter(as.character(repeat_type) %in% page_types) %>%
        mutate(repeat_type = factor(as.character(repeat_type), levels = page_types))

      n_rows <- ceiling(length(page_types) / ncol_plot)
      plot_height <- max(12, min(48, n_rows * 3.2))

      p_repeat_types <- ggplot(
        page_dt,
        aes(x = post_mean_bin, y = value, color = stat, group = stat)
      ) +
        geom_line(linewidth = 0.8) +
        geom_point(size = 1.3) +
        scale_x_continuous(breaks = 1:10) +
        scale_color_manual(values = c("Mean" = "#3a86ff", "Median" = "#fb5607")) +
        facet_wrap(~repeat_type, scales = "free_y", ncol = ncol_plot) +
        labs(
          title = paste0("Repeat Type Counts Across s_het Deciles (", window_label_map[[window_label]], " window)"),
          subtitle = paste0("Mean and median counts by repeat subtype (page ", page_idx, " of ", length(type_pages), ")"),
          x = "s_het post_mean decile",
          y = "Repeat count",
          color = "Statistic"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          legend.position = "top",
          strip.text = element_text(size = 9)
        )

      out_name <- if (length(type_pages) > 1L) {
        paste0("postmean_repeat_type_trends_mean_median_", window_label, "_p", page_idx, ".pdf")
      } else {
        paste0("postmean_repeat_type_trends_mean_median_", window_label, ".pdf")
      }

      ggsave(
        filename = file.path(cfg$output_dir, "plots", out_name),
        plot = p_repeat_types,
        width = 20,
        height = plot_height,
        limitsize = FALSE
      )
    }
  }
}

repeat_family_cols <- names(analysis_tbl) %>%
  str_subset("^repeat_family_.*_count_(100kb|250kb|500kb|1mb)$")

if (length(repeat_family_cols) > 0L) {
  repeat_family_postmean_stats <- analysis_tbl %>%
    select(post_mean_bin, all_of(repeat_family_cols)) %>%
    pivot_longer(
      cols = -post_mean_bin,
      names_to = "repeat_family_col",
      values_to = "value"
    ) %>%
    mutate(
      repeat_family = repeat_family_col %>%
        str_remove("^repeat_family_") %>%
        str_remove("_count_(100kb|250kb|500kb|1mb)$"),
      window = repeat_family_col %>%
        str_extract("(100kb|250kb|500kb|1mb)$"),
      window = factor(window, levels = c("100kb", "250kb", "500kb", "1mb"))
    ) %>%
    group_by(post_mean_bin, repeat_family, window) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      median_value = median(value, na.rm = TRUE),
      prop_nonzero = mean(value > 0, na.rm = TRUE),
      .groups = "drop"
    )

  fwrite(
    as.data.table(repeat_family_postmean_stats),
    file.path(cfg$output_dir, "postmean_bin_repeat_family_summary.tsv"),
    sep = "\t"
  )

  window_levels <- c("100kb", "250kb", "500kb", "1mb")
  window_label_map <- c(
    "100kb" = "100 kb",
    "250kb" = "250 kb",
    "500kb" = "500 kb",
    "1mb" = "1 Mb"
  )

  for (window_label in window_levels) {
    window_dt <- repeat_family_postmean_stats %>%
      filter(as.character(window) == window_label)

    if (nrow(window_dt) == 0L) {
      next
    }

    top_families <- window_dt %>%
      group_by(repeat_family) %>%
      summarise(mean_median = mean(median_value, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(mean_median)) %>%
      dplyr::slice_head(n = 25) %>%
      pull(repeat_family)

    family_heat_dt <- window_dt %>%
      filter(repeat_family %in% top_families) %>%
      mutate(repeat_family = factor(repeat_family, levels = top_families))

    p_repeat_family_heat <- ggplot(
      family_heat_dt,
      aes(x = factor(post_mean_bin), y = repeat_family, fill = median_value)
    ) +
      geom_tile() +
      scale_fill_viridis_c(option = "C") +
      labs(
        title = paste0("Top Repeat Families by s_het Decile (", window_label_map[[window_label]], " window)"),
        subtitle = "Median repeat-family counts by decile",
        x = "s_het post_mean decile",
        y = "Repeat family",
        fill = "Median count"
      ) +
      theme_minimal(base_size = 11)

    ggsave(
      filename = file.path(
        cfg$output_dir,
        "plots",
        paste0("postmean_repeat_family_heatmap_top25_", window_label, ".pdf")
      ),
      plot = p_repeat_family_heat,
      width = 13,
      height = 10
    )
  }
}

if (nrow(repeat_filt_observed_summary) > 0L) {
  filtered_trend_long <- repeat_filt_observed_summary %>%
    select(filter_set, window, post_mean_bin, observed_mean, observed_median) %>%
    pivot_longer(
      cols = c(observed_mean, observed_median),
      names_to = "stat",
      values_to = "value"
    ) %>%
    mutate(
      stat = recode(
        stat,
        observed_mean = "Mean",
        observed_median = "Median"
      ),
      filter_set = factor(
        filter_set,
        levels = repeat_filter_profiles$filter_set
      )
    )

  for (window_label in names(window_label_map)) {
    window_dt <- filtered_trend_long %>%
      filter(window == window_label)

    if (nrow(window_dt) == 0L) {
      next
    }

    p_repeat_filtered <- ggplot(
      window_dt,
      aes(x = post_mean_bin, y = value, color = stat, group = stat)
    ) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 1.6) +
      scale_x_continuous(breaks = 1:10) +
      scale_color_manual(values = c("Mean" = "#3a86ff", "Median" = "#ff6f00")) +
      facet_wrap(~filter_set, scales = "free_y", ncol = 3) +
      labs(
        title = paste0(
          "Filtered Repeat Burden Across s_het Deciles (",
          window_label_map[[window_label]],
          " window)"
        ),
        subtitle = "Literature-inspired repeat filters; mean and median counts per gene",
        x = "s_het post_mean decile",
        y = "Repeat count",
        color = "Statistic"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "top",
        strip.text = element_text(size = 10)
      )

    ggsave(
      filename = file.path(
        cfg$output_dir,
        "plots",
        paste0("postmean_repeat_filtered_trends_mean_median_", window_label, ".pdf")
      ),
      plot = p_repeat_filtered,
      width = 18,
      height = 12
    )
  }
}

if (nrow(repeat_filtered_enrichment) > 0L) {
  enrichment_plot_dt <- repeat_filtered_enrichment %>%
    mutate(
      filter_set = factor(
        filter_set,
        levels = repeat_filter_profiles$filter_set
      )
    )

  for (window_label in names(window_label_map)) {
    window_dt <- enrichment_plot_dt %>%
      filter(window == window_label)

    if (nrow(window_dt) == 0L) {
      next
    }

    p_obs_bg <- ggplot(window_dt, aes(x = post_mean_bin)) +
      geom_ribbon(
        aes(ymin = bg_median_q025, ymax = bg_median_q975),
        fill = "#d9d9d9",
        alpha = 0.7
      ) +
      geom_line(
        aes(y = bg_median_mean, color = "Background median"),
        linewidth = 0.9
      ) +
      geom_line(
        aes(y = observed_median, color = "Observed median"),
        linewidth = 1.0
      ) +
      geom_point(
        aes(y = observed_median, color = "Observed median"),
        size = 1.4
      ) +
      scale_x_continuous(breaks = 1:10) +
      scale_color_manual(
        values = c(
          "Observed median" = "#d1495b",
          "Background median" = "#2a9d8f"
        )
      ) +
      facet_wrap(~filter_set, scales = "free_y", ncol = 3) +
      labs(
        title = paste0(
          "Observed vs Random-Background Repeat Burden (",
          window_label_map[[window_label]],
          " window)"
        ),
        subtitle = "Ribbon = 95% interval from random insertion simulations",
        x = "s_het post_mean decile",
        y = "Median repeat count per gene",
        color = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "top",
        strip.text = element_text(size = 10)
      )

    ggsave(
      filename = file.path(
        cfg$output_dir,
        "plots",
        paste0("postmean_repeat_filtered_observed_vs_background_", window_label, ".pdf")
      ),
      plot = p_obs_bg,
      width = 18,
      height = 12
    )

    p_enrichment_ratio <- ggplot(
      window_dt,
      aes(x = post_mean_bin, y = enrichment_ratio_median)
    ) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey45") +
      geom_line(color = "#7b2cbf", linewidth = 0.9) +
      geom_point(color = "#7b2cbf", size = 1.4) +
      scale_x_continuous(breaks = 1:10) +
      facet_wrap(~filter_set, scales = "free_y", ncol = 3) +
      labs(
        title = paste0(
          "Repeat Enrichment Ratio vs Random Background (",
          window_label_map[[window_label]],
          " window)"
        ),
        subtitle = "Enrichment ratio = observed median / simulated background median",
        x = "s_het post_mean decile",
        y = "Median enrichment ratio"
      ) +
      theme_minimal(base_size = 12) +
      theme(strip.text = element_text(size = 10))

    ggsave(
      filename = file.path(
        cfg$output_dir,
        "plots",
        paste0("postmean_repeat_filtered_enrichment_ratio_", window_label, ".pdf")
      ),
      plot = p_enrichment_ratio,
      width = 18,
      height = 12
    )
  }
}

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

p_h2_box <- analysis_dt %>%
  as_tibble() %>%
  filter(is.finite(h2_GREML), !is.na(post_mean_bin)) %>%
  ggplot(aes(x = factor(post_mean_bin), y = h2_GREML)) +
  geom_boxplot(
    fill = "#8ecae6",
    color = "#1d3557",
    outlier.alpha = 0.35,
    outlier.size = 0.8
  ) +
  labs(
    title = "h2_GREML Distribution by s_het post_mean Decile",
    x = "s_het post_mean decile",
    y = "h2_GREML"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(cfg$output_dir, "plots", "regulatory_burden_vs_h2.pdf"),
  plot = p_h2_box,
  width = 9,
  height = 5.5
)

repeat_class_cols <- names(analysis_dt) %>%
  str_subset("^repeat_class_.*_count_100kb$")

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
