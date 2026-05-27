#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(stringr)
})

# --------------------------------------------------
# Settings
# --------------------------------------------------
runs_dir <- Sys.getenv(
  "RUNS_DIR",
  "/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260428_GE_GEUVADIS_v2/GeneExpression/runs"
)
run_suffix <- Sys.getenv("RUN_SUFFIX", "_peerauto_pmg0_npc5")

method_raw <- Sys.getenv("METHOD_RAW", "all_snps_tmm_raw")
method_irnt <- Sys.getenv("METHOD_IRNT", "all_snps_tmm_irnt")
force_tmm_h2 <- tolower(Sys.getenv("FORCE_TMM_H2", "true")) %in% c("1", "true", "yes", "y")

method_raw_parts <- str_match(method_raw, "^(all_snps|hm3_no_mhc)_(tmm|tpm)_(raw|irnt|inverse_normal)$")
method_irnt_parts <- str_match(method_irnt, "^(all_snps|hm3_no_mhc)_(tmm|tpm)_(raw|irnt|inverse_normal)$")

method_h2_raw <- method_raw
method_h2_irnt <- method_irnt
if (force_tmm_h2) {
  if (!is.na(method_raw_parts[1, 1])) {
    method_h2_raw <- paste(method_raw_parts[1, 2], "tmm", "raw", sep = "_")
  }
  if (!is.na(method_irnt_parts[1, 1])) {
    method_h2_irnt <- paste(method_irnt_parts[1, 2], "tmm", "irnt", sep = "_")
  }
}

shet_xlsx <- Sys.getenv("SHET_XLSX", "/gpfs/data/mostafavilab/shared_data/gene_information/s_het_info.xlsx")
shet_sheet <- Sys.getenv("SHET_SHEET", "Supplementary Table 1")
n_deciles <- suppressWarnings(as.integer(Sys.getenv("N_DECILES", "10")))
if (!is.finite(n_deciles) || n_deciles != 10L) stop("N_DECILES must be 10.")

analysis_root <- file.path(runs_dir, "_analysis", "shet_tmm_h2_v2")
plots_dir <- file.path(analysis_root, "plots")
tables_dir <- file.path(analysis_root, "tables")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

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

find_run_dir <- function(method, require_summary = TRUE) {
  all_dirs <- list.dirs(runs_dir, recursive = FALSE, full.names = TRUE)
  if (length(all_dirs) == 0L) stop("No run directories found under: ", runs_dir)

  hits <- lapply(method_candidates(method), function(m) {
    patt <- paste0("^", m, run_suffix, "([_].+)?$")
    all_dirs[str_detect(basename(all_dirs), patt)]
  }) %>%
    unlist() %>%
    unique()

  if (length(hits) == 0L) stop("No run directory found for method: ", method)

  valid <- tibble(run_dir = hits) %>%
    mutate(
      data_dir = file.path(run_dir, "results", "data"),
      summary_dir = file.path(run_dir, "results", "summary"),
      has_map = file.exists(file.path(data_dir, list.files(data_dir, pattern = "^pheno_.*\\.gene_index_map\\.txt$", full.names = FALSE)[1])),
      has_summary = file.exists(file.path(summary_dir, list.files(summary_dir, pattern = "^final_heritability_summary.*\\.tsv$", full.names = FALSE)[1]))
    ) %>%
    filter(has_map, if (require_summary) has_summary else TRUE)

  if (nrow(valid) == 0L) stop("Run directory exists but required files are missing for method: ", method)

  exact <- valid %>% filter(basename(run_dir) == paste0(method, run_suffix))
  if (nrow(exact) > 0L) return(exact$run_dir[[1]])
  valid$run_dir[[1]]
}

resolve_gene_ids_for_summary <- function(gene_raw, run_dir) {
  gene_chr <- as.character(gene_raw)
  idx_mask <- str_detect(str_trim(gene_chr), "^[0-9]+$")
  if (mean(idx_mask, na.rm = TRUE) < 0.8) return(normalize_gene_id(gene_chr))

  map_file <- list.files(file.path(run_dir, "results", "data"), pattern = "^pheno_.*\\.gene_index_map\\.txt$", full.names = TRUE)
  if (length(map_file) == 0L) return(normalize_gene_id(gene_chr))

  map_dt <- fread(map_file[1], sep = "\t")
  if (!all(c("gene_name", "mpheno_index") %in% names(map_dt))) return(normalize_gene_id(gene_chr))

  idx_tbl <- tibble(mpheno_index = suppressWarnings(as.integer(gene_chr)))
  mapped <- idx_tbl %>%
    left_join(map_dt %>% transmute(mpheno_index = as.integer(mpheno_index), gene_name = as.character(gene_name)), by = "mpheno_index")

  out <- gene_chr
  replace_mask <- idx_mask & !is.na(mapped$gene_name)
  out[replace_mask] <- mapped$gene_name[replace_mask]
  normalize_gene_id(out)
}

read_h2_bundle <- function(method_h2, ds_tbl, h2_col_name, label) {
  run_dir <- find_run_dir(method_h2, require_summary = TRUE)
  sum_file <- list.files(file.path(run_dir, "results", "summary"), pattern = "^final_heritability_summary.*\\.tsv$", full.names = TRUE)
  if (length(sum_file) == 0L) stop("Missing summary file in: ", run_dir)

  pass_tbl <- fread(sum_file[1]) %>%
    as_tibble() %>%
    transmute(
      Gene = resolve_gene_ids_for_summary(Gene, run_dir),
      Status = toupper(str_trim(as.character(Status))),
      h2 = suppressWarnings(as.numeric(h2_GREML))
    ) %>%
    filter(Status == "PASS", is.finite(h2))

  all_tbl <- inner_join(
    ds_tbl,
    pass_tbl %>% group_by(Gene) %>% summarise(!!h2_col_name := mean(h2, na.rm = TRUE), .groups = "drop"),
    by = "Gene"
  )

  if (nrow(all_tbl) == 0L) stop("No overlap between Shet and h2 for method: ", method_h2)

  epsilon <- min(pass_tbl$h2, na.rm = TRUE)
  ex_tbl <- pass_tbl %>%
    filter(h2 > epsilon) %>%
    group_by(Gene) %>%
    summarise(!!h2_col_name := mean(h2, na.rm = TRUE), .groups = "drop") %>%
    inner_join(ds_tbl, by = "Gene")

  if (nrow(ex_tbl) == 0L) stop("No rows left after epsilon filtering for method: ", method_h2)

  list(
    method_h2 = method_h2,
    run_dir = run_dir,
    epsilon = epsilon,
    n_pass_rows = nrow(pass_tbl),
    all_tbl = all_tbl,
    ex_tbl = ex_tbl,
    label = label
  )
}

# --------------------------------------------------
# Read Shet + h2
# --------------------------------------------------
if (!file.exists(shet_xlsx)) stop("SHET_XLSX not found: ", shet_xlsx)

ds <- read_xlsx(shet_xlsx, sheet = shet_sheet) %>%
  as_tibble() %>%
  transmute(
    Gene = normalize_gene_id(ensg),
    post_mean = suppressWarnings(as.numeric(post_mean))
  ) %>%
  filter(!is.na(Gene), Gene != "", is.finite(post_mean)) %>%
  group_by(Gene) %>%
  summarise(post_mean = mean(post_mean, na.rm = TRUE), .groups = "drop")

raw_bundle <- read_h2_bundle(method_h2_raw, ds, "h2_raw", "RAW")
irnt_bundle <- read_h2_bundle(method_h2_irnt, ds, "h2_irnt", "IRNT")

pair_all <- raw_bundle$all_tbl %>%
  select(Gene, post_mean, h2_raw) %>%
  inner_join(irnt_bundle$all_tbl %>% select(Gene, h2_irnt), by = "Gene")

pair_ex <- raw_bundle$ex_tbl %>%
  select(Gene, post_mean, h2_raw) %>%
  inner_join(irnt_bundle$ex_tbl %>% select(Gene, h2_irnt), by = "Gene") %>%
  mutate(decile = ntile(post_mean, n_deciles))

if (nrow(pair_ex) == 0L) stop("No shared genes between RAW and IRNT after epsilon filtering.")

summary_raw <- pair_ex %>%
  group_by(decile) %>%
  summarise(
    mean_h2 = mean(h2_raw, na.rm = TRUE),
    se_h2 = sd(h2_raw, na.rm = TRUE) / sqrt(n()),
    mean_post = mean(post_mean, na.rm = TRUE),
    n_genes = n(),
    .groups = "drop"
  ) %>%
  mutate(h2_type = "RAW")

summary_irnt <- pair_ex %>%
  group_by(decile) %>%
  summarise(
    mean_h2 = mean(h2_irnt, na.rm = TRUE),
    se_h2 = sd(h2_irnt, na.rm = TRUE) / sqrt(n()),
    mean_post = mean(post_mean, na.rm = TRUE),
    n_genes = n(),
    .groups = "drop"
  ) %>%
  mutate(h2_type = "IRNT")

summary_decile <- bind_rows(summary_raw, summary_irnt) %>%
  arrange(h2_type, decile)

# --------------------------------------------------
# Save tables/audit
# --------------------------------------------------
fwrite(as.data.table(raw_bundle$all_tbl), file.path(tables_dir, "raw_joined_h2_shet_all.tsv"), sep = "\t")
fwrite(as.data.table(raw_bundle$ex_tbl), file.path(tables_dir, "raw_joined_h2_shet_excluding_epsilon.tsv"), sep = "\t")
fwrite(as.data.table(irnt_bundle$all_tbl), file.path(tables_dir, "irnt_joined_h2_shet_all.tsv"), sep = "\t")
fwrite(as.data.table(irnt_bundle$ex_tbl), file.path(tables_dir, "irnt_joined_h2_shet_excluding_epsilon.tsv"), sep = "\t")
fwrite(as.data.table(pair_all), file.path(tables_dir, "raw_irnt_overlap_all.tsv"), sep = "\t")
fwrite(as.data.table(pair_ex), file.path(tables_dir, "raw_irnt_overlap_excluding_epsilon.tsv"), sep = "\t")
fwrite(as.data.table(summary_decile), file.path(tables_dir, "decile_summary_mean_se_raw_irnt.tsv"), sep = "\t")

audit <- tibble(
  metric = c(
    "force_tmm_h2",
    "method_h2_raw",
    "method_h2_irnt",
    "run_dir_raw",
    "run_dir_irnt",
    "epsilon_raw",
    "epsilon_irnt",
    "n_pass_rows_raw",
    "n_pass_rows_irnt",
    "n_pair_all",
    "n_pair_ex"
  ),
  value = c(
    as.character(force_tmm_h2),
    raw_bundle$method_h2,
    irnt_bundle$method_h2,
    raw_bundle$run_dir,
    irnt_bundle$run_dir,
    as.character(raw_bundle$epsilon),
    as.character(irnt_bundle$epsilon),
    as.character(raw_bundle$n_pass_rows),
    as.character(irnt_bundle$n_pass_rows),
    as.character(nrow(pair_all)),
    as.character(nrow(pair_ex))
  )
)
fwrite(as.data.table(audit), file.path(tables_dir, "run_audit.tsv"), sep = "\t")

# --------------------------------------------------
# Plots (both RAW + IRNT in all plots)
# --------------------------------------------------
p_decile <- ggplot(summary_decile, aes(x = decile, y = mean_h2, color = h2_type, group = h2_type)) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = mean_h2 - 1.96 * se_h2, ymax = mean_h2 + 1.96 * se_h2), width = 0.15) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    title = "Mean TMM h2 by Shet decile (RAW + IRNT)",
    subtitle = "Deciles recalculated after RAW/IRNT overlap merge; epsilon-filtered",
    x = "Shet decile (1-10)",
    y = "Mean h2 (GREML)",
    color = "Normalization"
  ) +
  theme_classic(base_size = 12)
ggsave(file.path(plots_dir, "plot1_tmm_h2_mean_se_by_decile_raw_irnt.png"), p_decile, width = 7.5, height = 4.8, dpi = 300)

p_post <- ggplot(summary_decile, aes(x = mean_post, y = mean_h2, color = h2_type, group = h2_type)) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = mean_h2 - 1.96 * se_h2, ymax = mean_h2 + 1.96 * se_h2), width = 0.002) +
  geom_line(linewidth = 0.9) +
  labs(
    title = "Mean TMM h2 by mean Shet post_mean (RAW + IRNT)",
    subtitle = "Summary computed by Shet decile after RAW/IRNT overlap merge",
    x = "Mean posterior mean (Shet)",
    y = "Mean h2 (GREML)",
    color = "Normalization"
  ) +
  theme_classic(base_size = 12)
ggsave(file.path(plots_dir, "plot2_tmm_h2_mean_se_by_mean_post_raw_irnt.png"), p_post, width = 7.5, height = 4.8, dpi = 300)

message("Saved v2 TMM h2 analysis (RAW+IRNT) to: ", analysis_root)
