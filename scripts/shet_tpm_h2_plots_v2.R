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
method_h2 <- Sys.getenv("METHOD_H2", "all_snps_tpm_raw")
if (!nzchar(method_h2)) stop("METHOD_H2 is empty.")

shet_xlsx <- Sys.getenv("SHET_XLSX", "/gpfs/data/mostafavilab/shared_data/gene_information/s_het_info.xlsx")
shet_sheet <- Sys.getenv("SHET_SHEET", "Supplementary Table 1")
n_deciles <- suppressWarnings(as.integer(Sys.getenv("N_DECILES", "10")))
if (!is.finite(n_deciles) || n_deciles != 10L) {
  stop("N_DECILES must be 10.")
}

analysis_root <- file.path(runs_dir, "_analysis", "shet_tpm_v2")
dir.create(analysis_root, recursive = TRUE, showWarnings = FALSE)

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

# --------------------------------------------------
# Read Shet Once
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
  summarise(post_mean = mean(post_mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(decile = ntile(post_mean, n_deciles))

run_dir <- find_run_dir(method_h2, require_summary = TRUE)
sum_file <- list.files(file.path(run_dir, "results", "summary"), pattern = "^final_heritability_summary.*\\.tsv$", full.names = TRUE)
if (length(sum_file) == 0L) stop("Missing summary file in: ", run_dir)

plots_dir <- file.path(analysis_root, "plots")
tables_dir <- file.path(analysis_root, "tables")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

dh <- fread(sum_file[1]) %>%
  as_tibble() %>%
  transmute(
    Gene = resolve_gene_ids_for_summary(Gene, run_dir),
    Status = toupper(str_trim(as.character(Status))),
    h2_GREML = suppressWarnings(as.numeric(h2_GREML))
  ) %>%
  filter(Status == "PASS", is.finite(h2_GREML))

dt <- inner_join(dh, ds, by = "Gene") %>%
  filter(is.finite(h2_GREML), is.finite(post_mean))
if (nrow(dt) == 0L) stop("No overlapping genes between h2 summary and Shet table for: ", method_h2)

epsilon <- min(dt$h2_GREML, na.rm = TRUE)
dt_ex <- dt %>% filter(h2_GREML > epsilon)
if (nrow(dt_ex) == 0L) stop("No rows left after epsilon filter for: ", method_h2)

d_summ <- dt_ex %>%
  group_by(decile) %>%
  summarise(
    mean_h2 = mean(h2_GREML, na.rm = TRUE),
    se_h2 = sd(h2_GREML, na.rm = TRUE) / sqrt(n()),
    mean_post = mean(post_mean, na.rm = TRUE),
    n_genes = n(),
    .groups = "drop"
  ) %>%
  arrange(decile)

fwrite(as.data.table(dt), file.path(tables_dir, "joined_h2_shet_all.tsv"), sep = "\t")
fwrite(as.data.table(dt_ex), file.path(tables_dir, "joined_h2_shet_excluding_epsilon.tsv"), sep = "\t")
fwrite(as.data.table(d_summ), file.path(tables_dir, "decile_summary_mean_se.tsv"), sep = "\t")

audit <- tibble(
  metric = c("method_h2", "run_dir", "n_overlap_all", "epsilon", "n_after_epsilon"),
  value = c(method_h2, run_dir, as.character(nrow(dt)), as.character(epsilon), as.character(nrow(dt_ex)))
)
fwrite(as.data.table(audit), file.path(tables_dir, "run_audit.tsv"), sep = "\t")

p_decile <- ggplot(d_summ, aes(x = decile, y = mean_h2)) +
  geom_point(color = "steelblue", size = 2.2) +
  geom_errorbar(
    aes(ymin = mean_h2 - 1.96 * se_h2, ymax = mean_h2 + 1.96 * se_h2),
    width = 0.15,
    color = "steelblue"
  ) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    title = paste0("Mean TPM h2 by Shet decile (epsilon-filtered): ", method_h2),
    subtitle = "Summary computed after filtering h2_GREML > min(h2_GREML)",
    x = "Shet decile (1-10)",
    y = "Mean h2 (GREML)"
  ) +
  theme_classic(base_size = 12)
ggsave(file.path(plots_dir, "plot1_tpm_h2_mean_se_by_decile.png"), p_decile, width = 7.5, height = 4.8, dpi = 300)

p_post <- ggplot(d_summ, aes(x = mean_post, y = mean_h2)) +
  geom_point(color = "steelblue", size = 2.2) +
  geom_errorbar(
    aes(ymin = mean_h2 - 1.96 * se_h2, ymax = mean_h2 + 1.96 * se_h2),
    width = 0.002,
    color = "steelblue"
  ) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  labs(
    title = paste0("Mean TPM h2 by mean Shet posterior mean: ", method_h2),
    subtitle = "Summary computed by Shet decile; GWAS layer omitted",
    x = "Mean posterior mean (Shet)",
    y = "Mean h2 (GREML)"
  ) +
  theme_classic(base_size = 12)
ggsave(file.path(plots_dir, "plot2_tpm_h2_mean_se_by_mean_post.png"), p_post, width = 7.5, height = 4.8, dpi = 300)

message("Saved v2 TPM-only Shet h2 analysis to: ", analysis_root)
