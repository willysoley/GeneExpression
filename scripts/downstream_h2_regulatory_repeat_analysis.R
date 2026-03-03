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
  tss_window_bp = 100000L,
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

sanitize_name <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("[^A-Za-z0-9]+", "_")
}

normalize_chr <- function(x) {
  y <- as.character(x)
  y <- str_replace(y, "^chr", "")
  y <- toupper(y)
  y[y == "23"] <- "X"
  y[y == "24"] <- "Y"
  y[y == "MT"] <- "M"
  paste0("chr", y)
}

normalize_repeat_class <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("\\?", "") %>%
    str_squish() %>%
    na_if("")
}

add_missing_numeric_cols <- function(df, cols, default_value = 0) {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) == 0L) {
    return(df)
  }

  default_expr <- rep(list(default_value), length(missing_cols))
  names(default_expr) <- missing_cols
  mutate(df, !!!default_expr)
}

is_chr_like <- function(x) {
  normalize_chr(x) %in% standard_chr
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
      mean(is_chr_like(x), na.rm = TRUE)
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
  chr <- normalize_chr(dt[[chr_col]])
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

expand_interval_window <- function(gr, flank_bp) {
  starts <- pmax(1L, start(gr) - flank_bp)
  ends <- end(gr) + flank_bp

  GRanges(
    seqnames = seqnames(gr),
    ranges = IRanges(start = starts, end = ends),
    strand = "*",
    gene_id_clean = mcols(gr)$gene_id_clean
  )
}

detect_shet_coord_columns <- function(dt) {
  chr_col <- intersect(
    c("chr", "chrom", "chromosome", "seqnames", "gene_chr", "gene_chrom"),
    names(dt)
  )[1]
  start_col <- intersect(
    c("start", "gene_start", "txStart", "tx_start", "geneStart", "begin"),
    names(dt)
  )[1]
  end_col <- intersect(
    c("end", "gene_end", "txEnd", "tx_end", "geneEnd", "stop"),
    names(dt)
  )[1]
  gene_name_col <- intersect(
    c("gene_name", "symbol", "hgnc_symbol", "gene_symbol"),
    names(dt)
  )[1]

  list(
    chr_col = chr_col,
    start_col = start_col,
    end_col = end_col,
    gene_name_col = gene_name_col
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
  tibble(
    gene_id_clean = mcols(window_gr)$gene_id_clean,
    !!col_name := as.integer(countOverlaps(window_gr, feature_gr, ignore.strand = TRUE))
  ) %>%
    as.data.table()
}

count_by_group <- function(window_gr, feature_gr, feature_group, prefix, suffix = "") {
  group_levels <- feature_group %>%
    unique() %>%
    sort()

  group_levels <- group_levels[!is.na(group_levels) & group_levels != ""]

  out <- tibble(gene_id_clean = mcols(window_gr)$gene_id_clean)
  if (length(group_levels) == 0L) {
    return(as.data.table(out))
  }

  group_map <- tibble(raw = as.character(group_levels)) %>%
    mutate(
      clean = gsub("[^A-Za-z0-9]+", "_", raw),
      clean = make.unique(clean)
    )

  class_tables <- group_map %>%
    split(.$raw) %>%
    lapply(function(x) {
      raw_name <- x$raw[[1]]
      clean_name <- x$clean[[1]]
      idx <- which(feature_group == raw_name)

      tibble(
        gene_id_clean = mcols(window_gr)$gene_id_clean,
        !!paste0(prefix, "_", clean_name, suffix) := as.integer(
          countOverlaps(window_gr, feature_gr[idx], ignore.strand = TRUE)
        )
      ) %>%
        as.data.table()
    })

  Reduce(
    f = function(left, right) merge(left, right, by = "gene_id_clean", all = TRUE),
    x = class_tables,
    init = as.data.table(out)
  )
}

resolve_roadmap_link_files <- function(cfg) {
  if (!is.null(cfg$roadmap_links_dir) && dir.exists(cfg$roadmap_links_dir)) {
    local_files <- list.files(
      cfg$roadmap_links_dir,
      pattern = cfg$roadmap_links_pattern,
      full.names = TRUE,
      recursive = TRUE
    )
    if (length(local_files) > 0L) {
      message("Using local Roadmap link files from: ", cfg$roadmap_links_dir)
      return(local_files)
    }
  }

  local_dir <- file.path(cfg$resource_dir, "Roadmap_links")
  dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)

  cached_files <- list.files(
    local_dir,
    pattern = cfg$roadmap_links_pattern,
    full.names = TRUE,
    recursive = TRUE
  )
  if (length(cached_files) > 0L) {
    message("Using cached Roadmap link files in: ", local_dir)
    return(cached_files)
  }

  index_urls <- unique(c(
    cfg$roadmap_links_urls,
    cfg$roadmap_links_url
  ))

  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(as.integer(old_timeout), as.integer(cfg$download_timeout_sec)))

  if (!is.null(cfg$roadmap_links_zip_url) && nzchar(cfg$roadmap_links_zip_url)) {
    zip_path <- file.path(local_dir, "RoadmapLinks.zip")
    unzip_dir <- file.path(local_dir, "RoadmapLinks_unzipped")
    dir.create(unzip_dir, recursive = TRUE, showWarnings = FALSE)

    zip_ok <- tryCatch(
      {
        safe_download(cfg$roadmap_links_zip_url, zip_path)
        TRUE
      },
      error = function(e) {
        message("Roadmap zip download failed: ", conditionMessage(e))
        FALSE
      }
    )

    if (zip_ok && file.exists(zip_path)) {
      unzip_ok <- tryCatch(
        {
          suppressWarnings(unzip(zip_path, exdir = unzip_dir))
          TRUE
        },
        error = function(e) {
          message("Roadmap zip extraction failed: ", conditionMessage(e))
          FALSE
        }
      )

      if (unzip_ok) {
        zip_files <- list.files(
          unzip_dir,
          pattern = cfg$roadmap_links_pattern,
          full.names = TRUE,
          recursive = TRUE
        )
        if (length(zip_files) > 0L) {
          message("Using Roadmap link files from downloaded ZIP.")
          return(zip_files)
        }
      }
    }
  }

  for (idx_url in index_urls) {
    message("Fetching Roadmap links index: ", idx_url)
    index_lines <- tryCatch(
      readLines(idx_url, warn = FALSE),
      error = function(e) {
        message("Roadmap index fetch failed for ", idx_url, ": ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(index_lines)) {
      next
    }

    href_matches <- str_match_all(index_lines, 'href="([^"]+)"')[[1]]
    file_names <- if (nrow(href_matches) == 0L) {
      character(0)
    } else {
      href_matches[, 2] %>%
        basename() %>%
        keep(~ str_detect(.x, cfg$roadmap_links_pattern)) %>%
        unique()
    }

    if (length(file_names) == 0L) {
      message("No matching Roadmap files found at index: ", idx_url)
      next
    }

    downloaded <- file_names %>%
      lapply(function(fname) {
        remote <- paste0(idx_url, fname)
        tryCatch(
          safe_download(
            url = remote,
            dest_path = file.path(local_dir, fname)
          ),
          error = function(e) {
            message("Roadmap file download failed for ", fname, ": ", conditionMessage(e))
            NULL
          }
        )
      }) %>%
      unlist()

    downloaded <- downloaded[!is.na(downloaded) & file.exists(downloaded)]
    if (length(downloaded) > 0L) {
      return(downloaded)
    }
  }

  stop(
    paste0(
      "Could not obtain Roadmap link files. ",
      "Set cfg$roadmap_links_dir to a local directory containing files matching pattern ",
      cfg$roadmap_links_pattern,
      ", or use cfg$enhancer_source = 'window_count'."
    )
  )
}

extract_biosample_id <- function(path) {
  fname <- basename(path)
  b <- str_extract(fname, "E\\d{3}")
  if (is.na(b)) {
    b <- fname
  }
  b
}

extract_roadmap_links <- function(path, gene_lookup_name_dt, gene_lookup_id_dt) {
  empty_links <- tibble(
    gene_id_clean = character(0),
    biosample = character(0),
    chr = character(0),
    start = numeric(0),
    end = numeric(0)
  )

  biosample <- extract_biosample_id(path)
  raw_dt <- fread(path, sep = "\t", header = FALSE, fill = TRUE, showProgress = FALSE)
  if (nrow(raw_dt) == 0L || ncol(raw_dt) < 4L) {
    return(as.data.table(empty_links))
  }

  chr_ratio <- mean(is_chr_like(raw_dt[[1]]), na.rm = TRUE)
  ensg_ratio <- mean(grepl("^ENSG", as.character(raw_dt[[4]])), na.rm = TRUE)

  if (is.finite(chr_ratio) && is.finite(ensg_ratio) && chr_ratio > 0.6 && ensg_ratio > 0.5) {
    link_dt <- tibble(
      chr = normalize_chr(raw_dt[[1]]),
      start = suppressWarnings(as.numeric(raw_dt[[2]])),
      end = suppressWarnings(as.numeric(raw_dt[[3]])),
      gene_id_clean = strip_gene_version(as.character(raw_dt[[4]]))
    ) %>%
      filter(
        chr %in% standard_chr,
        is.finite(start),
        is.finite(end),
        end > start,
        !is.na(gene_id_clean),
        gene_id_clean != ""
      ) %>%
      inner_join(as_tibble(gene_lookup_id_dt), by = "gene_id_clean") %>%
      mutate(
        biosample = biosample
      ) %>%
      select(gene_id_clean, biosample, chr, start, end)

    return(as.data.table(link_dt))
  }

  dt <- fread(path, sep = "\t", header = TRUE, fill = TRUE, showProgress = FALSE)
  names(dt) <- make.names(names(dt), unique = TRUE)

  gene_col <- intersect(c("targetGene", "target.gene", "TargetGene", "gene", "Gene", "ensg", "gene_id"), names(dt))[1]
  start_col <- intersect(c("start", "Start"), names(dt))[1]
  end_col <- intersect(c("end", "End"), names(dt))[1]
  if (is.na(gene_col) || is.na(start_col) || is.na(end_col)) {
    return(as.data.table(empty_links))
  }

  chr_col <- intersect(c("chr", "chrom", "chromosome", "X", "V1"), names(dt))[1]
  chr_raw <- if (!is.na(chr_col)) {
    dt[[chr_col]]
  } else {
    idx <- detect_chr_col(dt)
    if (is.na(idx)) {
      return(as.data.table(empty_links))
    }
    dt[[idx]]
  }

  gene_raw <- as.character(dt[[gene_col]])
  gene_is_ensg <- mean(grepl("^ENSG", gene_raw), na.rm = TRUE) > 0.5

  link_dt <- tibble(
    chr = normalize_chr(chr_raw),
    start = suppressWarnings(as.numeric(dt[[start_col]])),
    end = suppressWarnings(as.numeric(dt[[end_col]])),
    gene_raw = gene_raw
  ) %>%
    filter(
      chr %in% standard_chr,
      is.finite(start),
      is.finite(end),
      end > start,
      !is.na(gene_raw),
      gene_raw != ""
    )

  if (gene_is_ensg) {
    link_dt <- link_dt %>%
      mutate(gene_id_clean = strip_gene_version(gene_raw)) %>%
      inner_join(as_tibble(gene_lookup_id_dt), by = "gene_id_clean")
  } else {
    link_dt <- link_dt %>%
      mutate(gene_name = gene_raw) %>%
      inner_join(as_tibble(gene_lookup_name_dt), by = "gene_name")
  }

  link_dt %>%
    mutate(biosample = biosample) %>%
    select(gene_id_clean, biosample, chr, start, end) %>%
    as.data.table()
}

build_roadmap_enhancer_features <- function(roadmap_files, gene_tbl) {
  gene_tbl_tbl <- gene_tbl %>% as_tibble()

  gene_lookup_name_dt <- gene_tbl_tbl %>%
    filter(!is.na(gene_name), gene_name != "") %>%
    group_by(gene_name) %>%
    summarise(
      n_gene_ids = n_distinct(gene_id_clean),
      gene_id_clean = first(gene_id_clean),
      .groups = "drop"
    ) %>%
    filter(n_gene_ids == 1) %>%
    select(gene_name, gene_id_clean) %>%
    as.data.table()

  gene_lookup_id_dt <- gene_tbl_tbl %>%
    filter(!is.na(gene_id_clean), gene_id_clean != "") %>%
    distinct(gene_id_clean) %>%
    as.data.table()

  if (nrow(gene_lookup_id_dt) == 0L) {
    stop("No gene_id mapping available for Roadmap enhancer links.")
  }

  per_link <- roadmap_files %>%
    map_dfr(function(f) {
      message("Reading Roadmap links: ", basename(f))
      extract_roadmap_links(
        path = f,
        gene_lookup_name_dt = gene_lookup_name_dt,
        gene_lookup_id_dt = gene_lookup_id_dt
      )
    })

  out <- tibble(gene_id_clean = gene_tbl_tbl$gene_id_clean) %>%
    distinct(gene_id_clean)

  if (nrow(per_link) == 0L) {
    return(out %>%
      mutate(
        enh_link_active_biosample_n = 0,
        enh_link_mean_count_active = 0
      ) %>%
      as.data.table())
  }

  per_biosample <- per_link %>%
    distinct(gene_id_clean, biosample, chr, start, end, .keep_all = TRUE) %>%
    group_by(gene_id_clean, biosample) %>%
    summarise(
      enh_link_n = n(),
      .groups = "drop"
    )

  summary_dt <- per_biosample %>%
    group_by(gene_id_clean) %>%
    summarise(
      enh_link_active_biosample_n = n_distinct(biosample),
      enh_link_mean_count_active = mean(enh_link_n, na.rm = TRUE),
      .groups = "drop"
    )

  out %>%
    left_join(summary_dt, by = "gene_id_clean") %>%
    mutate(
      across(
        c(enh_link_active_biosample_n, enh_link_mean_count_active),
        ~ replace_na(.x, 0)
      )
    ) %>%
    as.data.table()
}

build_window_enhancer_features <- function(
  cfg,
  gene_win_small,
  gene_win_mid,
  gene_win_cis,
  tss_win_small,
  tss_win_mid,
  tss_win_cis,
  gene_body_gr
) {
  enhancer_path <- cfg$enhancer_bed
  if (is.null(enhancer_path) || enhancer_path == "") {
    enhancer_path <- safe_download(
      defaults$enhancer_url,
      file.path(cfg$resource_dir, basename(defaults$enhancer_url))
    )
  }

  enh_gr <- import_bed_like(enhancer_path, label = "enhancer annotations") %>%
    to_ucsc_style() %>%
    keep_standard_chr()

  feature_dt <- list(
    count_overlaps_dt(gene_win_small, enh_gr, "enh_count_100kb") %>% as_tibble(),
    count_overlaps_dt(gene_win_mid, enh_gr, "enh_count_500kb") %>% as_tibble(),
    count_overlaps_dt(gene_win_cis, enh_gr, "enh_count_1mb") %>% as_tibble(),
    count_overlaps_dt(tss_win_small, enh_gr, "enh_count_tss_100kb") %>% as_tibble(),
    count_overlaps_dt(tss_win_mid, enh_gr, "enh_count_tss_500kb") %>% as_tibble(),
    count_overlaps_dt(tss_win_cis, enh_gr, "enh_count_tss_1mb") %>% as_tibble(),
    tibble(
      gene_id_clean = mcols(gene_body_gr)$gene_id_clean,
      enh_overlap_gene_body_n = as.integer(countOverlaps(gene_body_gr, enh_gr, ignore.strand = TRUE))
    )
  ) %>%
    reduce(left_join, by = "gene_id_clean") %>%
    as.data.table()

  list(features = feature_dt, enhancer_gr = enh_gr)
}

build_log_terms <- function(dt, columns) {
  valid_cols <- columns %>%
    intersect(names(dt)) %>%
    keep(function(col) {
      x <- dt[[col]]
      ok <- is.finite(x)
      sum(ok) >= 5 && dplyr::n_distinct(x[ok]) >= 2
    })

  sprintf("scale(log1p(`%s`))", valid_cols)
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
  as_tibble() %>%
  filter(.data[[ancestry_col]] %in% cfg$eur_pops) %>%
  pull(.data[[run_col]]) %>%
  unique()

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
  as_tibble() %>%
  select(all_of(keep_cols))

eur_sample_cols <- setdiff(names(tpm_eur), gene_col)
if (length(eur_sample_cols) == 0L) {
  stop("No European sample columns matched TPM columns after cleaning suffixes.")
}

tpm_eur <- tpm_eur %>%
  mutate(across(all_of(eur_sample_cols), as.numeric)) %>%
  filter(if_any(all_of(eur_sample_cols), ~ .x > 0)) %>%
  mutate(
    gene_id_clean = strip_gene_version(.data[[gene_col]]),
    mean_tpm = rowMeans(select(., all_of(eur_sample_cols)), na.rm = TRUE)
  ) %>%
  as.data.table()

expressing_gene_ids <- tpm_eur$gene_id_clean %>% unique()

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
  as_tibble() %>%
  select(Gene, h2_GREML, SE_GREML, Pval_GREML) %>%
  mutate(gene_id_clean = strip_gene_version(Gene))

shet_dt <- read_xlsx(cfg$shet_xlsx, sheet = cfg$shet_sheet) %>%
  as_tibble()

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

tss_pos <- ifelse(gene_tbl$strand == "-", gene_tbl$end, gene_tbl$start)

tss_gr <- GRanges(
  seqnames = gene_tbl$chr,
  ranges = IRanges(start = tss_pos, end = tss_pos),
  strand = "*",
  gene_id_clean = gene_tbl$gene_id_clean
)

tss_win_small <- symm_window(tss_gr, cfg$tss_window_bp)
tss_win_mid <- symm_window(tss_gr, cfg$mid_window_bp)
tss_win_cis <- symm_window(tss_gr, cfg$cis_window_bp)

gene_body_gr <- GRanges(
  seqnames = gene_tbl$chr,
  ranges = IRanges(start = gene_tbl$start, end = gene_tbl$end),
  strand = "*",
  gene_id_clean = gene_tbl$gene_id_clean
)

gene_win_small <- expand_interval_window(gene_body_gr, cfg$tss_window_bp)
gene_win_mid <- expand_interval_window(gene_body_gr, cfg$mid_window_bp)
gene_win_cis <- expand_interval_window(gene_body_gr, cfg$cis_window_bp)

# ----------------------- 4) REGULATORY FEATURE COUNTS -------------------------
message("Step 5: build enhancer and open-chromatin features")

reg_features <- gene_tbl %>%
  select(gene_id_clean) %>%
  distinct(gene_id_clean)

enh_window_obj <- build_window_enhancer_features(
  cfg = cfg,
  gene_win_small = gene_win_small,
  gene_win_mid = gene_win_mid,
  gene_win_cis = gene_win_cis,
  tss_win_small = tss_win_small,
  tss_win_mid = tss_win_mid,
  tss_win_cis = tss_win_cis,
  gene_body_gr = gene_body_gr
)
enh_count_features <- enh_window_obj$features
enh_gr <- enh_window_obj$enhancer_gr

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
  count_overlaps_dt(gene_win_mid, open_gr, "open_count_500kb") %>% as_tibble(),
  count_overlaps_dt(gene_win_cis, open_gr, "open_count_1mb") %>% as_tibble(),
  count_overlaps_dt(tss_win_small, open_gr, "open_count_tss_100kb") %>% as_tibble(),
  count_overlaps_dt(tss_win_mid, open_gr, "open_count_tss_500kb") %>% as_tibble(),
  count_overlaps_dt(tss_win_cis, open_gr, "open_count_tss_1mb") %>% as_tibble(),
  tibble(
    gene_id_clean = mcols(gene_body_gr)$gene_id_clean,
    open_overlap_gene_body_n = as.integer(countOverlaps(gene_body_gr, open_gr, ignore.strand = TRUE))
  )
) %>%
  reduce(left_join, by = "gene_id_clean")

reg_features <- reg_features %>%
  left_join(open_features, by = "gene_id_clean")

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

names(rmsk)[1:17] <- c(
  "bin", "swScore", "milliDiv", "milliDel", "milliIns", "genoName", "genoStart",
  "genoEnd", "genoLeft", "strand", "repName", "repClass", "repFamily",
  "repStart", "repEnd", "repLeft", "id"
)

rmsk <- rmsk %>%
  as_tibble() %>%
  mutate(repClass = normalize_repeat_class(repClass)) %>%
  filter(
    genoName %in% paste0("chr", c(1:22, "X", "Y")),
    !is.na(repClass),
    genoEnd > genoStart
  )

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
  repeat_count_500kb = as.integer(countOverlaps(gene_win_mid, rep_gr, ignore.strand = TRUE)),
  repeat_count_1mb = as.integer(countOverlaps(gene_win_cis, rep_gr, ignore.strand = TRUE)),
  repeat_count_tss_100kb = as.integer(countOverlaps(tss_win_small, rep_gr, ignore.strand = TRUE)),
  repeat_count_tss_500kb = as.integer(countOverlaps(tss_win_mid, rep_gr, ignore.strand = TRUE)),
  repeat_count_tss_1mb = as.integer(countOverlaps(tss_win_cis, rep_gr, ignore.strand = TRUE)),
  repeat_overlap_gene_body_n = as.integer(countOverlaps(gene_body_gr, rep_gr, ignore.strand = TRUE))
)

repeat_class_dt <- count_by_group(
  window_gr = gene_win_small,
  feature_gr = rep_gr,
  feature_group = mcols(rep_gr)$repClass,
  prefix = "repeat_class",
  suffix = "_count_100kb"
) %>%
  as_tibble()

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
# repeat_<repClass>_count, sourced from repeat_class_<repClass>_count_100kb.
repeat_class_count_cols <- names(repeat_class_dt) %>%
  str_subset("^repeat_class_.*_count_100kb$")

if (length(repeat_class_count_cols) > 0L) {
  for (src_col in repeat_class_count_cols) {
    rep_class_name <- src_col %>%
      str_remove("^repeat_class_") %>%
      str_remove("_count_100kb$")
    alias_col <- paste0("repeat_", rep_class_name, "_count")

    if (!alias_col %in% names(repeat_class_dt)) {
      repeat_class_dt <- repeat_class_dt %>%
        mutate(!!alias_col := .data[[src_col]])
    }
  }
}

# --------------------------- 6) FINAL ANALYSIS TABLE --------------------------
analysis_dt <- merged %>%
  left_join(
    gene_tbl %>% select(gene_id_clean, gene_name, chr, start, end, gene_length),
    by = "gene_id_clean"
  ) %>%
  left_join(reg_features, by = "gene_id_clean") %>%
  left_join(as_tibble(repeat_features), by = "gene_id_clean") %>%
  left_join(repeat_class_dt, by = "gene_id_clean")

count_cols <- names(analysis_dt) %>%
  str_subset("^(enh_count_|enh_overlap_|open_count_|open_overlap_|repeat_count_|repeat_overlap_|repeat_class_|repeat_.*_count$)")

if (length(count_cols) > 0L) {
  analysis_dt <- analysis_dt %>%
    mutate(
      across(
        all_of(count_cols),
        ~ replace_na(as.numeric(.x), 0)
      )
    )
}

fwrite(as.data.table(analysis_dt), file.path(cfg$output_dir, "gene_level_features.tsv"), sep = "\t")

# ----------------------------- 7) SUMMARIES -----------------------------------
expected_feature_cols <- c(
  "enh_count_100kb",
  "enh_count_500kb",
  "enh_count_1mb",
  "enh_count_tss_100kb",
  "enh_count_tss_500kb",
  "enh_count_tss_1mb",
  "enh_overlap_gene_body_n",
  "open_count_100kb",
  "open_count_500kb",
  "open_count_1mb",
  "open_count_tss_100kb",
  "open_count_tss_500kb",
  "open_count_tss_1mb",
  "open_overlap_gene_body_n",
  "repeat_count_100kb",
  "repeat_count_500kb",
  "repeat_count_1mb",
  "repeat_count_tss_100kb",
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
    median_enh_500kb = median(enh_count_500kb, na.rm = TRUE),
    median_enh_1mb = median(enh_count_1mb, na.rm = TRUE),
    median_enh_overlap_gene_body_n = median(enh_overlap_gene_body_n, na.rm = TRUE),
    median_open_100kb = median(open_count_100kb, na.rm = TRUE),
    median_open_500kb = median(open_count_500kb, na.rm = TRUE),
    median_open_1mb = median(open_count_1mb, na.rm = TRUE),
    median_open_overlap_gene_body_n = median(open_overlap_gene_body_n, na.rm = TRUE),
    median_repeat_100kb = median(repeat_count_100kb, na.rm = TRUE),
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
  safe_spearman(analysis_dt$post_mean, analysis_dt$enh_count_500kb, "post_mean", "enh_count_500kb"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$enh_overlap_gene_body_n, "post_mean", "enh_overlap_gene_body_n"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$open_count_100kb, "post_mean", "open_count_100kb"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$open_count_500kb, "post_mean", "open_count_500kb"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$open_overlap_gene_body_n, "post_mean", "open_overlap_gene_body_n"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$repeat_count_100kb, "post_mean", "repeat_count_100kb"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$repeat_count_500kb, "post_mean", "repeat_count_500kb"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$repeat_overlap_gene_body_n, "post_mean", "repeat_overlap_gene_body_n"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$repeat_class_LINE_count_100kb, "post_mean", "repeat_class_LINE_count_100kb"),
  safe_spearman(analysis_dt$post_mean, analysis_dt$repeat_class_SINE_count_100kb, "post_mean", "repeat_class_SINE_count_100kb"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$enh_count_100kb, "h2_GREML", "enh_count_100kb"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$enh_count_500kb, "h2_GREML", "enh_count_500kb"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$enh_overlap_gene_body_n, "h2_GREML", "enh_overlap_gene_body_n"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$open_count_100kb, "h2_GREML", "open_count_100kb"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$open_count_500kb, "h2_GREML", "open_count_500kb"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$open_overlap_gene_body_n, "h2_GREML", "open_overlap_gene_body_n"),
  safe_spearman(analysis_dt$h2_GREML, analysis_dt$repeat_count_100kb, "h2_GREML", "repeat_count_100kb"),
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
  "enh_count_500kb",
  "enh_count_1mb",
  "enh_overlap_gene_body_n",
  "open_count_100kb",
  "open_count_500kb",
  "open_count_1mb",
  "open_overlap_gene_body_n",
  "repeat_count_100kb",
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
    slice_head(n = 6) %>%
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
  "enh_count_500kb",
  "enh_count_1mb",
  "enh_overlap_gene_body_n",
  "open_count_100kb",
  "open_count_500kb",
  "open_count_1mb",
  "open_overlap_gene_body_n",
  "repeat_count_100kb",
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

repeat_type_cols <- names(analysis_tbl) %>%
  str_subset("^repeat_[A-Za-z0-9_]+_count$") %>%
  setdiff(c(
    "repeat_count_100kb",
    "repeat_count_500kb",
    "repeat_count_1mb",
    "repeat_count_tss_100kb",
    "repeat_count_tss_500kb",
    "repeat_count_tss_1mb"
  ))

if (length(repeat_type_cols) > 0L) {
  repeat_type_postmean_stats <- analysis_tbl %>%
    select(post_mean_bin, all_of(repeat_type_cols)) %>%
    pivot_longer(
      cols = -post_mean_bin,
      names_to = "repeat_type",
      values_to = "value"
    ) %>%
    group_by(post_mean_bin, repeat_type) %>%
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
    select(post_mean_bin, repeat_type, mean_value, median_value) %>%
    pivot_longer(
      cols = c(mean_value, median_value),
      names_to = "stat",
      values_to = "value"
    ) %>%
    mutate(
      stat = recode(
        stat,
        mean_value = "Mean",
        median_value = "Median"
      )
    )

  p_repeat_types <- ggplot(
    repeat_type_trend_long,
    aes(x = post_mean_bin, y = value, color = stat, group = stat)
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.4) +
    scale_x_continuous(breaks = 1:10) +
    scale_color_manual(values = c("Mean" = "#3a86ff", "Median" = "#fb5607")) +
    facet_wrap(~repeat_type, scales = "free_y", ncol = 4) +
    labs(
      title = "Repeat Type Counts Across s_het Deciles",
      subtitle = "Mean and median counts by repeat subtype",
      x = "s_het post_mean decile",
      y = "Repeat count",
      color = "Statistic"
    ) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "top")

  ggsave(
    filename = file.path(cfg$output_dir, "plots", "postmean_repeat_type_trends_mean_median.pdf"),
    plot = p_repeat_types,
    width = 16,
    height = 10
  )
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
