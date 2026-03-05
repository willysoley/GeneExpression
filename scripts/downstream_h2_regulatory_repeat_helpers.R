#!/usr/bin/env Rscript

# Helper functions for downstream_h2_regulatory_repeat_analysis.R

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("Missing %s at: %s", label, path))
  }
}

sanity_check <- function(condition, pass_msg, fail_msg) {
  if (!isTRUE(condition)) {
    stop(fail_msg)
  }
  message("[Sanity] ", pass_msg)
}

sanity_rows <- function(df, label, min_rows = 1L) {
  n <- nrow(df)
  sanity_check(
    condition = is.finite(n) && n >= min_rows,
    pass_msg = paste0(label, ": n = ", format(n, big.mark = ",")),
    fail_msg = paste0(label, " has too few rows: n = ", n, " (expected >= ", min_rows, ").")
  )
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
    str_replace("^class[_\\s]+", "") %>%
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

ordered_feature_levels <- function(feature_names) {
  ordered_core <- c(
    "enh_count_100kb", "enh_count_250kb", "enh_count_500kb", "enh_count_1mb",
    "enh_count_tss_100kb", "enh_count_tss_250kb", "enh_count_tss_500kb", "enh_count_tss_1mb",
    "enh_overlap_gene_body_n",
    "open_count_100kb", "open_count_250kb", "open_count_500kb", "open_count_1mb",
    "open_count_tss_100kb", "open_count_tss_250kb", "open_count_tss_500kb", "open_count_tss_1mb",
    "open_overlap_gene_body_n",
    "repeat_count_100kb", "repeat_count_250kb", "repeat_count_500kb", "repeat_count_1mb",
    "repeat_count_tss_100kb", "repeat_count_tss_250kb", "repeat_count_tss_500kb", "repeat_count_tss_1mb",
    "repeat_overlap_gene_body_n"
  )

  core_levels <- ordered_core[ordered_core %in% feature_names]
  extra_levels <- setdiff(feature_names, ordered_core) %>% sort()
  c(core_levels, extra_levels)
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

detect_shet_decile_col <- function(dt, preferred_col = NULL) {
  if (!is.null(preferred_col) && nzchar(preferred_col) && preferred_col %in% names(dt)) {
    return(preferred_col)
  }

  candidate_names <- c(
    "post_mean_decile",
    "post_mean_bin",
    "s_het_decile",
    "shet_decile",
    "decile",
    "bin"
  )
  exact_candidates <- names(dt)[tolower(names(dt)) %in% tolower(candidate_names)]
  if (length(exact_candidates) > 0L) {
    return(exact_candidates[[1]])
  }

  is_valid_decile <- function(col_name) {
    x <- suppressWarnings(as.numeric(dt[[col_name]]))
    ok <- is.finite(x)
    if (sum(ok) < 10L) {
      return(FALSE)
    }

    u <- sort(unique(x[ok]))
    length(u) >= 5L &&
      length(u) <= 10L &&
      all(abs(u - round(u)) < 1e-8) &&
      min(u) >= 1 &&
      max(u) <= 10
  }

  pattern_cols <- names(dt)[str_detect(tolower(names(dt)), "decile|_bin$|bin$")]
  valid_pattern_cols <- pattern_cols[vapply(pattern_cols, is_valid_decile, logical(1))]
  if (length(valid_pattern_cols) > 0L) {
    return(valid_pattern_cols[[1]])
  }

  valid_any_cols <- names(dt)[vapply(names(dt), is_valid_decile, logical(1))]
  if (length(valid_any_cols) == 1L) {
    return(valid_any_cols[[1]])
  }

  NA_character_
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
  gene_win_quarter,
  gene_win_mid,
  gene_win_cis,
  tss_win_small,
  tss_win_quarter,
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
    count_overlaps_dt(gene_win_quarter, enh_gr, "enh_count_250kb") %>% as_tibble(),
    count_overlaps_dt(gene_win_mid, enh_gr, "enh_count_500kb") %>% as_tibble(),
    count_overlaps_dt(gene_win_cis, enh_gr, "enh_count_1mb") %>% as_tibble(),
    count_overlaps_dt(tss_win_small, enh_gr, "enh_count_tss_100kb") %>% as_tibble(),
    count_overlaps_dt(tss_win_quarter, enh_gr, "enh_count_tss_250kb") %>% as_tibble(),
    count_overlaps_dt(tss_win_mid, enh_gr, "enh_count_tss_500kb") %>% as_tibble(),
    count_overlaps_dt(tss_win_cis, enh_gr, "enh_count_tss_1mb") %>% as_tibble(),
    tibble(
      gene_id_clean = mcols(gene_body_gr)$gene_id_clean,
      enh_overlap_gene_body_n = as.integer(countOverlaps(gene_body_gr, enh_gr, ignore.strand = TRUE))
    )
  ) %>%
    purrr::reduce(left_join, by = "gene_id_clean") %>%
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

build_repeat_filter_profiles <- function() {
  # Literature-inspired defaults:
  # - Nat Genet rMEI methods commonly use LINE/SINE/SVA with milliDiv < 200.
  # - TE-focused analyses often keep interspersed TE classes and exclude low-complexity/simple repeats.
  # - Genome Research TEMR analyses commonly require minimum repeat size >= 100 bp.
  tibble::tribble(
    ~filter_set, ~include_class_regex, ~include_family_regex, ~exclude_class_regex, ~max_milliDiv, ~min_len_bp, ~source_note, ~source_url,
    "te_core_interspersed",
    "^(LINE|SINE|LTR|DNA|RC|Retroposon)$",
    NA_character_,
    "^(Simple_repeat|Low_complexity|Satellite|RNA|rRNA|scRNA|snRNA|srpRNA|tRNA)$",
    NA_real_,
    100,
    "TE-core interspersed classes; drop simple/low-complexity style classes",
    "https://genome.cshlp.org/content/early/2024/01/03/gr.278203.123.full.pdf",
    "LINE_young_milliDiv200",
    "^LINE$",
    NA_character_,
    NA_character_,
    200,
    100,
    "Nat Genet rMEI-style young/intermediate divergence threshold",
    "https://www.nature.com/articles/s41588-023-01390-2",
    "SINE_young_milliDiv200",
    "^SINE$",
    NA_character_,
    NA_character_,
    200,
    100,
    "Nat Genet rMEI-style young/intermediate divergence threshold",
    "https://www.nature.com/articles/s41588-023-01390-2",
    "SVA_young_milliDiv200",
    "^(SVA|Retroposon)$",
    "^SVA",
    NA_character_,
    200,
    100,
    "Nat Genet rMEI-style SVA-focused divergence threshold",
    "https://www.nature.com/articles/s41588-023-01390-2",
    "LTR_young_milliDiv200",
    "^LTR$",
    NA_character_,
    NA_character_,
    200,
    100,
    "Young LTR-focused set with same divergence scale",
    "https://www.nature.com/articles/s41588-023-01390-2",
    "DNA_young_milliDiv200",
    "^DNA$",
    NA_character_,
    NA_character_,
    200,
    100,
    "Young DNA transposon-focused set with same divergence scale",
    "https://www.nature.com/articles/s41588-023-01390-2"
  ) %>%
    mutate(
      filter_set = as.character(filter_set),
      include_class_regex = as.character(include_class_regex),
      include_family_regex = as.character(include_family_regex),
      exclude_class_regex = as.character(exclude_class_regex),
      max_milliDiv = as.numeric(max_milliDiv),
      min_len_bp = as.integer(min_len_bp),
      source_note = as.character(source_note),
      source_url = as.character(source_url)
    )
}

apply_repeat_filter_profile <- function(rmsk_tbl, profile_row) {
  out <- rmsk_tbl

  if (!is.na(profile_row$exclude_class_regex) && nzchar(profile_row$exclude_class_regex)) {
    out <- out %>%
      filter(!str_detect(repClass_norm, regex(profile_row$exclude_class_regex, ignore_case = TRUE)))
  }

  if (!is.na(profile_row$include_class_regex) && nzchar(profile_row$include_class_regex)) {
    out <- out %>%
      filter(str_detect(repClass_norm, regex(profile_row$include_class_regex, ignore_case = TRUE)))
  }

  if (!is.na(profile_row$include_family_regex) && nzchar(profile_row$include_family_regex)) {
    out <- out %>%
      filter(str_detect(repFamily_norm, regex(profile_row$include_family_regex, ignore_case = TRUE)))
  }

  if (!is.na(profile_row$max_milliDiv) && is.finite(profile_row$max_milliDiv)) {
    out <- out %>%
      filter(is.finite(milliDiv), milliDiv <= profile_row$max_milliDiv)
  }

  if (!is.na(profile_row$min_len_bp) && is.finite(profile_row$min_len_bp)) {
    out <- out %>%
      filter(is.finite(rep_len_bp), rep_len_bp >= profile_row$min_len_bp)
  }

  out
}

build_filtered_repeat_window_counts <- function(gene_windows, rep_gr, filter_set) {
  gene_ids <- mcols(gene_windows[[1]])$gene_id_clean
  out <- tibble(gene_id_clean = gene_ids)
  filter_clean <- sanitize_name(filter_set)

  for (window_label in names(gene_windows)) {
    col_name <- paste0("repeat_filt_", filter_clean, "_count_", window_label)
    out <- out %>%
      mutate(
        !!col_name := as.integer(
          countOverlaps(gene_windows[[window_label]], rep_gr, ignore.strand = TRUE)
        )
      )
  }

  out
}

load_chrom_sizes <- function(cfg, rmsk_tbl) {
  chrom_path <- cfg$chrom_info_tsv

  if (is.null(chrom_path) || chrom_path == "") {
    chrom_path <- tryCatch(
      safe_download(
        defaults$chrom_info_url,
        file.path(cfg$resource_dir, basename(defaults$chrom_info_url))
      ),
      error = function(e) {
        message("chromInfo download failed; falling back to chromosome sizes from rmsk.")
        NA_character_
      }
    )
  }

  if (!is.na(chrom_path) && file.exists(chrom_path)) {
    chrom_dt <- fread(chrom_path, sep = "\t", header = FALSE, fill = TRUE) %>%
      as_tibble() %>%
      transmute(
        chr = normalize_chr(V1),
        chrom_size = suppressWarnings(as.numeric(V2))
      ) %>%
      filter(
        chr %in% paste0("chr", c(1:22, "X", "Y")),
        is.finite(chrom_size),
        chrom_size > 0
      ) %>%
      distinct(chr, .keep_all = TRUE)

    if (nrow(chrom_dt) > 0L) {
      return(chrom_dt)
    }
  }

  rmsk_tbl %>%
    group_by(genoName) %>%
    summarise(chrom_size = max(genoEnd, na.rm = TRUE), .groups = "drop") %>%
    transmute(
      chr = normalize_chr(genoName),
      chrom_size = as.numeric(chrom_size)
    ) %>%
    filter(
      chr %in% paste0("chr", c(1:22, "X", "Y")),
      is.finite(chrom_size),
      chrom_size > 0
    ) %>%
    distinct(chr, .keep_all = TRUE)
}

simulate_poisson_background <- function(sim_input, n_iter, seed) {
  if (nrow(sim_input) == 0L || n_iter <= 0L) {
    return(tibble(
      iter = integer(0),
      post_mean_bin = integer(0),
      bg_mean = numeric(0),
      bg_median = numeric(0),
      bg_prop_nonzero = numeric(0)
    ))
  }

  set.seed(seed)
  sim_mat <- replicate(
    n = as.integer(n_iter),
    expr = rpois(nrow(sim_input), lambda = sim_input$lambda)
  )

  if (!is.matrix(sim_mat)) {
    sim_mat <- matrix(sim_mat, nrow = nrow(sim_input), ncol = as.integer(n_iter))
  }

  sim_dt <- as_tibble(sim_mat) %>%
    setNames(paste0("iter_", seq_len(ncol(sim_mat)))) %>%
    mutate(post_mean_bin = sim_input$post_mean_bin)

  sim_dt %>%
    pivot_longer(
      cols = starts_with("iter_"),
      names_to = "iter",
      values_to = "bg_count"
    ) %>%
    mutate(iter = as.integer(str_remove(iter, "^iter_"))) %>%
    group_by(iter, post_mean_bin) %>%
    summarise(
      bg_mean = mean(bg_count, na.rm = TRUE),
      bg_median = median(bg_count, na.rm = TRUE),
      bg_prop_nonzero = mean(bg_count > 0, na.rm = TRUE),
      .groups = "drop"
    )
}
