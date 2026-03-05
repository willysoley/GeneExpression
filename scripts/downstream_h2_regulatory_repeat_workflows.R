#!/usr/bin/env Rscript

# Workflow-level helper functions for downstream_h2_regulatory_repeat_analysis.R

run_repeat_background_analysis <- function(
  cfg,
  analysis_dt,
  repeat_filter_profiles,
  repeat_profile_tbls,
  repeat_profile_chr_stats,
  gene_window_list,
  window_width_col_map,
  rmsk,
  output_dir
) {
  window_label_map <- c(
    "100kb" = "100 kb",
    "250kb" = "250 kb",
    "500kb" = "500 kb",
    "1mb" = "1 Mb"
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
      file.path(output_dir, "postmean_repeat_filtered_observed_summary.tsv"),
      sep = "\t"
    )

    chrom_sizes_dt <- load_chrom_sizes(cfg, rmsk)
    sanity_rows(chrom_sizes_dt, "Chromosome size table", min_rows = 20L)
    fwrite(
      as.data.table(chrom_sizes_dt),
      file.path(output_dir, "chrom_sizes_used.tsv"),
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
        file.path(output_dir, "repeat_background_seed_info.tsv"),
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
      reference_scale <- if (
        is.finite(reference_fraction) &&
          reference_fraction > 0 &&
          is.finite(as.numeric(cfg$repeat_bg_repeat_fraction))
      ) {
        as.numeric(cfg$repeat_bg_repeat_fraction) / reference_fraction
      } else {
        0
      }

      sim_plan <- filter_bp_dt %>%
        mutate(
          reference_filter_set = reference_filter_set,
          target_fraction = pmin(
            0.95,
            pmax(0, observed_fraction * reference_scale)
          ),
          repeat_bg_method = bg_method,
          type_col = as.character(cfg$repeat_bg_type_col),
          composition_basis = as.character(cfg$repeat_bg_composition_basis)
        )

      fwrite(
        as.data.table(sim_plan),
        file.path(output_dir, "repeat_background_simulation_plan.tsv"),
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
          file.path(output_dir, "repeat_background_seed_map.tsv"),
          sep = "\t"
        )
      }

      if (length(bg_iter_tbls) > 0L) {
        repeat_background_iter <- bind_rows(bg_iter_tbls) %>%
          arrange(filter_set, window, seed_used, iter, post_mean_bin)

        fwrite(
          as.data.table(repeat_background_iter),
          file.path(output_dir, "repeat_filtered_background_iter.tsv"),
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
      file.path(output_dir, "postmean_repeat_filtered_background_summary.tsv"),
      sep = "\t"
    )
    fwrite(
      as.data.table(repeat_filtered_enrichment),
      file.path(output_dir, "postmean_repeat_filtered_enrichment.tsv"),
      sep = "\t"
    )
    sanity_rows(repeat_filtered_enrichment, "Filtered repeat enrichment summary", min_rows = 10L)
  }

  list(
    repeat_filt_observed_summary = repeat_filt_observed_summary,
    repeat_background_iter = repeat_background_iter,
    repeat_background_summary = repeat_background_summary,
    repeat_filtered_enrichment = repeat_filtered_enrichment,
    repeat_bg_seed_map = repeat_bg_seed_map,
    window_label_map = window_label_map,
    repeat_filter_map = repeat_filter_map
  )
}


generate_summary_plots <- function(
  cfg,
  analysis_dt,
  decile_summary,
  repeat_filt_observed_summary,
  repeat_filtered_enrichment,
  repeat_filter_profiles,
  output_dir
) {
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
      file.path(output_dir, "postmean_bin_feature_summary.tsv"),
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
      filename = file.path(output_dir, "plots", "postmean_all_feature_trends_mean_median.pdf"),
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
      filename = file.path(output_dir, "plots", "postmean_regulatory_repeat_trends_mean_median.pdf"),
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
      file.path(output_dir, "postmean_bin_repeat_type_summary.tsv"),
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
          filename = file.path(output_dir, "plots", out_name),
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
      file.path(output_dir, "postmean_bin_repeat_family_summary.tsv"),
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
          output_dir,
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

    window_label_map <- c(
      "100kb" = "100 kb",
      "250kb" = "250 kb",
      "500kb" = "500 kb",
      "1mb" = "1 Mb"
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
          output_dir,
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

    window_label_map <- c(
      "100kb" = "100 kb",
      "250kb" = "250 kb",
      "500kb" = "500 kb",
      "1mb" = "1 Mb"
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
          output_dir,
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
          output_dir,
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
    filename = file.path(output_dir, "plots", "median_h2_by_decile.pdf"),
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
    filename = file.path(output_dir, "plots", "regulatory_burden_vs_h2.pdf"),
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
      filename = file.path(output_dir, "plots", "repeat_class_heatmap.pdf"),
      plot = p_heat,
      width = 9,
      height = 6
    )
  }
}
