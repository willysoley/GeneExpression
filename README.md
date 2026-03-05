# GeneExpression: Downstream h2 / s_het Regulatory + Repeat Analysis

## Purpose
This repository runs downstream analyses linking:
- gene expression heritability (`h2_GREML`)
- evolutionary constraint (`s_het post_mean` and deciles)
- nearby regulatory burden (enhancer/open-chromatin counts)
- nearby repeat burden (total, class, family, and literature-filtered sets)

The main script is:
- `scripts/downstream_h2_regulatory_repeat_analysis.R`

It builds a gene-level feature matrix and summary plots/tables for decile trends, association models, repeat subtype behavior, and random-genomic-background enrichment.

## Repository Objects
### Core code
- `scripts/downstream_h2_regulatory_repeat_analysis.R`
- `scripts/downstream_h2_regulatory_repeat_helpers.R`

### Documentation
- `docs/how_to_run_downstream_analysis.md`
- `docs/downstream_analysis_methods.md`

## Inputs
### Required user inputs (set in `cfg`)
- `summary_tsv`: GREML heritability summary (`Gene`, `h2_GREML`, `SE_GREML`, `Pval_GREML`)
- `shet_xlsx`: s_het table (`ensg`, `post_mean`; coordinates optional)
- `tpm_tsv`: TPM matrix from Salmon/GEUVADIS processing
- `gene_gtf`: GENCODE hg38 GTF

### Optional user overrides
- `enhancer_bed`: enhancer BED
- `open_bed`: open-chromatin BED
- `repeat_rmsk`: UCSC `rmsk` file
- `chrom_info_tsv`: chromosome size table
- `roadmap_links_dir`: local Roadmap link files (optional)

### External defaults used if optional paths are `NULL`
- ENCODE SCREEN enhancer BED (`GRCh38-cCREs.ELS.bed`)
- ENCODE SCREEN open BED (`GRCh38-cCREs.bed`)
- UCSC `hg38/database/rmsk.txt.gz`
- UCSC `hg38/database/chromInfo.txt.gz`
- GEUVADIS SDRF URL

## Pipeline Steps
1. Load SDRF and keep European GEUVADIS runs.
2. Load TPM matrix and apply expression filter.
   - Default is GTEx-style TPM filter: `TPM >= 0.1` in `>=20%` of EU samples.
3. Merge heritability summary with s_het table and build `post_mean_bin` (1-10).
4. Build gene coordinates (S_het coords if available, otherwise GTF), then construct windows:
   - gene-centered: `100kb`, `250kb`, `500kb`, `1mb`
   - TSS-centered: `100kb`, `250kb`, `500kb`, `1mb`
5. Count enhancer/open-chromatin features in windows.
6. Load RepeatMasker annotations and compute:
   - total repeat counts
   - repeat class counts
   - repeat family counts
   - repeat subtype aliases (`repeat_<type>_count_<window>`)
7. Apply literature-style repeat filters (class/divergence/length).
8. Simulate random genomic background (default 100 iterations) with fixed seed strategy.
9. Generate summaries, models, and plots by `post_mean` decile.

## Key In-Memory Objects
- `merged`: h2 + s_het + expression-filtered gene table
- `gene_tbl`: gene metadata and window widths
- `reg_features`: enhancer/open-chromatin feature table
- `repeat_features`: total repeat burden features
- `repeat_class_dt`: repeat class burden features
- `repeat_family_dt`: repeat family burden features
- `analysis_dt`: final gene-level merged feature matrix
- `repeat_filt_observed_summary`: observed filtered-repeat decile summaries
- `repeat_background_iter`: per-iteration random-background summaries
- `repeat_filtered_enrichment`: observed vs background enrichment table

## Expectations and Sanity Checks
The script emits `[Sanity]` messages and stops early on failures. Checks include:
- input table row counts (SDRF/TPM/s_het/summary/rmsk)
- matched EU sample count and valid expression thresholds
- non-empty merged gene set and decile coverage
- non-empty enhancer/open/repeat annotations
- positive window widths
- final `analysis_dt` row count matching merged genes
- uniqueness of background simulation seeds

Typical expectation:
- thousands of genes retained after expression filter and merge
- 10 s_het bins present
- non-zero counts in at least some enhancer/open/repeat columns

## Outputs
All outputs are written under:
- `results/downstream_h2_regulatory_repeat/`

### Core tables
- `gene_level_features.tsv`
- `decile_feature_summary.tsv`
- `spearman_correlations.tsv`
- `model_lm_h2.tsv`
- `model_glm_h2sig.tsv`
- `model_lm_repeat_classes.tsv` (when predictors are valid)

### Expression-filter metadata
- `expression_filter_counts.tsv`
- `expression_filter_method.tsv`
- `non_zero_eur_gene_ids.txt`

### Repeat filtering/background metadata
- `repeat_filter_criteria.tsv`
- `repeat_filter_set_counts.tsv`
- `chrom_sizes_used.tsv`
- `repeat_background_seed_info.tsv`
- `repeat_background_seed_map.tsv`
- `repeat_filtered_background_iter.tsv`
- `postmean_repeat_filtered_observed_summary.tsv`
- `postmean_repeat_filtered_background_summary.tsv`
- `postmean_repeat_filtered_enrichment.tsv`

### Decile feature summaries
- `postmean_bin_feature_summary.tsv`
- `postmean_bin_repeat_type_summary.tsv`
- `postmean_bin_repeat_family_summary.tsv`

### Plots (`results/downstream_h2_regulatory_repeat/plots/`)
- Global trends: h2, regulatory burden, repeat burden
- Repeat type trends per window (`100kb`, `250kb`, `500kb`, `1mb`)
- Repeat family heatmaps per window
- Filtered-repeat trends per window
- Observed-vs-background filtered-repeat plots per window
- Background enrichment-ratio plots per window

## Reproducibility
Set these in `cfg`:
- `repeat_bg_seed`: base seed for background simulations
- `repeat_bg_n_iter`: number of random iterations

Exact seeds used are recorded in:
- `repeat_background_seed_map.tsv`

Seed strategy:
- `seed_used = repeat_bg_seed + simulation_index_in_filter_window_loop`

If config and inputs are unchanged, reruns produce identical background outputs.

## Run
From repository root:

```bash
Rscript scripts/downstream_h2_regulatory_repeat_analysis.R
```

## Notes
- Current expression filtering implements the TPM component of GTEx-style filtering.  
  (GTEx eQTL prep also uses a read-count threshold when counts are available.)
- Windows are gene-centered for primary burden analyses; TSS-centered features are also exported for comparison.
