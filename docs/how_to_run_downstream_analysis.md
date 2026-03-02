# Run Guide: h2 / s_het Regulatory + Repeat Pipeline

## 1) Edit config paths

Open:
- `scripts/downstream_h2_regulatory_repeat_analysis.R`

Update at minimum:
- `cfg$summary_tsv`
- `cfg$shet_xlsx`
- `cfg$tpm_tsv`
- `cfg$gene_gtf`

Optional overrides:
- `cfg$enhancer_source` (`"roadmap_links"` or `"window_count"`)
- `cfg$enhancer_fallback_source` (`"window_count"` or `"none"`)
- `cfg$roadmap_links_dir` (local folder of Roadmap enhancer-gene link files)
- `cfg$roadmap_links_url` (remote index URL if downloading)
- `cfg$roadmap_links_urls` (ordered list of fallback index URLs)
- `cfg$roadmap_links_zip_url` (Zenodo ZIP URL containing Roadmap links)
- `cfg$roadmap_links_pattern` (regex for selecting Roadmap files)
- `cfg$download_timeout_sec` (network timeout for downloads/index fetch)
- `cfg$enhancer_bed` (cell-type-specific enhancer BED)
- `cfg$open_bed` (cell-type-specific ATAC/DNase BED)
- `cfg$repeat_rmsk` (local UCSC `rmsk.txt(.gz)`)

If optional paths are `NULL`, the script downloads defaults:
- ENCODE SCREEN `GRCh38-cCREs.ELS.bed`
- ENCODE SCREEN `GRCh38-cCREs.bed`
- UCSC `hg38/database/rmsk.txt.gz`
- Roadmap link files from ErnstLab index and/or Zenodo ZIP (if `enhancer_source = "roadmap_links"` and no local dir is provided)

Timeout handling:
- If Roadmap index/download fails, the script can automatically fall back to `window_count` mode when `cfg$enhancer_fallback_source = "window_count"`.
- To force strict paper-mode only, set `cfg$enhancer_fallback_source = "none"`.
- Recommended for HPC reproducibility: pre-download/extract Roadmap link files and set `cfg$roadmap_links_dir` to that local path.

## 2) Run

```bash
Rscript scripts/downstream_h2_regulatory_repeat_analysis.R
```

## 3) Main outputs

In `results/downstream_h2_regulatory_repeat/`:
- `non_zero_eur_gene_ids.txt`
- `expression_filter_counts.tsv`
- `gene_level_features.tsv`
- `decile_feature_summary.tsv`
- `spearman_correlations.tsv`
- `model_lm_h2.tsv`
- `model_glm_h2sig.tsv`
- `model_lm_repeat_classes.tsv` (if repeat classes present)
- `plots/median_h2_by_decile.pdf`
- `plots/regulatory_burden_vs_h2.pdf`
- `plots/repeat_class_heatmap.pdf`

Interpretation note:
- `repeat_count_100kb`, `repeat_count_1mb`, and `repeat_class_*` are counts of repeat intervals within fixed TSS-centered windows, not indicators of whether the gene body overlaps repeats.
- `enh_link_active_biosample_n` and `enh_link_mean_total_bp_active` are paper-style enhancer features from linked enhancer-gene maps (Mostafavi et al. 2023 / Liu et al. links), not fixed-window enhancer overlaps.
- `repeat_class_LINE_*` and `repeat_class_SINE_*` provide subtype-specific repeat burden features.
