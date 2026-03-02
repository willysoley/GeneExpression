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
- `cfg$enhancer_bed` (cell-type-specific enhancer BED)
- `cfg$open_bed` (cell-type-specific ATAC/DNase BED)
- `cfg$repeat_rmsk` (local UCSC `rmsk.txt(.gz)`)

If optional paths are `NULL`, the script downloads defaults:
- ENCODE SCREEN `GRCh38-cCREs.ELS.bed`
- ENCODE SCREEN `GRCh38-cCREs.bed`
- UCSC `hg38/database/rmsk.txt.gz`

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

