# Run Guide: GTEx-Style Phenotype Preparation (TMM + INT)

## Purpose
Prepare gene expression phenotypes for heritability/QTL analyses using a GTEx-style workflow:
1. filter genes by expression thresholds
2. TMM-normalize read counts
3. inverse-normal transform each gene across samples

The script is:
- `scripts/prepare_gtex_like_expression_phenotype.R`

## Inputs
Set these in `cfg`:
- `tpm_tsv`: TPM matrix (`gene x sample`)
- `counts_tsv`: raw read count matrix (`gene x sample`)
- `sdrf_url`: GEUVADIS SDRF (for selecting EUR runs)
- `eur_pops`: ancestry labels to keep

Optional:
- `sample_id_map_tsv`: 2-column mapping (`sample_id`, `participant_id`) if phenotype columns need renaming to genotype IDs.

## Filtering modes
- `gtex_tpm_only`:
  - keep genes with `TPM >= 0.1` in `>= 20%` samples
- `gtex_tpm_and_counts`:
  - keep genes with `TPM >= 0.1` in `>= 20%` samples
  - and raw counts `>= 6` in `>= 20%` samples

## Run
From repository root:

```bash
Rscript scripts/prepare_gtex_like_expression_phenotype.R
```

## Outputs
Under `results/gtex_like_expression_phenotype/`:
- `filtered_gene_ids.txt`
- `phenotype_sample_ids.txt`
- `tmm_cpm_gene_by_sample.tsv.gz`
- `phenotype_int_tmm_gene_by_sample.tsv.gz`
- `phenotype_preparation_summary.tsv`

## Notes
- TMM is computed on raw counts using `edgeR::calcNormFactors`.
- Inverse-normal transform is done per gene across samples.
- If `sample_id_map_tsv` is provided, output sample columns are renamed to participant IDs.
