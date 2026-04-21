# Run Guide: GREML Nextflow Pipeline

## Overview
This workflow runs:
1. EUR sample filtering from GEUVADIS SDRF
2. genotype GRM + PCA generation
3. phenotype preparation using mandatory GTEx-style gene filtering:
   - `TPM >= 0.1` in `>= 20%` samples
   - `raw counts >= 6` in `>= 20%` samples
4. TMM normalization on counts
5. inverse-normal transform per gene
6. PEER + PCA covariate assembly
7. per-gene GCTA GREML and HEreg
8. summary table export

Core files:
- `nf/main.nf`
- `nf/nextflow.config`
- `nf/bin/prepare_phenotypes.R`
- `run_greml.sh`

## Configure
Edit `nf/nextflow.config`:
- `params.tpm_file`
- `params.counts_file`
- `params.plink_prefix`
- `params.gcta_path`
- optional thresholds: `gtex_tpm_threshold`, `gtex_count_threshold`, `gtex_sample_frac_threshold`

## Submit
From repository root:

```bash
sbatch run_greml.sh
```

Optional parameter overrides at submit-time:

```bash
sbatch --export=ALL run_greml.sh --peer_nk 45 --outdir results_custom
```

## Outputs
Under `results/`:
- `data/final.phenotypes.tsv`
- `data/final.qcovar`
- `data/gene_index_map.txt`
- `data/filtered_gene_ids.txt`
- `pca/geno_pca.eigenvec`
- `pca/genotype_grm.*`
- `summary/final_heritability_summary.tsv`

Logs:
- `nextflow_logs/nextflow.log`
- `nextflow_logs/nextflow_driver.out`
- `nextflow_logs/nextflow_driver.err`
