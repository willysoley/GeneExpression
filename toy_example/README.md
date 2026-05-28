# GEUVADIS GREML Toy Example

This folder is a small smoke test for the latest GEUVADIS GREML Nextflow
workflow in `nf/main.nf`.

It reuses the existing small genotype/SDRF fixture under
`/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260106_0_test_GREML`
and generates a tiny local TPM/count matrix so the pipeline can be exercised
without the full production expression matrices.

Files:

- `numbered_scripts/0_make_smoke_test_inputs.sh` creates the synthetic TPM and count matrices under `smoke_test_data/`.
- `numbered_scripts/1_run_smoke_test_greml.sh` runs the Nextflow workflow locally against the toy inputs.
- `numbered_scripts/2_check_smoke_test_outputs.sh` checks the phenotype, covariate, and summary outputs.
- `toy_smoke_test.config` overrides the production paths and lowers PEER to a tiny value for this smoke test.

Run the toy example in order:

```bash
cd /gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260428_GE_GEUVADIS_v2/GeneExpression
bash toy_example/numbered_scripts/0_make_smoke_test_inputs.sh
bash toy_example/numbered_scripts/1_run_smoke_test_greml.sh
bash toy_example/numbered_scripts/2_check_smoke_test_outputs.sh
```

The smoke test is intentionally small. It checks that the workflow can:

- map synthetic expression onto the existing European sample set,
- build genotype PCs and a GRM from the existing small PLINK fixture,
- prepare phenotypes and covariates,
- and finish the per-gene GREML summary stage.
