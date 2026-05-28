# GEUVADIS GREML Toy Example

This folder is a small smoke test for the latest GEUVADIS GREML Nextflow
workflow in `nf/main.nf`.

It reuses the existing small genotype/SDRF fixture under
`/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260106_0_test_GREML`
and generates a tiny local TPM/count matrix so the pipeline can be exercised
without the full production expression matrices.

The toy example is launched through the single wrapper script
`toy_example/run_toy_example.sh`.

If you need to change paths or parameters for a local smoke test, edit
`toy_smoke_test.config` first. Then run:

```bash
bash toy_example/run_toy_example.sh
```

The wrapper submits the regular `run_greml.sh` driver with the toy config and
the Slurm resources needed for the small local smoke test.

The toy data under `toy_example/smoke_test_data/` are already committed, so the
user does not need to run the helper scripts unless they want to regenerate or
debug the toy inputs.

Files:

- `toy_smoke_test.config` overrides the production paths and lowers PEER to a tiny value for this smoke test.
- `run_toy_example.sh` is the one-command launch wrapper.
- `numbered_scripts/` contains the input-generator and checker helpers used during development.

The smoke test is intentionally small. It checks that the workflow can:

- map synthetic expression onto the existing European sample set,
- build genotype PCs and a GRM from the existing small PLINK fixture,
- prepare phenotypes and covariates,
- and finish the per-gene GREML summary stage.
