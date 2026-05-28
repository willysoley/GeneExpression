# GEUVADIS GREML Toy Scripts

Run these scripts in order:

1. `0_make_smoke_test_inputs.sh` creates the synthetic TPM and count tables and copies the fixture SDRF into the toy folder.
2. `1_run_smoke_test_greml.sh` launches the Nextflow workflow locally.
3. `2_check_smoke_test_outputs.sh` verifies that the phenotype and heritability summary files exist.

This is a smoke test for workflow wiring and file schemas, not a scientific GREML analysis.
