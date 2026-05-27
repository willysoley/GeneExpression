#!/usr/bin/env Rscript

# Fixed entrypoint: RAW + IRNT analysis
# Uses the shared pipeline in shet_tpm_h2_plots.R with pinned norm2 settings.

Sys.setenv(
  METHOD_NORM2 = "all_snps_tmm_irnt",
  METHOD_H2_TPM_NORM2 = "all_snps_tpm_irnt",
  ANALYSIS_LABEL = "shet_tmm_h2_plots_v1_irnt"
)

source("/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260428_GE_GEUVADIS_v2/GeneExpression/scripts/shet_tpm_h2_plots.R", local = FALSE)
