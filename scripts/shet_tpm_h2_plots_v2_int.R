#!/usr/bin/env Rscript

# Fixed entrypoint: RAW + INT (inverse normal) analysis
# Uses the shared pipeline in shet_tpm_h2_plots.R with pinned norm2 settings.

Sys.setenv(
  METHOD_NORM2 = "all_snps_tmm_inverse_normal",
  METHOD_H2_TPM_NORM2 = "all_snps_tpm_inverse_normal",
  ANALYSIS_LABEL = "shet_tmm_h2_plots_v2_int"
)

source("/Users/sl8085/Documents/MostafaviLab/Git/GeneExpression/scripts/shet_tpm_h2_plots.R", local = FALSE)
