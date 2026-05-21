# GEUVADIS Gene Expression GREML Workflow Handoff

Last updated: 2026-05-21

This document is a knowledge-transfer guide for the GEUVADIS gene expression
heritability workflow in this repository. It is written for the next person who
needs to understand, rerun, audit, or extend the analysis after the original
workflow owner leaves the project.

The workflow estimates gene-level SNP heritability of RNA expression phenotypes
using GEUVADIS RNA-seq data, 1000 Genomes genotype data, GTEx-inspired expression
processing, and GCTA GREML. The repository also contains downstream diagnostic
analyses to compare expression choices, normalization choices, SNP sets, and
expression distribution properties.

## 1. Scientific Question

The central question is:

How much of gene expression variation across European GEUVADIS individuals can
be explained by genotype-derived relatedness, and how sensitive are those
heritability estimates to expression quantification, normalization, SNP set, and
expression distribution shape?

Operationally, each gene is treated as a quantitative phenotype. GCTA GREML is
used to estimate the fraction of phenotypic variance explained by a genomic
relationship matrix:

```text
h2_GREML = V(G) / Vp
```

where `V(G)` is the additive genetic variance component estimated from the GRM
and `Vp` is total phenotypic variance after the model accounts for covariates.

The workflow asks this question under multiple preprocessing choices:

- SNP set: all available SNPs in the genotype prefix versus HapMap3/no-MHC SNPs.
- Expression source: TPM-derived expression versus TMM-normalized count-derived expression.
- Phenotype normalization: raw expression scale versus inverse-normal transformed expression.
- Covariates: genotype PCs and PEER factors.

The broader scientific motivation is that expression heritability estimates can
be sensitive to the phenotype scale. Low expression, zero inflation, skewed
expression distributions, and different normalization choices can all change how
GREML partitions variance.

## 2. Previous Research And What This Workflow Mimics

### GEUVADIS and 1000 Genomes

GEUVADIS generated RNA-seq data for lymphoblastoid cell lines from 1000 Genomes
individuals. The public GEUVADIS collection contains mRNA and small RNA
sequencing from over 460 samples in CEU, FIN, GBR, TSI, and YRI. This workflow
focuses on the European-ancestry groups used in many GEUVADIS analyses:

- British: GBR
- Finnish: FIN
- Tuscan: TSI
- Utah residents with Northern and Western European ancestry: CEU

The project uses GEUVADIS expression data because it has both population-scale
RNA-seq and matched 1000 Genomes genotype data. That makes it a good test case
for expression heritability: for the same individuals, we can construct gene
expression phenotypes and genotype-derived GRMs.

### Hernandez et al.-style expression heritability framing

Hernandez et al. studied gene expression heritability in GEUVADIS and emphasized
how variant frequency classes, rare variants, and expression architecture affect
estimated cis heritability. In their GEUVADIS European analysis, they started
from European individuals and removed cryptically related samples using 1000
Genomes relatedness/identity-by-state information, ending with an unrelated EUR
set.

This workflow mimics the Hernandez-style sample logic by:

- Starting from GEUVADIS EUR metadata.
- Mapping RNA-seq run IDs to 1000 Genomes genotype IDs.
- Optionally intersecting the mapped EUR sample set with `unrelated_keep_file`.
- Using that sample set consistently for GRM, PCA, phenotype prep, and GREML.

Important distinction:

```text
Hernandez-style pruning = external 1000G unrelated sample keep list
not necessarily GCTA --grm-cutoff 0.025
```

The workflow now exposes this as:

```text
params.unrelated_keep_file
```

### GTEx-style expression phenotype processing

The workflow intentionally follows GTEx-like expression processing choices:

- Keep genes with TPM >= 0.1 in at least 20% of samples.
- Keep genes with counts >= 6 in at least 20% of samples.
- Use TMM normalization for counts.
- Use inverse-normal transformation per gene for normalized expression phenotypes.
- Use genotype PCs and PEER factors as covariates.
- Use a sample-size-based rule for PEER factors:
  `15, 30, 45, 60` for increasing sample-size bins.

GTEx used these types of steps to make expression phenotypes more robust for
linear modeling and eQTL mapping. This workflow reuses that logic before GREML.

### TMM normalization rationale

TMM, the trimmed mean of M-values method from Robinson and Oshlack, was developed
because simple library-size scaling can be biased when RNA composition differs
between samples. If a sample has a few highly expressed genes, those genes take
up sequencing "real estate" and can make all other genes appear lower by simple
proportional scaling. TMM estimates sample-specific scaling factors after
trimming extreme log-fold changes and intensity values, producing effective
library sizes for cross-sample comparison.

In this pipeline, when `expression_source=tmm`, counts are converted to
TMM-normalized CPM:

```text
DGEList(counts)
calcNormFactors(method = "TMM")
cpm(normalized.lib.sizes = TRUE, log = FALSE, prior.count = 0)
```

### GCTA GREML rationale

GCTA GREML uses a genomic relationship matrix as the covariance structure for
random additive genetic effects. If genetically more similar individuals also
have more similar expression for a gene, GREML attributes some expression
variance to genotype. The output `V(G)/Vp` is interpreted as the SNP heritability
captured by the SNPs used to build the GRM.

In this workflow, the GRM is built once per run setting and the same GRM is used
across genes. Each gene is then selected from a multi-phenotype file using
GCTA's `--mpheno` option.

## 3. Repository Layout

Primary files:

```text
README.md
run_greml.sh
nf/main.nf
nf/nextflow.config
nf/bin/prepare_phenotypes.R
scripts/
docs/
```

Important documentation already in the repo:

```text
docs/how_to_run_greml_nextflow.md
docs/how_to_prepare_gtex_like_phenotypes.md
docs/how_to_run_downstream_analysis.md
docs/downstream_analysis_methods.md
```

Important analysis scripts:

```text
scripts/plot_greml_h2_summary_boxplot.R
scripts/pairwise_h2_comparisons.R
scripts/tmm_expression_h2_deciles.R
scripts/tmm_raw_skewness_analysis.R
scripts/tpm_tmm_peer_correlations.R
scripts/downstream_h2_regulatory_repeat_analysis.R
```

## 4. Main HPC Paths

Project root on HPC:

```text
/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260428_GE_GEUVADIS_v2/GeneExpression
```

Run directory root:

```text
/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260428_GE_GEUVADIS_v2/GeneExpression/runs
```

Nextflow scratch/work root:

```text
/gpfs/scratch/sl8085/nextflow_work/greml_production
```

RNA-seq matrices from the previous Salmon workflow:

```text
/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260125_salmon_nextflow/results/matrices/gene_tpm.tsv
/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260125_salmon_nextflow/results/matrices/gene_counts.tsv
```

Default genotype prefix:

```text
/gpfs/data/mostafavilab/shared_data/LDSC_data/1000G_EUR_Phase3_plink/1000G.EUR.QC.22
```

This resolves to:

```text
1000G.EUR.QC.22.bed
1000G.EUR.QC.22.bim
1000G.EUR.QC.22.fam
```

The current prefix ends with `.22`. A successor should verify whether this is
intentionally chromosome 22 only or whether the intended production GRM should be
built from genome-wide merged chromosomes. This is one of the most important
handoff checks.

GCTA executable:

```text
/gpfs/data/mostafavilab/programs/GCTA/gcta-1.95.0-linux-kernel-3-x86_64/squashfs-root/AppRun
```

HM3/no-MHC SNP list:

```text
/gpfs/data/mostafavilab/shared_data/LDSC_data/hm3_no_mhc_snps.txt
```

GEUVADIS SDRF metadata:

```text
https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt
```

## 5. Upstream RNA-seq Quantification

The previous Nextflow workflow quantified GEUVADIS RNA-seq with Salmon using
GRCh37/hg19 Ensembl 75 transcriptome references:

```text
Homo_sapiens.GRCh37.75.gtf.gz
Homo_sapiens.GRCh37.75.cdna.all.fa.gz
```

It did the following:

1. Build a Salmon transcriptome index from the cDNA fasta.
2. Build a `tx2gene.tsv` map from the GTF.
3. Run `salmon quant` per GEUVADIS FASTQ using:

```text
--validateMappings
--seqBias
--gcBias
--geneMap tx2gene.tsv
```

4. Merge per-sample `quant.genes.sf` into:

```text
gene_counts.tsv
gene_tpm.tsv
```

The Salmon workflow itself produced TPM and estimated counts. It did not create
TMM. TMM is calculated later during phenotype preparation from the gene-level
count matrix.

## 6. Main Nextflow Workflow Logic

The main workflow is in:

```text
nf/main.nf
```

The main configuration is in:

```text
nf/nextflow.config
```

The Slurm driver is:

```text
run_greml.sh
```

### Step 1: Slurm driver and fanout

The usual entry point is:

```bash
sbatch run_greml.sh
```

If no CLI arguments are supplied, `run_greml.sh` defaults to fanout mode and
submits 12 child jobs:

```text
2 SNP sets x 2 expression sources x 3 normalization modes = 12 runs
```

Default fanout dimensions:

```text
SNP sets: hm3_no_mhc, all_snps
Expression sources: tpm, tmm
Normalization types: irnt, inverse_normal, raw
```

Current child run labels look like:

```text
all_snps_tmm_raw_peerauto_pmg0_npc5
all_snps_tmm_irnt_peerauto_pmg0_npc5
hm3_no_mhc_tpm_inverse_normal_peerauto_pmg0_npc5
```

The script prints a combo hash, but the current run directory name does not
include that hash.

### Step 2: FILTER_EUROPEANS

This process:

1. Reads GEUVADIS SDRF metadata.
2. Detects the population column by searching for `ancestry category`.
3. Detects the RNA run column by searching for `ENA_RUN`.
4. Keeps European groups:

```text
British
Finnish
Tuscan
Utah
```

5. Maps GEUVADIS samples to PLINK `.fam` IIDs.
6. Optionally intersects with `unrelated_keep_file`.
7. Emits:

```text
eur_keep.txt
eur_map.tsv
```

`eur_keep.txt` is the FID/IID keep list used by GCTA. `eur_map.tsv` maps RNA run
IDs to genotype FID/IID and is used by phenotype preparation.

### Step 3: GENERATE_PCA

This process builds the GRM:

```bash
gcta --bfile <plink prefix> \
  --keep eur_keep.txt \
  <optional --extract hm3_no_mhc_snps.txt> \
  --make-grm \
  --out genotype_grm \
  --thread-num <cpus>
```

Then it runs PCA from the GRM:

```bash
gcta --grm genotype_grm \
  --pca 5 \
  --out geno_pca \
  --thread-num <cpus>
```

Outputs are published to:

```text
results/pca/<runLabel>/
```

Key outputs:

```text
geno_pca.eigenvec
genotype_grm.grm.bin
genotype_grm.grm.N.bin
genotype_grm.grm.id
genotype_grm.log
```

Scientific rationale:

- The GRM defines expected genetic covariance among individuals.
- PCA captures broad ancestry structure and is included as a fixed covariate.
- HM3/no-MHC reduces complexity from high-LD and difficult regions and improves
  comparability with common SNP reference analyses.

### Step 4: PREPARE_PHENOTYPES

This process runs:

```text
nf/bin/prepare_phenotypes.R
```

It receives:

```text
TPM matrix
count matrix
eur_map.tsv
geno_pca.eigenvec
PLINK .fam
phenotype prefix
PEER settings
GTEx-style filter thresholds
expression source
normalization type
```

The script performs:

1. Sample mapping and ordering.
2. Gene ID cleanup.
3. Duplicate gene handling.
4. Mandatory expression filter.
5. TMM or TPM phenotype-source selection.
6. RAW or inverse-normal phenotype normalization.
7. PEER factor estimation.
8. PCA and PEER covariate merge.
9. GCTA-compatible phenotype and qcovar file writing.

Outputs are published to:

```text
results/data/
```

Important files:

```text
pheno_<runLabel>.phenotypes.tsv
pheno_<runLabel>.qcovar
pheno_<runLabel>.gene_index_map.txt
pheno_<runLabel>.filtered_gene_ids.txt
pheno_<runLabel>.phenotypes.tsv.gz
```

The phenotype and qcovar files intentionally have no header because GCTA expects
that format. The gene index map has a header.

### Step 5: ESTIMATE_HERITABILITY

This process runs once per gene. The workflow reads the gene index map, then
submits one `ESTIMATE_HERITABILITY` task per gene.

GREML command:

```bash
gcta --grm genotype_grm \
  --pheno pheno_<runLabel>.phenotypes.tsv \
  --mpheno <gene index> \
  --qcovar pheno_<runLabel>.qcovar \
  --reml \
  --out <gene>_greml \
  --thread-num 1
```

The output `.hsq` file contains:

```text
V(G)
V(e)
Vp
V(G)/Vp
Pval
```

The pipeline extracts:

```text
h2_GREML
SE_GREML
Pval_GREML
```

HE regression is implemented but disabled by default:

```text
run_he = false
```

This was disabled because earlier HE runs failed due an Intel MKL shared-library
issue in the GCTA AppRun environment. GREML is the primary output.

### Step 6: SUMMARIZE_RESULTS

This process concatenates one `.stats` file per gene into:

```text
results/summary/final_heritability_summary_<runLabel>.tsv
```

Columns:

```text
Gene
Status
Index
h2_GREML
SE_GREML
Pval_GREML
h2_HE
SE_HE
```

Each gene gets `PASS` or `FAIL`. Per-gene failures do not necessarily fail the
whole workflow because the workflow writes a `.stats` row and exits the gene job
cleanly.

## 7. Parameter Table And Scientific Meaning

### Sample and genotype parameters

```text
eur_pops = ["British", "Finnish", "Tuscan", "Utah"]
```

Meaning:

- Keeps the four European GEUVADIS populations.
- Excludes YRI because the primary analysis is EUR-only.
- Mimics GEUVADIS/Hernandez-style European analysis.

```text
unrelated_keep_file = ""
```

Meaning:

- Optional one-column IID or two-column FID/IID file.
- When supplied, it restricts mapped GEUVADIS EUR samples to an external
  unrelated sample set.
- This is the preferred way to mimic Hernandez-style cryptic-relatedness pruning.

```text
use_hm3_no_hla = true
```

Meaning:

- Uses HapMap3 SNPs excluding the MHC/HLA region when true.
- Rationale: HM3 SNPs are a common well-imputed/common-variant reference set;
  excluding MHC reduces extreme LD and complex-region behavior.

```text
n_pcs = 5
```

Meaning:

- Uses the top five genotype PCs as covariates.
- Rationale: genotype PCs capture broad ancestry structure even within EUR.
- This is slightly different from GTEx V7, which listed top three genotype PCs,
  but five PCs is a conservative choice for a smaller population-genetic analysis.

### Expression filtering parameters

```text
gtex_tpm_threshold = 0.1
gtex_count_threshold = 6
gtex_sample_frac_threshold = 0.2
```

Meaning:

- A gene must have TPM >= 0.1 in at least 20% of samples.
- The same gene must also have count >= 6 in at least 20% of samples.
- Both conditions must be true.

Scientific rationale:

- TPM threshold removes near-background transcription.
- Count threshold removes genes without enough read evidence.
- The 20% rule prevents a gene expressed in only a few outlier individuals from
  entering the GREML analysis.
- This follows GTEx-style eQTL expression filtering.

### Expression-source parameters

```text
expression_source = "tmm"
```

Allowed:

```text
tmm
tpm
```

If `tmm`, the pipeline starts from counts, calculates TMM factors with edgeR, and
uses TMM-normalized CPM.

If `tpm`, the pipeline uses Salmon-derived TPM directly.

Scientific difference:

- TPM corrects within sample for gene/transcript length and sequencing depth.
- TMM uses raw count composition across samples to estimate effective library
  size scaling factors.
- TPM is useful for abundance; TMM is often preferred for cross-sample count
  comparisons because it addresses RNA composition bias.

### Normalization parameters

```text
normalization_type = "irnt"
```

Allowed:

```text
raw
irnt
inverse_normal
```

`raw`:

- Uses TMM CPM or TPM values directly.
- Most biologically interpretable scale.
- More vulnerable to zero inflation and skewness.

`irnt`:

- Uses `log2(x + 1)`.
- Applies per-gene inverse-normal transform with Blom offset `3/8`.
- Makes each gene approximately standard-normal across samples.
- Less interpretable in expression units but better aligned with linear-model assumptions.

`inverse_normal`:

- Similar to `irnt`, but uses rank offset `0.5`.
- Kept mainly for comparison.

### PEER parameters

```text
peer_nk = "auto"
peer_max_genes = 0
```

GTEx-style auto rule:

```text
N < 150: 15 PEER factors
150 <= N < 250: 30 PEER factors
250 <= N < 350: 45 PEER factors
N >= 350: 60 PEER factors
```

Scientific rationale:

- PEER estimates hidden factors that capture technical and biological unwanted
  variation, such as batch effects or broad latent expression patterns.
- GTEx chose the sample-size rule based on optimizing eGene discovery.
- `peer_max_genes = 0` means use all filtered genes for PEER.
- Setting `peer_max_genes > 0` uses top-variance genes only and is a speed option.

### GCTA parameters

```text
run_he = false
```

Meaning:

- Only GREML is run by default.
- Haseman-Elston regression is optional and currently skipped.

```text
gcta_path = /gpfs/data/mostafavilab/programs/GCTA/...
```

Meaning:

- Exact GCTA executable used by the workflow.

## 8. How To Run The Workflow

Go to the project root on HPC:

```bash
cd /gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260428_GE_GEUVADIS_v2/GeneExpression
```

Submit the default 12-run fanout:

```bash
sbatch run_greml.sh
```

Submit a single run:

```bash
GREML_FANOUT=0 sbatch run_greml.sh \
  --use_hm3_no_hla false \
  --expression_source tmm \
  --normalization_type raw
```

Run with Hernandez-style unrelated keep list:

```bash
GREML_FANOUT=0 sbatch run_greml.sh \
  --use_hm3_no_hla false \
  --expression_source tmm \
  --normalization_type raw \
  --unrelated_keep_file /path/to/geuvadis_1000g_unrelated_eur.keep
```

Submit 12 fanout runs with unrelated keep file:

```bash
GREML_FANOUT=1 GREML_FANOUT_INTERVAL_SEC=0 sbatch run_greml.sh \
  --unrelated_keep_file /path/to/geuvadis_1000g_unrelated_eur.keep
```

Disable resume:

```bash
GREML_RESUME=0 sbatch run_greml.sh ...
```

Use a different run root:

```bash
GREML_RUN_DIR=/path/to/output_root sbatch run_greml.sh ...
```

Use a different scratch root:

```bash
NXF_WORK_ROOT=/gpfs/scratch/<user>/nextflow_work/greml_production sbatch run_greml.sh ...
```

## 9. Expected Run Directory Layout

A run directory looks like:

```text
runs/all_snps_tmm_raw_peerauto_pmg0_npc5/
  nextflow_logs/
  results/
    data/
    pca/
    summary/
```

Data folder:

```text
results/data/pheno_all_snps_tmm_raw.phenotypes.tsv
results/data/pheno_all_snps_tmm_raw.phenotypes.tsv.gz
results/data/pheno_all_snps_tmm_raw.qcovar
results/data/pheno_all_snps_tmm_raw.gene_index_map.txt
results/data/pheno_all_snps_tmm_raw.filtered_gene_ids.txt
```

PCA/GRM folder:

```text
results/pca/all_snps_tmm_raw/genotype_grm.grm.bin
results/pca/all_snps_tmm_raw/genotype_grm.grm.N.bin
results/pca/all_snps_tmm_raw/genotype_grm.grm.id
results/pca/all_snps_tmm_raw/geno_pca.eigenvec
```

Summary folder:

```text
results/summary/final_heritability_summary_all_snps_tmm_raw.tsv
```

## 10. Quality Checks After Running

Check final summary exists:

```bash
ls runs/*_peerauto_pmg0_npc5/results/summary/final_heritability_summary*.tsv
```

Check PASS/FAIL counts:

```bash
awk -F'\t' 'NR>1 {n[$2]++} END{for (k in n) print k,n[k]}' \
  runs/all_snps_tmm_raw_peerauto_pmg0_npc5/results/summary/final_heritability_summary*.tsv
```

Check sample count in GRM:

```bash
wc -l runs/all_snps_tmm_raw_peerauto_pmg0_npc5/results/pca/all_snps_tmm_raw/genotype_grm.grm.id
```

Check samples against PLINK `.fam`:

```bash
PLINK_PREFIX=/gpfs/data/mostafavilab/shared_data/LDSC_data/1000G_EUR_Phase3_plink/1000G.EUR.QC.22
RUN_DIR=/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260428_GE_GEUVADIS_v2/GeneExpression/runs/all_snps_tmm_raw_peerauto_pmg0_npc5
comm -23 \
  <(awk '{print $2}' "$RUN_DIR/results/pca/all_snps_tmm_raw/genotype_grm.grm.id" | sort -u) \
  <(awk '{print $2}' "${PLINK_PREFIX}.fam" | sort -u)
```

No output means all used samples are in the genotype panel.

## 11. Downstream Analysis Scripts And Plots

All diagnostic analyses write under:

```text
runs/_analysis/
```

### Global h2 summary boxplot

Script:

```text
scripts/plot_greml_h2_summary_boxplot.R
```

Outputs:

```text
runs/_analysis/plot_greml_h2_summary_boxplot/tables/summary_of_summaries.tsv
runs/_analysis/plot_greml_h2_summary_boxplot/plots/greml_h2_boxplot_by_setting.png
```

Scientific rationale:

- Shows distribution of PASS `h2_GREML` estimates across all settings.
- Helps identify which preprocessing choices inflate, deflate, or destabilize h2.
- Mean labels are useful because h2 distributions can be skewed.

### Pairwise h2 comparisons

Script:

```text
scripts/pairwise_h2_comparisons.R
```

Main plot categories:

```text
SNP set: ALL vs HM3/no-MHC
Expression source: TMM vs TPM
Normalization: RAW vs IRNT
```

Scientific rationale:

- Each scatter plot compares gene-by-gene h2 estimates between two settings.
- Identity line shows perfect agreement.
- Color by expression level tests whether low expression drives disagreement.
- Color by Wald Z tests whether discrepancies are concentrated among weak or
  strong h2 estimates.
- Discrepancy counts quantify genes where one method gives high h2 and another
  gives near-zero h2.

### Expression decile h2 analysis

Script:

```text
scripts/tmm_expression_h2_deciles.R
```

Outputs:

```text
runs/_analysis/tmm_expression_h2_deciles/plots/tmm_expression_h2_deciles_boxplot.png
runs/_analysis/tmm_expression_h2_deciles/plots/median_tpm_expression_h2_deciles_boxplot.png
runs/_analysis/tmm_expression_h2_deciles/plots/expression_skewness_by_decile_boxplot.png
```

Scientific rationale:

- Tests whether gene-level h2 is higher for more highly expressed genes.
- Uses TMM RAW mean expression and TPM RAW median expression as abundance bins.
- Filters very low expression genes to avoid bins dominated by near-zero signals.
- Adds skewness by expression bin to check whether distribution shape changes
  across abundance strata.

Default thresholds:

```text
MIN_LOG2_MEAN_TMM_PLUS1 = 1.5
MIN_LOG2_MEDIAN_TPM_PLUS1 = 0
```

### TMM raw skewness analysis

Script:

```text
scripts/tmm_raw_skewness_analysis.R
```

Outputs include:

```text
tmm_raw_skewness_hist_density_<METHOD>.png
tmm_raw_skewness_mean_vs_skew_<METHOD>.png
tmm_raw_skewness_selected_gene_distributions_<METHOD>.png
tmm_vs_tpm_h2_discrepancy_expr_distributions_<METHOD>.png
all_raw_h2_vs_hm3_raw_h2_tmm_colored_by_log2meanexpr_cap*.png
all_raw_h2_vs_hm3_raw_h2_tmm_colored_by_skewness_cap*.png
```

Scientific rationale:

- Skewness detects whether raw expression phenotypes are dominated by long right
  tails or zero inflation.
- Selected gene histograms give concrete examples for lab discussion.
- H2 discrepancy examples test whether TPM/TMM disagreement is visible at the
  expression distribution level.
- ALL-vs-HM3 colored scatter tests whether SNP-set sensitivity is linked to
  expression abundance or skewness.

### TPM vs TMM phenotype correlations

Script:

```text
scripts/tpm_tmm_peer_correlations.R
```

Outputs:

```text
tpm_tmm_peer_gene_correlations.tsv
tpm_tmm_peer_correlation_summary.tsv
tpm_tmm_peer_correlation_boxplot.png
tpm_tmm_peer_mean_scatter.png
```

Scientific rationale:

- Compares TPM-derived and TMM-derived phenotype values directly.
- High per-gene correlation means expression source choice is unlikely to change
  that gene's phenotype ranking across individuals.
- Low correlation identifies genes where TPM/TMM choice may alter h2.

### Regulatory/repeat downstream analysis

Script:

```text
scripts/downstream_h2_regulatory_repeat_analysis.R
```

This is a separate biological follow-up pipeline. It links h2 to:

- `s_het post_mean`
- gene regulatory burden
- enhancer/open chromatin overlap
- repeat classes/families/subtypes
- randomized repeat background

Scientific rationale:

- Tests whether gene expression heritability varies with evolutionary constraint.
- Tests whether regulatory architecture or repeat burden differs across gene bins.
- Observed-vs-background repeat analysis asks whether repeat patterns exceed
  random genomic expectation.

## 12. Interpretation Guide For Main Plots

### h2 boxplot by setting

Use this first. It answers:

- Which pipeline choices produce higher h2?
- Are distributions broad, compressed, or outlier-heavy?
- Are PASS counts similar across methods?

Do not interpret differences without checking PASS/FAIL counts.

### Pairwise h2 scatter plots

Use these to diagnose gene-level stability.

Interpretation:

- Points on diagonal: stable gene h2 across methods.
- Points above diagonal: y-axis method estimates higher h2.
- Points below diagonal: x-axis method estimates higher h2.
- Low-expression colored outliers: likely expression-scale or low-information issue.
- High-Z outliers: more likely biologically meaningful or at least statistically stable.

### Expression decile h2 plots

Use these to ask whether expression abundance drives h2.

Interpretation:

- Increasing h2 across deciles: low-expression genes may suppress h2 estimates.
- Flat h2 across deciles: h2 estimates are robust to abundance after filtering.
- Mean greater than median: h2 distribution is right-skewed in that bin.

### Skewness plots

Use these to ask whether raw expression distributions are suitable for GREML.

Interpretation:

- High positive skew: a few individuals have much higher expression.
- High zero fraction: raw-scale h2 may be unstable.
- IRNT should reduce sensitivity to skew because it rank-normalizes per gene.

### TPM/TMM correlation plots

Use these to ask whether expression source changes phenotype rankings.

Interpretation:

- High correlation: TPM and TMM mostly preserve individual-level ordering.
- Low correlation: TPM/TMM choice can change the phenotype and therefore h2.

## 13. Known Caveats

1. The default genotype prefix ends with `.22`. Confirm whether the analysis is
   intentionally chromosome 22 only or whether genome-wide genotype files should
   be merged or supplied.

2. `params.grm_prefix` exists but the code currently uses hard-coded
   `genotype_grm`.

3. `run_he=true` may fail unless the GCTA AppRun MKL runtime issue is fixed.
   The current workflow intentionally sets `run_he=false`.

4. Full reuse mode skips sample filtering and phenotype preparation. If
   `unrelated_keep_file` is set while reusing GRM and phenotype files, the keep
   file does not change those reused artifacts.

5. Phenotype-only reuse with a newly generated GRM can mismatch normalization
   and sample set. The workflow warns but does not stop.

6. Per-gene GREML failures do not necessarily fail the whole workflow. Always
   inspect `Status` in the summary table.

7. Several downstream scripts have hard-coded HPC paths. Before rerunning on a
   new project copy, inspect the first 20-40 lines of each R script.

8. Plot files are generated on the HPC run directory under `runs/_analysis/`;
   they are not stored in this local repo by default.

## 14. Minimal Handoff Checklist

Before a new analyst reruns the workflow:

1. Confirm the genotype prefix is correct and genome-wide if needed.
2. Confirm the GEUVADIS TPM/count matrices exist.
3. Confirm `gcta_path` works on compute nodes.
4. Confirm the `r_peer` environment has `data.table`, `edgeR`, `peer`, and ideally `matrixStats`.
5. Decide whether to provide `unrelated_keep_file`.
6. Decide whether to run all 12 settings or a single setting.
7. Keep `run_he=false` unless the MKL issue is solved.
8. Run a single small/default job first if changing input paths.
9. Check sample count in `genotype_grm.grm.id`.
10. Check PASS/FAIL counts in the final summary.
11. Run diagnostic plots under `scripts/`.
12. Archive both `nextflow_logs/` and `results/`.

## 15. Source List

GEUVADIS/1000G:

- International Genome Sample Resource GEUVADIS data collection:
  https://www.internationalgenome.org/data-portal/data-collections/geuvadis.html
- GEUVADIS project summary: RNA-seq from 1000 Genomes samples across CEU, FIN,
  GBR, TSI, and YRI.

Salmon:

- Patro et al. Salmon provides fast and bias-aware quantification of transcript
  expression. Nature Methods 2017.
  https://www.nature.com/articles/nmeth.4197

TMM:

- Robinson and Oshlack. A scaling normalization method for differential
  expression analysis of RNA-seq data. Genome Biology 2010.
  https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25

GTEx methods:

- GTEx Analysis Methods V7 PDF.
  https://storage.googleapis.com/gtex-public-data/Portal_Analysis_Methods_v7_09052017.pdf

GCTA/GREML:

- Yang et al. GCTA: A Tool for Genome-wide Complex Trait Analysis. American
  Journal of Human Genetics 2011.
  https://www.sciencedirect.com/science/article/pii/S0002929710005987
- GCTA documentation and notes on GRM and relatedness pruning:
  https://yanglab.westlake.edu.cn/software/gcta/

GEUVADIS expression heritability:

- Hernandez et al. Ultrarare variants drive substantial cis heritability of human
  gene expression. Nature Genetics 2019.
  https://www.nature.com/articles/s41588-019-0487-7
- Dahl et al. On Negative Heritability and Negative Estimates of Heritability.
  Genetics 2020.
  https://academic.oup.com/genetics/article/215/2/343/5930411

