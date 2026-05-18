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

## Normalization Stack (Immediately After RNA-seq Download)
Starting point for this section:
- `gene_tpm.tsv` from Salmon quantification aggregation
- `gene_counts.tsv` (estimated raw counts on the same gene/sample grid)

The normalization chain below follows `nf/bin/prepare_phenotypes.R` and runs before GREML.

### Layer 0: Salmon TPM concept (what TPM means here)
Inputs used:
- TPM is treated as a within-sample relative abundance scale derived from Salmon-estimated counts and effective lengths.

Conceptual TPM definition:
```text
rate(g,s) = count(g,s) / effective_length(g,s)
TPM(g,s)  = 1e6 * rate(g,s) / sum_g' rate(g',s)
```
- The `1e6` factor makes the sample-wide TPM sum equal to one million.

Rationale:
- Corrects for feature length and sequencing depth at the abundance-estimate level.
- Useful as a stable abundance gate (and optionally as direct phenotype source with `expression_source=tpm`).

### Layer 1: Mandatory GTEx-style expression filter (no bypass)
Code path:
- Applied to both TPM and counts before downstream normalization.

Exact rule:
```text
N = number of mapped EUR samples
min_samples = ceiling(sample_frac_threshold * N)

keep_gene(g) if:
  sum_s [TPM(g,s)   >= tpm_threshold]   >= min_samples
  AND
  sum_s [COUNT(g,s) >= count_threshold] >= min_samples
```

Default numeric parameters (`nf/nextflow.config`):
- `gtex_tpm_threshold = 0.1`
- `gtex_count_threshold = 6`
- `gtex_sample_frac_threshold = 0.2`
- Therefore `min_samples = ceiling(0.2 * N)` (e.g., if `N=373`, `min_samples=75`).

Rationale:
- `TPM >= 0.1` removes near-background abundance.
- `counts >= 6` enforces minimal read support.
- `>=20%` presence avoids retaining genes expressed in only a few outliers.

### Layer 2: TMM factors and effective library sizes (`expression_source=tmm`)
Code path:
- `dge <- DGEList(counts = counts_filt)`
- `dge <- calcNormFactors(dge, method = "TMM")`

Interpretation:
- For sample `s`, edgeR computes a TMM scale factor `f_s`.
- Effective library size is `L*_s = L_s * f_s`, where `L_s` is original library size.
- edgeR rescales factors symmetrically (product of factors is 1), so only relative between-sample scaling changes.

Rationale:
- Reduces composition bias so sample-level depth/composition differences do not dominate cross-sample comparisons.

### Layer 3: CPM conversion on effective library sizes
Code path:
- `cpm(dge, normalized.lib.sizes = TRUE, log = FALSE, prior.count = 0)`

Formula:
```text
CPM(g,s) = 1e6 * COUNT(g,s) / L*_s
```

Important numeric settings:
- `normalized.lib.sizes = TRUE` uses `L*_s` (not raw `L_s`).
- `log = FALSE` keeps linear CPM scale.
- `prior.count = 0` adds no pseudocount at this step.

Rationale:
- Puts all samples on a comparable per-million scale after TMM correction.

### Layer 4: log transform
Applied when `normalization_type` is `irnt` or `inverse_normal`:
```text
log_expr(g,s) = log2(expr_base(g,s) + 1)
```

Numeric choice:
- Pseudocount is fixed at `+1` before log2.

Rationale:
- Compresses heavy right tails and keeps zeros finite (`log2(1)=0`) before rank-based normalization.

### Layer 5: Rank-based inverse-normal transform per gene across samples
Applied per gene (row-wise across samples), after `log2(x+1)`.

General formula implemented:
```text
z(g,s) = qnorm((rank(g,s) - c) / (n_g - 2c + 1))
```
- `rank` uses average ranks for ties.
- `n_g` is number of finite values for gene `g`.
- If `n_g == 1`, transformed value is set to `0`.
- If `n_g == 0`, values remain `NA`.

Offsets by mode:
- `normalization_type = "irnt"`: `c = 3/8` (Blom-style offset)
- `normalization_type = "inverse_normal"`: `c = 0.5` (rankit-style offset)
- `normalization_type = "raw"`: skip this layer entirely

Rationale:
- Produces approximately Gaussian marginals per gene, improving robustness for linear-model assumptions downstream.

### Practical caveats
- `expression_source="tpm"` skips TMM and CPM layers entirely; downstream `log2+INT` then runs on TPM values.
- Mandatory filtering still uses both TPM and counts even when expression source is TPM.
- `min_samples` changes with cohort size because of `ceiling(0.2 * N)`.
- `prior.count=0` in CPM means zeros remain zero until the explicit `log2(x+1)` step.
- Gene IDs are version-stripped and duplicates are collapsed before filtering:
  - TPM duplicates are averaged.
  - Count duplicates are summed.
- Non-finite matrix entries are set to 0 before filtering; negative counts are clamped to 0.

## Detailed Expansion: Steps 6-10
This section expands the requested steps in the current implementation:
6) GRM construction, 7) PCA, 8) phenotype prep QC, 9) mandatory GTEx-style filtering, 10) expression source (TPM vs TMM).

### Step 6: GRM construction
What it does:
- Builds the sample-by-sample genomic relationship matrix (GRM) on the EUR sample set that passed `FILTER_EUROPEANS`.

Exact commands/logic:
```bash
${params.gcta_path} --bfile ${bed.baseName} \
  --keep ${keep_list} \
  ${hm3_extract_arg} \
  --make-grm --out genotype_grm --thread-num ${task.cpus}
```
- `hm3_extract_arg` is empty unless `use_hm3_no_hla=true`, where it becomes:
  `--extract ${params.hm3_no_hla_snplist}`
- Output basename is fixed to `genotype_grm` and GCTA writes:
  - `genotype_grm.grm.bin`
  - `genotype_grm.grm.N.bin`
  - `genotype_grm.grm.id`

Key defaults:
- `use_hm3_no_hla = true`
- `hm3_no_hla_snplist = /gpfs/data/mostafavilab/shared_data/LDSC_data/hm3_no_mhc_snps.txt`
- Process resources (`slurm` profile): `cpus=4`, `memory=16 GB`, `time=2h`

Scientific rationale:
- GREML models additive polygenic effects using a GRM; this matrix is the core random-effect covariance structure.
- Restricting to HM3/no-MHC by default reduces extreme LD/complex-region behavior and stabilizes cross-run comparability.

### Step 7: PCA
What it does:
- Computes genotype PCs from the just-built GRM and carries the top PCs into covariates.

Exact commands/logic:
```bash
${params.gcta_path} --grm genotype_grm --pca ${params.n_pcs} --out geno_pca --thread-num ${task.cpus}
```
- Produces `geno_pca.eigenvec`.
- During phenotype assembly, covariate PCs are selected as:
  `PC1..PCmin(n_pcs, available_PCs)`.

Key defaults:
- `n_pcs = 5`
- Same resource block as GRM (`cpus=4`, `memory=16 GB`, `time=2h`)

Scientific rationale:
- Genotype PCs absorb broad ancestry/population-structure signal and reduce confounding in GREML/eQTL-style modeling.

### Step 8: phenotype prep QC
What it does:
- Enforces script integrity, runtime diagnostics, and output existence checks before downstream per-gene heritability jobs run.

Exact commands/logic:
- The pipeline launches phenotype prep with full argument propagation and log capture:
```bash
Rscript ${params.prepare_pheno_script} \
  "${tpm}" "${counts}" "${id_map}" "${pca_file}" "${fam_file}" \
  "${pheno_prefix}" "${num_peer}" \
  "${params.gtex_tpm_threshold}" "${params.gtex_count_threshold}" \
  "${params.gtex_sample_frac_threshold}" "${params.n_pcs}" \
  "${params.expression_source}" "${params.normalization_type}" "${peer_max_genes}" \
  2>&1 | tee prepare_phenotypes.log
```
- Runtime hard checks in `PREPARE_PHENOTYPES`:
  - required diagnostic markers must appear in `prepare_phenotypes.log`
  - explicit propagation checks for `Normalization mode: ...` and `Expression source: ...`
  - required outputs must exist and be non-empty:
    - `<prefix>.phenotypes.tsv`
    - `<prefix>.qcovar`
    - `<prefix>.gene_index_map.txt`
    - `<prefix>.filtered_gene_ids.txt`
- Pre-run script integrity checks (`run_greml.sh` and `nf/main.nf`):
  - marker presence preflight
  - SHA-256 and mtime logging
  - workflow rejects a non-canonical/stale `prepare_phenotypes.R` copy (SHA mismatch)

Key defaults:
- `prepare_pheno_script` resolves to canonical `nf/bin/prepare_phenotypes.R`
- Auto phenotype prefix (if empty): `pheno_<hm3_no_mhc|all_snps>_<tmm|tpm>_<irnt|inverse_normal|raw>`

Scientific rationale:
- Fail-fast QC prevents silent parameter drift, stale-script reuse, and empty/partial phenotype artifacts that would invalidate high-throughput GREML runs.

### Step 9: mandatory GTEx-style filtering
What it does:
- Keeps only genes with both sufficient abundance (TPM) and sufficient read support (counts) in enough samples.

Exact commands/logic:
- In `nf/bin/prepare_phenotypes.R`:
```r
N <- ncol(tpm_mat)
min_samples <- ceiling(sample_frac_threshold * N)

tpm_keep   <- rowSums(tpm_mat    >= tpm_threshold)   >= min_samples
count_keep <- rowSums(counts_mat >= count_threshold) >= min_samples
keep <- tpm_keep & count_keep
```
- Filtering is mandatory (no bypass flag).
- If zero genes pass, the script stops with:
  `No genes passed mandatory TPM+count filtering.`

Key defaults:
- `gtex_tpm_threshold = 0.1`
- `gtex_count_threshold = 6`
- `gtex_sample_frac_threshold = 0.2` (i.e., `>=20%` of mapped EUR samples)

Scientific rationale:
- Dual-threshold filtering removes low-information genes that inflate noise, while preserving genes consistently expressed with read support across the cohort.

### Step 10: expression source (TPM vs TMM)
What it does:
- Chooses the gene-by-sample matrix used as phenotype input before optional rank-based normalization.

Exact commands/logic:
- Branch in `nf/bin/prepare_phenotypes.R`:
```r
expr_base_by_gene <- if (expression_source == "tpm") {
  tpm_filt
} else {
  dge <- DGEList(counts = counts_filt)
  dge <- calcNormFactors(dge, method = "TMM")
  cpm(dge, normalized.lib.sizes = TRUE, log = FALSE, prior.count = 0)
}
```
- Then:
  - `normalization_type=raw`: use `expr_base_by_gene` directly
  - `normalization_type=irnt` or `inverse_normal`: `log2(x+1)` then per-gene inverse-normal transform

Key defaults:
- `expression_source = "tmm"` (default)
- Allowed values: `tmm`, `tpm`
- `normalization_type = "irnt"` (legacy aliases normalized internally: `ukb_irnt -> irnt`, `tmm_only -> raw`)

Scientific rationale:
- `tmm`: starts from raw counts and corrects compositional/library-size bias, usually preferred for count-derived cross-sample comparability.
- `tpm`: keeps length-normalized abundance scale and can be useful for sensitivity analyses against count-normalization assumptions.
- Keeping both options in one pipeline enables controlled method comparisons without changing downstream GREML machinery.

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

Default launch behavior (important for Slurm sweep submissions):
- `run_greml.sh` now isolates each parameter combo into its own launch dir:
  `GREML_RUN_DIR/runs/<combo-label>_<combo-hash>/`
- This avoids Nextflow session-history lock collisions on shared
  `launchDir/.nextflow/history` when many `sbatch` jobs start together.

Optional parameter overrides at submit-time:

```bash
sbatch --export=ALL run_greml.sh --peer_nk 45 --outdir results_custom
```

If you submit from outside the repository root:

```bash
sbatch --export=ALL,GREML_PROJECT_ROOT=/absolute/path/to/GeneExpression \
  /absolute/path/to/GeneExpression/run_greml.sh
```

Set an explicit run root (recommended for organized sweeps):

```bash
sbatch --export=ALL,GREML_PROJECT_ROOT=/absolute/path/to/GeneExpression,GREML_RUN_DIR=/abs/path/greml_runs \
  /absolute/path/to/GeneExpression/run_greml.sh --peer_nk 45
```

Optional compatibility mode (legacy shared launch dir, not recommended for high parallel submit):

```bash
sbatch --export=ALL,GREML_ISOLATE_BY_COMBO=0 run_greml.sh
```

## Outputs
Under each isolated run directory:
- `GREML_RUN_DIR/runs/<combo-label>_<combo-hash>/results/`

Result files:
- `data/final.phenotypes.tsv`
- `data/final.qcovar`
- `data/gene_index_map.txt`
- `data/filtered_gene_ids.txt`
- `pca/geno_pca.eigenvec`
- `pca/genotype_grm.*`
- `summary/final_heritability_summary.tsv`

Logs:
- `nextflow_logs/nextflow_<slurm-jobid|manual>.log`
- Slurm driver logs are written from `#SBATCH -o/-e` settings.
