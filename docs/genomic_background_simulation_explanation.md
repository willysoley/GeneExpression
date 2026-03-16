# Genomic Background Simulation for Repeat-Window Analysis

This note explains how the genomic background simulation in this repository works.
It covers:

- the biological/statistical concept
- the exact algorithm used in the current code
- how the `old_milliDiv200plus` repeat sets are handled
- where each step lives in the codebase

## 1. Why simulate a genomic background?

The basic question is:

"If repeat elements were placed across the genome at random, while keeping some
important properties of the chosen repeat subset, how many repeat elements would
we expect near genes in each s_het decile?"

The observed analysis counts real repeat elements around genes. The background
simulation creates a null model. That null model lets us ask whether high- or
low-s_het genes have more repeat burden than expected from a randomized genome.

In this repository, the background is not a permutation of genes. It is a
randomized repeat genome:

- genes and gene windows stay fixed
- repeat annotations are replaced by simulated repeat intervals
- the simulated intervals preserve broad properties of the chosen repeat subset
- the exact genomic positions are broken

## 2. What is preserved and what is randomized?

### Preserved

- chromosome sizes
- the selected repeat subset (for example `SINE_old_milliDiv200plus`)
- the repeat-length distribution of that subset
- the chromosome distribution of that subset
- the gene windows used in the real analysis (`100kb`, `250kb`, `500kb`, `1mb`)
- the gene-to-decile assignments (`post_mean_bin`)

### Randomized

- the exact genomic coordinates of repeats
- local clustering and hotspot structure
- the original neighborhood around each element
- exact locus-level overlap patterns

So the null model answers:

"Given the overall burden, lengths, and chromosome occupancy of this repeat
subset, what would repeat-window counts look like if placements were random?"

## 3. Current configuration used by the script

The main configuration lives in:

- `scripts/downstream_h2_regulatory_repeat_analysis.R:56-63`

The key parameters are:

- `repeat_bg_method = "explicit_genome"`
- `repeat_bg_repeat_fraction = 0.5`
- `repeat_bg_type_col = "repClass_norm"`
- `repeat_bg_composition_basis = "bp"`
- `repeat_bg_n_iter = 100`
- `repeat_bg_n_cores = 6`
- `repeat_bg_seed = 20260305`
- `repeat_bg_use_cache = TRUE`

Interpretation:

- `explicit_genome` means the code explicitly simulates repeat intervals across
  chromosomes, rather than using a Poisson approximation.
- `repeat_bg_type_col = "repClass_norm"` means the simulated composition is
  organized by repeat class.
- `repeat_bg_composition_basis = "bp"` means class weights are based on total
  base pairs, not raw element counts.
- `repeat_bg_n_iter = 100` means 100 randomized genomes are generated per
  filter set.

## 4. Step-by-step algorithm

### Step 1. Build repeat filter sets

Filter definitions are created in:

- `scripts/downstream_h2_regulatory_repeat_helpers.R:731-840`

Examples include:

- `te_core_interspersed`
- `SINE_young_milliDiv200`
- `SINE_old_milliDiv200plus`
- `LINE_old_milliDiv200plus`

Each filter set specifies:

- which repeat class to include
- optional family restrictions
- optional `max_milliDiv`
- optional `min_milliDiv`
- minimum repeat length

### Step 2. Apply the chosen filter to the UCSC RepeatMasker table

Filtering is done in:

- `scripts/downstream_h2_regulatory_repeat_helpers.R:843-876`

This function applies, in order:

- excluded classes
- included classes
- included families
- `milliDiv <= max_milliDiv` if present
- `milliDiv >= min_milliDiv` if present
- `rep_len_bp >= min_len_bp`

For example, `SINE_old_milliDiv200plus` means:

- keep only `repClass_norm == "SINE"`
- keep only `milliDiv >= 200`
- keep only `rep_len_bp >= 100`

### Step 3. Measure how much of the genome the filtered set occupies

This happens in:

- `scripts/downstream_h2_regulatory_repeat_workflows.R:218-235`

For each filter set, the code computes:

- `observed_bp = sum(rep_len_bp)`
- `observed_fraction = observed_bp / genome_bp`

where `genome_bp` is the sum of chromosome sizes from `chromInfo`.

### Step 4. Convert observed burden into a simulation target fraction

This happens in:

- `scripts/downstream_h2_regulatory_repeat_workflows.R:237-277`

The code first chooses a reference filter set:

- by default `te_core_interspersed`

Then it rescales each filter set using:

`reference_scale = repeat_bg_repeat_fraction / reference_fraction`

and:

`target_fraction = observed_fraction * reference_scale`

with a maximum cap of `0.95`.

This is important:

the simulated background for a filter set is not simply "simulate the exact
observed number of elements." Instead, it is a scaled random-genome burden
derived from the observed fraction of that filter set relative to the reference
TE set.

### Step 5. Build a simulation catalog from the filtered subset

This is done in:

- `scripts/downstream_h2_regulatory_repeat_helpers.R:1156-1234`

For the filtered repeat subset, the code records:

- `type_stats`
- `chrom_prob`
- `length_pools`

`type_stats` contains, for each `repeat_type`:

- `obs_n`
- `obs_bp`
- `mean_len_bp`
- `type_prob`

Because `repeat_bg_type_col = "repClass_norm"`, `repeat_type` is the normalized
repeat class.

Because `repeat_bg_composition_basis = "bp"`, the class weights are:

`type_prob = class_total_bp / total_bp_across_types`

The chromosome probabilities are also based on base pairs:

`chrom_prob = type_bp_on_chr / total_type_bp`

The length pools are the empirical repeat lengths observed in the filtered set.

### Step 6. Simulate a randomized repeat genome

This is done in:

- `scripts/downstream_h2_regulatory_repeat_helpers.R:1237-1322`

For one iteration:

1. compute `target_total_bp = genome_bp * repeat_fraction`
2. split that target bp across repeat types using `type_prob`
3. for each type:
   - sample lengths from the empirical `length_pools`
   - sample chromosomes using `chrom_prob`
   - sample random start positions uniformly within chromosome bounds
4. combine all simulated intervals into one simulated repeat table

The critical point is that the simulation preserves:

- how long repeats tend to be
- which chromosomes they tend to occupy
- how much bp each repeat type contributes

but it removes the original loci.

### Step 7. Count simulated repeats near genes

This is done in:

- `scripts/downstream_h2_regulatory_repeat_workflows.R:283-390`

For each simulated genome:

- construct `GRanges` for the simulated repeats
- count overlaps with each gene window:
  - `100kb`
  - `250kb`
  - `500kb`
  - `1mb`
- summarize those counts by `post_mean_bin`

The summarized metrics are:

- `bg_mean`
- `bg_median`
- `bg_prop_nonzero`

### Step 8. Aggregate across 100 randomized genomes

This is done in:

- `scripts/downstream_h2_regulatory_repeat_workflows.R:523-607`

For each filter set, window, and `post_mean_bin`, the code computes:

- mean background mean
- 2.5th and 97.5th percentiles of the background mean
- mean background median
- 2.5th and 97.5th percentiles of the background median
- mean background nonzero proportion
- 2.5th and 97.5th percentiles of that nonzero proportion

It then compares observed versus background and writes enrichment tables.

## 5. How `old_milliDiv200plus` genomic background is calculated

Take `SINE_old_milliDiv200plus` as the example.

### Biological interpretation

`milliDiv` is the RepeatMasker divergence from the consensus sequence. Larger
values indicate more diverged, and therefore operationally "older," copies.

So:

- `SINE_young_milliDiv200` means `milliDiv <= 200`
- `SINE_old_milliDiv200plus` means `milliDiv >= 200`

### Code-level interpretation

For `SINE_old_milliDiv200plus`, the background uses only the old SINE subset:

1. filter `rmsk` to old SINEs
2. take the empirical old-SINE length distribution
3. take the empirical chromosome occupancy of old SINEs
4. scale the total simulated bp using the reference-set scaling rule
5. place those old-SINE-like intervals at random coordinates
6. count how many fall inside each gene window

Because `repeat_bg_type_col = "repClass_norm"`, this filter set usually has one
effective repeat type: `SINE`.

That means the simulation for `SINE_old_milliDiv200plus` preserves:

- old-SINE lengths
- old-SINE chromosome distribution
- old-SINE total burden after scaling

but not:

- the true old-SINE loci
- local insertion hotspots
- exact family-level structure within SINE unless you change `type_col`

## 6. Pseudocode for the current explicit-genome null

```r
for each filter_set:
  filtered_tbl <- apply_repeat_filter_profile(rmsk, filter_set)
  observed_fraction <- sum(filtered_tbl$rep_len_bp) / genome_bp
  target_fraction <- observed_fraction * reference_scale

  catalog <- build_repeat_sim_catalog(
    filtered_tbl,
    chrom_sizes_dt,
    type_col = "repClass_norm",
    composition_basis = "bp"
  )

  for iter in 1:n_iter:
    sim_dt <- simulate_repeat_genome_dt(
      catalog = catalog,
      chrom_sizes_dt = chrom_sizes_dt,
      repeat_fraction = target_fraction,
      seed = deterministic_seed
    )

    for each window in c("100kb", "250kb", "500kb", "1mb"):
      bg_counts <- countOverlaps(gene_window, sim_gr)
      summarize bg_counts by post_mean_bin
```

## 7. What output files are produced?

The simulation machinery writes:

- `repeat_background_simulation_plan.tsv`
- `repeat_background_seed_info.tsv`
- `repeat_background_seed_map.tsv`
- `repeat_filtered_background_iter.tsv`
- `postmean_repeat_filtered_background_summary.tsv`
- `postmean_repeat_filtered_enrichment.tsv`

These are generated from:

- `scripts/downstream_h2_regulatory_repeat_workflows.R:193-277`
- `scripts/downstream_h2_regulatory_repeat_workflows.R:505-607`

## 8. Important caveat

The current null is a broad genomic-background null, not an exact matched null.

It preserves broad repeat properties, but it does not preserve:

- family-level insertion preference unless you simulate by family
- GC content
- mappability
- local chromatin context
- gene-density context
- exact repeat count per chromosome or per gene

So this background is appropriate for asking:

"Is repeat burden near genes stronger than expected under random placement of
repeat-like elements with the same broad composition?"

It is not the same as asking:

"Is this gene enriched relative to an exactly matched local genomic null?"

## 9. If you want a stricter null later

There are three natural variants you could switch to:

- exact observed count instead of scaled bp fraction
- simulate by `repFamily_norm` instead of `repClass_norm`
- chromosome-stratified exact-count simulation per filter set

Those would produce a more restrictive background, especially for old vs young
subsets.
