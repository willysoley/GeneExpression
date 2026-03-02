# Downstream h2 / s_het Regulatory + Repeat Analysis Notes

This note documents the data sources and analysis design used in:
`scripts/downstream_h2_regulatory_repeat_analysis.R`.

## Why these resources

1. **ENCODE SCREEN cCREs (hg38)**
- Widely used registry of candidate cis-regulatory elements built from DNase accessibility and histone marks.
- Used here for:
  - enhancer counts near genes (`GRCh38-cCREs.ELS.bed`)
  - open-chromatin proxy counts near genes (`GRCh38-cCREs.bed`, all cCREs)

2. **UCSC RepeatMasker (`rmsk`)**
- Standard genome-wide repeat annotation used in many repeat/transposon analyses.
- Used here to compute repeat burden near gene TSS by total count and class-specific count.
- Clarification: repeat features are defined as **counts of repeat intervals inside fixed windows around each gene TSS** (for example, +/-100 kb and +/-1 Mb), not whether the gene body itself overlaps a repeat.

3. **GEUVADIS SDRF + TPM matrix**
- European run filtering from SDRF, then keep genes with TPM > 0 in at least one European sample.

4. **Cis-window conventions**
- The script computes both +/-100 kb (local promoter-proximal) and +/-1 Mb (typical cis-regulatory window) feature counts.

## Method choices mirrored from common practice

1. **Cis region around TSS**
- +/-1 Mb around TSS is a common cis window for expression-genetic analyses (GTEx-style cis definition).

2. **Enhancer/open-chromatin as gene-level burden**
- Count regulatory elements overlapping fixed windows around each gene TSS, then test association with `h2_GREML` and `Pval_GREML < 0.05`.
- This simple burden model is commonly used as a first-pass before more complex enhancer-gene linking frameworks (e.g., ABC maps).

3. **Repeat classes**
- Aggregate repeat counts by `repClass` (LINE, SINE, LTR, DNA, etc.) within TSS windows and model class-specific associations with heritability.
- The analysis does **not** test a binary gene-body overlap with repeats.

## Primary sources / links

1. ENCODE SCREEN download portal (cCRE BED files, GRCh38):
- https://screen.wenglab.org/downloads

2. ENCODE cCRE reference paper:
- Moore et al., *Nature* 2020 (Expanded encyclopaedias of DNA elements in the human and mouse genomes)
- https://www.nature.com/articles/s41586-020-2493-4

3. Enhancer-gene linking at scale (ABC framework):
- Nasser et al., *Nature* 2021
- https://www.nature.com/articles/s41586-021-03446-x

4. GTEx cis-window example (variants within 1 Mb of TSS):
- GTEx Consortium, *Nature* 2017
- https://pmc.ncbi.nlm.nih.gov/articles/PMC5776756/

5. UCSC Genome Browser downloads (hg38 database directory):
- https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/

6. UCSC `rmsk` table description:
- https://hgdownload.cse.ucsc.edu/goldenpath/help/bigRmsk.html

7. Roadmap Epigenomics (alternative enhancer states via ChromHMM 15-state):
- https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html
- https://www.nature.com/articles/nature14248

8. FANTOM5 enhancer atlas (alternative enhancer resource):
- https://fantom.gsc.riken.jp/5/
- https://www.nature.com/articles/nature12787

## Optional extensions

1. Replace default open-chromatin proxy (`GRCh38-cCREs.bed`) with **cell-type-specific ATAC/DNase narrowPeak BED** files and rerun the same counting/model code.
2. Add distance-weighted burden instead of raw counts (e.g., weight by `1/(distance_to_TSS + 1)`).
3. Swap proximity-based mapping with precomputed enhancer-gene links (ABC or SCREEN links) for a stricter assignment model.
