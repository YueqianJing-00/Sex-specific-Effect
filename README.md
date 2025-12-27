# Sex-specific Effect

This repository provides an end-to-end pipeline. Below is a toy example using **UK Biobank phenotype `Z87`** as the input trait. You can replace `Z87` with any other supported trait code.

## 1) Download primary summary statistics

Download the GWAS summary statistics from the Neale Lab and place them under:

```text
Primary_Summary_Statistics/
```

For the toy example (`Z87`), the expected directory structure is:

```text
Primary_Summary_Statistics/Z87/
  ├── Z87.gwas.imputed_v3.both_sexes.tsv.bgz
  ├── Z87.gwas.imputed_v3.female.tsv.bgz
  └── Z87.gwas.imputed_v3.male.tsv.bgz
```

## 2) Clean and organize summary statistics

Run:

```bash
python clean_sumstat.py Z87
```

This gives you the cleaned summary statistics with column SNP, beta, se, tstat, P, A1 (minor allele), A2 (reference allele), etc.

```text
Primary_Summary_Statistics/Z87/
  ├── Z87.sumstats.both_sexes.tsv
  ├── Z87.sumstats.male.tsv
  └── Z87.sumstats.female.tsv
```
## 3) Pairwise genetic correlation and trait grouping

Pairwise genetic correlations among 733 traits are computed using LDSC (filtering trait pairs with |r_g| ≥ 0.5):

```bash
python ldsc.py \
  --ref-ld-chr ./ldsc/eur_w_ld_chr/ \
  --w-ld-chr ./ldsc/eur_w_ld_chr/ \
  --rg Z87.ldsc.both_sexes.sumstats.gz,../Z85/Z85.ldsc.both_sexes.sumstats.gz \
  --out ../../Genetic_Correlations/Z87-Z85
```

Based on the resulting genetic correlation matrix, the 733 traits are grouped into 270 trait groups (saved under `traitgroups/`).

For each trait group, we construct group-level summary statistics by taking, for each variant, the minimum P-value across all traits within the group.

## 4）Clumping loci within a trait group (PLINK 1.9)

Assuming `Z87` belongs to trait group `1`, we perform LD clumping on the trait-group summary statistics using PLINK 1.9:

```bash
plink --bfile 1000G \
  --clump Trait_Groups/1/1.sumstats.tsv \
  --clump-kb 500 \
  --clump-p1 5e-8 \
  --clump-p2 1 \
  --clump-r2 0.1 \
  --out Trait_Groups/1/1
```


## 5) Group SNPs by clumped loci and extract sex-specific signals

This step uses the PLINK clumping results to group SNPs into loci, then outputs locus-annotated tables for each trait in the selected trait group.

Run (replace `<TRAITGROUP_NUM>` with an integer, e.g. `1`):

```bash
python primary_traitgroup_sse.py Z87
```

**Outputs (per trait, saved under `Primary_Summary_Statistics/<TRAIT>/`)**
- `<TRAIT>.clumped_all.tsv`: all SNPs grouped by locus
- `<TRAIT>.clumped_SSE.tsv`: sex-specific SNPs (genome-wide significant in one sex only); includes `p_diff_adj` (BH-FDR) when available

> Note: the script currently uses hard-coded absolute paths (e.g., `/gpfs/...`). Update those paths to match your environment before running.

