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

