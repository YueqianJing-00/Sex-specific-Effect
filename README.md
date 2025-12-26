# Sex-specific-Effect

We will present a toy example for using the whole pipline, using Z87 as an example, users can replace it with any trait.

First download Summary Statistics from Neale lab, save it in the folder Primary_Summary_Statistics, and we have Primary_Summary_Statistics/Z87/Z87.gwas.imputed_v3.both_sexes.tsv.bgz, etc.

Then use the following code to clean and organize the summary statistics:
```
python clean_sumstat.py Z87
```
