# covid-crispr-eqtl: inflation analysis

Quantifying the inflation in COVID-19 GWAS data for the GTEx lung eVariants for the top-ranked genes from the CRISPR screen

## Data

* Top-ranked genes in the CRISPR screen - includes top 500 genes ranked by three methods (RRA, RIGER, SB score) in low MOI experiment from [Daniloski et al. 2021, Cell](https://doi.org/10.1016/j.cell.2020.10.030), Figure S1C
* COVID-19 GWAS summary statistics by the [COVID-19 Host Genetics Initiative](https://www.covid19hg.org/) (round 5)
* cis-eQTLs in Lung [GTEx v8](https://gtexportal.org/home/)

## Analysis

```bash
summary_covid_gwas_round5.Rmd 500
```
