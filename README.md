# covid-crispr-eqtl

Integration of genome-scale CRISPR loss-of-function screens and eQTLs in diverse cell types and tissues to pinpoint genes underlying COVID-19 risk in the 3p21.31 locus

This repository contains code to conduct the analyses presented in the article "_Integrative approach identifies_ SLC6A20 _and_ CXCR6 _as putative causal genes for the COVID-19 GWAS signal in the 3p21.31 locus_" (Kasela et al., 2021), see the preprint in [medRxiv](https://doi.org/10.1101/2021.04.09.21255184)).

## Overview of this repository

* [crispr_genes](https://github.com/ainenLab/covid-crispr-eqtl/tree/master/crispr_genes) - script to explore the top-ranked genes from the CRISPR screen from [Daniloski et al. 2021, Cell](https://doi.org/10.1016/j.cell.2020.10.030)
* [covid_gwas](https://github.com/LappalainenLab/covid-crispr-eqtl/tree/master/covid_gwas) - scripts to perform inflation of COVID-19 GWAS signal for variants that are cis-eQTLs in Lung
* [eqtl_coloc](https://github.com/LappalainenLab/covid-crispr-eqtl/tree/master/eqtl_coloc) - scripts to perform colocalization analysis of molQTLs and COVID-19 GWAS signal in the the 3p21.31 locus
* [ldsc](https://github.com/LappalainenLab/covid-crispr-eqtl/tree/master/ldsc) - scripts to run the LDSC pipeline
