# covid-crispr-eqtl: colocalization analysis

## Data

* Using eQTLs and sQTLs from [GTEx v8](https://gtexportal.org/home/)
* Using molQTLs from the [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/)
* Using COVID-19 GWAS summary statistics by the [COVID-19 Host Genetics Initiative](https://www.covid19hg.org/) (round 5)

## Prepare data for analyses

Run colocalization analysis for

* B2_ALL_leave_23andme
* C2_ALL_leave_23andme

COVID-19 HGI Freeze 5 lead variant for severity is rs10490770 - 3:45823240:T:C

```bash
# Prepare metadata file - metadata_freeze5.txt
## code phenotype   population  n_cases n_controls  type    file
## B2_ALL_leave_23andme Hospitalized covid vs. population   All 12888   1295966 cc  COVID19_HGI_B2_ALL_leave_23andme_20210107.txt.gz
## C2_ALL_leave_23andme Covid vs. population    All 36590   1668938 cc  COVID19_HGI_C2_ALL_leave_23andme_20210107.txt.gz
```

Get summary stats for the chr3 locus for the selected GWASes:

```bash
get_covid19_chr3_sumstats.R
```

Get summary stats for the chr3 locus for molQTLs from the eQTL Catalogue:

```bash
METHOD="ge" #"ge" "microarray" "tx" "txrev"
get_eqtl_catalogue_chr3_sumstats.R ${METHOD}
```

## Colocalization analysis with [coloc](https://github.com/chr1swallace/coloc)

Run coloc assuming one causal variant per locus:

* For eQTLs, sQTLs and tuQTLs

```bash
GWAS="B2_ALL_leave_23andme" #C2_ALL_leave_23andme
QTL="GTEx_v8" #eQTL_catalogue GTEx_v8_sqtl eQTL_catalogue_tx eQTL_catalogue_txrev)

coloc_covid19_gwas_eqtls.R ${QTL} ${GWAS}
```

Run coloc-cond/mask:

* For eQTLs, sQTLs and tuQTLs

```bash
GWAS="B2_ALL_leave_23andme" #C2_ALL_leave_23andme
QTL="GTEx_v8 eQTL_catalogue" #GTEx_v8_sqtl eQTL_catalogue_tx eQTL_catalogue_txrev

coloc_covid19_gwas_eqtls.condmask.R ${QTL} ${GWAS}
```

Summary of colocalization analysis:

```bash
# B2 analysis and eQTLs
summary_coloc_chr3_genes.Rmd

# B2 analysis and other molQTLs
summary_coloc_chr3_genes_other_molqtls.Rmd
```

## Colocalization analysis with [Joint Likelihoood Mapping (JLIM)](https://github.com/cotsapaslab/jlim)

Primary trait summary statistics file:

* _Named by TraitName.CHR.STARTBP.ENDBP.tx_
* _Space-delimited and has CHR, BP, and SNP. It also has to carry one of STAT, T, or P. If P is specified, the two-sided P-value will be transformed into Z-statistic. STAT or T will be approximated as Z-statistic._

```bash
GWAS='B2_ALL_eur' #C2_ALL_eur
CHR='3'
STARTBP='45323240'
ENDBP='46323240'

cat covid_gwas/freeze5/chr3_coloc/input/${GWAS}_chr3_locus.txt | cut -f 1,2,5,9 | grep -v "CHR" | sed 's/\t/ /g' | awk 'BEGIN{OFS=" ";FS=" "} {if(NR==1){print("CHR","BP","SNP","P")}; print}' > covid_gwas/freeze5/chr3_coloc/input/${GWAS}.${CHR}.${STARTBP}.${ENDBP}.txt
```

Specify the indexSNP file:

```bash
# covid_gwas/freeze5/chr3_coloc/input/chr3-indexSNP.tsv
# CHR     SNP     BP      STARTBP ENDBP
# 3       3:45823240:T:C      45823240        45323240        46323240
```

Generate the reference LD file:

* _Reference LD files are provided one for each interval. It is a tab-delimited file without header. The file name is specified as locus.CHR.STARTBP.ENDBP.txt.gz. Each row is a marker, and it contains the following columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, and is followed by two alleles for each individual._

```bash
module load tabix
module load vcftools/0.1.14
# export PERL5LIB=/nfs/sw/vcftools/vcftools-0.1.14/perl/
PATH_1KG='~/lab/data/1kg'

# List of individuals in the 1000G files that only includes unrelated individuals
tabix -h ${PATH_1KG}/phase3_GRCh38/ALL.chr3.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz 3:1-2 | grep -w "#CHROM" | sed 's/\t/\n/g' > covid_gwas/freeze5/chr3_coloc/input/ld0/header_1000g.txt
# Pull out these from the 1000G populations file
cat ${PATH_1KG}/20130606_sample_info_sheet1.txt | grep -f covid_gwas/freeze5/chr3_coloc/input/ld0/header_1000g.txt >> covid_gwas/freeze5/chr3_coloc/input/ld0/indiv_unrelated_1000g.txt
# List of individuals to include from superpopulations EUR - CEU,TSI,FIN,GBR,IBS - n = 522
cat covid_gwas/freeze5/chr3_coloc/input/ld0/indiv_unrelated_1000g.txt | grep -E 'CEU|TSI|FIN|GBR|IBS' | cut -f 1 > covid_gwas/freeze5/chr3_coloc/input/ld0/indiv_unrelated_1000g_eur.txt
# paste to the get_ref_ld_eur_1000g.pl script
sed 'H;1h;$!d;x;y/\n/,/' covid_gwas/freeze5/chr3_coloc/input/ld0/indiv_unrelated_1000g_eur.txt

# Prepare LD file - 27,208 variants in the region
get_ref_ld_eur_1000g.pl ${PATH_1KG}/phase3_GRCh38/ covid_gwas/freeze5/chr3_coloc/input/chr3-indexSNP.tsv covid_gwas/freeze5/chr3_coloc/input/ld0/

# Add SNP IDs as chr:pos:ref:alt
zcat covid_gwas/freeze5/chr3_coloc/input/ld0/locus.3.45323240.46323240.txt.gz | awk 'BEGIN{OFS="\t";FS="\t"} {new_var=$1":"$2":"$4":"$5; $3=new_var; print}' > covid_gwas/freeze5/chr3_coloc/input/ld0/locus.3.45323240.46323240.txt
rm covid_gwas/freeze5/chr3_coloc/input/ld0/locus.3.45323240.46323240.txt.gz
gzip covid_gwas/freeze5/chr3_coloc/input/ld0/locus.3.45323240.46323240.txt
```

Run JLIM:

* eQTLs from GTEx and the eQTL Catalogue

```bash
DIR='~/projects/covid_crispr'
QTL='eQTL_catalogue' # 'GTEx_v8' 'eQTL_catalogue'
if [ ${QTL} == 'GTEx_v8' ]; then
    REF_DB='GTEx.v8.EUR'
    DIR_QTL='~/lab/data/gtex/v8/eqtl/GTEx_Analysis_v8_EUR_eQTL_all_associations'
    declare -a TISSUE=($(ls ${DIR_QTL}))
    echo ${#TISSUE[@]}
elif [ ${QTL} == 'eQTL_catalogue' ]; then
    REF_DB='eQTLCatalogue'
    # Use "ge" or "microarray" in the DIR_QTL
    DIR_QTL=${DIR}/covid_gwas/freeze5/chr3_coloc/input/eQTL_catalogue/microarray
    declare -a TISSUE=($(ls ${DIR_QTL} | sed 's/.chr3_locus.txt.gz//g'))
    echo ${#TISSUE[@]}
else
   echo "QTL can be GTEx_v8 or eQTLCatalogue"
fi
GWAS='B2_ALL_eur' #'C2_ALL_eur'
CHR='3'
BP='45823240'
STARTBP='45323240'
ENDBP='46323240'

for TISS in ${TISSUE[@]}; do
    echo ${TISS}
    if [ ${QTL} == 'GTEx_v8' ]; then
        SECTR_FILE=${DIR_QTL}/${TISS}/${TISS}.v8.EUR.allpairs.chr${CHR}.parquet
    else
        SECTR_FILE=${DIR_QTL}/${TISS}.chr3_locus.txt.gz
    fi
    # Run JLIM
    ~/programs/jlim/run_jlim_mod.sh --maintr-file ${DIR}/covid_gwas/freeze5/chr3_coloc/input/${GWAS}.${CHR}.${STARTBP}.${ENDBP}.txt \
        --sectr-file ${SECTR_FILE} \
        --ref-ld ${DIR}/covid_gwas/freeze5/chr3_coloc/input/ld0/locus.${CHR}.${STARTBP}.${ENDBP}.txt.gz \
        --index-snp ${CHR}:${BP} \
        --output-file ${DIR}/covid_gwas/freeze5/chr3_coloc/result_jlim/${QTL}/${GWAS}_${QTL}_${TISS}.out \
        --sectr-ref-db ${REF_DB} > ${DIR}/log_jlim_${GWAS}_${QTL}_${TISS}.out 2> ${DIR}/log_jlim_${GWAS}_${QTL}_${TISS}.error &
done
```

Compare JLIM and coloc results for eQTLs:

```bash
compare_with_coloc_jlim.Rmd
```
