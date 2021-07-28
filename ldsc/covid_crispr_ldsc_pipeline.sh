# We jointly modeled the annotation corresponding to the top-ranked genes and the 96 annotations in the "baseline model" (baselineLD_v2.2.tgz, https://alkesgroup.broadinstitute.org/LDSCORE/GRCh38/).
# We have used regression weights (weights.tgz, https://alkesgroup.broadinstitute.org/LDSCORE/GRCh38/) from HapMap3 SNPs, excluding the HLA region and restricted to SNPs with minor allele frequency (MAF) > 5% (plink_files.tgz, https://alkesgroup.broadinstitute.org/LDSCORE/GRCh38/).
# COVID19-hg GWAS meta-analyses round 5 summary statistics were downloaded from: https://www.covid19hg.org/results/r5/
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Step 1: Making annotation using the top-ranked genes identified in the CRISPR screen (n = 890 protein-coding genes across the three ranking methods) and using all protein-coding
# genes as background to add a column of 0 or 1, this column serves as annotation. Here we are using a window size of 100KB.
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash
#SBATCH --array=1-22

# Activate conda environment with specific version of python and dependencies for the LDSC package to work.
source activate ldsc
chr=$SLURM_ARRAY_TASK_ID
file_name=top_500_genes_SB_RIGER_RRA_NBC_BL_100KB
root_path="/ldsc/LDSCORE_GRCh38/plink_files"
annot_path="/anno"

# Make directory only if it doesn't exist already.
if [ ! -d $annot_path/${file_name}/ ]; then mkdir $annot_path/${file_name}/; fi

python $HOME/ldsc/make_annot.py \
--gene-set-file $annot_path/${file_name}.txt \
--gene-coord-file $annot_path/gene_annot_pcg_ensemblID.txt \
--windowsize 100000 \
--bimfile $root_path/1000G.EUR.hg38.${chr}.bim \
--annot-file $annot_path/${file_name}/${file_name}.chr${chr}.txt

# Remove duplicate SNPs.
awk '!x[$0]++' $annot_path/${file_name}/${file_name}.chr${chr}.txt | \
awk -v col_name="$file_name" '{OFS="\t"; if (NR==1) {print $1,$2,$3,$4,col_name} else {print $1,$2,$3,$4,$5} }' | \
gzip > $annot_path/${file_name}/${file_name}.chr${chr}.annot.gz

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Step 2: Making LD scores based on the annotations we created in Step 1.
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash
#SBATCH --array=1-22

# Activate conda environment with specific version of python and dependencies for the LDSC package to work.
source activate ldsc
chr=$SLURM_ARRAY_TASK_ID
file_name=top_500_genes_SB_RIGER_RRA_NBC_BL_100KB
root_path="/ldsc/LDSCORE_GRCh38"
annot_path="/anno/${file_name}"

echo $chr
echo "Calculating ldscore for chr${chr}"
python $HOME/ldsc/ldsc.py \
--l2 \
--bfile $root_path/plink_files/1000G.EUR.hg38.${chr} \
--ld-wind-cm 1 \
--annot $annot_path/${file_name}.chr${chr}.annot.gz \
--out $annot_path/${file_name}.chr${chr} \
--print-snps $root_path/weights/allsnps.list

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Step 3: GWAS summary statistics are usually not in the format that LDSC understands. munge_sumstats.py script converts summary statistics into the LDSC format.
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash

# Activate conda environment with specific version of python and dependencies for the LDSC package to work.
source activate ldsc
gwas_path="/COVID19_HGI_gwas/20210107"
ldscore_path="/ldsc/LDSCORE_GRCh38/weights/"
c2=${gwas_path}/COVID19_HGI_C2_ALL_eur_leave_23andme_20210107
b2=${gwas_path}/COVID19_HGI_B2_ALL_eur_leave_23andme_20210107
a2=${gwas_path}/COVID19_HGI_A2_ALL_eur_leave_23andme_20210107

for gwas_name in $c2 $b2 $a2
do
  # Rename column names of the file for munge_sumstats.py to work.
  sed -i -e '1s/#CHR/CHR/' -e '1s/all_inv_var_meta_beta/Beta/' -e '1s/all_inv_var_meta_sebeta/se/' $gwas_name.txt
  sed -i -e '1s/all_inv_var_meta_p/P/' -e '1s/all_meta_sample_N/N/' $gwas_name.txt
  sed -i -e '1s/all_meta_AF/EAF/' $gwas_name.txt

  echo "Formatting ${gwas_name}.txt"
  python $HOME/ldsc/munge_sumstats.py \
  --sumstats ${gwas_name}.txt \
  --out ${gwas_name}.Formatted \
  --merge-alleles ${ldscore_path}/w_hm3.noMHC.snplist \
  --a1 ALT \
  --a2 REF \
  --snp rsid \
  --chunksize 500000
done

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Step 4: Regressing GWAS summary statistics on LD scores generated from the baseline model and our gene list.
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash

# Activate conda environment with specific version of python and dependencies for the LDSC package to work.
source activate ldsc
root_path="/ldsc/LDSCORE_GRCh38"
gwas_path="/COVID19_HGI_gwas/20210107"
output_path="/anno/top_500_genes_SB_RIGER_RRA_NBC_BL_100KB"
file_name=top_500_genes_SB_RIGER_RRA_NBC_BL_100KB
c2=${gwas_path}/COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.Formatted.sumstats.gz
b2=${gwas_path}/COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.Formatted.sumstats.gz
a2=${gwas_path}/COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.Formatted.sumstats.gz

# Make results directory only if it doesn't exist already.
if [ ! -d $output_path/results/ ]; then mkdir $output_path/results/; fi

for gwas_name in $c2 $b2 $a2
do
  echo "Analysing ${gwas_name}"
  python $HOME/ldsc/ldsc.py \
  --h2 ${gwas_name} \
  --ref-ld-chr ${root_path}/baselineLD_v2.2/baselineLD.,${output_path}/${file_name}.chr \
  --w-ld-chr ${root_path}/weights/weights.hm3_noMHC.\
  --overlap-annot \
  --frqfile-chr ${root_path}/plink_files/1000G.EUR.hg38.\
  --out $output_path/results/${gwas_name}
  --print-coefficients
done