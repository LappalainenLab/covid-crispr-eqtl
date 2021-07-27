#!/usr/bin/env Rscript
#-------------------------------------------
# Get molQTL sum stats for the chr3 locus
# from the eQTL Catalogue
#-------------------------------------------

library(here)
library(seqminer)
library(data.table)

method <- commandArgs(trailingOnly = TRUE) # ge, microarray, tx, txrev

# Genes of interest in the chr3 locus
gene_df <- data.frame("gene_id" = c("ENSG00000144791", "ENSG00000211456", "ENSG00000163817", "ENSG00000163818", "ENSG00000173585", "ENSG00000163820", "ENSG00000172215", "ENSG00000173578"),
                      "gene_name" = c("LIMD1", "SACM1L", "SLC6A20", "LZTFL1", "CCR9", "FYCO1", "CXCR6", "XCR1"),
                      stringsAsFactors = F)

# Path to eQTL Catalogue data
tabix_paths <- read.table("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
tabix_paths <- tabix_paths[tabix_paths$quant_method %in% c("ge", "microarray", "tx", "txrev"),]

# eQTL summary stats for GTEx tissues by GTEx, but include tx and txrev events from the eQTL Catalogue
tabix_paths <- tabix_paths[!(tabix_paths$study %in% "GTEx" & tabix_paths$quant_method %in% c("ge", "exon")),]
table(tabix_paths$quant_method)

# Extract column names from first file
##column_names <- fread(tabix_paths$ftp_path[1], header = F, sep = "\t", stringsAsFactors = F, nrows = 1, data.table = F)
##column_names <- as.character(column_names[1,])
column_names <- c("molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", "ma_samples", "maf", "pvalue", "beta", "se", "type", "ac", "an", "r2", "molecular_trait_object_id", "gene_id", "median_tpm", "rsid")

# Get data for the chr3 locus genes for the selected method
cat("Get data for quant_method", method, fill = TRUE)
tabix_paths <- tabix_paths[tabix_paths$quant_method %in% method, ]
table(tabix_paths$quant_method)

# COVID-19 HGI Freeze 5 lead variant - rs10490770 - chr3:45823240:T:C
lead_gwas <- "3:45823240:T:C"
gwas_chr <- as.numeric(unlist(strsplit(lead_gwas, ":"))[1])
gwas_pos <- as.numeric(unlist(strsplit(lead_gwas, ":"))[2])
cat("----------------", fill = TRUE)

for (i in 1:nrow(tabix_paths)) {
  ftp_path =  tabix_paths$ftp_path[i]
  cat(ftp_path =  tabix_paths$ftp_path[i], fill = T)

  # Get eQTL data
  eqtl_range <- paste0(gwas_chr, ":",  gwas_pos - 1.5e6, "-", gwas_pos + 1.5e6)
  eqtl_data <- tabix.read.table(tabixFile = ftp_path, tabixRange = eqtl_range, stringsAsFactors = FALSE)

  stopifnot(ncol(eqtl_data) == length(column_names))
  colnames(eqtl_data) <- column_names

  # Only selected genes
  eqtl_data <- eqtl_data[eqtl_data$gene_id %in% gene_df$gene_id, ]

  # For tx and txrevise events, choosing one of the molecular_trait_id to represent the gene
  if (method %in% c("tx", "txrev")) {
    f <- sub("all.tsv", "permuted.tsv", ftp_path)
    if (grepl("GEUVADIS", tabix_paths$ftp_path[i])) {
      f <- sub("_tx_|_txrev_", ".", f)
    }
    permuted <- fread(f, header = T, sep = "\t", stringsAsFactors = F, data.table = F)
    permuted$gene_id <- sapply(permuted$molecular_trait_object_id, function(x) unlist(strsplit(x, "[.]"))[1])
    permuted <- permuted[permuted$gene_id %in% gene_df$gene_id,]
    stopifnot(permuted$molecular_trait_id %in% eqtl_data$molecular_trait_id)
    eqtl_data <- eqtl_data[eqtl_data$molecular_trait_id %in% permuted$molecular_trait_id, ]
    rm(permuted)
  }

  # Save
  fwrite(eqtl_data, file = here("covid_gwas", "freeze5", "chr3_coloc", "input", "eQTL_catalogue", method, paste0(tabix_paths$study[i], "_", tabix_paths$quant_method[i], "_", tabix_paths$qtl_group[i], ".chr3_locus.txt.gz")), col.names = T, row.names = F, sep = "\t", quote = F)

}

cat("----------------", fill = TRUE)
cat("Done!")
