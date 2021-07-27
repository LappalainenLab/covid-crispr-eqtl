#!/usr/bin/env Rscript
#-------------------------------------------------------
# Colocalization analysis for COVID-19 GWAS (freeze 5)
# and eQTLs from the eQTL Catalogue or GTEx v8
#   * Using coloc-cond/mask
#-------------------------------------------------------

# Idea: run coloc between COVID-19 GWAS and eQTLs, if there's an eVariant within 100kb of the GWAS variant with nominal p-value < 10-4
# Run coloc-condmask in a 2Mb region centered on the lead GWAS variant

# Use conditioning on GTEx v8 data
# Masking on eQTL Catalogue data, LD data from 1000G EUR population
# Method = 'single' for GWAS

library(here)
library(data.table)
library(seqminer)
library(dplyr)
library(ggplot2)
library(patchwork)
library(coloc)

# Source functions fro coloc-cond/mask
source(here("eqtl_coloc", "functions_run_coloc.R"))
source(here("eqtl_coloc", "functions_plot_locus.R")
source(here("eqtl_coloc", "functions_get_genomatrix.R")
source(here("eqtl_coloc", "functions_coloc_signals_adj.R")

# ggplot2 theme
theme_set(theme_classic() +
            theme(plot.title = element_text(size = 14),
                  plot.subtitle = element_text(size = 12),
                  axis.title = element_text(size = 12),
                  axis.title.y = element_text(vjust = 2),
                  axis.title.x = element_text(vjust = -0.75),
                  axis.text = element_text(size = 11),
                  strip.text = element_text(size = 11),
                  legend.title = element_text(size = 12),
                  legend.text = element_text(size = 11)))

# eQTL Catalogue or GTEx
args <- commandArgs(trailingOnly = TRUE)
eqtls <- args[1] # "GTEx_v8", "GTEx_v8_sqtl" or "eQTL_catalogue", "eQTL_catalogue_tx", "eQTL_catalogue_txrev"
print(eqtls)
gwas_analysis <- args[2] # "B2_ALL_leave_23andme", "C2_ALL_leave_23andme"
print(gwas_analysis)

# coloc-cond/mask
coloc <- "condmask"
print(coloc)
if (eqtls %in% "GTEx_v8") {
  QTL_COLOC <- 'beta,cond'
} else {
  QTL_COLOC <- 'beta,mask'
}
GWAS_COLOC <- "beta,single"
COLOC_MODE <- "iterative"

# Running coloc in a 2Mb region centered at the lead GWAS variant (+/- 1Mb from the SNP)
window <- 1e+06
pval_th <- 1e-4
distance_th <- 100000

# COVID-19 GWAS data (round 5) ----------------------
gwas_data <- fread(here("covid_gwas", "freeze5", "chr3_coloc", "input", paste0(gwas_analysis, "_chr3_locus.txt")), header = T, sep = "\t", stringsAsFactors = F, data.table = F)
gwas_meta <- read.table(here("covid_gwas", "freeze5", "chr3_coloc", "metadata_freeze5.txt"), header = T, sep = "\t", stringsAsFactors = F)
stopifnot(gwas_analysis %in% gwas_meta$code)
gwas_meta <- gwas_meta %>%
  filter(code %in% gwas_analysis)

# COVID-19 HGI Freeze 5 lead variant - rs10490770 - chr3:45823240:T:C
lead_gwas <- "3:45823240:T:C"
#lead_gwas <- gwas_data[which.min(gwas_data$all_inv_var_meta_p), "SNP"]
cat("Lead GWAS variant -", lead_gwas, fill = TRUE)
lead_gwas_id <- gsub(":", "_", paste0("chr", lead_gwas))
gwas_chr <- as.numeric(unlist(strsplit(lead_gwas, ":"))[1])
gwas_pos <- as.numeric(unlist(strsplit(lead_gwas, ":"))[2])

# sample size and type, center on the lead GWAS variant
gwas_data <- gwas_data %>%
  filter(POS > (gwas_pos - window) & POS < (gwas_pos + window)) %>%
  mutate(gwas_type = gwas_meta$type,
         gwas_n = gwas_meta$n_cases + gwas_meta$n_controls,
         gwas_s = gwas_meta$n_cases / gwas_n)
rm(gwas_meta)

# Genes of interest in that region
if (grepl("GTEx_v8", eqtls)) {
  gene_df <- data.frame("gene_id" = c("ENSG00000144791.9", "ENSG00000211456.10", "ENSG00000163817.15", "ENSG00000163818.16", "ENSG00000173585.15", "ENSG00000163820.14", "ENSG00000172215.5", "ENSG00000173578.7"),
                        "gene_name" = c("LIMD1", "SACM1L", "SLC6A20", "LZTFL1", "CCR9", "FYCO1", "CXCR6", "XCR1"),
                        stringsAsFactors = F)
} else {
  gene_df <- data.frame("gene_id" = c("ENSG00000144791", "ENSG00000211456", "ENSG00000163817", "ENSG00000163818", "ENSG00000173585", "ENSG00000163820", "ENSG00000172215", "ENSG00000173578"),
                        "gene_name" = c("LIMD1", "SACM1L", "SLC6A20", "LZTFL1", "CCR9", "FYCO1", "CXCR6", "XCR1"),
                        stringsAsFactors = F)
}
print(gene_df)

# Get eQTL summary stats from GTEx v8 or eQTL Catalogue --------
## Tutorial and functios to work with eQTL Catalogue data from here: http://htmlpreview.github.io/?https://github.com/kauralasoo/eQTL-Catalogue-resources/blob/master/scripts/tabix_use_case.html
if (grepl("GTEx_v8", eqtls)) {
  path_to_files <- ifelse(eqtls %in% "GTEx_v8", "~/gtex_v8/eqtl/GTEx_Analysis_v8_eQTL_all_associations_indexed/results/",
                          ifelse(eqtls %in% "GTEx_v8_sqtl", "~/gtex_v8/sqtl/GTEx_Analysis_v8_sQTL_all_associations_indexed/results/", NA))
  files <- list.files(path = path_to_files)
  files <- files[!grepl("tbi", files)]
  tabix_paths <- data.frame("study" = eqtls,
                            "qtl_group" = as.character(sapply(files, function(x) unlist(strsplit(x, "[.]"))[1])),
                            "ftp_path" = paste0(path_to_files, files),
                            "quant_method" = ifelse(grepl("sqtl", eqtls), "splicing", "ge"),
                            stringsAsFactors = F)
  # Sample size in GTEx
  sample_annot <- read.table("~/lab/data/gtex/v8/sample_annotations/tissue_sample_count.txt", header = F, sep = "\t", stringsAsFactors = F)
  stopifnot(tabix_paths$qtl_group %in% sample_annot$V1)
  tabix_paths$eqtl_n <- sample_annot[match(tabix_paths$qtl_group, sample_annot$V1), "V2"]
  rm(sample_annot)

} else {
  tabix_paths <- read.table("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")

  # eQTL summary stats for GTEx tissues by GTEx, but include tx and txrev events from the eQTL Catalogue
  if (eqtls %in% "eQTL_catalogue") {
    tabix_paths <- tabix_paths %>%
      filter(!study %in% "GTEx",
             quant_method %in% c("ge", "microarray"))
  } else {
    tabix_paths <- tabix_paths %>%
      filter(if (eqtls %in% "eQTL_catalogue_tx") quant_method %in% "tx" else quant_method %in% "txrev")
  }

  ## Extract column names from first file
  ###column_names <- fread(tabix_paths$ftp_path[1], header = F, sep = "\t", stringsAsFactors = F, nrows = 1, data.table = F)
  ###column_names <- as.character(column_names[1,])
  #column_names <- c("molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", "ma_samples", "maf", "pvalue", "beta", "se", "type", "ac", "an", "r2", "molecular_trait_object_id", "gene_id", "median_tpm", "rsid")

  # Pulled out data from the chr3 locus
  tabix_paths$ftp_path <- here("covid_gwas", "freeze5", "chr3_coloc", "input", "eQTL_catalogue", tabix_paths$quant_method, paste0(tabix_paths$study, "_", tabix_paths$quant_method, "_", tabix_paths$qtl_group, ".chr3_locus.txt.gz"))

}

res <- data.frame("gwas" = rep(gwas_analysis, ifelse(grepl("txrev", eqtls), 3000, 1500)), #don't know exactly how many rows there might be in total
                  "lead_gwas_variant" = NA,
                  "lead_gwas_variant_rsid" = NA,
                  "lead_gwas_meta_pval" = NA,
                  "lead_gwas_meta_beta" = NA,
                  "lead_gwas_meta_se" = NA,
                  "eqtl_gene_id" = NA,
                  "eqtl_gene_name" = NA,
                  "eqtl_molecular_trait_id" = NA,
                  "eqtl_data" = NA,
                  "eqtl_n" = NA,
                  "lead_eqtl_variant_larger_window" = NA,
                  "lead_eqtl_pval_larger_window" = NA,
                  "lead_eqtl_variant" = NA,
                  "lead_eqtl_variant_rsid" = NA,
                  "lead_eqtl_pval" = NA,
                  "lead_eqtl_beta" = NA,
                  "lead_eqtl_se" = NA,
                  "lead_gwas_variant_eqtl_pval" = NA,
                  #"coloc_nsnps" = NA,
                  #"coloc_pp0" = NA,
                  #"coloc_pp1" = NA,
                  #"coloc_pp2" = NA,
                  #"coloc_pp3" = NA,
                  #"coloc_pp4" = NA,
                  stringsAsFactors = F)

coloc_result <- list()
k <- 1
for (i in 1:nrow(tabix_paths)) {
  ftp_path =  tabix_paths$ftp_path[i]
  cat(ftp_path =  tabix_paths$ftp_path[i], fill = T)

  if (grepl("eQTL_catalogue", eqtls)) {
    eqtl_data <- fread(ftp_path, header = T, sep = "\t", stringsAsFactors = F, data.table = F)

    #stopifnot(ncol(eqtl_data) == length(column_names))
    #colnames(eqtl_data) <- column_names

  } else {
    eqtl_range <- paste0("chr", gwas_chr, ":",  gwas_pos - 1.5e6, "-", gwas_pos + 1.5e6)
    eqtl_data <- tabix.read.table(tabixFile = ftp_path, tabixRange = eqtl_range, stringsAsFactors = FALSE)

    # Colnames similar to eQTL Catalogue files
    idx <- sapply(c("pval_nominal", "slope", "slope_se", "chr", "pos"), function(x) which(colnames(eqtl_data) %in% x))
    colnames(eqtl_data)[idx] <- c("pvalue", "beta", "se", "chromosome", "position")
    eqtl_data$chromosome <- sub("chr", "", eqtl_data$chromosome)

    # Add ref and alt allele
    eqtl_data$ref <- sapply(eqtl_data$variant_id, function(x) unlist(strsplit(x, "_"))[3])
    eqtl_data$alt <- sapply(eqtl_data$variant_id, function(x) unlist(strsplit(x, "_"))[4])

    # Remove "_b38" at the end of variant
    eqtl_data$variant <- sapply(eqtl_data$variant_id, function(x) sub("_b38", "", x))

    # For sQTLs need to run coloc per clu cluster per gene_id, choosing one of the clusters only to represent the group
    if (eqtls %in% "GTEx_v8_sqtl") {
      sgenes <- read.table(paste0("~/gtex_v8/sqtl/GTEx_Analysis_v8_sQTL/", tabix_paths$qtl_group[i], ".v8.sgenes.txt.gz"), header = T, sep = "\t", stringsAsFactors = F)
      sgenes <- sgenes[sgenes$group_id %in% gene_df$gene_id, ]
      stopifnot(sgenes$gene_id %in% eqtl_data$gene_id)
      eqtl_data <- eqtl_data[eqtl_data$gene_id %in% sgenes$gene_id, ]
      eqtl_data$molecular_trait_id <- eqtl_data$gene_id
      eqtl_data$gene_id <- sapply(eqtl_data$molecular_trait_id, function(x) unlist(strsplit(x, ":"))[5])
      rm(sgenes)
    }
  }

  # Double-check
  stopifnot(substring(eqtl_data$gene_id[1], 1, 4) == "ENSG")

  # Genes of interest
  eqtl_data <- eqtl_data %>%
    filter(gene_id %in% gene_df$gene_id) %>%
    mutate(gene_name = gene_df[match(gene_id, gene_df$gene_id), "gene_name"])

  if (tabix_paths$quant_method[i] %in% "ge") {
    eqtl_data <- eqtl_data %>%
      mutate(trait_id = gene_name) %>%
      mutate(trait_id_f = factor(trait_id, levels = gene_df$gene_name))
  } else {
    comb <- sapply(gene_df$gene_name, function(x) {
      trait_id <- sort(unique(eqtl_data[eqtl_data$gene_name %in% x, "molecular_trait_id"]))
      if (length(trait_id) > 0) {
        return(paste0(x, "-", trait_id))
      } else {
        return(x)
      }
    })
    comb <- as.character(unlist(comb))
    eqtl_data <- eqtl_data %>%
      mutate(trait_id = paste0(gene_name, "-", molecular_trait_id)) %>%
      mutate(trait_id_f = factor(trait_id, levels = comb))
  }

  trait_id <- levels(eqtl_data$trait_id_f)
  # Run gene by gene (gene-molecular trait if microarray)
  # Run coloc if min p-value in the coloc-region is < 10-4
  for (j in 1:length(trait_id)) {
    trait_id_lookup <- trait_id[j]
    gene_name <- unlist(strsplit(trait_id_lookup, "-"))[1]
    molecular_trait_id <- unlist(strsplit(trait_id_lookup, "-"))[2]

    if (trait_id_lookup %in% eqtl_data$trait_id) {
      cat("Gene ", trait_id_lookup, fill = T)

      if (grepl("eQTL_catalogue", eqtls)) {
        # Select one gene, remove rsid duplicates and multi-allelic variants
        sel <- eqtl_data %>%
          filter(trait_id %in% trait_id_lookup) %>%
          select(-rsid) %>%
          distinct() %>% #rsid duplicates
          mutate(id = paste(chromosome, position, sep = ":")) %>%
          group_by(id) %>%
          mutate(row_count = n()) %>%
          ungroup() %>%
          filter(row_count == 1) %>% #Multialllics
          mutate(SNP = paste(chromosome, position, ref, alt, sep = ":")) %>%
          data.frame()
      } else {
        # Select one gene, remove rsid duplicates and multi-allelic variants
        sel <- eqtl_data %>%
          filter(trait_id %in% trait_id_lookup) %>%
          mutate(id = paste(chromosome, position, sep = ":")) %>%
          group_by(id) %>%
          mutate(row_count = n()) %>%
          ungroup() %>%
          filter(row_count == 1) %>% #Multialllics
          mutate(SNP = paste(chromosome, position, ref, alt, sep = ":")) %>%
          data.frame()
      }

      # merge data
      if (sum(sel$SNP %in% gwas_data$SNP) > 100) {
        mdata <- inner_join(x = gwas_data, y = sel, by = "SNP")

        # Using beta + se for coloc - missing values not allowed
        if (sum(is.na(mdata$se)) > 0) {
          cat("Warning: need to exclude", sum(is.na(mdata$se)), "rows with missing eQTL se", fill = TRUE)
          mdata <- mdata %>%
            filter(!is.na(se))
        }

        #res$lead_gwas_variant[k] <- mdata[which.min(mdata$all_inv_var_meta_p), "SNP"]
        #res$lead_gwas_pval[k] <- mdata[which.min(mdata$all_inv_var_meta_p), "all_inv_var_meta_p"]
        res$lead_gwas_variant[k] <- lead_gwas
        res$lead_gwas_variant_rsid[k] <- mdata[mdata$SNP %in% lead_gwas, "rsid"]
        res$lead_gwas_meta_pval[k] <- mdata[mdata$SNP %in% lead_gwas, "all_inv_var_meta_p"]
        res$lead_gwas_meta_beta[k] <- mdata[mdata$SNP %in% lead_gwas, "all_inv_var_meta_beta"]
        res$lead_gwas_meta_se[k] <- mdata[mdata$SNP %in% lead_gwas, "all_inv_var_meta_sebeta"]
        res$eqtl_gene_id[k] <- gene_df[match(gene_name, gene_df$gene_name), "gene_id"]
        res$eqtl_gene_name[k] <- gene_name
        res$eqtl_molecular_trait_id[k] <- ifelse(eqtls %in% "GTEx_v8_sqtl", unique(sel$molecular_trait_id), molecular_trait_id)
        res$eqtl_data[k] <- paste(tabix_paths[i, c("study", "quant_method", "qtl_group")], collapse = "_")
        res$eqtl_n[k] <- ifelse(grepl("eQTL_catalogue", eqtls), tabix_paths$sample_size[i], tabix_paths$eqtl_n[i])
        res$lead_eqtl_variant_larger_window[k] <- sel[which.min(sel$pvalue), "SNP"]
        res$lead_eqtl_pval_larger_window[k] <- sel[which.min(sel$pvalue), "pvalue"]
        lead_eqtl <-  mdata[which.min(mdata$pvalue), "SNP"]
        res$lead_eqtl_variant[k] <- lead_eqtl
        res$lead_eqtl_variant_rsid[k] <- mdata[mdata$SNP %in% lead_eqtl, "rsid"]
        res$lead_eqtl_pval[k] <- mdata[mdata$SNP %in% lead_eqtl, "pvalue"]
        res$lead_eqtl_beta[k] <- mdata[mdata$SNP %in% lead_eqtl, "beta"]
        res$lead_eqtl_se[k] <- mdata[mdata$SNP %in% lead_eqtl, "se"]
        #res$lead_gwas_variant_eqtl_pval[k] <- mdata[which.min(mdata$all_inv_var_meta_p), "pvalue"]
        res$lead_gwas_variant_eqtl_pval[k] <- mdata[mdata$SNP %in% lead_gwas, "pvalue"]

        # Check if there's eQTL variants with P-value < 10-4 within 100kb of the GWAS variant
        min_p_check <- sel %>%
            filter(position > gwas_pos - distance_th & position < gwas_pos + distance_th,
                   pvalue < pval_th)

        # eQTL variants with P-value < 10-4 within 100kb of the GWAS variant
        if (nrow(min_p_check) > 0) {
          # coloc-condmask
          mdata$phenotype_id <- mdata$trait_id
          mdata$id <- sapply(mdata$variant, function(x) paste(unlist(strsplit(x, "_"))[1:2], collapse = "_"))

          mdata$gwas_pval <- mdata$all_inv_var_meta_p
          mdata$gwas_beta <- mdata$all_inv_var_meta_beta
          mdata$gwas_se <- mdata$all_inv_var_meta_sebeta

          mdata$qtl_pval <- mdata$pvalue
          mdata$qtl_beta <- mdata$beta
          mdata$qtl_se <- mdata$se
          mdata$qtl_maf <- mdata$maf
          mdata$qtl_n <- ifelse(grepl("eQTL_catalogue", eqtls), tabix_paths$sample_size[i], tabix_paths$eqtl_n[i])

          mdata$chr <- mdata$chromosome
          mdata$pos <- mdata$position

          # Coloc mode
          qtl_coloc_input <- unlist(strsplit(QTL_COLOC, ","))[1]
          qtl_coloc_mode <- unlist(strsplit(QTL_COLOC, ","))[2]
          gwas_coloc_input <- unlist(strsplit(GWAS_COLOC, ","))[1]
          gwas_coloc_mode <- unlist(strsplit(GWAS_COLOC, ","))[2]

          # LD data
          if (grepl("GTEx_v8", eqtls)) {
            gtex_cov <- read.table(paste0("~/gtex_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/", tabix_paths$qtl_group[i], ".v8.covariates.txt"), header = T, sep = "\t", stringsAsFactors = F, row.names = 1, nrows = 5)
            INDIV <- colnames(gtex_cov)
            GENOFILE <- "~/gtex_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz"
          } else {
            mdata$variant_id <- paste0("chr", mdata$chromosome, "_", mdata$position, "_", mdata$ref, "_", mdata$alt, "_b38")
            mdata$tss_distance <- NA

            pop_1kg_eur <- read.table(here("covid_gwas", "freeze5", "chr3_coloc", "input", "ld0", "indiv_unrelated_1000g_eur.txt"), header = F, sep = "\t", stringsAsFactors = F)
            cat("For LD using EUR population from 1000G", fill = T)
            INDIV <- pop_1kg_eur$V1
            GENOFILE <- paste0("~/lab/data/1kg/phase3_GRCh38/ALL.chr", gwas_chr, ".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz")
          }

          # paste0("fig_locusplot_", eqtls, "_condmask")
          coloc_result[[k]] <- run_coloc(phenotype = trait_id_lookup, data = mdata,
                                         qtl_coloc_input, qtl_coloc_mode,
                                         gwas_coloc_input, gwas_coloc_mode,
                                         genofile = GENOFILE, indiv = INDIV,
                                         locuscompare = TRUE,
                                         locuscompare_fig = here("covid_gwas", "freeze5", "chr3_coloc", "result_condmask", "fig_locuscompare", paste0("locuscompare.", gwas_analysis, ".", trait_id_lookup, "_", res$eqtl_data[k])),
                                         locuszoom = TRUE,
                                         locuszoom_fig = here("covid_gwas", "freeze5", "chr3_coloc", "result_condmask", "fig_locuszoom", paste0("locuszoom.", gwas_analysis, ".", trait_id_lookup, "_", res$eqtl_data[k])),
                                         plot_main = paste0(gwas_analysis, " and ", trait_id_lookup, "\n", tabix_paths$study[i], " - ", tabix_paths$qtl_group[i], " (", tabix_paths$quant_method[i], ")"),
                                         ld_show = TRUE, ld_variant = "qtl",
                                         p1 = 1e-4, p2 = 1e-4, p12 = 5e-6,
                                         pthr = 1e-04,
                                         coloc_mode = COLOC_MODE,
                                         sensitivity = FALSE,
                                         qtl_trait = eqtls,
                                         gwas_trait = gwas_analysis)

          coloc_result[[k]] <- cbind((coloc_result[[k]])[, 1:13], res[k,])

        } else {
          cat("No eQTL variants with P-value < 10-4 within 100kb of the GWAS variant", fill = T)
        }
      } else {
        cat("Less than 100 SNPs", fill = T)
      }
    } else {
      cat("Gene ", trait_id_lookup, " not found", fill = T)
      res$eqtl_gene_id[k] <- gene_df[match(gene_name, gene_df$gene_name), "gene_id"]
      res$eqtl_gene_name[k] <- gene_name
      res$eqtl_molecular_trait_id[k] <- molecular_trait_id
      res$eqtl_data[k] <- paste(tabix_paths[i, c("study", "quant_method", "qtl_group")], collapse = "_")
    }
    k <- k + 1
    cat("--", fill = T)
  }
  cat("-------------------------", fill = T)
}

res <- do.call(rbind, coloc_result)

# Modify column names for sQTL and txrev studies
if (grepl("sqtl", eqtls)) {
  colnames(res) <- gsub("eqtl", "sqtl", colnames(res))
}
if (grepl("tx", eqtls)) {
  colnames(res) <- gsub("eqtl", "tuqtl", colnames(res))
}

write.table(res, file = here("covid_gwas", "freeze5", "chr3_coloc", "result_condmask", paste0("summary_coloc_condmask.", gwas_analysis, ".", eqtls, ".txt")), col.names = T, row.names = F, sep = "\t", quote = F)

cat("Done!")
