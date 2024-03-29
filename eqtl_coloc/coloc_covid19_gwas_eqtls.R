#!/usr/bin/env Rscript
#-------------------------------------------------------
# Colocalization analysis for COVID-19 GWAS (freeze 5)
# and molQTLs from the eQTL Catalogue or GTEx v8
#-------------------------------------------------------

# Idea: run coloc between COVID-19 GWAS and molQTLs, if there's an eVariant within 100kb of the GWAS variant with nominal p-value < 10-4
# Run coloc in a 1Mb region centered on the lead GWAS variant, if more than 100 variants in the coloc-region

library(here)
library(data.table)
library(seqminer)
library(dplyr)
library(ggplot2)
library(patchwork)
library(coloc)

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

# Running coloc in a 1Mb region centered at the lead GWAS variant (+/- 500kb from the SNP)
window <- 500000
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
                  "coloc_nsnps" = NA,
                  "coloc_pp0" = NA,
                  "coloc_pp1" = NA,
                  "coloc_pp2" = NA,
                  "coloc_pp3" = NA,
                  "coloc_pp4" = NA,
                  stringsAsFactors = F)

k <- 1
for (i in 1:nrow(tabix_paths)) {
  ftp_path =  tabix_paths$ftp_path[i]
  cat(ftp_path =  tabix_paths$ftp_path[i], fill = T)

  if (grepl("eQTL_catalogue", eqtls)) {
    eqtl_data <- fread(ftp_path, header = T, sep = "\t", stringsAsFactors = F, data.table = F)

    #stopifnot(ncol(eqtl_data) == length(column_names))
    #colnames(eqtl_data) <- column_names

    ## For txrevise events, choosing one of the molecular_trait_id to represent the gene
    #if (eqtls %in% c("eQTL_catalogue_tx", "eQTL_catalogue_txrev")) {
    #  f <- sub("all.tsv", "permuted.tsv", ftp_path)
    #  if (grepl("GEUVADIS", tabix_paths$ftp_path[i])) {
    #    f <- sub("_tx_|_txrev_", ".", f)
    #  }
    #  permuted <- data.table::fread(f, header = T, sep = "\t", stringsAsFactors = F, data.table = F)
    #  permuted$gene_id <- sapply(permuted$molecular_trait_object_id, function(x) unlist(strsplit(x, "[.]"))[1])
    #  permuted <- permuted[permuted$gene_id %in% gene_df$gene_id,]
    #  stopifnot(permuted$molecular_trait_id %in% eqtl_data$molecular_trait_id)
    #  eqtl_data <- eqtl_data[eqtl_data$molecular_trait_id %in% permuted$molecular_trait_id, ]
    #  rm(permuted)
    #}

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

  if (tabix_paths$quant_method[i] %in% c("ge", "splicing")) {
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

  # Figure
  x_points <- seq(round(gwas_pos - 1e6, -6), gwas_pos + 1e6, 1e6/2) # define x-axis breaks
  top_variant <- eqtl_data %>%
    filter(variant %in% lead_gwas_id)

  # individual y_max => need to create blank data
  min_p <- tapply(eqtl_data$pvalue, eqtl_data$trait_id, min)
  min_p <- min_p[levels(eqtl_data$trait_id_f)]
  blank_data <- data.frame(position = gwas_pos,
                           pvalue = c(ifelse(min_p > 10**-6 | is.na(min_p) , 10**-6, min_p), rep(1, length(min_p))),
                           trait_id_f = factor(rep(levels(eqtl_data$trait_id_f), 2), levels = levels(eqtl_data$trait_id_f)),
                           stringsAsFactors = F)

  # ncol and width
  ncol <- ifelse(length(levels(eqtl_data$trait_id_f)) %% nrow(gene_df) == 0,
                 length(levels(eqtl_data$trait_id_f)) / nrow(gene_df), trunc(length(levels(eqtl_data$trait_id_f)) / nrow(gene_df)) + 1)
  width <- ifelse(ncol > 1, 3*ncol*0.9, 3)

  g <- ggplot(eqtl_data, aes(x = position, y = -log10(pvalue))) +
    labs(title = tabix_paths$study[i],
         subtitle = paste0(tabix_paths$qtl_group[i], " (", tabix_paths$quant_method[i], ")"),
         x = paste0("Position on chr", gwas_chr, " (Mb)"),
         y = bquote(-log[10] ~ "(eQTL P-value)")) +
    geom_point(col = "grey20") +
    geom_hline(yintercept = -log10(1e-4), col = "indianred", lty = 2) +
    geom_point(data = top_variant, col = "mediumorchid1") +
    geom_blank(data = blank_data) +
    scale_x_continuous(breaks = x_points, labels = x_points/1e6, expand = expansion(mult = c(0.01, 0.01), add = c(0, 0))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2), add = c(0, 0))) +
    facet_wrap(~trait_id_f, drop = FALSE, nrow = nrow(gene_df), scales = "free_y")
  ggsave(here("covid_gwas", "freeze5", "chr3_coloc", "fig_cis_region", paste0(paste(tabix_paths[i, c("study", "quant_method", "qtl_group")], collapse = "_"), ".png")), g, width = width, height = 12)

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
          res_coloc <- coloc.abf(dataset1 = list(beta = mdata$all_inv_var_meta_beta,
                                                 varbeta = mdata$all_inv_var_meta_sebeta**2,
                                                 type = mdata$gwas_type[1],
                                                 s = mdata$gwas_s[1],
                                                 N = mdata$gwas_n[1]),
                                 dataset2 = list(beta = mdata$beta,
                                                 varbeta = mdata$se**2,
                                                 type = "quant",
                                                 MAF = mdata$maf,
                                                 N = ifelse(grepl("eQTL_catalogue", eqtls), tabix_paths$sample_size[i], tabix_paths$eqtl_n[i])),
                                 p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)

          res$coloc_nsnps[k] <- res_coloc$summary["nsnps"]
          res$coloc_pp0[k] <- res_coloc$summary["PP.H0.abf"]
          res$coloc_pp1[k] <- res_coloc$summary["PP.H1.abf"]
          res$coloc_pp2[k] <- res_coloc$summary["PP.H2.abf"]
          res$coloc_pp3[k] <- res_coloc$summary["PP.H3.abf"]
          res$coloc_pp4[k] <- res_coloc$summary["PP.H4.abf"]

          if (res_coloc$summary["PP.H4.abf"] > 0.25) {
            qtl_type <- ifelse(grepl("sqtl", eqtls), "sQTL",
                               ifelse(grepl("tx", eqtls), "tuQTL", "eQTL"))
            # Locuscompare and locuszoom figures
            lc <- ggplot(data = mdata, aes(x = -log10(all_inv_var_meta_p), y = -log10(pvalue))) +
              labs(#title = paste0(gwas_analysis, " - ", trait_id_lookup),
                   #subtitle = paste0(tabix_paths$study[i], " - ", tabix_paths$qtl_group[i], " (", tabix_paths$quant_method[i], ")"),
                   x = bquote(-log[10] ~ "(GWAS P-value)"),
                   y = bquote(-log[10] ~ "(" ~ .(qtl_type) ~ " P-value)")) +
              geom_point(pch = 21)

            x_points <- seq(round(gwas_pos - window, -5), gwas_pos + window, 2e5)
            x_points <- x_points[x_points > gwas_pos - window & x_points < gwas_pos + window]
            top_variant <- data.frame(POS = gwas_pos,
                                      pval = as.numeric(mdata[mdata$SNP %in% lead_gwas, c("pvalue", "all_inv_var_meta_p")]),
                                      group = c(qtl_type, "GWAS"),
                                      stringsAsFactors = F)

            if (qtl_type %in% "eQTL") {
              mdata_pivot <- mdata %>%
                mutate(GWAS = all_inv_var_meta_p,
                       eQTL = pvalue) %>%
                tidyr::pivot_longer(cols = c("GWAS", "eQTL"), names_to = "group", values_to = "pval")
            } else if (qtl_type %in% "sQTL") {
              mdata_pivot <- mdata %>%
                mutate(GWAS = all_inv_var_meta_p,
                       sQTL = pvalue) %>%
                tidyr::pivot_longer(cols = c("GWAS", "sQTL"), names_to = "group", values_to = "pval")
            } else {
              mdata_pivot <- mdata %>%
                mutate(GWAS = all_inv_var_meta_p,
                       tuQTL = pvalue) %>%
                tidyr::pivot_longer(cols = c("GWAS", "tuQTL"), names_to = "group", values_to = "pval")
            }

            lz <- ggplot(data = mdata_pivot, aes(x = POS, y = -log10(pval))) +
              labs(#title = tabix_paths$study[i],
                  #subtitle = paste0(tabix_paths$qtl_group[i], " (", tabix_paths$quant_method[i], ")"),
                  #subtitle = paste0("PP3 = ", signif(res_coloc$summary["PP.H3.abf"], 2),
                  #                  ", PP4 = ", signif(res_coloc$summary["PP.H4.abf"], 2)),
                  x = paste0("Position on chr", gwas_chr, " (Mb)"),
                  y = bquote(-log[10] ~ "(P-value)")) +
              geom_point(pch = 21) +
              geom_point(data = top_variant, col = "purple", size = 1.5) +
              facet_grid(group ~ ., scales = "free_y") +
              scale_x_continuous(breaks = x_points, labels = x_points/1e6, expand = expansion(mult = c(0.01, 0.01), add = c(0, 0))) +
              theme(plot.subtitle = element_text(size = 11))

            fig <- lc + lz +
              plot_annotation(tag_levels = 'A',
                              title = paste0(gwas_analysis, " - ", qtl_type, "s for ", trait_id_lookup," from ", tabix_paths$study[i], " - ", tabix_paths$qtl_group[i], " (", tabix_paths$quant_method[i], ")"),
                              subtitle = paste0("PP3 = ", signif(res_coloc$summary["PP.H3.abf"], 2), ", PP4 = ", signif(res_coloc$summary["PP.H4.abf"], 2)))
            ggsave(here("covid_gwas", "freeze5", "chr3_coloc", paste0("fig_locusplot_", eqtls), paste0("locusplots.", gwas_analysis, ".", trait_id_lookup, "_", res$eqtl_data[k], ".pp4_", round(res_coloc$summary["PP.H4.abf"], 2), ".pdf")), plot = fig, width = 8, height = 4.5)
          } else {
            cat("PP4 < 0.25", fill = T)
          }
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

res <- res %>%
  filter(!is.na(eqtl_data))

# Modify column names for sQTL and txrev studies
if (grepl("sqtl", eqtls)) {
  colnames(res) <- gsub("eqtl", "sqtl", colnames(res))
}
if (grepl("tx", eqtls)) {
  colnames(res) <- gsub("eqtl", "tuqtl", colnames(res))
}

write.table(res, file = here("covid_gwas", "freeze5", "chr3_coloc", "result", paste0("summary_coloc.", gwas_analysis, ".", eqtls, ".txt")), col.names = T, row.names = F, sep = "\t", quote = F)

cat("Done!")
