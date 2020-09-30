#!/usr/bin/env Rscript
#-------------------------------------------------------
# Colocalization analysis for COVID-19 GWAS
# and eQTLs from the eQTL Catalogue or GTEx v8
#     * Using coloc-standard or coloc-cond/mask
#-------------------------------------------------------

# Idea: run coloc between COVID-19 GWAS and eQTLs, if there's a eVariant within 100kb of the GWAS variant with nominal p-value < 10-4
# Run coloc-standard in a 1Mb region centered on the lead GWAS variant (+/-500 kb from the lead GWAS variant)
# Run coloc-condmask in a 2Mb region centered on the lead GWAS variant (+/-1 Mb from the lead GWAS variant)

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
source(here("eqtl_coloc", "functions_plot_locus.R"))
source(here("eqtl_coloc", "functions_get_genomatrix.R"))
source(here("eqtl_coloc", "functions_coloc_signals_adj.R"))

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

args <- commandArgs(trailingOnly = TRUE)

# Which eQTL sorice to use: "GTEx_v8" or "eQTL_catalogue"
eqtls <- args[1] # 
print(eqtls)

# Which coloc to use: "standard" or "condmask"
coloc <- args[2]
print(coloc)

# coloc-condmask
if (coloc %in% "condmask") {
  if (eqtls %in% "GTEx_v8") {
    QTL_COLOC <- 'beta,cond'
  } else {
    QTL_COLOC <- 'pvalues,mask'
  }
  GWAS_COLOC <- "beta,single"
  COLOC_MODE <- "iterative"
}

# Running coloc in a 1Mb- or 2Mb-region centered at the lead GWAS variant (+/- 500kb or +/- 1Mb from the SNP)
window <- ifelse(coloc %in% "condmask", 1e+06, 5e+05)
pval_th <- 1e-4
distance_th <- 1e+05

# COVID-19 GWAS data: ANA_C2_V2 (round 3) ----------------------
# Lead variant: 3:45867022:C:G
gwas_analysis <- "ANA_C2_V2"
gwas_file <- list.files(path = here("data"), pattern = gwas_analysis)
gwas_file <- gwas_file[!grepl(".tbi", gwas_file)]
gwas_chr <- 3
gwas_pos <- 45867022
range <- paste0(gwas_chr, ":",  gwas_pos - window, "-", gwas_pos + window)
gwas_data <- tabix.read.table(tabixFile = here("data", gwas_file), tabixRange = range, stringsAsFactors = FALSE)

# sample size and type from COVID-19 HGI webpage
gwas_data$gwas_type <- "cc"
gwas_data$gwas_n <- 6696 + 1073072
gwas_data$gwas_s <- 6696 / gwas_data$gwas_n

# lead gwas variant
gwas_data[which.min(gwas_data$all_inv_var_meta_p), "SNP"]
lead_gwas <- "3:45867022:C:G"
lead_gwas_id <- "chr3_45867022_C_G"

# Genes of interest in that region
if (eqtls %in% "GTEx_v8") {
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
if (eqtls %in% "GTEx_v8") {
  files <- list.files("~/gtex_v8/eqtl/GTEx_Analysis_v8_eQTL_all_associations_indexed/results/")
  files <- files[!grepl("tbi", files)]
  tabix_paths <- data.frame("study" = "GTEx_v8",
                            "qtl_group" = as.character(sapply(files, function(x) unlist(strsplit(x, "[.]"))[1])),
                            "ftp_path" = paste0("~/gtex_v8/eqtl/GTEx_Analysis_v8_eQTL_all_associations_indexed/results/", files),
                            "quant_method" = "ge",
                            stringsAsFactors = F)
  # Sample size in GTEx
  sample_annot <- read.table("~/lab/data/gtex/v8/sample_annotations/tissue_sample_count.txt", header = F, sep = "\t", stringsAsFactors = F)
  stopifnot(tabix_paths$qtl_group %in% sample_annot$V1)
  tabix_paths$eqtl_n <- sample_annot[match(tabix_paths$qtl_group, sample_annot$V1), "V2"]
  rm(sample_annot)
} else {
  ## Tutorial and functios to work with eQTL Catalogue data from here: http://htmlpreview.github.io/?https://github.com/kauralasoo/eQTL-Catalogue-resources/blob/master/scripts/tabix_use_case.html
  tabix_paths <- read.table("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  tabix_paths <- tabix_paths %>%
    filter(quant_method %in% c("ge", "microarray"))

  # Extract column names from first file
  ##column_names <- fread(tabix_paths$ftp_path[1], header = F, sep = "\t", stringsAsFactors = F, nrows = 1, data.table = F)
  ##column_names <- as.character(column_names[1,])
  column_names <- c("variant", "type", "rsid", "ref", "r2", "pvalue", "position", "molecular_trait_object_id", "molecular_trait_id", "median_tpm", "maf", "gene_id", "chromosome", "beta", "an", "alt", "ac")
}

res <- data.frame("gwas" = rep(gwas_analysis, 1000), #don't know exactly how many rows there might be in total
                  "lead_gwas_variant" = NA,
                  "lead_gwas_pval" = NA,
                  "eqtl_gene_id" = NA,
                  "eqtl_gene_name" = NA,
                  "eqtl_molecular_trait_id" = NA,
                  "eqtl_data" = NA,
                  "eqtl_n" = NA,
                  "lead_eqtl_variant_larger_window" = NA,
                  "lead_eqtl_pval_larger_window" = NA,
                  "lead_eqtl_variant" = NA,
                  "lead_eqtl_pval" = NA,
                  "lead_gwas_variant_eqtl_pval" = NA,
                  "coloc_nsnps" = NA,
                  "coloc_pp0" = NA,
                  "coloc_pp1" = NA,
                  "coloc_pp2" = NA,
                  "coloc_pp3" = NA,
                  "coloc_pp4" = NA,
                  stringsAsFactors = F)

if (coloc %in% "condmask") {
  coloc_result <- list()
}

k <- 1
for (i in 1:nrow(tabix_paths)) {
  ftp_path =  tabix_paths$ftp_path[i]
  cat(ftp_path =  tabix_paths$ftp_path[i], fill = T)

  # Get eQTL data
  eqtl_range <- ifelse(eqtls %in% "GTEx_v8",
                       paste0("chr", gwas_chr, ":",  gwas_pos - 1.5e6, "-", gwas_pos + 1.5e6),
                       paste0(gwas_chr, ":",  gwas_pos - 1.5e6, "-", gwas_pos + 1.5e6))
  eqtl_data <- tabix.read.table(tabixFile = ftp_path, tabixRange = eqtl_range, stringsAsFactors = FALSE)

  if (eqtls %in% "eQTL_catalogue") {
    if (ncol(eqtl_data) == 19) {
      # Some file have 19 columns?
      cat("Warning: 19 columns instead of 17", fill = T)
      colnames(eqtl_data) <- c("molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", "ma_samples", "an", "ac", "maf", "pvalue", "beta", "se", "molecular_trait_object_id", "gene_id", "median_tpm", "r2", "type", "rsid")
    } else {
      stopifnot(ncol(ftp_path) == length(column_names))
      colnames(eqtl_data) <- column_names
    }
  } else {
    # Make colnames similar to eQTL Catalogue files
    idx <- sapply(c("pval_nominal", "slope", "slope_se", "chr", "pos"), function(x) which(colnames(eqtl_data) %in% x))
    colnames(eqtl_data)[idx] <- c("pvalue", "beta", "se", "chromosome", "position")
    eqtl_data$chromosome <- sub("chr", "", eqtl_data$chromosome)

    # Add ref and alt allele
    eqtl_data$ref <- sapply(eqtl_data$variant_id, function(x) unlist(strsplit(x, "_"))[3])
    eqtl_data$alt <- sapply(eqtl_data$variant_id, function(x) unlist(strsplit(x, "_"))[4])

    # Remove "_b38" at the end of variant
    eqtl_data$variant <- sapply(eqtl_data$variant_id, function(x) sub("_b38", "", x))
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

  # Figure of cis-region of the genes
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
  ggsave(here("covid_gwas", "chr3_coloc", "fig_cis_region", paste0(paste(tabix_paths[i, c("study", "quant_method", "qtl_group")], collapse = "_"), ".png")), g,
           width = width, height = 12)

  trait_id <- levels(eqtl_data$trait_id_f)
  # Run gene by gene (gene-molecular trait if microarray)
  # Run coloc if min p-value in the coloc-region is < 10-4
  for (j in 1:length(trait_id)) {
    trait_id_lookup <- trait_id[j]
    gene_name <- unlist(strsplit(trait_id_lookup, "-"))[1]
    molecular_trait_id <- unlist(strsplit(trait_id_lookup, "-"))[2]

    if (trait_id_lookup %in% eqtl_data$trait_id) {
      cat("Gene ", trait_id_lookup, fill = T)

      if (eqtls %in% "eQTL_catalogue") {
        # Select one gene, remove rsid duplicates and multi-allelic variants
        sel <- eqtl_data %>%
          filter(trait_id %in% trait_id_lookup) %>%
          select(-rsid) %>%
          distinct() %>% #rsid duplicates
          mutate(id = paste(chromosome, position, sep = ":")) %>%
          group_by(id) %>%
          mutate(row_count = n()) %>%
          ungroup() %>%
          filter(row_count == 1) %>% # multi-allelic
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

        res$lead_gwas_variant[k] <- mdata[which.min(mdata$all_inv_var_meta_p), "SNP"]
        res$lead_gwas_pval[k] <- mdata[which.min(mdata$all_inv_var_meta_p), "all_inv_var_meta_p"]
        res$eqtl_gene_id[k] <- gene_df[match(gene_name, gene_df$gene_name), "gene_id"]
        res$eqtl_gene_name[k] <- gene_name
        res$eqtl_molecular_trait_id[k] <- molecular_trait_id
        res$eqtl_data[k] <- paste(tabix_paths[i, c("study", "quant_method", "qtl_group")], collapse = "_")
        res$eqtl_n[k] <- ifelse(eqtls %in% "eQTL_catalogue", mdata$an[1]/2, tabix_paths$eqtl_n[i])
        res$lead_eqtl_variant_larger_window[k] <- sel[which.min(sel$pvalue), "SNP"]
        res$lead_eqtl_pval_larger_window[k] <- sel[which.min(sel$pvalue), "pvalue"]
        res$lead_eqtl_variant[k] <- mdata[which.min(mdata$pvalue), "SNP"]
        res$lead_eqtl_pval[k] <- mdata[which.min(mdata$pvalue), "pvalue"]
        res$lead_gwas_variant_eqtl_pval[k] <- mdata[which.min(mdata$all_inv_var_meta_p), "pvalue"]

        # Check if there's eQTL variants with P-value < 10-4 within 100kb of the GWAS variant
        min_p_check <- sel %>%
            filter(position > gwas_pos - distance_th & position < gwas_pos + distance_th,
                   pvalue < pval_th)

        # eQTL variants with P-value < 10-4 within 100kb of the GWAS variant
        if (nrow(min_p_check) > 0) {
          if (coloc %in% "standard") {
            # coloc-standard
            res_coloc <- coloc.abf(dataset1 = list(beta = mdata$all_inv_var_meta_beta,
                                                   varbeta = mdata$all_inv_var_meta_sebeta**2,
                                                   type = mdata$gwas_type[1],
                                                   s = mdata$gwas_s[1],
                                                   N = mdata$gwas_n[1]),
                                   dataset2 = list(pvalues = mdata$pvalue,
                                                   type = "quant",
                                                   MAF = mdata$maf,
                                                   N = ifelse(eqtls %in% "eQTL_catalogue", mdata$an[1]/2, tabix_paths$eqtl_n[i])),
                                   p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)

            res$coloc_nsnps[k] <- res_coloc$summary["nsnps"]
            res$coloc_pp0[k] <- res_coloc$summary["PP.H0.abf"]
            res$coloc_pp1[k] <- res_coloc$summary["PP.H1.abf"]
            res$coloc_pp2[k] <- res_coloc$summary["PP.H2.abf"]
            res$coloc_pp3[k] <- res_coloc$summary["PP.H3.abf"]
            res$coloc_pp4[k] <- res_coloc$summary["PP.H4.abf"]

            if (res_coloc$summary["PP.H4.abf"] > 0) {
              # Locuscompare and locuszoom figures
              lc <- ggplot(data = mdata, aes(x = -log10(all_inv_var_meta_p), y = -log10(pvalue))) +
                labs(title = paste0(gwas_analysis, " - ", trait_id_lookup),
                     subtitle = paste0(tabix_paths$study[i], " - ", tabix_paths$qtl_group[i], " (", tabix_paths$quant_method[i], ")"),
                     x = bquote(-log[10] ~ "(GWAS P-value)"),
                     y = bquote(-log[10] ~ "(eQTL P-value)")) +
                geom_point(pch = 21)

              x_points <- seq(round(gwas_pos - window, -5), gwas_pos + window, 2e5)
              x_points <- x_points[x_points > gwas_pos - window & x_points < gwas_pos + window]
              top_variant <- data.frame(POS = gwas_pos,
                                        pval = as.numeric(mdata[mdata$SNP %in% lead_gwas, c("pvalue", "all_inv_var_meta_p")]),
                                        group = c("eQTL", "GWAS"),
                                        stringsAsFactors = F)

              lz <- mdata %>%
                mutate(GWAS = all_inv_var_meta_p,
                       eQTL = pvalue) %>%
                tidyr::pivot_longer(cols = c("GWAS", "eQTL"), names_to = "group", values_to = "pval") %>%
                ggplot(aes(x = POS, y = -log10(pval))) +
                labs(#title = tabix_paths$study[i],
                    #subtitle = paste0(tabix_paths$qtl_group[i], " (", tabix_paths$quant_method[i], ")"),
                    subtitle = paste0("PP3 = ", signif(res_coloc$summary["PP.H3.abf"], 2),
                                      ", PP4 = ", signif(res_coloc$summary["PP.H4.abf"], 2)),
                    x = paste0("Position on chr", gwas_chr, " (Mb)"),
                    y = bquote(-log[10] ~ "(P-value)")) +
                geom_point(pch = 21) +
                geom_point(data = top_variant, col = "purple", size = 1.5) +
                facet_grid(group ~ ., scales = "free_y") +
                scale_x_continuous(breaks = x_points, labels = x_points/1e6, expand = expansion(mult = c(0.01, 0.01), add = c(0, 0))) +
                theme(plot.subtitle = element_text(size = 11))

              fig <- lc / lz +
                plot_layout(heights = c(1, 2)) +
                plot_annotation(tag_levels = 'A')
              ggsave(here("eqtl_coloc", paste0("fig_locusplot_", eqtls), paste0("locusplots.", gwas_analysis, ".", trait_id_lookup, "_", res$eqtl_data[k], ".pp4_", round(res_coloc$summary["PP.H4.abf"], 2), ".pdf")),
                     plot = fig, width = 4, height = 9)
            } else {
              cat("PP4 < 0.25", fill = T)
            }
          } else {
            # coloc-condmask
            mdata$phenotype_id <- mdata$trait_id
            mdata$id <- sapply(mdata$variant, function(x) paste(unlist(strsplit(x, "_"))[1:2], collapse = "_"))

            mdata$gwas_pval <- mdata$all_inv_var_meta_p
            mdata$gwas_beta <- mdata$all_inv_var_meta_beta
            mdata$gwas_se <- mdata$all_inv_var_meta_sebeta

            mdata$qtl_pval <- mdata$pvalue
            mdata$qtl_maf <- mdata$maf
            mdata$qtl_n <- ifelse(eqtls %in% "eQTL_catalogue", mdata$an[1]/2, tabix_paths$eqtl_n[i])

            mdata$chr <- mdata$chromosome
            mdata$pos <- mdata$position

            # Coloc mode
            qtl_coloc_input <- unlist(strsplit(QTL_COLOC, ","))[1]
            qtl_coloc_mode <- unlist(strsplit(QTL_COLOC, ","))[2]
            gwas_coloc_input <- unlist(strsplit(GWAS_COLOC, ","))[1]
            gwas_coloc_mode <- unlist(strsplit(GWAS_COLOC, ","))[2]

            # LD data
            if (eqtls %in% "GTEx_v8") {
              mdata$qtl_beta <- mdata$beta
              mdata$qtl_se <- mdata$se
              mdata$qtl_sdY <- coloc:::sdY.est(mdata$qtl_se**2, mdata$qtl_maf, mdata$qtl_n)

              gtex_cov <- read.table(paste0("~/gtex_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/", tabix_paths$qtl_group[i], ".v8.covariates.txt"), header = T, sep = "\t", stringsAsFactors = F, row.names = 1, nrows = 5)
              INDIV <- colnames(gtex_cov)
              GENOFILE <- "~/gtex_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz"
            } else {
              mdata$variant_id <- paste0("chr", mdata$chromosome, "_", mdata$position, "_", mdata$ref, "_", mdata$alt, "_b38")

              pop_1kg <- read.table("~/data/1kg/gazal_et_al_2019.table_S4.filtered_unrelated_outbred.txt", header = T, sep = "\t", stringsAsFactors = F)
              cat("For LD using CEU population from 1000G", fill = T)
              INDIV <- pop_1kg[pop_1kg$POP %in% "CEU", "IID"]
              GENOFILE <- paste0("~/data/1kg/phase3_GRCh38/ALL.chr", gwas_chr, ".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz")
            }

            coloc_result[[k]] <- run_coloc(phenotype = trait_id_lookup, data = mdata,
                                            qtl_coloc_input, qtl_coloc_mode,
                                            gwas_coloc_input, gwas_coloc_mode,
                                            genofile = GENOFILE, indiv = INDIV,
                                            locuscompare = TRUE,
                                            locuscompare_fig = here("eqtl_coloc", paste0("fig_locusplot_", eqtls, "_condmask"), paste0("locuscompare.", gwas_analysis, ".", trait_id_lookup, "_", res$eqtl_data[k])),
                                            locuszoom = TRUE,
                                            locuszoom_fig = here("eqtl_coloc", paste0("fig_locusplot_", eqtls, "_condmask"), paste0("locuszoom.", gwas_analysis, ".", trait_id_lookup, "_", res$eqtl_data[k])),
                                            plot_main = paste0(gwas_analysis, " and ", trait_id_lookup, "\n", tabix_paths$study[i], " - ", tabix_paths$qtl_group[i], " (", tabix_paths$quant_method[i], ")"),
                                            ld_show = TRUE, ld_variant = "qtl",
                                            p1 = 1e-4, p2 = 1e-4, p12 = 5e-6,
                                            pthr = 1e-04,
                                            coloc_mode = COLOC_MODE,
                                            sensitivity = FALSE,
                                            qtl_trait = "eQTL",
                                            gwas_trait = "GWAS")

            coloc_result[[k]] <- cbind(coloc_result[[k]], res[k,])
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

if (coloc %in% "condmask") {
  res <- do.call(rbind, coloc_result)
} else {
  res <- res %>%
    filter(!is.na(eqtl_data))
}

write.table(res, file = here("eqtl_coloc", paste0("summary_coloc_", eqtls, "_", coloc, ".txt")), col.names = T, row.names = F, sep = "\t", quote = F)


