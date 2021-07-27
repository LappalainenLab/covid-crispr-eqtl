# -------------------------------------------------------------------------------------------------------------
# Functions for running coloc analysis w/ possibility to make locuscompare, locuszoom, and sensitivity plots
# Author: Silva Kasela
# -------------------------------------------------------------------------------------------------------------

# library(coloc)
# library(ggplot2)

# source(here("eqtl_coloc", "functions_coloc_signals_adj.R"))
# source(here("eqtl_coloc", "functions_get_genomatrix.R"))
# source(here("eqtl_coloc", "functions_plot_locus.R"))

## Function for running coloc and making locuscompare figure and/or locuszoom

run_coloc <- function(phenotype, data,
                      qtl_coloc_input, qtl_coloc_mode,
                      gwas_coloc_input, gwas_coloc_mode,
                      genofile, covfile = NULL, indiv = NULL,
                      locuscompare = FALSE, locuscompare_fig = "locuscompare",
                      locuszoom = FALSE, locuszoom_fig = "locuszoom",
                      plot_main = "",
                      ld_show = FALSE, ld_variant = "qtl",
                      p1, p2, p12,
                      maxhits = 3, pthr = 1e-05, coloc_mode = "allbutone",
                      sensitivity = FALSE, sensitivity_fig = "sensitivity", sensitivity_file = NULL,
                      sensitivity_p12_vec = c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4),
                      qtl_trait = "qtl", gwas_trait = "gwas")
{
  ## `phenotype` - ID of the molecualr phenotype, e.g., gene ID or probe ID
  ## `qtl_coloc_input` or `gwas_coloc_input` - specifies input data for coloc, "beta" or "pvalues"
  ## `qtl_coloc_mode` or `gwas_coloc_mode` - specifies the mode of coloc, "single, "cond", "mask"
  ## `genofile` - /path/to/geno.vcf.gz or NA, needed for coloc-cond/mask or showing LD on locusplots
  ## `covfile` - /path/to/covariates.txt used in molQTL to get the list of individuals to use for calculating LD
  ## `indiv` - list of individuals to use for calculating LD
  ## `locuscompare` or `locuszoom` - logical, to make figure or not
  ## `locuscompare_fig` or `locuszoom_fig` - path/to/prefix of the file name that is "paste0(locuscompare_fig, ".pp4_", round(result["PP.H4.abf"], 2), ".pdf")"
  ## `plot_main` - title for the locuscompare and locuszoom plots
  ## `ld_show` - logical, to color the points on locuscompare and locuszoom plots according to LD or not
  ## `ld_variant` - "qtl" or "gwas", variant with the lowest p-value in the QTL or GWAS dataset to use for calculating LD, respectively
  ## `p1`, `p2`, `p12` - priors for coloc
  ## `sensitivity` - logical, to run sensitivity analysis or not
  ## `sensitivity_fig` - path/to/prefix of the file name that is "paste0(sensitivity_fig, ".pp4_", round(result["PP.H4.abf"], 2), ".pdf")"
  ## `sensitivity_p12_vec` - vector of different p12 priors to be used in sensitivity analysis
  ## `qtl_trait` or `gwas_trait` - specification of the QTL or GWAS trait in the output file

  select <- data[data$phenotype_id %in% phenotype,]
  # Remove duplicate positions -----
  dupl <- select$id[duplicated(select$id)]
  message(paste0(phenotype, ": remove duplicated positions - ", length(dupl)))
  select <- select[!select$id %in% dupl, ]
  message(paste0(phenotype, " and ", nrow(select), " SNPs in the cis-region"))
  cat(paste0(phenotype, " and ", nrow(select), " SNPs in the cis-region"), fill = T)

  ## Genotype data ----
  ### Might need to exclude monomorphic variants from `select`
  chr <- select$chr[1]
  pos_start <- min(select$pos)
  pos_end <- max(select$pos)

  if (ld_show | qtl_coloc_mode %in% c("cond", "mask") | gwas_coloc_mode %in% c("cond", "mask")) {
    if (!is.na(genofile)) {
      #range <- paste0("chr", chr, ":", pos_start, "-", pos_end)
      range <- paste0(chr, ":", pos_start, "-", pos_end)
      cat("Genotype data from ", genofile, fill = T)
      geno <- get_geno_matrix(genofile = genofile, range = range, variant_id = select$variant_id, covfile = covfile, indiv = indiv)

      if (ncol(geno) != nrow(select)) {
        select <- select[select$variant_id %in% colnames(geno), ]
      }

      # delete variants that have no variability - monomorphic
      unique_genogroups <- apply(geno, 2, function(x){
        length(unique(x))
      })
      if (sum(unique_genogroups == 1) > 0) {
        cat("Excluding monomorphic variants - ", sum(unique_genogroups == 1), fill = T)
        geno <- geno[, unique_genogroups > 1]
        select <- select[!select$variant_id %in% names(unique_genogroups[unique_genogroups == 1]), ]
      }
    } else {
      stop("GENOFILE not specified!")
    }
  }

  # Some annotations -------
  gwas_lead_snv <- ifelse("gwas_rsid" %in% colnames(select), select[which.min(select$gwas_pval), "gwas_rsid"], select[which.min(select$gwas_pval), "id"])
  gwas_lead_snv_pval <- select[which.min(select$gwas_pval), "gwas_pval"]
  gwas_lead_snv_pos <- select[which.min(select$gwas_pval), "pos"]
  if ("tss_distance" %in% colnames(select)) {
    pos_pheno <- select$pos[1] - select$tss_distance[1] - 1
  } else {
    pos_pheno <- NA
  }
  qtl_lead_snv <- ifelse("qtl_rsid" %in% colnames(select), select[which.min(select$qtl_pval), "qtl_rsid"], select[which.min(select$qtl_pval), "variant_id"])
  qtl_lead_snv_pval <- select[which.min(select$qtl_pval), "qtl_pval"]
  qtl_lead_snv_pos <- select[which.min(select$qtl_pval), "pos"]
  gwas_hit_nearest_log10p_5 <- min(abs(select[select$gwas_pval < 1e-5, "pos"] - as.numeric(qtl_lead_snv_pos)))
  if (!is.na(genofile)) {
    r2_lead_snvs <- (cor(geno[, select[c(which.min(select$qtl_pval), which.min(select$gwas_pval)), "variant_id"]])[1, 2])**2
    cat("r2 between the lead QTL and lead GWAS variant: ", r2_lead_snvs, fill = T)
  } else {
    r2_lead_snvs <- NULL
  }

  # Input data for coloc -----
  ## NB! coloc input needs to be V not se! - using se^2 in run_coloc() function
  ## QTL dataset - beta + se
  if (qtl_coloc_input == "beta") {
    if ("qtl_sdY" %in% colnames(data)) {
      dat_2 <- list(snp = select$id, beta = select$qtl_beta, varbeta = select$qtl_se^2, type = "quant", sdY = select$qtl_sdY[1], N = select$qtl_n[1], MAF = select$qtl_maf, method = qtl_coloc_mode)
    } else {
    dat_2 <- list(snp = select$id, beta = select$qtl_beta, varbeta = select$qtl_se^2, type = "quant", N = select$qtl_n[1], MAF = select$qtl_maf, method = qtl_coloc_mode)
    dat_2$sdY <- coloc:::sdY.est(dat_2$varbeta, dat_2$MAF, dat_2$N)
    }
  } else {
    dat_2 <- list(snp = select$id, pvalues = select$qtl_pval, type = "quant", N = select$qtl_n[1], MAF = select$qtl_maf, method = qtl_coloc_mode)
  }
  ## GWAS dataset - use beta + se, if available ------
  if (gwas_coloc_input == "beta") {
    if (select$gwas_type[1] == "quant") {
      dat_1 <- list(snp = select$id, beta = select$gwas_beta, varbeta = select$gwas_se^2, type = "quant", N = select$gwas_n[1], MAF = select$gwas_maf, method = gwas_coloc_mode)
      dat_1$sdY <- coloc:::sdY.est(dat_1$varbeta, dat_1$MAF, dat_1$N)
    } else {
      dat_1 <- list(snp = select$id, beta = select$gwas_beta, varbeta = select$gwas_se^2, type = "cc", N = select$gwas_n[1], s = select$gwas_s[1], method = gwas_coloc_mode)
    }
  } else if (gwas_coloc_input == "pvalues") {
    if (select$gwas_type[1] == "quant") {
      dat_1 <- list(snp = select$id, pvalues = select$gwas_pval, type = "quant", N = select$gwas_n[1], MAF = select$gwas_maf, method = gwas_coloc_mode)
    } else {
      dat_1 <- list(snp = select$id, pvalues = select$gwas_pval, type = "cc", N = select$gwas_n[1], s = select$gwas_s[1], MAF = select$gwas_maf, method = gwas_coloc_mode)
    }
  } else {
    stop("Coloc mode can be either beta or pvalues")
  }

  # Run coloc ------
  if (qtl_coloc_mode == "single" & gwas_coloc_mode == "single") {
    ## Coloc single using coloc.abf()
    cat("Running coloc.abf", fill = TRUE)
    cat("df1:", names(dat_1), "=", dat_1[["method"]], fill = T)
    cat("df2:", names(dat_2), "=", dat_2[["method"]],  fill = T)
    result <- tryCatch(coloc.abf(dataset1 = dat_1,
                                 dataset2 = dat_2,
                                 p1 = p1, p2 = p2, p12 = p12)$summary,
                       error = function(e) {
                         print(e)
                         return(c("nsnps" = phenotype, "PP.H0.abf" = NA, "PP.H1.abf" = NA, "PP.H2.abf" = NA, "PP.H3.abf" = NA, "PP.H4.abf" = NA))
                       })
    result <- data.frame(t(result), stringsAsFactors = F)

    # sensitivity analysis if required and if PP4 > 0.5
    if (sensitivity & result["PP.H4.abf"] > 0.5 & !is.na(result["PP.H4.abf"])) {
      cat("Plotting sensitivity figure", fill = T)
      sens <- lapply(sensitivity_p12_vec, function(x) {
        s <- coloc.abf(dataset1 = dat_1, dataset2 = dat_2, p1 = p1, p2 = p2, p12 = x)$summary
        s["p12"] <- x
        return(s)
      })
      sens <- do.call(rbind, sens)
      sens <- data.frame(sens, stringsAsFactors = F)
      sens$phenotype_id <- phenotype
      # figure
      plot_sens <- tidyr::pivot_longer(data = sens, cols = 2:6,
                                  names_to = "PP.abf", names_pattern = "PP\\.(..)\\.abf",
                                  values_to = "PP", values_ptypes = list(PP = numeric()))
      g <- ggplot(data = plot_sens, aes(x = p12, y = PP, col = PP.abf, group = PP.abf)) +
        labs(title = plot_main,
             subtitle = paste0(phenotype, " locus"),
             x  = expression(p[12]),
             y = "Posterior probability",
             col = "") +
        geom_point() +
        geom_line() +
        geom_hline(yintercept = 0.5, lty = 2, col = "grey75") +
        scale_x_log10() +
        scale_y_continuous(limits = c(0,1)) +
        scale_color_brewer(type = "qual", palette = 6) +
        theme_classic() +
        theme(plot.title = element_text(size = 14),
              plot.subtitle = element_text(size = 13),
              axis.title = element_text(size = 13),
              axis.title.y = element_text(vjust = 2),
              axis.title.x = element_text(vjust = -0.75),
              axis.text = element_text(size = 12),
              strip.text = element_text(size = 12),
              legend.title = element_text(size = 13),
              legend.text = element_text(size = 12),
              legend.position = "top")
      ggsave(filename = paste0(sensitivity_fig, ".", qtl_lead_snv, ".pp4_", round(result["PP.H4.abf"], 2), ".pdf"), plot = g, width = 5.2, height = 4.2)
      rm(plot_sens)

      # save sensitivity file
      if (!is.null(sensitivity_file)) {
        cat("Writing out sensitivity file", fill = T)
        write.table(sens, file = sensitivity_file, col.names = T, row.names = F, sep = "\t", quote = F)
      }
    }
  } else {
    # Get genotype correlations (ie r, not r^2)
    ld_r <- cor(geno)
    stopifnot(colnames(ld_r) == select$variant_id)
    colnames(ld_r) <- rownames(ld_r) <- select$id
    ## Coloc with conditioning or masking
    ## Using modified coloc.signals function to allow to finemap.signals with mode = mask using pvalues
    cat(paste0("Running coloc.signals in ", coloc_mode, " mode"), fill = TRUE)
    cat("df1:", names(dat_1), "=", dat_1[["method"]], fill = T)
    cat("df2:", names(dat_2), "=", dat_2[["method"]], fill = T)

    # Independent signals in teh GWAS dataset
    finemapped_signals <- finemap.signals_adj(D = dat_2, LD = ld_r, method = dat_2[["method"]],
                                              maxhits = maxhits, r2thr = 0.01, pthr = pthr)
    cat("Independent signals in dat_2", fill = T)
    print(finemapped_signals)

    # Run coloc-condmask
    result <- coloc.signals_adj(dataset1 = dat_1,
                                dataset2 = dat_2,
                                LD = ld_r,
                                mode = coloc_mode,
                                p1 = p1, p2 = p2, p12 = p12,
                                maxhits = maxhits, r2thr = 0.01, pthr = pthr)
    # Sensitivity figure
    if (sensitivity) { # sum(result$summary$PP.H4.abf > 0.5) > 0
      n <- nrow(result$summary)
      for (j in 1:min(n, 3)) {
        pdf(paste0(sensitivity_fig, ".loci_", j, ".pp4_", round(result$summary[j, "PP.H4.abf"], 2), ".pdf"), width = 8, height = 5)
        tryCatch(sensitivity(result, "H4 > 0.5", row = j), error = function(e) plot(1:10))
        dev.off()
      }
    }
    if (locuscompare | locuszoom) {
      result_plot <- result$results %>% data.frame()
    }
    result <- result$summary
    result <- data.frame(result, stringsAsFactors = F)

    # Same order as in finemapped_signals - result file from coloc-condmask is sorted by hit2 and hit1
    if (length(names(finemapped_signals)) > 0) {
      stopifnot(names(finemapped_signals) %in% result$hit2)
      result <- result[match(names(finemapped_signals), result$hit2), ]
    }

    # No locusplot if only one indep signals with PP4 < 0.5 OR no loci has PP4 > 0.5
    if ( (sum(result[, "PP.H4.abf"] < 0.5) == 1 & nrow(result) == 1) | (sum(result[, "PP.H4.abf"] > 0.5) == 0) ) {
      locuscompare <- FALSE
      locuszoom <- FALSE
    }
  }

  # Make locuscompare and/or locuszoom plot ------
  if (locuscompare | locuszoom) {
    if (sum(!is.na(result[, "PP.H4.abf"])) > 0) {
      ## make figure if the loci would pass post-processing filters
      # if (gwas_hit_nearest_log10p_5 < 250000) {
        # loop over independent signals in the data
        for (j in 1:min(nrow(result), 3)) {
          if (qtl_coloc_mode %in% c("cond", "mask") | gwas_coloc_mode %in% c("cond", "mask")) {
            plot_data <- select
            plot_data <- plot_data[match(result_plot$snp, plot_data$id), ]
            # gwas is df1 and qtl is df2
            idx <- which(colnames(result_plot) == paste0("z.df1.row", j) | colnames(result_plot) == paste0("z.df1"))
            plot_data$gwas_pval <- pnorm(-abs(result_plot[,idx]))*2
            idx <- which(colnames(result_plot) == paste0("z.df2.row", j) | colnames(result_plot) == paste0("z.df2"))
            plot_data$qtl_pval <- pnorm(-abs(result_plot[,idx]))*2
            file_name_end <- paste0(".loci_condmasked_", j, ".pdf")
          } else {
            # Locusplot figures if using using coloc-standard
            plot_data <- select
            file_name_end <- ".pdf"
          }

          if (ld_show) { # locuscompare plot with colors
            if (ld_variant == "gwas") { # lead GWAS variant for showing LD
              ld_variant_id <- plot_data$variant_id[which.min(plot_data$gwas_pval)]
              ld_variant_group <- "gwas"
            } else {
              ld_variant_id <- plot_data$variant_id[which.min(plot_data$qtl_pval)]
              ld_variant_group <- "qtl"
            }
          } else {
            ld_variant_id = NULL
            ld_pos_start = NULL
            ld_pos_end = NULL
            geno = FALSE
          }

          g <- locusplot_ld(data = plot_data,
                            gwas_label = bquote(-log[10] ~ "(GWAS p-value)"),
                            qtl_label = bquote(-log[10] ~ "(eQTL p-value)"),
                            main = plot_main,
                            pp4 = result[j, "PP.H4.abf"],
                            highlight_index_gwas = which.min(plot_data$gwas_pval),
                            highlight_index_qtl = which.min(plot_data$qtl_pval),
                            highlight_index_gwas.id = plot_data[which.min(plot_data$gwas_pval), "id"],
                            highlight_index_qtl.id = plot_data[which.min(plot_data$qtl_pval), "id"],
                            highlight_index_gwas.label = plot_data[which.min(plot_data$gwas_pval), "id"],
                            highlight_index_qtl.label = plot_data[which.min(plot_data$qtl_pval), "id"],
                            ld_variant_id = ld_variant_id,
                            ld_variant_group = ld_variant_group,
                            geno = geno,
                            locuscompare = locuscompare,
                            locuszoom = locuszoom)

          if ("locuscompare" %in% names(g)) {
            idx <- which(names(g) == "locuscompare")
            ggsave(filename = paste0(locuscompare_fig, ".", qtl_lead_snv, ".pp4_", round(result[j, "PP.H4.abf"], 2), file_name_end), plot = g[[idx]], width = 4.5, height = 4.35)
          }
          if ("locuszoom" %in% names(g)) {
            idx <- which(names(g) == "locuszoom")
            ggsave(filename = paste0(locuszoom_fig, ".", qtl_lead_snv, ".pp4_", round(result[j, "PP.H4.abf"], 2), file_name_end), plot = g[[idx]], width = 5, height = 5)
          }
        }
      # }
    }
  }

  out <- cbind(result, "qtl" = qtl_trait, "phenotype_id" = phenotype,
                 "qtl_lead_snv" = qtl_lead_snv, "qtl_lead_snv_pval" = qtl_lead_snv_pval, "qtl_lead_snv_pos" = qtl_lead_snv_pos,
                 "gwas" = gwas_trait, "gwas_lead_snv" = gwas_lead_snv, "gwas_lead_snv_pval" = gwas_lead_snv_pval, "gwas_lead_snv_pos" = gwas_lead_snv_pos, "gwas_hit_nearest_log10p_5" = gwas_hit_nearest_log10p_5,
                 "chr" = paste0("chr", chr), pos_phenotype = pos_pheno, region_start = pos_start, region_end = pos_end)

  return(out)
}
