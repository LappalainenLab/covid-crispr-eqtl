# -------------------------------------------------------------
# Function: Get (mean imputed) genotype matrix using seqminer
# Author: Silva Kasela
# -------------------------------------------------------------

# library(seqminer)

get_geno_matrix <- function(genofile, range, variant_id, covfile = NULL, indiv = NULL) {
  # Arguments -----
  # genofile - "/path/to/biallelic.vcf.gz"
  # range - chrXX:start-end
  # variant_id - list of variants
  # cov - path to covariates matrix used in QTL mapping to get the list of individuals to filter genotype matrix
  # indiv - list of individuals to include from the genotype matrix
  # ----------------

  geno <- tabix.read.table(tabixFile = genofile, tabixRange = range, stringsAsFactors = FALSE)
  # Try starting range with "chr"
  if (nrow(geno) == 0) {
    cat("Warning: no genotypes found. Adding chr to range.", fill = T)
    range <- paste0("chr", range)
    geno <- tabix.read.table(tabixFile = genofile, tabixRange = range, stringsAsFactors = FALSE)
  }
  # Set ID if missing
  if (all(geno$ID %in% ".")) {
    geno$ID <- paste0("chr", geno$CHROM, "_", geno$POS, "_", geno$REF, "_", geno$ALT, "_b38")
  }
  # Can all variants be found from the genotype file
  if (!all(variant_id %in% geno$ID)) {
    cat("Warning: could not get genotype for all the variants - ", sum(!variant_id %in% geno$ID), "missing", fill = T)
    variant_id <- variant_id[variant_id %in% geno$ID]
  }
  rownames(geno) <- geno$ID
  geno <- geno[match(variant_id, geno$ID), -c(1:9)]

  # Filter geno file based on list of individuals
  if (!is.null(covfile) | !is.null(indiv)) {
    if (!is.null(covfile)) {
      cov <- read.table(covfile, header = T, sep = "\t", stringsAsFactors = F, row.names = 1, nrows = 1)
      indiv <- colnames(cov)
      if (length(indiv) == 1) {
        # try transposing the matrix
        cov <- read.table(covfile, header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
        indiv <- rownames(cov)
      }
      rm(cov)
    }
    stopifnot(indiv %in% colnames(geno))
    geno <- geno[, indiv]
  }

  ind <- colnames(geno)
  geno <- apply(geno, 1, function(x) {
    out <- ifelse(x %in% c("0|0", "0/0"), 0,
                  ifelse(x %in% c("0|1", "1|0", "0/1", "1/0"), 1,
                         ifelse(x %in% c("1|1", "1/1"), 2, NA)))
  })
  rownames(geno) <- ind

  # Impute missing genotypes
  n_missing <- apply(geno, 2, function(x) sum(is.na(x)))
  if (sum(n_missing) > 0) {
    geno <- impute_mean(geno, variant_id = names(n_missing[n_missing > 0]))
  }

  return(geno)
}

impute_mean <- function(geno, variant_id) {
  ## Impute missing genotypes NA to mean
  mean <- apply(geno, 2, function(x) mean(x[!is.na(x)]))
  # Imputed genotypes back to 0, 1, 2
  mean <- ifelse(mean < 0.5, 0, ifelse(mean < 1.5, 1, 2))
  for (j in 1:length(variant_id)) {
    geno[is.na(geno[,variant_id[j]]), variant_id[j]] <- mean[variant_id[j]]
  }
  cat("Imputed at least 1 sample in ", length(variant_id), " sites", fill = T)
  return(geno)
}
