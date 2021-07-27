#!/usr/bin/env Rscript
#-------------------------------------------
# Get COVID-19 summary stats for the chr3 locus
#-------------------------------------------

library(here)
library(data.table)

# Read in COVID-19 metadata ----------------
meta <- read.table(here("covid_gwas", "freeze5", "chr3_coloc", "metadata_freeze5.txt"), header = T, sep = "\t", stringsAsFactors = F)

# Get summary stats for the chr3 locus ---------
# Around the lead variant from Freeze3 - Lead variant: 3:45867022:C:G
chr <- 3
start <- 45867022 - 1000000
end <- 45867022 + 1000000

pdf(here("covid_gwas", "freeze5", "chr3_coloc", "input", "fig_covid19_gwas_chr3_locus.pdf"), width = 6, height = 4)
for (i in 1:nrow(meta)) {
  cat(meta$code[i], fill = TRUE)
  data <- fread(paste0("https://storage.googleapis.com/covid19-hg-public/20201215/results/20210107/", meta$file[i]), header = T, sep = "\t", stringsAsFactors = FALSE, data.table = F)
  data <- data[data$`#CHR` %in% chr,]
  data <- data[data$POS > start & data$POS < end, ]
  print(head(data[order(data$all_inv_var_meta_p), ]))
  write.table(data, file = here("covid_gwas", "freeze5", "chr3_coloc", "input", paste0(meta$code[i], "_chr3_locus.txt")), col.names = T, row.names = F, sep = "\t", quote = F)

  # Figure
  data$het_p_flag <- ifelse(data$all_inv_var_het_p < 0.001, 1, 0)
  plot(data$POS, -log10(data$all_inv_var_meta_p), col = ifelse(data$het_p_flag %in% 1, "green", "black"),
       main = meta$code[i], xlab = paste0("Position on chr", chr), ylab = "-log10(all_inv_var_meta_p)")
  legend("topleft", legend = c("het_p < 0.001", "het_p >= 0.001"), col = c("green", "black"), pch = c(1, 1), bty = "n")

  cat("----------------", fill = TRUE)
}

dev.off()

cat("Done!", fill = TRUE)
