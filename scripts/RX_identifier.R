#!/usr/bin/env Rscript

# RX_identifier for Canis lupus familiaris (Dog)
# Usage: Rscript RX_identifier.R <sample_prefix>
# Based on: "Accurate sex identification of ancient elephant and other animal remains using low-coverage DNA shotgun sequencing data" by de Flamingh et 2020 in G3: Genes, Genomes, Genetics. https://doi.org/10.1534/g3.119.400833

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Error: No sample prefix provided. Usage: Rscript RX_identifier.R <sample_prefix>")
}
PREFIX <- args[1]

# 1. Load data (Dog: 37 autosomes + X chromosome = 38 rows)
idx_file <- paste0(PREFIX, '.idxstats')
if (!file.exists(idx_file)) {
  stop(paste("Error: File", idx_file, "not found."))
}

# Samtools idxstats: [Chr_name] [Chr_len] [Mapped_reads] [Unmapped_reads]
idxstats <- read.table(idx_file, header=FALSE, nrows=38)
colnames(idxstats) <- c("chr", "len", "mapped", "unmapped")

# 2. Total counts
total_ref_len <- sum(as.numeric(idxstats$len))
total_mapped_reads <- sum(as.numeric(idxstats$mapped))

# 3. Calculate Normalized Ratios (Rt) for all 38 chromosomes
# Rt = (Mapped_reads / Total_Mapped) / (Chr_len / Total_Ref_Len)
idxstats$Rt <- (idxstats$mapped / total_mapped_reads) / (idxstats$len / total_ref_len)

# 4. Calculate Rx Ratio (X coverage relative to each autosome)
# X is row 38, Autosomes are rows 1-37
Rx_X <- idxstats$Rt[38]
autosomes_Rt <- idxstats$Rt[1:37]
tot_ratios <- Rx_X / autosomes_Rt

# 5. Statistics
Rx <- mean(tot_ratios)
confinterval <- 1.96 * (sd(tot_ratios) / sqrt(37))
CI1 <- Rx - confinterval
CI2 <- Rx + confinterval

# 6. Output Results
cat("--------------------------------------------\n")
cat("Sample ID:", PREFIX, "\n")
cat("Rx Ratio: ", round(Rx, 4), "\n")
cat("95% CI:   [", round(CI1, 4), ",", round(CI2, 4), "]\n")
cat("--------------------------------------------\n")

if (CI1 > 0.8) {
  result <- "Female (XX)"
} else if (CI2 < 0.6) {
  result <- "Male (XY)"
} else if (Rx > 0.6 && Rx < 0.8) {
  result <- "Consistent with XY, but low confidence"
} else {
  result <- "Inconclusive / Could not be assigned"
}

cat("Result:   ", result, "\n")
