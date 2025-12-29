#!/usr/bin/env Rscript

# Usage: Rscript scripts/visualize_length.R <input_file> <output_pdf> <sample_id>

# 1. Load Dependencies
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(data.table)
})

# 2. Parse Command Line Arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript visualize_length.R <input_file> <output_pdf> <sample_id>")
}

data_file   <- args[1]
output_file <- args[2]
sample_name <- args[3]

# 3. Data Ingestion
# Using fread for high-performance reading
if (!file.exists(data_file)) stop("Input file not found!")
data_raw <- fread(data_file, header = FALSE, col.names = "length")

# 4. Data Sampling for Performance
# For ancient DNA length plots, 1 million points are sufficient for a smooth distribution
set.seed(123)
max_obs <- 1e6 
if (nrow(data_raw) > max_obs) {
  message(paste("Sampling 1,000,000 observations for visualization from", nrow(data_raw)))
  data_raw <- data_raw[sample(1:.N, size = max_obs)]
}

# 5. Visualization
# Standard aDNA length distribution typically ranges from 30bp to 150bp
p <- ggplot(data_raw, aes(x = length)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_density(aes(y = ..count..), color = "darkred", size = 0.8) +
  labs(
    title = paste("Fragment Length Distribution:", sample_name),
    subtitle = paste("Total Reads Analyzed:", format(nrow(data_raw), big.mark=",")),
    x = "Read Length (bp)",
    y = "Frequency"
  ) +
  scale_x_continuous(limits = c(20, 160), breaks = seq(20, 160, 20)) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    panel.grid.major.y = element_line(color = "grey90")
  )

# 6. Save Output
ggsave(output_file, plot = p, width = 8, height = 6, device = "pdf")
message(paste("Plot saved to:", output_file))
