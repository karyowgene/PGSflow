#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Input: Raw PRS scores, reference panel (e.g., 1000 Genomes PRS)
# Output: Z-scores and percentiles
library(data.table)

prs_raw <- fread(args[1])
ref_data <- fread(args[2])  # Reference population PRS distribution

# Calculate Z-scores
mean_ref <- mean(ref_data$PRS)
sd_ref <- sd(ref_data$PRS)
prs_raw[, PRS_Z := (PRS - mean_ref) / sd_ref]

# Calculate percentile
prs_raw[, PRS_Percentile := round(ecdf(ref_data$PRS)(PRS) * 100, by = .I]

# Risk classification (configurable thresholds)
risk_thresholds <- c(High = 90, Moderate = 75, Low = 0)
prs_raw[, Risk := fcase(
  PRS_Percentile >= risk_thresholds["High"], "High",
  PRS_Percentile >= risk_thresholds["Moderate"], "Moderate",
  default = "Low"
)]

# Save output
fwrite(prs_raw, file = args[3], sep = "\t")

# Optional: Covariate adjustment (BMI/smoking)
if (length(args) > 3) {
  covariates <- fread(args[4])
  merged_data <- merge(prs_raw, covariates, by = "IID")
  model <- lm(PRS ~ BMI + Smoking, data = merged_data)
  merged_data[, Adjusted_PRS := residuals(model)]
  fwrite(merged_data, file = "output/prs_scores/adjusted_prs.tsv", sep = "\t")
}
