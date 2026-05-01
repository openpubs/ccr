library(openxlsx)

############################################################
# Load bootstrap results and application results
############################################################

load("app_boot_weighted.RData")
load("app_ace_weighted.RData")

############################################################
# Covariate names in the order of the model
############################################################

covariate_names <- c(
  "Female vs. Male",
  "High vs. Below High",
  "Attend college vs. Below High",
  "Graduate from college vs. Below High",
  "15-25 vs. <15",
  "25-35 vs. <15",
  "35-50 vs. <15",
  ">50 vs. <15",
  "25-34 vs. 18-24",
  "35-44 vs. 18-24",
  "45-54 vs. 18-24",
  "55-64 vs. 18-24",
  ">65 vs. 18-24",
  "Asian vs. AIAN",
  "African American vs. AIAN",
  "Hispanic vs. AIAN",
  "Other vs. AIAN",
  "White vs. AIAN"
)

colnames(boot_other) <- covariate_names
names(estimated_other) <- covariate_names

############################################################
# Bootstrap CI for CCR index
############################################################

beta_result <- data.frame(
  Variable = "Adverse Childhood Experience Index",
  Estimate = round(estimated_beta, 3),
  CI_lower = round(quantile(boot_beta, probs = 0.025, na.rm = TRUE), 3),
  CI_upper = round(quantile(boot_beta, probs = 0.975, na.rm = TRUE), 3)
)

############################################################
# Bootstrap CI for other covariates
############################################################

covariate_result <- data.frame(
  Variable = covariate_names,
  Estimate = round(estimated_other, 3),
  CI_lower = round(apply(boot_other, 2, quantile, probs = 0.025, na.rm = TRUE), 3),
  CI_upper = round(apply(boot_other, 2, quantile, probs = 0.975, na.rm = TRUE), 3)
)

############################################################
# Final numeric output
############################################################

app_CCR_numbers <- rbind(beta_result, covariate_result)

############################################################
# Save numeric output
############################################################

write.xlsx(
  app_CCR_numbers,
  file = "app_CCR_logOR_95CI_numbers.xlsx",
  rowNames = FALSE,
  overwrite = TRUE
)

############################################################
# Print numeric output
############################################################

app_CCR_numbers