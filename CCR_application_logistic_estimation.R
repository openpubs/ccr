#####  Define log partial likelihood ####

library(tidyverse)
library(glmnet)
library(writexl)
library(openxlsx)

load("data_1923_ungrp.RData")

##### Parameters #####

wqs_penalty = 1e-1
p = 17
np = 18

initial_weight = c(3, rep(1/p, p))

##### Outcome #####

Y = data_selected_cc_ad$o_depressive_disorder2
table(Y)

##### Survey Weight #####

W = data_selected_cc_ad$s_weight
W = W / mean(W)   # stabilize optimization

##### Design Matrix #####

data.XZ <- data_selected_cc_ad %>%
  select(-matches("^(o_|s_)"))

XZ_only = as.matrix(data.XZ)

names(data.XZ)

##### Weighted Log-Likelihood #####

log_partial_likelihood <- function(theta, X, Y, W) {
  
  beta0 <- theta[1]
  beta <- theta[2]
  w <- theta[3:(p+2)]
  beta.other <- matrix(theta[(p+3):length(theta)], ncol=1)
  
  X.1 = X[,1:p]
  X.2 = X[,(p+1):dim(X)[2]]
  
  q_sums <- rowSums(t(t(X.1) * w))
  
  eta.qz <- beta0 + beta * q_sums + X.2 %*% beta.other
  
  log.est = Y * eta.qz - log(1 + exp(eta.qz))
  
  return_function <- - sum(W * log.est) / sum(W)
  
  return(return_function)
}

##### Objective Function #####

objective_function <- function(theta, X, Y, W, lambda_penalty) {
  
  w <- theta[3:(p+2)]
  
  penalty <- lambda_penalty * (sum(w) - 1)^2
  
  return(log_partial_likelihood(theta, X, Y, W) + penalty)
}

##### Initial Values from GLM #####

lmfit.0 <- glm(
  o_depressive_disorder2 ~ z_sex_female + 
    z_education_high + z_education_attend_univ + z_education_univ + 
    z_income_15_25k + z_income_25_35k + z_income_35_50k + z_income_gt50k + 
    z_age_25_34 + z_age_35_44 + z_age_45_54 + z_age_55_64 + z_age_65plus + 
    z_race_asian + z_race_black + z_race_hispanic + z_race_other + 
    z_race_white,
  data = data_selected_cc_ad,
  family = binomial()
)

initial_theta = c(
  coef(lmfit.0)[1],
  initial_weight,
  coef(lmfit.0)[2:length(coef(lmfit.0))]
)

##### Check Likelihood #####

log_partial_likelihood(initial_theta, XZ_only, Y, W)

##### Optimization #####

result <- nlminb(
  start = initial_theta,
  objective = function(theta)
    objective_function(theta, XZ_only, Y, W, lambda_penalty = wqs_penalty),
  lower = c(-5, -5, rep(0, p), rep(-5, np)),
  upper = c(5, 5, rep(1, p), rep(5, np)),
  control = list(abs.tol = 1e-8)
)

##### Results #####

estimated_intercept <- result$par[1]
estimated_beta <- result$par[2]
estimated_weights <- result$par[3:(p+2)]

variable_names <- c(
  "x_violence_between_parents_once",
  "x_violence_between_parents_multiple",
  "x_parental_separation_yes",
  "x_caregiver_physical_violence_once",
  "x_caregiver_physical_violence_multiple",
  "x_caregiver_verbal_violence_once",
  "x_caregiver_verbal_violence_multiple",
  "x_coerced_sexual_touching_once",
  "x_coerced_sexual_touching_multiple",
  "x_sexually_touched_once",
  "x_sexually_touched_multiple",
  "x_forced_to_have_sex_once",
  "x_forced_to_have_sex_multiple",
  "x_household_depression_suicidality_yes",
  "x_household_alcohol_abuse_yes",
  "x_household_drug_abuse_yes",
  "x_household_prison_history_yes"
)

names(estimated_weights) <- variable_names

##### Save Weights #####

df_weights <- data.frame(Weight = as.numeric(estimated_weights)) * 100
rownames(df_weights) <- names(estimated_weights)

write.xlsx(
  df_weights,
  file = "estimated_weights_survey_weighted.xlsx",
  rowNames = TRUE
)

##### Other Covariates #####

estimated_other <- result$par[(p+3):length(initial_theta)]

##### Weight Sum #####

weight_sum = sum(estimated_weights)

##### Domain Weights #####

domain <- data.frame(
  "Violence between parents" = estimated_weights[1] + estimated_weights[2],
  "Parental Separation" = estimated_weights[3],
  "Caregiver Physical Violence" = estimated_weights[4] + estimated_weights[5],
  "Caregiver Verbal Violence" = estimated_weights[6] + estimated_weights[7],
  "Coerced Sexual Touching" = estimated_weights[8] + estimated_weights[9],
  "Sexually Touched" = estimated_weights[10] + estimated_weights[11],
  "Forced to Have Sex" = estimated_weights[12] + estimated_weights[13],
  "Household Depression/Suicidality" = estimated_weights[14],
  "Household Alcohol Abuse" = estimated_weights[15],
  "Household Drug Abuse" = estimated_weights[16],
  "Household Prison History" = estimated_weights[17],
  check.names = FALSE
)

##### Selection Threshold #####

threshold = 1/11

select.weight <- as.data.frame(
  lapply(domain, function(x) ifelse(x > threshold, 1, 0)),
  check.names = FALSE
)

##### Save Workspace #####

save.image("app_ace_weighted.RData")