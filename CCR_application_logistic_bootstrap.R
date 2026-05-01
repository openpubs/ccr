library(dplyr)
library(boot)

load("data_1923_ungrp.RData")

# data_selected_cc_ad = data_selected_cc_ad[1:1000,]

############################################################
# Setup
############################################################
wqs_penalty <- 1e-1
p <- 17
np <- 18
initial_weight <- c(1, rep(1 / p, p))  # beta, and weights

############################################################
# Outcome, survey weight, and design matrix
############################################################
Y_all <- data_selected_cc_ad$o_depressive_disorder2

W_all <- data_selected_cc_ad$s_weight
W_all <- W_all / mean(W_all)   # stabilize optimization

data.XZ_all <- data_selected_cc_ad %>%
  select(-matches("^(o_|s_)"))

XZ_only_all <- as.matrix(data.XZ_all)

############################################################
# Weighted log-likelihood function
############################################################
log_partial_likelihood <- function(theta, X, Y, W) {
  beta0 <- theta[1]
  beta <- theta[2]
  w <- theta[3:(p + 2)]
  beta.other <- matrix(theta[(p + 3):length(theta)], ncol = 1)
  
  X.1 <- X[, 1:p, drop = FALSE]
  X.2 <- X[, (p + 1):ncol(X), drop = FALSE]
  
  q_sums <- rowSums(t(t(X.1) * w))
  eta.qz <- beta0 + beta * q_sums + X.2 %*% beta.other
  log.est <- Y * eta.qz - log(1 + exp(eta.qz))
  
  return(-sum(W * log.est) / sum(W))   # weighted negative log-likelihood
}

############################################################
# Penalized objective
############################################################
objective_function <- function(theta, X, Y, W, lambda_penalty) {
  w <- theta[3:(p + 2)] 
  penalty <- lambda_penalty * (sum(w) - 1)^2
  return(log_partial_likelihood(theta, X, Y, W) + penalty)
}

############################################################
# Initial GLM for starting values
############################################################
lmfit.0 <- glm(
  o_depressive_disorder2 ~ z_sex_female + 
    z_education_high + z_education_attend_univ + z_education_univ + 
    z_income_15_25k + z_income_25_35k + z_income_35_50k + z_income_gt50k + 
    z_age_25_34 + z_age_35_44 + z_age_45_54 + z_age_55_64 + z_age_65plus + 
    z_race_asian + z_race_black + z_race_hispanic + z_race_other + 
    z_race_white,
  data = data_selected_cc_ad,
  family = binomial(),
  weights = s_weight
)

initial_theta <- c(coef(lmfit.0)[1], initial_weight, coef(lmfit.0)[-1])

############################################################
# Bootstrap function
############################################################
bootstrap_wqs <- function(B = 1000) {
  set.seed(123456789)
  n <- nrow(data_selected_cc_ad)
  
  results <- vector("list", B)
  
  for (b in 1:B) {
    # Resample with replacement
    idx <- sample(1:n, size = n, replace = TRUE)
    
    Y <- Y_all[idx]
    XZ_only <- XZ_only_all[idx, , drop = FALSE]
    W <- W_all[idx]
    W <- W / mean(W)   # re-stabilize within bootstrap sample
    
    # Run optimization
    fit <- tryCatch({
      nlminb(
        start = initial_theta,
        objective = function(theta) {
          objective_function(
            theta,
            XZ_only,
            Y,
            W,
            lambda_penalty = wqs_penalty
          )
        },
        lower = c(-5, -5, rep(0, p), rep(-5, np)),
        upper = c(5, 5, rep(1, p), rep(5, np)),
        control = list(abs.tol = 1e-8)
      )
    }, error = function(e) return(NULL))
    
    if (!is.null(fit)) {
      results[[b]] <- list(
        beta0 = fit$par[1],
        beta = fit$par[2],
        weights = fit$par[3:(p + 2)],
        beta_other = fit$par[(p + 3):length(fit$par)]
      )
    }
  }
  
  # Filter out failed fits
  results <- results[!sapply(results, is.null)]
  return(results)
}

############################################################
# Run bootstrap
############################################################
B <- 1000  
boot_results <- bootstrap_wqs(B = B)

############################################################
# Convert results to matrix for summary
############################################################
boot_weights <- do.call(rbind, lapply(boot_results, function(x) x$weights))
boot_beta <- sapply(boot_results, function(x) x$beta)
boot_other <- do.call(rbind, lapply(boot_results, function(x) x$beta_other))

############################################################
# Variable names
############################################################
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

colnames(boot_weights) <- variable_names

# Summary (e.g., mean weights)
# colMeans(boot_weights)

save.image("app_boot_weighted.RData")