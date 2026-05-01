library(MASS)
library(pROC)
library(openxlsx)

# This script reproduces the unadjusted simulation study for the CCR method.
# It includes both Gaussian and binary outcomes under three correlation settings.

n_sim <- 1000
n <- 1000
p_o <- 10
p <- 15
K <- 5
rho_values <- c(0.2, 0.5, 0.8)

# The first five exposures have three levels and are represented by two dummy variables each.
# The last five exposures have two levels and are represented by one dummy variable each.
beta_true <- c(c(2, 2), c(2, 1), c(2, 0), c(0, 0), c(0, 0),
               2, 0, 0, 0, 0)
beta0_gaussian <- -3
beta0_binary <- -6

ccr_threshold <- 1 / p_o

# Aggregate dummy-specific CCR weights to the exposure level.
# This aggregation is used for identifying key exposures.
to_exposure <- function(w) {
  c(
    sum(w[1:2]), sum(w[3:4]), sum(w[5:6]),
    sum(w[7:8]), sum(w[9:10]),
    w[11], w[12], w[13], w[14], w[15]
  )
}

calc_bias_rmse <- function(est_mat, true_vec) {
  diff_mat <- sweep(est_mat, 2, true_vec, "-")
  list(
    bias = colMeans(diff_mat, na.rm = TRUE),
    rmse = sqrt(colMeans(diff_mat^2, na.rm = TRUE))
  )
}

round_df <- function(df, digits = 3) {
  is_num <- sapply(df, is.numeric)
  df[is_num] <- lapply(df[is_num], round, digits = digits)
  df
}

# Numerically stable computation of log(1 + exp(x)) for logistic likelihood.
safe_log1pexp <- function(x) {
  ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}

# Generate categorical exposures from latent multivariate normal variables,
# then create dummy variables using the same coding scheme as in the manuscript.
make_dummy_matrix <- function(latent) {
  X1_5 <- apply(latent[, 1:5], 2, function(x)
    cut(x, quantile(x, c(0, 0.5, 0.8, 1)),
        labels = 0:2, include.lowest = TRUE))
  
  X6_10 <- apply(latent[, 6:10], 2, function(x)
    cut(x, quantile(x, c(0, 0.6, 1)),
        labels = 0:1, include.lowest = TRUE))
  
  X <- data.frame(cbind(X1_5, X6_10))
  X[] <- lapply(X, function(x) as.numeric(as.character(x)))
  colnames(X) <- paste0("Z", 1:10)
  
  X_factor <- data.frame(
    lapply(X, function(x) factor(x, levels = sort(unique(x))))
  )
  colnames(X_factor) <- paste0("Z", 1:10)
  
  X_dummy <- model.matrix(~ ., data = X_factor)[, -1, drop = FALSE]
  colnames(X_dummy) <- paste0("X", 1:ncol(X_dummy))
  
  if (ncol(X_dummy) != p) {
    stop(sprintf("X_dummy has %d columns, but p = %d", ncol(X_dummy), p))
  }
  
  X_dummy
}

# Fit CCR for a Gaussian outcome by minimizing the mean squared error
# with a penalty enforcing the sum-to-one constraint on the weights.
fit_ccr_gaussian <- function(X, Y, penalty = 1) {
  objective <- function(theta) {
    b0 <- theta[1]
    b <- theta[2]
    w <- theta[3:(p + 2)]
    pred <- b0 + b * (X %*% w)
    mean((Y - pred)^2) + penalty * (sum(w) - 1)^2
  }
  
  init <- c(mean(Y), 1, rep(1 / p, p))
  nlminb(
    start = init,
    objective = objective,
    lower = c(-100, -100, rep(0, p)),
    upper = c(100, 100, rep(1, p))
  )
}

# Fit CCR for a binary outcome by minimizing the negative logistic log-likelihood
# with the same sum-to-one weight penalty.
fit_ccr_binary <- function(X, Y, penalty = 0.1) {
  objective <- function(theta) {
    b0 <- theta[1]
    b <- theta[2]
    w <- theta[3:(p + 2)]
    eta <- as.numeric(b0 + b * (X %*% w))
    neg_loglik <- -mean(Y * eta - safe_log1pexp(eta))
    neg_loglik + penalty * (sum(w) - 1)^2
  }
  
  init_prob <- min(max(mean(Y), 1e-4), 1 - 1e-4)
  init <- c(qlogis(init_prob), 1, rep(1 / p, p))
  nlminb(
    start = init,
    objective = objective,
    lower = c(-100, -100, rep(0, p)),
    upper = c(100, 100, rep(1, p))
  )
}

# Run one simulation setting defined by outcome type and exposure correlation rho.
run_one_setting <- function(outcome = c("gaussian", "binary"), rho) {
  outcome <- match.arg(outcome)
  
  select_ccr <- rep(0, p_o)
  pred_metric <- rep(NA_real_, n_sim)
  logloss <- rep(NA_real_, n_sim)
  
  beta_ccr_mat <- matrix(NA_real_, n_sim, p)
  zero_prop_ccr_mat <- matrix(NA_real_, n_sim, p)
  
  warning_log <- data.frame(
    sim = integer(),
    fold = integer(),
    message = character(),
    stringsAsFactors = FALSE
  )
  
  log_warning <- function(sim, fold, w) {
    warning_log <<- rbind(
      warning_log,
      data.frame(
        sim = sim,
        fold = fold,
        message = conditionMessage(w),
        stringsAsFactors = FALSE
      )
    )
  }
  
  for (sim in 1:n_sim) {
    set.seed(sim)
    
    # Generate correlated categorical exposures through latent Gaussian variables.
    corr_matrix <- matrix(rho, nrow = p_o, ncol = p_o)
    diag(corr_matrix) <- 1
    latent <- mvrnorm(n = n, mu = rep(0, p_o), Sigma = corr_matrix)
    X_dummy <- make_dummy_matrix(latent)
    
    # Generate the outcome under the corresponding CCR data-generating model.
    if (outcome == "gaussian") {
      Y <- as.numeric(beta0_gaussian + X_dummy %*% beta_true + rnorm(n))
      fit <- withCallingHandlers(
        fit_ccr_gaussian(X_dummy, Y),
        warning = function(w) {
          log_warning(sim, NA, w)
          invokeRestart("muffleWarning")
        }
      )
    } else {
      eta <- as.numeric(beta0_binary + X_dummy %*% beta_true)
      prob <- 1 / (1 + exp(-eta))
      Y <- rbinom(n, 1, prob)
      fit <- withCallingHandlers(
        fit_ccr_binary(X_dummy, Y),
        warning = function(w) {
          log_warning(sim, NA, w)
          invokeRestart("muffleWarning")
        }
      )
    }
    
    # Store CCR estimates from the full sample fit.
    b_hat <- fit$par[2]
    w_hat <- fit$par[3:(p + 2)]
    
    beta_ccr_mat[sim, ] <- b_hat * w_hat
    zero_prop_ccr_mat[sim, ] <- as.numeric(w_hat < 1e-4)
    
    # Identify key exposures using the aggregated-weight threshold 1 / number of exposures.
    exposure_w <- to_exposure(w_hat)
    select_ccr <- select_ccr + (exposure_w > ccr_threshold)
    
    # Evaluate predictive performance using K-fold cross-validation.
    folds <- sample(rep(1:K, length.out = n))
    pred_ccr <- rep(NA_real_, n)
    
    for (k in 1:K) {
      test <- which(folds == k)
      train <- setdiff(1:n, test)
      
      X_train <- X_dummy[train, , drop = FALSE]
      X_test <- X_dummy[test, , drop = FALSE]
      Y_train <- Y[train]
      
      if (outcome == "gaussian") {
        fit_cv <- withCallingHandlers(
          fit_ccr_gaussian(X_train, Y_train),
          warning = function(w) {
            log_warning(sim, k, w)
            invokeRestart("muffleWarning")
          }
        )
        pred_ccr[test] <- as.numeric(
          fit_cv$par[1] + fit_cv$par[2] * (X_test %*% fit_cv$par[3:(p + 2)])
        )
      } else {
        if (length(unique(Y_train)) < 2) next
        fit_cv <- withCallingHandlers(
          fit_ccr_binary(X_train, Y_train),
          warning = function(w) {
            log_warning(sim, k, w)
            invokeRestart("muffleWarning")
          }
        )
        eta_test <- as.numeric(
          fit_cv$par[1] + fit_cv$par[2] * (X_test %*% fit_cv$par[3:(p + 2)])
        )
        pred_ccr[test] <- 1 / (1 + exp(-eta_test))
      }
    }
    
    valid <- !is.na(pred_ccr)
    
    if (outcome == "gaussian") {
      pred_metric[sim] <- mean((Y[valid] - pred_ccr[valid])^2)
    } else if (sum(valid) > 0 && length(unique(Y[valid])) == 2) {
      roc_obj <- roc(Y[valid], pred_ccr[valid], quiet = TRUE)
      pred_metric[sim] <- as.numeric(auc(roc_obj))
      logloss[sim] <- -mean(
        Y[valid] * log(pred_ccr[valid] + 1e-8) +
          (1 - Y[valid]) * log(1 - pred_ccr[valid] + 1e-8)
      )
    }
    
    if (sim %% 50 == 0) {
      message("Completed ", outcome, ", rho = ", rho, ", simulation ", sim)
    }
  }
  
  # Summarize identification, coefficient estimation, prediction, and zero-weight frequency.
  beta_res <- calc_bias_rmse(beta_ccr_mat, beta_true)
  zero_ccr <- colMeans(zero_prop_ccr_mat, na.rm = TRUE)
  
  list(
    outcome = outcome,
    rho = rho,
    selection = select_ccr / n_sim,
    prediction = pred_metric,
    logloss = logloss,
    beta = beta_res,
    zero = zero_ccr,
    warning = warning_log
  )
}

make_tables <- function(res) {
  prediction_name <- ifelse(res$outcome == "gaussian", "MSE", "AUC")
  prediction_value <- mean(res$prediction, na.rm = TRUE)
  
  tables <- list(
    selection = data.frame(
      Exposure = paste0("V", 1:p_o),
      CCR = res$selection
    ),
    prediction = data.frame(
      Outcome = res$outcome,
      Rho = res$rho,
      Metric = prediction_name,
      CCR = prediction_value
    ),
    beta_bias = data.frame(
      Beta = paste0("beta", 1:p),
      CCR = res$beta$bias
    ),
    beta_rmse = data.frame(
      Beta = paste0("beta", 1:p),
      CCR = res$beta$rmse
    ),
    zero = data.frame(
      Weight = paste0("w", 1:p),
      CCR = res$zero
    ),
    warning_detail = res$warning
  )
  
  if (res$outcome == "binary") {
    tables$logloss <- data.frame(
      Outcome = res$outcome,
      Rho = res$rho,
      Metric = "LogLoss",
      CCR = mean(res$logloss, na.rm = TRUE)
    )
  }
  
  lapply(tables, round_df, digits = 3)
}

# Run all unadjusted simulation settings and export the results.
all_results <- list()
wb <- createWorkbook()
num_style <- createStyle(numFmt = "0.000")

for (outcome in c("gaussian", "binary")) {
  for (rho in rho_values) {
    res <- run_one_setting(outcome = outcome, rho = rho)
    all_results[[paste(outcome, rho, sep = "_rho")]] <- res
    
    tables <- make_tables(res)
    outcome_label <- ifelse(outcome == "gaussian", "G", "B")
    rho_label <- gsub("\\.", "", as.character(rho))
    prefix <- paste0(outcome_label, "_r", rho_label, "_")
    
    for (table_name in names(tables)) {
      sheet_name <- substr(paste0(prefix, table_name), 1, 31)
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet_name, tables[[table_name]])
      
      if (ncol(tables[[table_name]]) >= 2 && nrow(tables[[table_name]]) >= 1) {
        addStyle(
          wb, sheet_name, style = num_style,
          rows = 2:(nrow(tables[[table_name]]) + 1),
          cols = 2:ncol(tables[[table_name]]),
          gridExpand = TRUE, stack = TRUE
        )
      }
    }
  }
}

saveWorkbook(wb, "simulation_unadjusted_results.xlsx", overwrite = TRUE)

