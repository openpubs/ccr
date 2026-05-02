library(MASS)
library(boot)

set.seed(123)
n_sim <- 1000
n <- 1000
p <- 2 * 5 + 5
p.o <- 10
cci_penalty <- 0.5
initial_weight <- c(0, rep(1 / p, p))
beta.0 <- -5

beta_true <- c(c(2, 2), c(2, 1), c(2, 0), c(0, 0), c(0, 0),
               2, 0, 0, 0, 0) * 0

type1_errors <- 0

corr_matrix <- matrix(0.8, nrow = p.o, ncol = p.o) # change to 0.5 and 0.2 for results under different correlation result
diag(corr_matrix) <- 1

for (sim in 1:n_sim) {
  latent_data <- mvrnorm(n = n, mu = rep(0, p.o), Sigma = corr_matrix)
  
  X1_to_X5 <- apply(latent_data[, 1:5], 2, function(x) {
    cut(x, breaks = quantile(x, probs = c(0, 0.5, 0.8, 1)),
        labels = 0:2, include.lowest = TRUE)
  })
  
  X6_to_X10 <- apply(latent_data[, 6:10], 2, function(x) {
    cut(x, breaks = quantile(x, probs = c(0, 0.6, 1)),
        labels = 0:1, include.lowest = TRUE)
  })
  
  X <- data.frame(cbind(X1_to_X5, X6_to_X10))
  X[] <- lapply(X, function(col) as.numeric(as.character(col)))
  X_factor <- data.frame(lapply(X, function(x) factor(x, levels = sort(unique(x)))))
  X_dummy <- model.matrix(~ . - 1, data = X_factor)[, -1]
  colnames(X_dummy) <- paste0("X", 1:ncol(X_dummy))
  
  Z1 <- rbinom(n, 1, 0.4)
  Z2 <- rnorm(n, 40, 5)
  Z <- cbind(Z1, Z2)
  
  linp <- beta.0 + X_dummy %*% beta_true - 1.5 * Z1 + 0.1 * Z2
  sigma_noise <- 1
  Y <- as.numeric(linp + rnorm(n, mean = 0, sd = sigma_noise))
  
  data.Fx <- data.frame(X_dummy, Z, Y)
  
  bootstrap_function <- function(data, indices) {
    boot_sample <- data[indices, ]
    boot_Y <- boot_sample$Y
    boot_XZ <- as.matrix(boot_sample[, 1:17])
    boot_Z <- as.matrix(boot_sample[, 16:17])
    
    log_partial_likelihood <- function(theta, X, Y) {
      beta0 <- theta[1]
      beta <- theta[2]
      w <- theta[3:(p + 2)]
      beta.other <- matrix(theta[(p + 3):length(theta)], ncol = 1)
      q_sums <- rowSums(t(t(X[, 1:p]) * w))
      eta <- beta0 + beta * q_sums + boot_Z %*% beta.other
      return(mean((Y - eta)^2))
    }
    
    objective_function <- function(theta, X, Y, lambda_penalty) {
      w <- theta[3:(p + 2)]
      penalty <- lambda_penalty * (sum(w) - 1)^2
      return(log_partial_likelihood(theta, X, Y) + penalty)
    }
    
    lmfit.0 <- glm(Y ~ Z1 + Z2, data = boot_sample, family = gaussian())
    initial_theta <- c(coef(lmfit.0)[1], initial_weight,
                       coef(lmfit.0)[2:length(coef(lmfit.0))])
    
    result <- nlminb(
      start = initial_theta,
      objective = function(theta) objective_function(theta, boot_XZ, boot_Y, cci_penalty),
      lower = c(-100, -100, rep(0, p), -100, -100),
      upper = c(100, 100, rep(1, p), 100, 100),
      control = list(abs.tol = 1e-8)
    )
    
    return(result$par)
  }
  
  bootstrap_results <- boot(data = data.Fx,
                            statistic = bootstrap_function,
                            R = 500)
  beta.estimate <- bootstrap_results$t[, 2]
  beta.ci <- c(quantile(beta.estimate, probs = 0.025),
               quantile(beta.estimate, probs = 0.975))
  
  type1_errors <- ifelse(beta.ci[1] > 0 | beta.ci[2] < 0, 1, 0) + type1_errors
}

save.image("test_gau_1000.RData")

library(MASS)
library(boot)

set.seed(123)
n_sim <- 1000
n <- 1000
p <- 2 * 5 + 5
p.o <- 10
cci_penalty <- 0.5
initial_weight <- c(0, rep(1 / p, p))
beta.0 <- -5

beta_true <- c(c(2, 2), c(2, 1), c(2, 0), c(0, 0), c(0, 0),
               2, 0, 0, 0, 0) * 0

type1_errors <- 0

corr_matrix <- matrix(0.8, nrow = p.o, ncol = p.o) # change to 0.5 and 0.2 for results under different correlation result
diag(corr_matrix) <- 1

for (sim in 1:n_sim) {
  latent_data <- mvrnorm(n = n, mu = rep(0, p.o), Sigma = corr_matrix)
  
  X1_to_X5 <- apply(latent_data[, 1:5], 2, function(x) {
    cut(x, breaks = quantile(x, probs = c(0, 0.5, 0.8, 1)),
        labels = 0:2, include.lowest = TRUE)
  })
  
  X6_to_X10 <- apply(latent_data[, 6:10], 2, function(x) {
    cut(x, breaks = quantile(x, probs = c(0, 0.6, 1)),
        labels = 0:1, include.lowest = TRUE)
  })
  
  X <- data.frame(cbind(X1_to_X5, X6_to_X10))
  X[] <- lapply(X, function(col) as.numeric(as.character(col)))
  X_factor <- data.frame(lapply(X, function(x) factor(x, levels = sort(unique(x)))))
  X_dummy <- model.matrix(~ . - 1, data = X_factor)[, -1]
  colnames(X_dummy) <- paste0("X", 1:ncol(X_dummy))
  
  Z1 <- rbinom(n, 1, 0.4)
  Z2 <- rnorm(n, 40, 5)
  Z <- cbind(Z1, Z2)
  
  linp <- beta.0 + X_dummy %*% beta_true - 1.5 * Z1 + 0.1 * Z2
  pY <- 1 / (1 + exp(-linp))
  Y <- rbinom(n, 1, pY)
  
  data.Fx <- data.frame(X_dummy, Z, Y)
  
  bootstrap_function <- function(data, indices) {
    boot_sample <- data[indices, ]
    boot_Y <- boot_sample$Y
    boot_XZ <- as.matrix(boot_sample[, 1:17])
    boot_Z <- as.matrix(boot_sample[, 16:17])
    
    log_partial_likelihood <- function(theta, X, Y) {
      beta0 <- theta[1]
      beta <- theta[2]
      w <- theta[3:(p + 2)]
      beta.other <- matrix(theta[(p + 3):length(theta)], ncol = 1)
      q_sums <- rowSums(t(t(X[, 1:p]) * w))
      eta <- beta0 + beta * q_sums + boot_Z %*% beta.other
      log.est <- Y * eta - log(1 + exp(eta))
      return(-mean(log.est))
    }
    
    objective_function <- function(theta, X, Y, lambda_penalty) {
      w <- theta[3:(p + 2)]
      penalty <- lambda_penalty * (sum(w) - 1)^2
      return(log_partial_likelihood(theta, X, Y) + penalty)
    }
    
    lmfit.0 <- glm(Y ~ Z1 + Z2, data = boot_sample, family = binomial())
    initial_theta <- c(coef(lmfit.0)[1], initial_weight,
                       coef(lmfit.0)[2:length(coef(lmfit.0))])
    
    result <- nlminb(
      start = initial_theta,
      objective = function(theta) objective_function(theta, boot_XZ, boot_Y, cci_penalty),
      lower = c(-100, -100, rep(0, p), -100, -100),
      upper = c(100, 100, rep(1, p), 100, 100),
      control = list(abs.tol = 1e-8)
    )
    
    return(result$par)
  }
  
  bootstrap_results <- boot(data = data.Fx,
                            statistic = bootstrap_function,
                            R = 500)
  beta.estimate <- bootstrap_results$t[, 2]
  beta.ci <- c(quantile(beta.estimate, probs = 0.025),
               quantile(beta.estimate, probs = 0.975))
  
  type1_errors <- ifelse(beta.ci[1] > 0 | beta.ci[2] < 0, 1, 0) + type1_errors
}

save.image("test_logit_1000.RData")
