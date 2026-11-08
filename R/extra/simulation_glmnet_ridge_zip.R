# Clear environment
rm(list = ls())

# Libraries
library(MASS)
library(pscl)
library(openxlsx)
library(glmnet)
library(dplyr)

# Generate data
generate_X <- function(n, p, rho) {
  Z <- matrix(rnorm(n * p), nrow = n)
  X <- sqrt(1 - rho^2) * Z + rho * matrix(rnorm(n), nrow = n, ncol = p)
  scale(X)
}

generate_Y_zip <- function(X, beta, logit_intercept) {
  eta <- X %*% beta
  mu <- exp(pmin(eta, 13.8))
  pi0 <- plogis(logit_intercept)
  Y_raw <- rpois(nrow(X), mu)
  is_zero <- rbinom(nrow(X), 1, pi0)
  Y <- ifelse(is_zero == 1, 0, Y_raw)
  as.integer(Y)
}

# Ridge ZIP log-likelihood
ridge_zip_ll <- function(params, X, Y, lambda) {
  p <- ncol(X)
  beta <- params[1:p]
  gamma <- params[p + 1]
  eta <- X %*% beta
  mu <- exp(pmin(eta, 13.8))
  pi0 <- plogis(gamma)
  loglik <- ifelse(Y == 0,
                   log(pmax(pi0 + (1 - pi0) * dpois(0, mu), 1e-10)),
                   log(pmax((1 - pi0) * dpois(Y, mu), 1e-10)))
  -sum(loglik) + lambda * sum(beta^2)
}

# Parameters
n_vals <- c(100,200,300)
p_vals <- c(10)
rho_vals <- c(0.2, 0.90, 0.99)
logit_intercepts <- c(0, 1, 2)
nsim <- 10
results <- data.frame()

# Main loop
for (logit_intercept in logit_intercepts) {
  for (n in n_vals) {
    for (p in p_vals) {
      for (rho in rho_vals) {
        cat("Running: n =", n, ", p =", p, ", rho =", rho,
            ", intercept =", logit_intercept, "\n")
        
        mse_mat <- matrix(NA, nrow = nsim, ncol = 5)
        colnames(mse_mat) <- c("ZIP", "ZINB", "HurdlePoisson",
                               "RidgeHurdlePoisson", "RidgeZIP")
        
        for (sim in 1:nsim) {
          set.seed(sim)
          X <- generate_X(n, p, rho)
          XtX <- t(X) %*% X
          eig <- eigen(XtX)
          beta_true <- eig$vectors %*% diag(eig$values)
          beta_true <- beta_true / sqrt(sum(beta_true^2))
          Y <- generate_Y_zip(X, beta_true, logit_intercept)
          df <- data.frame(Y = Y, X)
          
          # ZIP
          fit1 <- try(zeroinfl(Y ~ . | 1, data = df, dist = "poisson"), silent = TRUE)
          if (!inherits(fit1, "try-error")) {
            b <- coef(fit1)[grep("^count_", names(coef(fit1)))]
            if (length(b) >= p) mse_mat[sim, "ZIP"] <- mean((b[1:p] - beta_true)^2)
          }
          
          # ZINB
          fit2 <- try(zeroinfl(Y ~ . | 1, data = df, dist = "negbin"), silent = TRUE)
          if (!inherits(fit2, "try-error")) {
            b <- coef(fit2)[grep("^count_", names(coef(fit2)))]
            if (length(b) >= p) mse_mat[sim, "ZINB"] <- mean((b[1:p] - beta_true)^2)
          }
          
          # Hurdle Poisson
          try({
            logit_fit <- glm((Y > 0) ~ ., data = df, family = binomial)
            count_fit <- glm(Y ~ ., data = df[Y > 0, ], family = poisson)
            b <- coef(count_fit)[-1]
            if (length(b) == p) mse_mat[sim, "HurdlePoisson"] <- mean((b - beta_true)^2)
          }, silent = TRUE)
          
          # Ridge Hurdle Poisson
          try({
            y_pos <- Y > 0
            if (sum(y_pos) > p + 1) {
              lambda_log <- cv.glmnet(X, y_pos, alpha = 0, family = "binomial")$lambda.min
              lambda_pois <- cv.glmnet(X[y_pos, ], Y[y_pos], alpha = 0, family = "poisson")$lambda.min
              
              fit_log <- glmnet(X, y_pos, alpha = 0, family = "binomial", lambda = lambda_log)
              fit_pois <- glmnet(X[y_pos, ], Y[y_pos], alpha = 0, family = "poisson", lambda = lambda_pois)
              beta_pois <- as.vector(coef(fit_pois))[-1]
              
              mse_mat[sim, "RidgeHurdlePoisson"] <- mean((beta_pois - beta_true)^2)
            }
          }, silent = TRUE)
          
          # Ridge ZIP
          try({
            lambda_zip <- cv.glmnet(X, Y, alpha = 0, family = "poisson")$lambda.min
            par_init <- c(rep(0, p), 0)
            ridge_zip_fit <- optim(par_init, ridge_zip_ll, X = X, Y = Y, lambda = lambda_zip,
                                   method = "BFGS", control = list(maxit = 1000))
            if (ridge_zip_fit$convergence == 0) {
              b <- ridge_zip_fit$par[1:p]
              mse_mat[sim, "RidgeZIP"] <- mean((b - beta_true)^2)
            }
          }, silent = TRUE)
        }
        
        # Store averaged MSEs
        results <- rbind(results, data.frame(
          InterceptLogit = logit_intercept,
          n = n, p = p, rho = rho,
          ZIP = mean(mse_mat[, "ZIP"], na.rm = TRUE),
          ZINB = mean(mse_mat[, "ZINB"], na.rm = TRUE),
          HurdlePoisson = mean(mse_mat[, "HurdlePoisson"], na.rm = TRUE),
          RidgeHurdlePoisson = mean(mse_mat[, "RidgeHurdlePoisson"], na.rm = TRUE),
          RidgeZIP = mean(mse_mat[, "RidgeZIP"], na.rm = TRUE)
        ))
      }
    }
  }
}

# Save result
write.xlsx(results, "C:\\Users\\hnayem\\OneDrive - Florida International University\\Documents\\Private\\Phd work\\Folder for Simulation work\\ZI VS Multico\\Results\\zi_sim_full_ridge_models_revised456.xlsx")
