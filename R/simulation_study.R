# Clear environment
rm(list = ls())

# Libraries
library(MASS)
library(pscl)
library(openxlsx)
library(dplyr)

# Generate predictors
generate_X <- function(n, p, rho) {
  Z <- matrix(rnorm(n * p), nrow = n)
  X <- sqrt(1 - rho^2) * Z + rho * matrix(rnorm(n), nrow = n, ncol = p)
  scale(X)
}

# Generate zero-inflated Negative Binomial outcomes
generate_Y_zinb <- function(X, beta, logit_intercept, theta) {
  eta <- X %*% beta
  mu <- exp(pmin(eta, 13.8))
  pi0 <- plogis(logit_intercept)
  Y_raw <- rnbinom(nrow(X), mu = mu, size = theta)
  is_zero <- rbinom(nrow(X), 1, pi0)
  Y <- ifelse(is_zero == 1, 0, Y_raw)
  as.integer(Y)
}

# Ridge Poisson estimation
ridge_poisson_weights <- function(X, y, lambda) {
  p <- ncol(X)
  XtX <- t(X) %*% X
  XtY <- t(X) %*% y
  beta_ridge <- solve(XtX + lambda * diag(p), XtY)
  return(beta_ridge)
}

# Ridge Logistic regression via IRLS
ridge_logistic_weights <- function(X, y, lambda, max_iter = 25) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)
  for (i in 1:max_iter) {
    eta <- X %*% beta
    mu <- plogis(eta)
    W <- diag(as.vector(mu * (1 - mu)), n, n)
    z <- eta + (y - mu) / pmax(mu * (1 - mu), 1e-5)
    XtWX <- t(X) %*% W %*% X + lambda * diag(p)
    XtWz <- t(X) %*% W %*% z
    beta <- solve(XtWX, XtWz)
  }
  return(beta)
}

# Cross-validated lambda for Poisson Ridge
cv_ridge_poisson <- function(X, y, lambdas = 10^seq(-4, 2, length.out = 10), k = 5) {
  n <- nrow(X)
  folds <- sample(rep(1:k, length.out = n))
  mse_mat <- matrix(NA, nrow = k, ncol = length(lambdas))
  
  for (i in 1:k) {
    test_idx <- which(folds == i)
    train_idx <- setdiff(1:n, test_idx)
    X_train <- X[train_idx, ]
    y_train <- y[train_idx]
    X_test <- X[test_idx, ]
    y_test <- y[test_idx]
    
    for (j in seq_along(lambdas)) {
      beta <- ridge_poisson_weights(X_train, y_train, lambdas[j])
      y_pred <- exp(X_test %*% beta)
      mse_mat[i, j] <- mean((y_test - y_pred)^2)
    }
  }
  lambda_best <- lambdas[which.min(colMeans(mse_mat))]
  return(lambda_best)
}

# Cross-validated lambda for Logistic Ridge
cv_ridge_logistic <- function(X, y, lambdas = 10^seq(-4, 2, length.out = 10), k = 5) {
  n <- nrow(X)
  folds <- sample(rep(1:k, length.out = n))
  logloss_mat <- matrix(NA, nrow = k, ncol = length(lambdas))
  
  for (i in 1:k) {
    test_idx <- which(folds == i)
    train_idx <- setdiff(1:n, test_idx)
    X_train <- X[train_idx, ]
    y_train <- y[train_idx]
    X_test <- X[test_idx, ]
    y_test <- y[test_idx]
    
    for (j in seq_along(lambdas)) {
      beta <- ridge_logistic_weights(X_train, y_train, lambdas[j])
      p_pred <- plogis(X_test %*% beta)
      logloss_mat[i, j] <- -mean(y_test * log(pmax(p_pred, 1e-6)) + (1 - y_test) * log(pmax(1 - p_pred, 1e-6)))
    }
  }
  lambda_best <- lambdas[which.min(colMeans(logloss_mat))]
  return(lambda_best)
}

# Parameters
n_vals <- c(100,200,500)
p_vals <- c(10,20)
rho_vals <- c(0.8,0.90,0.95,0.99)
logit_intercepts <- c(0, 1)
over_dis_vals <- c(1, 5)
nsim <- 10
results <- data.frame()

# Main loop
for (logit_intercept in logit_intercepts) {
  for (n in n_vals) {
    for (p in p_vals) {
      for (rho in rho_vals) {
        for (theta in over_dis_vals) {
          cat("Running: n =", n, ", p =", p, ", rho =", rho,
              ", intercept =", logit_intercept, ", theta =", theta, "\n")
          
          mse_mat <- matrix(NA, nrow = nsim, ncol = 5)
          colnames(mse_mat) <- c("ZIP", "ZINB", "HurdlePoisson", "RidgeHurdlePoisson", "RidgeZIP")
          
          for (sim in 1:nsim) {
            set.seed(sim)
            X <- generate_X(n, p, rho)
            XtX <- t(X) %*% X
            eig <- eigen(XtX)
            beta_true <- eig$vectors %*% diag(eig$values)
            beta_true <- beta_true / sqrt(sum(beta_true^2))
            gamma_true <- eig$vectors %*% diag(eig$values)
            gamma_true <- gamma_true / sqrt(sum(gamma_true^2))
            
            Y <- generate_Y_zinb(X, beta_true, logit_intercept, theta)
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
            
            # Ridge projection function
            ridge_projection <- function(X, beta_mle, lambda) {
              p <- ncol(X)
              XtX <- t(X) %*% X
              ridge_mat <- solve(XtX + lambda * diag(p)) %*% XtX
              beta_ridge <- ridge_mat %*% beta_mle
              return(as.vector(beta_ridge))
            }
            
            # Ridge Hurdle Poisson
            try({
              y_pos <- Y > 0
              if (sum(y_pos) > p + 1) {
                glm_logit <- glm(y_pos ~ X, family = binomial)
                lambda_log <- cv_ridge_logistic(X, y_pos)
                beta_logit_mle <- coef(glm_logit)[-1]
                beta_ridge_logit <- ridge_projection(X, beta_logit_mle, lambda_log)
                
                glm_pois <- glm(Y[y_pos] ~ X[y_pos, ], family = poisson)
                lambda_pois <- cv_ridge_poisson(X[y_pos, ], Y[y_pos])
                beta_pois_mle <- coef(glm_pois)[-1]
                beta_ridge_pois <- ridge_projection(X[y_pos, ], beta_pois_mle, lambda_pois)
                
                mse_logit <- mean((beta_ridge_logit - gamma_true)^2)
                mse_pois  <- mean((beta_ridge_pois - beta_true)^2)
                mse_mat[sim, "RidgeHurdlePoisson"] <- (mse_logit + mse_pois) / 2
              }
            }, silent = TRUE)
            
            # Ridge ZIP
            try({
              glm_logit <- glm((Y == 0) ~ X, family = binomial)
              lambda_log <- cv_ridge_logistic(X, Y == 0)
              beta_logit_mle <- coef(glm_logit)[-1]
              beta_ridge_logit <- ridge_projection(X, beta_logit_mle, lambda_log)
              
              glm_zip <- glm(Y ~ X, family = poisson)
              lambda_zip <- cv_ridge_poisson(X, Y)
              beta_zip_mle <- coef(glm_zip)[-1]
              beta_ridge_zip <- ridge_projection(X, beta_zip_mle, lambda_zip)
              
              mse_logit <- mean((beta_ridge_logit - gamma_true)^2)
              mse_pois  <- mean((beta_ridge_zip - beta_true)^2)
              mse_mat[sim, "RidgeZIP"] <- (mse_logit + mse_pois) / 2
            }, silent = TRUE)
          }
          
          results <- rbind(results, data.frame(
            InterceptLogit = logit_intercept,
            n = n, p = p, rho = rho, theta = theta,
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
}

# Save results
write.xlsx(results, "C:\\Users\\hnayem\\OneDrive - Florida International University\\Documents\\Private\\Phd work\\Folder for Simulation work\\ZI VS Multico\\Results\\zi_sim_full_ridge_models_revised_CV_with_overdispersion23456.xlsx")
