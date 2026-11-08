# Clear environment
rm(list = ls())

# Libraries
library(pscl)
library(dplyr)
library(openxlsx)

# Load data
data("bioChemists")
df <- bioChemists

# Convert categorical to binary
df$fem <- as.numeric(df$fem == "Women")
df$mar <- as.numeric(df$mar == "Married")

############


# Correlation matrix
cor_matrix <- cor(X)

# Correlation plot
library(corrplot)
corrplot::corrplot(cor_matrix, method = "number", addCoef.col = "black", 
                   tl.col = "black", number.cex = 1.5, type = "lower")


# Response variable
y <- df$art

# Mean and variance
mean_y <- mean(y)
var_y <- var(y)

# Overdispersion ratio
overdispersion_ratio <- var_y / mean_y
overdispersion_ratio

# Fit Poisson model
poisson_fit <- glm(art ~ fem + mar + kid5 + phd + ment, family = poisson, data = df)

# Predicted lambda
lambda_hat <- predict(poisson_fit, type = "response")

# Expected number of zeros
expected_zeros <- sum(dpois(0, lambda_hat))
expected_zeros



sum(df$art == 0)  # Number of observed zeros


#   Condition Number

# Design matrix (excluding the intercept)
XX <- model.matrix(poisson_fit)[, -1]  # Remove intercept column

# Condition number
cond_number <- kappa(XX)
cond_number



###########

# Design matrix
X_raw <- as.matrix(df[, c("fem", "mar", "kid5", "phd", "ment")])
X <- scale(X_raw)
Y <- df$art
n <- nrow(X)
p <- ncol(X)

# Eigen-based true beta
XtX <- t(X) %*% X
eig <- eigen(XtX)
beta_true <- eig$vectors %*% diag(eig$values)
beta_true <- beta_true / sqrt(sum(beta_true^2))  # normalize

# Ridge projection helper
ridge_projection <- function(X, beta_mle, lambda) {
  p <- ncol(X)
  XtX <- t(X) %*% X
  ridge_mat <- solve(XtX + lambda * diag(p)) %*% XtX
  beta_ridge <- ridge_mat %*% beta_mle
  return(as.vector(beta_ridge))
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
      logloss_mat[i, j] <- -mean(y_test * log(pmax(p_pred, 1e-6)) + 
                                   (1 - y_test) * log(pmax(1 - p_pred, 1e-6)))
    }
  }
  lambda_best <- lambdas[which.min(colMeans(logloss_mat))]
  return(lambda_best)
}

# Model results storage
mse_results <- data.frame(Method = character(), MSE = numeric())

# ZIP
fit_zip <- zeroinfl(Y ~ X | 1, dist = "poisson")
beta_zip <- coef(fit_zip)[grep("^count_", names(coef(fit_zip)))]
mse_zip <- mean((beta_zip - beta_true)^2)
mse_results <- rbind(mse_results, data.frame(Method = "ZIP", MSE = mse_zip))

# ZINB
fit_zinb <- zeroinfl(Y ~ X | 1, dist = "negbin")
beta_zinb <- coef(fit_zinb)[grep("^count_", names(coef(fit_zinb)))]
mse_zinb <- mean((beta_zinb - beta_true)^2)
mse_results <- rbind(mse_results, data.frame(Method = "ZINB", MSE = mse_zinb))

# Hurdle Poisson
y_pos <- Y > 0
fit_hurdle <- glm(Y ~ X, family = poisson, subset = y_pos)
beta_hurdle <- coef(fit_hurdle)[-1]
mse_hurdle <- mean((beta_hurdle - beta_true)^2)
mse_results <- rbind(mse_results, data.frame(Method = "HurdlePoisson", MSE = mse_hurdle))

# Ridge ZIP
lambda_log <- cv_ridge_logistic(X, Y == 0)
lambda_pois <- cv_ridge_poisson(X, Y)

glm_logit <- glm((Y == 0) ~ X, family = binomial)
beta_logit_mle <- coef(glm_logit)[-1]
beta_ridge_logit <- ridge_projection(X, beta_logit_mle, lambda_log)

glm_pois <- glm(Y ~ X, family = poisson)
beta_pois_mle <- coef(glm_pois)[-1]
beta_ridge_pois <- ridge_projection(X, beta_pois_mle, lambda_pois)

mse_ridge_zip <- mean((beta_ridge_pois - beta_true)^2)
mse_results <- rbind(mse_results, data.frame(Method = "RidgeZIP", MSE = mse_ridge_zip))

# Ridge Hurdle Poisson
lambda_log_h <- cv_ridge_logistic(X, y_pos)
lambda_pois_h <- cv_ridge_poisson(X[y_pos, ], Y[y_pos])

glm_logit_h <- glm((Y > 0) ~ X, family = binomial)
beta_logit_h_mle <- coef(glm_logit_h)[-1]
beta_ridge_logit_h <- ridge_projection(X, beta_logit_h_mle, lambda_log_h)

glm_pois_h <- glm(Y ~ X, family = poisson, subset = y_pos)
beta_pois_h_mle <- coef(glm_pois_h)[-1]
beta_ridge_pois_h <- ridge_projection(X[y_pos, ], beta_pois_h_mle, lambda_pois_h)

mse_ridge_hurdle <- mean((beta_ridge_pois_h - beta_true)^2)
mse_results <- rbind(mse_results, data.frame(Method = "RidgeHurdlePoisson", MSE = mse_ridge_hurdle))

# Show results
print(mse_results)

# Save
write.xlsx(mse_results, "C:/Users/hnayem/OneDrive - Florida International University/Documents/Private/Phd work/Folder for Simulation work/ZI VS Multico/Results/bioChemists_real_MSE_results_ridge_zi.xlsx")
