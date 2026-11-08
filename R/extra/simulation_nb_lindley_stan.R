# Clear environment
rm(list = ls())

# Load libraries
library(MASS)
library(pscl)
library(glmnet)
library(openxlsx)
library(rstan)
library(dplyr)

# Stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Data generation functions
generate_X <- function(n, p, rho) {
  Z <- matrix(rnorm(n * p), nrow = n)
  X <- sqrt(1 - rho) * Z + rho * matrix(rnorm(n), nrow = n, ncol = p)
  return(scale(X))
}

generate_Y_zip <- function(X, beta, zero_prob) {
  eta <- X %*% beta
  lambda <- exp(pmin(eta, 13.8))  # exp(13.8) ≈ 1e6
  
  is_zero <- rbinom(nrow(X), 1, zero_prob)
  Y <- ifelse(is_zero == 1, 0, rpois(nrow(X), lambda))
  Y <- as.integer(Y)
  stopifnot(all(Y >= 0), all(Y == floor(Y)), all(!is.na(Y)))
  return(Y)
}

# Stan model code
stan_code <- "
functions {
  real lindley_mixture_lpdf(real epsilon, real theta) {
    return log_mix(
      1 / (1 + theta),
      gamma_lpdf(epsilon | 1, theta),
      gamma_lpdf(epsilon | 2, theta)
    );
  }
}
data {
  int<lower=0> N;
  int<lower=0> K;
  matrix[N, K] X;
  int<lower=0> Y[N];
}
parameters {
  real alpha;
  vector[K] beta;
  real<lower=0> theta;
  real<lower=0> phi;
  vector<lower=0>[N] epsilon;
}
model {
  vector[N] eta = X * beta + alpha;
  vector[N] mu = exp(eta) .* epsilon;

  alpha ~ normal(0, 5);
  beta ~ normal(0, 5);
  theta ~ gamma(2, 1);
  phi ~ gamma(2, 0.5);

  for (n in 1:N) {
    target += lindley_mixture_lpdf(epsilon[n] | theta);
    Y[n] ~ neg_binomial_2(mu[n], phi);
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] mu;
  for (n in 1:N) {
    mu[n] = exp(dot_product(X[n], beta) + alpha) * epsilon[n];
    log_lik[n] = neg_binomial_2_lpmf(Y[n] | mu[n], phi);
  }
}
"
writeLines(stan_code, con = "lindley_nb_latent.stan")

# Compile Stan model
stan_model <- stan_model("lindley_nb_latent.stan")

# Parameter grid
n_vals <- c(100, 200, 500)
p_vals <- c(20, 30)
rho_vals <- c(0.80, 0.90, 0.95)
zero_probs <- c(0.3,0.50,.80)
nsim <- 1  # simulations per setting

results <- data.frame()

for (zero_prob in zero_probs) {
  for (n in n_vals) {
    for (p in p_vals) {
      for (rho in rho_vals) {
        cat("Running: n =", n, ", p =", p, ", rho =", rho, ", zero_prob =", zero_prob, "\n")
        
        mse_mat <- matrix(0, nrow = nsim, ncol = 6)
        colnames(mse_mat) <- c("ZIP", "ZINB", "NBLindley", "Ridge", "Lasso", "ElasticNet")
        
        for (sim in 1:nsim) {
          X <- generate_X(n, p, rho)
          XtX <- t(X) %*% X
          eigen_decomp <- eigen(XtX)
          eigenvectors <- eigen_decomp$vectors
          beta_true <- eigenvectors[, 1]  # Use first principal direction
          
          Y <- generate_Y_zip(X, beta_true, zero_prob)
          

          
          # ZIP
          zip_fit <- try(zeroinfl(Y ~ X - 1 | 1, dist = "poisson"), silent = TRUE)
          if (inherits(zip_fit, "try-error") || any(is.na(coef(zip_fit)))) {
            mse_mat[sim, "ZIP"] <- NA
          } else {
            zip_beta <- coef(zip_fit)[2:(p + 1)]
            mse_mat[sim, "ZIP"] <- mean((zip_beta - beta_true)^2)
          }
          
          # ZINB
          zinb_fit <- try(zeroinfl(Y ~ X - 1 | 1, dist = "negbin"), silent = TRUE)
          if (inherits(zinb_fit, "try-error") || any(is.na(coef(zinb_fit)))) {
            mse_mat[sim, "ZINB"] <- NA
          } else {
            zinb_beta <- coef(zinb_fit)[2:(p + 1)]
            mse_mat[sim, "ZINB"] <- mean((zinb_beta - beta_true)^2)
          }
          
          
          # Ridge
          cv_ridge <- cv.glmnet(X, Y, alpha = 0, family = "poisson")
          ridge_beta <- as.numeric(coef(cv_ridge))[-1]
          mse_mat[sim, "Ridge"] <- mean((ridge_beta - beta_true)^2)
          
          # Lasso
          cv_lasso <- cv.glmnet(X, Y, alpha = 1, family = "poisson")
          lasso_beta <- as.numeric(coef(cv_lasso))[-1]
          mse_mat[sim, "Lasso"] <- mean((lasso_beta - beta_true)^2)
          
          # Elastic Net
          cv_enet <- cv.glmnet(X, Y, alpha = 0.5, family = "poisson")
          enet_beta <- as.numeric(coef(cv_enet))[-1]
          mse_mat[sim, "ElasticNet"] <- mean((enet_beta - beta_true)^2)
          
          # NB Lindley (Stan)
          stan_data <- list(N = n, K = p, X = X, Y = Y)
          fit <- sampling(stan_model, data = stan_data,
                          iter = 10000, warmup = 500, chains = 2,
                          seed = 123,
                          control = list(max_treedepth = 12, adapt_delta = 0.95))
          
          # Extract posterior samples for beta
          posterior_draws <- extract(fit, pars = "beta")$beta
          
          # Compute posterior mean (or median) of beta
          beta_hat_nbl <- colMeans(posterior_draws)  # or use apply(posterior_draws, 2, median)
          
          # Compute MSE between estimated beta and true beta
          mse_mat[sim, "NBLindley"] <- mean((beta_hat_nbl - beta_true)^2)
          
        }
        
        mean_mses <- colMeans(mse_mat)
        results <- rbind(results, data.frame(
          ZeroProb = zero_prob,
          n = n, p = p, rho = rho,
          ZIP = mean_mses["ZIP"],
          ZINB = mean_mses["ZINB"],
          NBLindley = mean_mses["NBLindley"],
          Ridge = mean_mses["Ridge"],
          Lasso = mean_mses["Lasso"],
          ElasticNet = mean_mses["ElasticNet"]
        ))
        
      }
    }
  }
}

# Save results
# write.xlsx(results, "C:\\Users\\hnayem\\OneDrive - Florida International University\\Documents\\Private\\Phd work\\Folder for Simulation work\\ZI VS Multico\\Results\\simulation_zero_inflated_results_all_probs.xlsx")
print(results)



write.xlsx(results, "C:\\Users\\hnayem\\OneDrive - Florida International University\\Documents\\Private\\Phd work\\Folder for Simulation work\\ZI VS Multico\\Results\\ziro_probs_52925_final.xlsx")
