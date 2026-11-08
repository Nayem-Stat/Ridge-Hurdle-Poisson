# Clear environment
rm(list = ls())

# Load libraries
library(MASS)
library(pscl) # For standard ZIP and ZINB
library(openxlsx)
library(brms) # For Ridge and Hurdle models
library(dplyr)
library(glmnet) # Although not directly used in the brms section, kept as it was in original script

# Stan options
options(mc.cores = parallel::detectCores())

# Data generation functions
generate_X <- function(n, p, rho) {
  Z <- matrix(rnorm(n * p), nrow = n)
  X <- sqrt(1 - rho^2) * Z + rho * matrix(rnorm(n), nrow = n, ncol = p)
  return(scale(X))
}

generate_Y_zip <- function(X, beta, zero_prob) {
  eta <- X %*% beta
  lambda <- exp(pmin(eta, 13.8))  # cap to avoid overflow
  is_zero <- rbinom(nrow(X), 1, zero_prob)
  Y <- ifelse(is_zero == 1, 0, rpois(nrow(X), lambda))
  Y <- as.integer(Y)
  stopifnot(all(Y >= 0), all(Y == floor(Y)), all(!is.na(Y)))
  return(Y)
}

# Simulation settings
n_vals <- c(100, 200)
p_vals <- c(5, 10)
rho_vals <- c(0.95, 0.99)
zero_probs <- c(0.3, 0.5, 0.8)
nsim <- 100 # Increased for more robust results, consider 100 or more for real simulations

results <- data.frame()

# Additional helper to check convergence
brms_converged <- function(fit) {
  summ <- try(summary(fit), silent = TRUE)
  if (inherits(summ, "try-error")) return(FALSE)
  rhat_vals <- summ$fixed[, "Rhat"]
  ess_vals <- summ$fixed[, "Bulk_ESS"]
  return(all(rhat_vals < 1.1, na.rm = TRUE) && all(ess_vals > 100, na.rm = TRUE))
}

for (zero_prob in zero_probs) {
  for (n in n_vals) {
    for (p in p_vals) {
      for (rho in rho_vals) {
        cat("Running: n =", n, ", p =", p, ", rho =", rho, ", zero_prob =", zero_prob, "\n")
        
        mse_mat <- matrix(NA, nrow = nsim, ncol = 8)
        colnames(mse_mat) <- c("ZIP", "ZINB", "RidgeZIP", "RidgeZINB",
                               "HurdlePoisson", "HurdleNB", "RidgeHurdlePoisson", "RidgeHurdleNB")
        
        for (sim in 1:nsim) {
          X <- generate_X(n, p, rho)
          XtX <- t(X) %*% X
          beta_true <- eigen(XtX)$vectors[, 1]
          Y <- generate_Y_zip(X, beta_true, zero_prob)
          data_brms <- data.frame(Y = Y, X = X, ID = 1:n)
          
          # ZIP
          zip_fit <- try(zeroinfl(Y ~ X - 1 | 1, data = data.frame(Y = Y, X = X), dist = "poisson"), silent = TRUE)
          if (!inherits(zip_fit, "try-error")) {
            zip_beta <- coef(zip_fit)[2:(p + 1)]
            if (!any(is.na(zip_beta))) mse_mat[sim, "ZIP"] <- mean((zip_beta - beta_true)^2)
          }
          
          # ZINB
          zinb_fit <- try(zeroinfl(Y ~ X - 1 | 1, data = data.frame(Y = Y, X = X), dist = "negbin"), silent = TRUE)
          if (!inherits(zinb_fit, "try-error")) {
            zinb_beta <- coef(zinb_fit)[2:(p + 1)]
            if (!any(is.na(zinb_beta))) mse_mat[sim, "ZINB"] <- mean((zinb_beta - beta_true)^2)
          }
          
          model_list <- list(
            RidgeZIP = list(family = zero_inflated_poisson(), prior = prior(normal(0, 0.1), class = "b")),
            RidgeZINB = list(family = zero_inflated_negbinomial(), prior = prior(normal(0, 0.1), class = "b")),
            HurdlePoisson = list(family = hurdle_poisson(), prior = NULL),
            HurdleNB = list(family = hurdle_negbinomial(), prior = NULL),
            RidgeHurdlePoisson = list(family = hurdle_poisson(), prior = prior(normal(0, 0.1), class = "b")),
            RidgeHurdleNB = list(family = hurdle_negbinomial(), prior = prior(normal(0, 0.1), class = "b"))
          )
          
          for (model_name in names(model_list)) {
            family <- model_list[[model_name]]$family
            prior_ <- model_list[[model_name]]$prior
            fit <- try(
              brm(
                formula = bf(Y ~ 0 + X),
                data = data_brms,
                family = family,
                prior = prior_,
                chains = 2,
                iter = 4000,
                warmup = 2000,
                refresh = 0,
                control = list(adapt_delta = 0.999, max_treedepth = 15),
                backend = "cmdstanr",  # Use cmdstanr backend
                threads = threading(2),  # Optional speed-up
                silent = FALSE
              ),
              silent = TRUE
            )
            
        
            
            if (!inherits(fit, "try-error") && brms_converged(fit)) {
              beta_draws <- as_draws_matrix(fit, variable = "b_X")
              if (!anyNA(beta_draws)) {
                beta_est <- colMeans(beta_draws)
                mse_mat[sim, model_name] <- mean((beta_est - beta_true)^2)
              }
            } else {
              cat("Warning: Model", model_name, "failed or didn't converge at n =", n, ", p =", p, ", rho =", rho, "\n")
            }
          }
        }
        
        # Aggregate results
        mean_mses <- colMeans(mse_mat, na.rm = TRUE)
        results <- rbind(results, data.frame(
          ZeroProb = zero_prob,
          n = n, p = p, rho = rho,
          ZIP = mean_mses["ZIP"],
          ZINB = mean_mses["ZINB"],
          RidgeZIP = mean_mses["RidgeZIP"],
          RidgeZINB = mean_mses["RidgeZINB"],
          HurdlePoisson = mean_mses["HurdlePoisson"],
          HurdleNB = mean_mses["HurdleNB"],
          RidgeHurdlePoisson = mean_mses["RidgeHurdlePoisson"],
          RidgeHurdleNB = mean_mses["RidgeHurdleNB"]
        ))
      }
    }
  }
}

# Save results
write.xlsx(results, "C:\\Users\\hnayem\\OneDrive - Florida International University\\Documents\\Private\\Phd work\\Folder for Simulation work\\ZI VS Multico\\Results\\zi_sim_with_ridge_zip_zinb_hurdle_cmp1234.xlsx")
