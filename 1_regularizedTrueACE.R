true_ace_est <- function(outcome_model_star_coef,
                         outcome_model_star_coef_t = 0,
                         covmodels_star_coef_A,
                         covmodels_star_coef,
                         K,
                         a = 1,
                         N = 1e6, seed = 2024) {
  true_haz <- c()
  set.seed(seed)

  L <- matrix(rbinom(N * nrow(covmodels_star_coef), 1, 0.5), nrow = N)
  A <- rep(a, N)
  Y <- plogis(outcome_model_star_coef[1] +
                cbind(A, L) %*% matrix(outcome_model_star_coef[-1]))
  true_haz <- c(true_haz, mean(Y))

  for (i in 1:(K - 1)) {
    L <- matrix(rbinom(N * nrow(covmodels_star_coef), 1,
                       plogis(L %*% covmodels_star_coef +
                                A %*% t(covmodels_star_coef_A))),
                nrow = N)
    A <- rep(a, N)
    Y <- plogis(outcome_model_star_coef_t * i +
                  outcome_model_star_coef[1] +
                  cbind(A, L) %*% matrix(outcome_model_star_coef[-1]))
    true_haz <- c(true_haz, mean(Y))
  }

  1 - cumprod(1 - true_haz)
}

true_risk1 <- true_ace_est(outcome_model_star_coef     = outcome_model_star_coef,
                           outcome_model_star_coef_t   = outcome_model_star_coef_t,
                           covmodels_star_coef_A       = covmodels_star_coef_A,
                           covmodels_star_coef         = covmodels_star_coef,
                           K = K, a = 1)
true_risk0 <- true_ace_est(outcome_model_star_coef     = outcome_model_star_coef,
                           outcome_model_star_coef_t   = outcome_model_star_coef_t,
                           covmodels_star_coef_A       = covmodels_star_coef_A,
                           covmodels_star_coef         = covmodels_star_coef,
                           K = K, a = 0)
