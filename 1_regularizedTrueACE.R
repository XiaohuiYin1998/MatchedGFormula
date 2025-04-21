true_ace_est <- function(outcome_model_star_coef, 
                         outcome_model_star_coef_t = 0,
                         covmodels_star_coef_A,
                         covmodels_star_coef,
                         K,
                         a = 1, 
                         N = 1e6, seed = 2024){
  true_haz <- c()
  set.seed(seed)
  ## t=0
  L <- matrix(rbinom(N * nrow(covmodels_star_coef), 1, 0.5), nrow = N)
  A <- rep(a,N)
  Y <- plogis(outcome_model_star_coef[1] + cbind(A,L) %*% matrix(outcome_model_star_coef[-1] ))
  true_haz <- c(true_haz, mean(Y))
  
  ## t=1 to K-1
  for(i in 1:(K-1)){
    L <- matrix(rbinom(N * nrow(covmodels_star_coef), 1, 
                       plogis(L %*% covmodels_star_coef +  A %*% t(covmodels_star_coef_A))), nrow = N)
    A <- rep(a,N)
    Y <- plogis(outcome_model_star_coef_t * i + outcome_model_star_coef[1] + cbind(A,L) %*% matrix(outcome_model_star_coef[-1] ))
    true_haz <- c(true_haz, mean(Y))
  }
  
  true_risks <- 1-cumprod(1-true_haz)
  return(true_risks)
}
true_risk1 <- true_ace_est(outcome_model_star_coef = outcome_model_star_coef, 
                           outcome_model_star_coef_t = outcome_model_star_coef_t, 
                           covmodels_star_coef_A = covmodels_star_coef_A, 
                           covmodels_star_coef = covmodels_star_coef, 
                           K=K, a = 1)
true_risk0 <- true_ace_est(outcome_model_star_coef = outcome_model_star_coef, 
                           outcome_model_star_coef_t = outcome_model_star_coef_t, 
                           covmodels_star_coef_A = covmodels_star_coef_A, 
                           covmodels_star_coef = covmodels_star_coef, 
                           K=K, a = 0)

################################################################################
###### EXAMPLE
################################################################################
# true_haz <- c()
# true_Y_prob = c() ## another calculation
# N <- 1e6
# set.seed(2024)
# 
# ## t=0
# L <- matrix(rbinom(N * nrow(covmodels_star_coef), 1, 0.5), nrow = N)
# A <- rep(1,N)
# Y <- plogis(outcome_model_star_coef[1] + cbind(A,L) %*% matrix(outcome_model_star_coef[-1] ))
# true_haz <- c(true_haz, mean(Y))
# true_Y_prob = cbind(true_Y_prob, Y)
# 
# ## t=1 to K-1
# for(i in 1:(K-1)){
#   L <- matrix(rbinom(N * nrow(covmodels_star_coef), 1, 
#                      plogis(L %*% covmodels_star_coef +  A %*% t(covmodels_star_coef_A))), nrow = N)
#   A <- rep(1,N)
#   Y <- plogis(outcome_model_star_coef[1] + cbind(A,L) %*% matrix(outcome_model_star_coef[-1] ))
#   true_haz <- c(true_haz, mean(Y))
#   true_Y_prob = cbind(true_Y_prob, Y)
# }
# 
# true_Y_prob = cbind(true_Y_prob[,1], 
#                     1-(1-true_Y_prob[,1])*(1-true_Y_prob[,2]), 
#                     1-(1-true_Y_prob[,1])*(1-true_Y_prob[,2])*(1-true_Y_prob[,3]),
#                     1-(1-true_Y_prob[,1])*(1-true_Y_prob[,2])*(1-true_Y_prob[,3])*(1-true_Y_prob[,4]),
#                     1-(1-true_Y_prob[,1])*(1-true_Y_prob[,2])*(1-true_Y_prob[,3])*(1-true_Y_prob[,4])*(1-true_Y_prob[,5]),
#                     1-(1-true_Y_prob[,1])*(1-true_Y_prob[,2])*(1-true_Y_prob[,3])*(1-true_Y_prob[,4])*(1-true_Y_prob[,5])*(1-true_Y_prob[,6]))
# 
# colMeans(true_Y_prob)
# 
# true_risks <- c(true_haz[1],
#                 1-(1-true_haz[1])*(1-true_haz[2]),
#                 1-(1-true_haz[1])*(1-true_haz[2])*(1-true_haz[3]),
#                 1-(1-true_haz[1])*(1-true_haz[2])*(1-true_haz[3])*(1-true_haz[4]),
#                 1-(1-true_haz[1])*(1-true_haz[2])*(1-true_haz[3])*(1-true_haz[4])*(1-true_haz[5]),
#                 1-(1-true_haz[1])*(1-true_haz[2])*(1-true_haz[3])*(1-true_haz[4])*(1-true_haz[5])*(1-true_haz[6]))
# 
# 
# true_haz <- c()
# true_Y_prob0 = c() ## another calculation
# N <- 1e6
# set.seed(2024)
# 
# ## t=0
# L <- matrix(rbinom(N * nrow(covmodels_star_coef), 1, 0.5), nrow = N)
# A <- rep(0,N)
# Y <- plogis(outcome_model_star_coef[1] + cbind(A,L) %*% matrix(outcome_model_star_coef[-1] ))
# true_haz <- c(true_haz, mean(Y))
# true_Y_prob0 = cbind(true_Y_prob0, Y)
# 
# ## t=1 to K-1
# for(i in 1:(K-1)){
#   L <- matrix(rbinom(N * nrow(covmodels_star_coef), 1, 
#                      plogis(L %*% covmodels_star_coef + A %*% t(covmodels_star_coef_A))), nrow = N)
#   A <- rep(0,N)
#   Y <- plogis(outcome_model_star_coef[1] + cbind(A,L) %*% matrix(outcome_model_star_coef[-1] ))
#   true_haz <- c(true_haz, mean(Y))
#   true_Y_prob0 = cbind(true_Y_prob0, Y)
# }
# 
# true_Y_prob0 = cbind(true_Y_prob0[,1], 
#                      1-(1-true_Y_prob0[,1])*(1-true_Y_prob0[,2]), 
#                      1-(1-true_Y_prob0[,1])*(1-true_Y_prob0[,2])*(1-true_Y_prob0[,3]),
#                      1-(1-true_Y_prob0[,1])*(1-true_Y_prob0[,2])*(1-true_Y_prob0[,3])*(1-true_Y_prob0[,4]),
#                      1-(1-true_Y_prob0[,1])*(1-true_Y_prob0[,2])*(1-true_Y_prob0[,3])*(1-true_Y_prob0[,4])*(1-true_Y_prob0[,5]),
#                      1-(1-true_Y_prob0[,1])*(1-true_Y_prob0[,2])*(1-true_Y_prob0[,3])*(1-true_Y_prob0[,4])*(1-true_Y_prob0[,5])*(1-true_Y_prob0[,6]))
# 
# colMeans(true_Y_prob0)
# 
# true_risks_notreat <- c(true_haz[1],
#                         1-(1-true_haz[1])*(1-true_haz[2]),
#                         1-(1-true_haz[1])*(1-true_haz[2])*(1-true_haz[3]),
#                         1-(1-true_haz[1])*(1-true_haz[2])*(1-true_haz[3])*(1-true_haz[4]),
#                         1-(1-true_haz[1])*(1-true_haz[2])*(1-true_haz[3])*(1-true_haz[4])*(1-true_haz[5]),
#                         1-(1-true_haz[1])*(1-true_haz[2])*(1-true_haz[3])*(1-true_haz[4])*(1-true_haz[5])*(1-true_haz[6]))
# 
# 
