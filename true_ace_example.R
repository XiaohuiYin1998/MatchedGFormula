p <- length(covnames)
set.seed(2024)
covmodels_star_coef_A <- rnorm(p)
covmodels_star_coef <- matrix(rnorm(p * p), nrow = p)
covmodels_star_coef <- matrix(0, nrow = p, ncol = p)
intervention_model_star_coef <- c(0, 1, rnorm(p))
censor_model_star_coef <- c(-3, 1, rnorm(p))
outcome_model_star_coef <- c(-6, 2, rnorm(p))


## -6; -4; -3; 
true_haz <- c()
N <- 1e6
set.seed(2024)
## t=0
L <- matrix(rbinom(N * nrow(covmodels_star_coef), 1, 0.5), nrow = N)
# A <- rbinom(N, 1, plogis(intervention_model_star_coef[1] + L %*% matrix(intervention_model_star_coef[-(1:2)])))
A <- rep(1,N)
Y <- plogis(outcome_model_star_coef[1] + cbind(A,L) %*% matrix(outcome_model_star_coef[-1] ))
true_haz <- c(true_haz, mean(Y))
### another calculation
true_Y_prob = c()
true_Y_prob = cbind(true_Y_prob, Y)

## t=1
for(i in 1:(K-1)){
  L <- matrix(rbinom(N * nrow(covmodels_star_coef), 1, 
                     plogis(L %*% covmodels_star_coef +  A %*% t(covmodels_star_coef_A))), nrow = N)
  # A <- rbinom(N, 1, plogis(intervention_model_star_coef[1] + intervention_model_star_coef[2] * A + L %*% matrix(intervention_model_star_coef[-(1:2)])))
  A <- rep(1,N)
  Y <- plogis(outcome_model_star_coef[1] + cbind(A,L) %*% matrix(outcome_model_star_coef[-1] ))
  true_haz <- c(true_haz, mean(Y))
  ### add the following
  true_Y_prob = cbind(true_Y_prob, Y)
}

true_Y_prob = cbind(true_Y_prob[,1], 
                    1-(1-true_Y_prob[,1])*(1-true_Y_prob[,2]), 
                    1-(1-true_Y_prob[,1])*(1-true_Y_prob[,2])*(1-true_Y_prob[,3]),
                    1-(1-true_Y_prob[,1])*(1-true_Y_prob[,2])*(1-true_Y_prob[,3])*(1-true_Y_prob[,4]),
                    1-(1-true_Y_prob[,1])*(1-true_Y_prob[,2])*(1-true_Y_prob[,3])*(1-true_Y_prob[,4])*(1-true_Y_prob[,5]),
                    1-(1-true_Y_prob[,1])*(1-true_Y_prob[,2])*(1-true_Y_prob[,3])*(1-true_Y_prob[,4])*(1-true_Y_prob[,5])*(1-true_Y_prob[,6]))

colMeans(true_Y_prob)

true_risks <- c(true_haz[1],
                1-(1-true_haz[1])*(1-true_haz[2]),
                1-(1-true_haz[1])*(1-true_haz[2])*(1-true_haz[3]),
                1-(1-true_haz[1])*(1-true_haz[2])*(1-true_haz[3])*(1-true_haz[4]),
                1-(1-true_haz[1])*(1-true_haz[2])*(1-true_haz[3])*(1-true_haz[4])*(1-true_haz[5]),
                1-(1-true_haz[1])*(1-true_haz[2])*(1-true_haz[3])*(1-true_haz[4])*(1-true_haz[5])*(1-true_haz[6]))


true_haz <- c()
N <- 1e6
set.seed(2024)
## t=0
L <- matrix(rbinom(N * nrow(covmodels_star_coef), 1, 0.5), nrow = N)
# A <- rbinom(N, 1, plogis(intervention_model_star_coef[1] + L %*% matrix(intervention_model_star_coef[-(1:2)])))
A <- rep(0,N)
Y <- plogis(outcome_model_star_coef[1] + cbind(A,L) %*% matrix(outcome_model_star_coef[-1] ))
true_haz <- c(true_haz, mean(Y))
### another calculation
true_Y_prob0 = c()
true_Y_prob0 = cbind(true_Y_prob0, Y)

## t=1
for(i in 1:(K-1)){
  L <- matrix(rbinom(N * nrow(covmodels_star_coef), 1, 
                     plogis(L %*% covmodels_star_coef + A %*% t(covmodels_star_coef_A))), nrow = N)
  # A <- rbinom(N, 1, plogis(intervention_model_star_coef[1] + intervention_model_star_coef[2] * A + L %*% matrix(intervention_model_star_coef[-(1:2)])))
  A <- rep(0,N)
  Y <- plogis(outcome_model_star_coef[1] + cbind(A,L) %*% matrix(outcome_model_star_coef[-1] ))
  true_haz <- c(true_haz, mean(Y))
  ### add the following
  true_Y_prob0 = cbind(true_Y_prob0, Y)
}

true_Y_prob0 = cbind(true_Y_prob0[,1], 
                     1-(1-true_Y_prob0[,1])*(1-true_Y_prob0[,2]), 
                     1-(1-true_Y_prob0[,1])*(1-true_Y_prob0[,2])*(1-true_Y_prob0[,3]),
                     1-(1-true_Y_prob0[,1])*(1-true_Y_prob0[,2])*(1-true_Y_prob0[,3])*(1-true_Y_prob0[,4]),
                     1-(1-true_Y_prob0[,1])*(1-true_Y_prob0[,2])*(1-true_Y_prob0[,3])*(1-true_Y_prob0[,4])*(1-true_Y_prob0[,5]),
                     1-(1-true_Y_prob0[,1])*(1-true_Y_prob0[,2])*(1-true_Y_prob0[,3])*(1-true_Y_prob0[,4])*(1-true_Y_prob0[,5])*(1-true_Y_prob0[,6]))

colMeans(true_Y_prob0)

true_risks_notreat <- c(true_haz[1],
                        1-(1-true_haz[1])*(1-true_haz[2]),
                        1-(1-true_haz[1])*(1-true_haz[2])*(1-true_haz[3]),
                        1-(1-true_haz[1])*(1-true_haz[2])*(1-true_haz[3])*(1-true_haz[4]),
                        1-(1-true_haz[1])*(1-true_haz[2])*(1-true_haz[3])*(1-true_haz[4])*(1-true_haz[5]),
                        1-(1-true_haz[1])*(1-true_haz[2])*(1-true_haz[3])*(1-true_haz[4])*(1-true_haz[5])*(1-true_haz[6]))


# rst <- lapply(1:100, function(j){
#   readRDS(file = sprintf('./results/repli_%d.RDS', j))
# })
# df_treat <- lapply(1:100, function(j){
#   rst[[j]]$df_treat
# })
# df_treat <- do.call('rbind', df_treat)
# df_notreat <- lapply(1:100, function(j){
#   rst[[j]]$df_notreat
# })
# df_notreat <- do.call('rbind', df_notreat)
# time_treat <- lapply(1:100, function(j){
#   rst[[j]]$time_treat
# })
# time_treat <- do.call('rbind', time_treat)
# time_notreat <- lapply(1:100, function(j){
#   rst[[j]]$time_notreat
# })
# time_notreat <- do.call('rbind', time_notreat)
# 
# saveRDS(list(df_treat = df_treat, df_notreat = df_notreat,
#      time_treat = time_treat, time_notreat = time_notreat), file = 'gform_sim.RDS')
# 
# 
# 
# 

true_risk_calc <- function(N=1e6, K=6, a=1, 
                           covmodels_star_coef_A, 
                           intervention_model_star_coef, 
                           censor_model_star_coef, 
                           outcome_model_star_coef){
  set.seed(2024)
  true_haz <- c()
  N <- 1e6
  set.seed(2024)
  ## t=0
  L <- matrix(rbinom(N * nrow(covmodels_star_coef), 1, 0.5), nrow = N)
  A <- rep(a,N)
  Y <- plogis(outcome_model_star_coef[1] + cbind(A,L) %*% matrix(outcome_model_star_coef[-1] ))
  true_haz <- c(true_haz, mean(Y))
  
  ## t=1
  for(i in 1:(K-1)){
    L <- matrix(rbinom(N * nrow(covmodels_star_coef), 1, 
                       plogis(L %*% covmodels_star_coef +  A * covmodels_star_coef_A)), nrow = N)
    # A <- rbinom(N, 1, plogis(intervention_model_star_coef[1] + intervention_model_star_coef[2] * A + L %*% matrix(intervention_model_star_coef[-(1:2)])))
    A <- rep(a,N)
    Y <- plogis(outcome_model_star_coef[1] + cbind(A,L) %*% matrix(outcome_model_star_coef[-1] ))
    true_haz <- c(true_haz, mean(Y))
  }
  
  
  true_risks <- c(true_haz[1],
                  1-(1-true_haz[1])*(1-true_haz[2]),
                  1-(1-true_haz[1])*(1-true_haz[2])*(1-true_haz[3]),
                  1-(1-true_haz[1])*(1-true_haz[2])*(1-true_haz[3])*(1-true_haz[4]),
                  1-(1-true_haz[1])*(1-true_haz[2])*(1-true_haz[3])*(1-true_haz[4])*(1-true_haz[5]),
                  1-(1-true_haz[1])*(1-true_haz[2])*(1-true_haz[3])*(1-true_haz[4])*(1-true_haz[5])*(1-true_haz[6]))
  true_risks
}
