library(data.table)
outcome_name <- 'Y'
censor_name <- 'C'
intervention_name = 'combined_social'
mhds <- c('alch', 'drug', 'sch', 'mdd', 'bipd', 'ptsd', 'anx')
cmbds <- c('mi', 'chf', 'pvd', 'cevd', 'dementia', 'cpd', 'rheumd', 
           'pud', 'mld', 'diab', 'diabwc', 'hp', 'rend', 'canc', 
           'msld', 'metacanc', 'aids')
covnames = c(mhds, cmbds)[1:10] # set 10 covariates

# ## Sparse setups
# p <- length(covnames); p0 <- 3 ## set the first 3 covariates as active covariates
# set.seed(2024)
# covmodels_star_coef_A <- rnorm(p)
# outcome_model_star_coef_t <- 0
# covmodels_star_coef <- matrix(0, nrow = p, ncol = p)
# covmodels_star_coef[1:p0, 1:p0] <- rnorm(p0*p0)
# intervention_model_star_coef <- c(0, 1, c(rnorm(p0), rep(0, p - p0)))

## General setups
p <- length(covnames)
set.seed(2024)
covmodels_star_coef_A <- rnorm(p)
outcome_model_star_coef_t <- 0
covmodels_star_coef <- matrix(rnorm(p*p), nrow = p, ncol = p)
intervention_model_star_coef <- c(0, 1, rnorm(p))

## Setting 1
# if(setup == 1){
#   censor_model_star_coef <- c(-1, 1, c(rnorm(p0), rep(0, p - p0)))
#   outcome_model_star_coef <- c(-6, 2, c(rnorm(p0), rep(0, p - p0)))
# }
# # Setting 2
# setup = 2
# if(setup == 2){
#   # censor_model_star_coef <- c(-6, 1, c(rnorm(p0), rep(0, p - p0)))
#   # outcome_model_star_coef <- c(-8, 2, c(rnorm(p0), rep(0, p - p0)))
#   censor_model_star_coef <- c(-6, 1, rnorm(p))
#   outcome_model_star_coef <- c(-8, 2, rnorm(p))
# }
# # ## Setting 3
# setup = 3
# if(setup == 3){
#   # censor_model_star_coef <- c(-6, 1, c(rnorm(p0), rep(0, p - p0)))
#   # outcome_model_star_coef <- c(-4, 2, c(rnorm(p0), rep(0, p - p0)))
#   censor_model_star_coef <- c(-6, 1, rnorm(p))
#   outcome_model_star_coef <- c(-4, 2, rnorm(p))
# }
# ## Setting 4
# if(setup == 4){
#   censor_model_star_coef <- c(-1, 1, c(rnorm(p0), rep(0, p - p0)))
#   outcome_model_star_coef <- c(-8, 2, c(rnorm(p0), rep(0, p - p0)))
# }
# Setting 5
setup = 5
if(setup == 5){
  censor_model_star_coef <- c(-6, 1, rnorm(p))
  outcome_model_star_coef <- c(-6, 2, rnorm(p))
}

K <- 6
va_simu_data <- function(i, K, 
                         covmodels_star_coef_A,
                         covmodels_star_coef,
                         censor_model_star_coef,
                         outcome_model_star_coef_t = 0,
                         outcome_model_star_coef){
  id <- as.numeric(i)
  p <- nrow(covmodels_star_coef)
  # Data at time 0
  Lbase <- c(sample(x = c('F', 'M'), size = 1, prob = rep(1, 2), replace = TRUE),
             sample(c('30_39', '40_49', '50_59', '60_69', '70_79', '79_'), size = 1, prob = rep(1,6), replace = TRUE),
             sample(c('Asian', 'Black', 'Hawaiian', 'Native_Indians', 'Unknown'), size = 1, prob = c(5,2,1,1,4), replace = TRUE),
             sample(c('Single', 'Married', 'Divorced', 'Unknow'), size = 1, prob = c(2,3,1,4), replace = TRUE))
  names(Lbase) <- c('Gender', 'Age', 'Race', 'Marital')
  L <- rbinom(nrow(covmodels_star_coef), 1, 0.5)
  A <- rbinom(1, 1, plogis(intervention_model_star_coef[1] + sum(intervention_model_star_coef[-(1:2)] * L)))
  C <- rbinom(1, 1, plogis(sum(censor_model_star_coef[1] + sum(censor_model_star_coef[-1] * c(A, L)))))
  Y <- rbinom(1, 1, plogis(outcome_model_star_coef[1] + sum(outcome_model_star_coef[-1] * c(A, L))))
  temp <- data.frame(id = id, t0 = 0, t(L), A, C, Y)
  
  if (temp$C==1){
    temp$Y <- NA
  }else if (temp$Y != 1){
    for (j in 2:K){
      t0 <- j-1
      ### , change to + ?
      Lstar <- rbinom(nrow(covmodels_star_coef), 1, 
                      plogis(covmodels_star_coef_A * temp$A[t0]+ 
                               covmodels_star_coef %*% matrix(as.numeric(temp[t0, paste0('X', 1:p)]))))
      ### change temp[t0, paste0('X', 1:p)] to Lstar
      Astar <- rbinom(1, 1, plogis(intervention_model_star_coef[1] + intervention_model_star_coef[2] * temp$A[t0] + sum(intervention_model_star_coef[-(1:2)] * Lstar)))
      ### change temp[t0, paste0('X', 1:p)] to Lstar
      Cstar <- rbinom(1, 1, plogis(censor_model_star_coef[1] + censor_model_star_coef[2] * Astar + sum(censor_model_star_coef[-(1:2)] * Lstar)))
      if (Cstar == 1){
        Ystar <- NA
        temp <- rbind(temp, c(id, t0, Lstar, Astar, Cstar, Ystar))
        break
      }
      else{
        Ystar <- rbinom(1, 1, plogis(
          outcome_model_star_coef[1] + 
            outcome_model_star_coef_t * t0 + 
            outcome_model_star_coef[2] * Astar + sum(outcome_model_star_coef[-(1:2)] * temp[t0, paste0('X', 1:p)])))
      }
      temp <- rbind(temp, c(id, t0, Lstar, Astar, Cstar, Ystar))
      if(Ystar == 1){break}
    }
  }
  data.frame('id' = temp$id, 't0' = temp$t0, t(Lbase), temp[, -(1:2)])
}