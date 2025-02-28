source('gformula_revised_latest.R')
library(data.table)
outcome_name <- 'Y'
censor_name <- 'C'
intervention_name = 'combined_social'
mhds <- c('alch', 'drug', 'sch', 'mdd', 'bipd', 'ptsd', 'anx')
cmbds <- c('mi', 'chf', 'pvd', 'cevd', 'dementia', 'cpd', 'rheumd', 
           'pud', 'mld', 'diab', 'diabwc', 'hp', 'rend', 'canc', 
           'msld', 'metacanc', 'aids')
covnames = c(mhds, cmbds)
p <- length(covnames)
set.seed(2024)
covmodels_star_coef_A <- rnorm(p)
covmodels_star_coef <- matrix(rnorm(p * p), nrow = p)
### difference
covmodels_star_coef <- matrix(0, nrow = p, ncol = p)
intervention_model_star_coef <- c(0, 1, rnorm(p))
censor_model_star_coef <- c(-1, 1, rnorm(p))
outcome_model_star_coef <- c(-6, 2, rnorm(p))

va_simu_data <- function(i, K, 
                         covmodels_star_coef,
                         censor_model_star_coef,
                         outcome_model_star_coef){
  id <- as.numeric(i)
  p <- nrow(covmodels_star_coef)
  # Data at time 0
  Lbase <- c(sample(x = c('F', 'M'), size = 1, prob = rep(1, 2), replace = TRUE),
             sample(c('30-39', '40-49', '50-59', '60-69', '70-79', '>79'), size = 1, prob = rep(1,6), replace = TRUE),
             sample(c('Asian', 'Black', 'Hawaiian', 'Native Indians', 'Unknown'), size = 1, prob = c(5,2,1,1,4), replace = TRUE),
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
        Ystar <- rbinom(1, 1, plogis(outcome_model_star_coef[1] + outcome_model_star_coef[2] * Astar + sum(outcome_model_star_coef[-(1:2)] * temp[t0, paste0('X', 1:p)])))
      }
      temp <- rbind(temp, c(id, t0, Lstar, Astar, Cstar, Ystar))
      if(Ystar == 1){break}
    }
  }
  temp
  data.frame('id'=temp$id, 't0'=temp$t0, t(Lbase), temp[, -(1:2)])
}

set.seed(2000)
myseeds <- sample(1:1e5, size = 100)
simuID <- Sys.getenv("SLURM_ARRAY_TASK_ID")
simuID <- as.integer(simuID)
if(is.na(simuID)){simuID <- 1}
myseed <- myseeds[simuID]
set.seed(myseed)

n <- 10000; K <- 6
df <- lapply(as.list(1:n), FUN=function(ind){
  va_simu_data(ind, K=K,covmodels_star_coef,
               censor_model_star_coef,
               outcome_model_star_coef)
})
dffull <- do.call('rbind', df)
colnames(dffull)[grep('X', colnames(dffull))] <- covnames
colnames(dffull)[colnames(dffull) == 'A'] <- intervention_name

dffull <- as.data.table(dffull)


time_name = 't0'
id_name = 'id'
base_covnames <- c('Gender', 'Age', 'Race', 'Marital')

ymodel <- as.formula(sprintf('Y~%s', paste0(c(base_covnames, covnames, intervention_name), collapse = '+')))
censor_model <- as.formula(sprintf('C~%s', paste0(c(covnames, intervention_name), collapse = '+')))
covtypes <- rep('binary', length(covnames) + 1)
### the last model for treatment doesn't include previous A?
covmodels <- lapply(seq(c(covnames, intervention_name)), function(j){
  as.formula(sprintf('%s~%s', c(covnames, intervention_name)[j], 
                     paste0(paste0('lag1_', c(covnames, intervention_name)[-j]), collapse = '+')))
})
### include previous treatment in the model for treatment
covmodels[[25]] = as.formula(sprintf('%s+%s',covmodels[25],paste0('lag1_',intervention_name)))
cov_mintimes <- rep(1, length(covnames) + 1)
histvars <- c(covnames, intervention_name)
histvals <- 1 # which lag variables are needed; or histvals <- 1:2

t0 <- Sys.time()
rst1 <- 
  gform_noniter_complete(dffull, K = K, 
                         time_name = 't0', id_name = 'id',
                         outcome_name = outcome_name, ymodel = ymodel,
                         outcome_mintime = 0, 
                         censor_name = censor_name, censor_model = censor_model, 
                         censor_mintime = 0, 
                         intervention_name = intervention_name,
                         intervention = rep(1, K),
                         covnames = covnames, covtypes = covtypes, 
                         covmodels = covmodels, 
                         base_covnames = base_covnames,
                         cov_mintimes = cov_mintimes,
                         histvars = histvars, 
                         histvals = histvals,
                         seed = myseed)
t1 <- Sys.time()
rst2 <- 
  gform_noniter_match(dffull, K = K, J = 10,
                      outcome_name = outcome_name, ymodel = ymodel,
                      outcome_mintime = 0, 
                      censor_name = censor_name, censor_model = censor_model, 
                      censor_mintime = 0, 
                      intervention_name = intervention_name,
                      intervention = rep(1, K),
                      covnames = covnames, covtypes = covtypes, 
                      covmodels = covmodels, 
                      base_covnames = base_covnames,
                      cov_mintimes = cov_mintimes,
                      histvars = histvars, 
                      histvals = histvals,
                      seed = myseed)
t2 <- Sys.time()
rst3 <- 
  gform_iter_complete(dffull, K = K, 
                      id_name = 'id',
                      outcome_name = outcome_name, ymodel = ymodel,
                      outcome_mintime = 0, 
                      intervention_name = intervention_name,
                      intervention = rep(1, K),
                      cov_mintimes = cov_mintimes,
                      histvars = histvars, 
                      histvals = histvals,
                      seed = myseed)
t3 <- Sys.time()
rst4 <- 
  gform_iter_match(dffull, K = K, Js = c(5,5,10,10,20,20),
                   outcome_name = outcome_name, ymodel = ymodel,
                   outcome_mintime = 0, 
                   base_covnames = base_covnames,
                   intervention_name = intervention_name,
                   intervention = rep(1, K),
                   histvars = histvars, 
                   histvals = histvals,
                   seed = myseed)
t4 <- Sys.time()
# rst_pack <- gformula(obs_data = dffull, id = 'id',
#                      time_name = 't0', covnames = c(covnames, intervention_name),
#                      outcome_name = 'Y', outcome_type = 'survival',
#                      censor_name = censor_name, censor_model = censor_model,
#                      covtypes = covtypes,
#                      covparams = list(covmodels = covmodels),
#                      basecovs = base_covnames,
#                      ymodel = ymodel,
#                      histories = c(lagged), histvars = list(histvars),
#                      # intvars = list(intervention_name, intervention_name),
#                      # interventions = list(list(c(static, rep(0, K))),
#                      #                      list(c(static, rep(1, K)))),
#                      # int_descript = c('Never treat', 'Always treat'),
#                      intvars = list(intervention_name),
#                      interventions = list(list(c(static, rep(1, K)))),
#                      int_descript = c('Always treat'),
#                      model_fits = TRUE,
#                      seed = myseed)
# t5 <- Sys.time()
# rst5 <- list(); rst5$risks <- rst_pack$result$`g-form risk`[rst_pack$result$Interv.==1]
# rst5$risks


df_treat <- rbind(data.frame(Time = 1:K, Risk = rst1$risks, Method = 'noniter_complete'), 
                  data.frame(Time = 1:K, Risk = rst2$risks, Method = 'noniter_match'),
                  data.frame(Time = 1:K, Risk = rst3$risks, Method = 'iter_complete'),
                  data.frame(Time = 1:K, Risk = rst4$risks, Method = 'iter_match'))
# data.frame(Time = 1:K, Risk = rst5$risks, Method = 'noniter_complete_package'))
time_treat <- rbind(data.frame(Time = t1-t0, Method = 'noniter_complete'), 
                    data.frame(Time = t2-t1, Method = 'noniter_match'),
                    data.frame(Time = t3-t2, Method = 'iter_complete'),
                    data.frame(Time = t4-t3, Method = 'iter_match'))
# data.frame(Time = t5-t4, Method = 'noniter_complete_package'))

t0 <- Sys.time()
rst1 <- 
  gform_noniter_complete(dffull, K = K, 
                         time_name = 't0', id_name = 'id',
                         outcome_name = outcome_name, ymodel = ymodel,
                         outcome_mintime = 0, 
                         censor_name = censor_name, censor_model = censor_model, 
                         censor_mintime = 0, 
                         intervention_name = intervention_name,
                         intervention = rep(0, K),
                         covnames = covnames, covtypes = covtypes, 
                         covmodels = covmodels, 
                         base_covnames = base_covnames,
                         cov_mintimes = cov_mintimes,
                         histvars = histvars, 
                         histvals = histvals,
                         seed = myseed)
t1 <- Sys.time()
rst2 <- 
  gform_noniter_match(dffull, K = K, J = 5,
                      outcome_name = outcome_name, ymodel = ymodel,
                      outcome_mintime = 0, 
                      censor_name = censor_name, censor_model = censor_model, 
                      censor_mintime = 0, 
                      intervention_name = intervention_name,
                      intervention = rep(0, K),
                      covnames = covnames, covtypes = covtypes, 
                      covmodels = covmodels, 
                      base_covnames = base_covnames,
                      cov_mintimes = cov_mintimes,
                      histvars = histvars, 
                      histvals = histvals,
                      seed = myseed)
t2 <- Sys.time()
rst3 <- 
  gform_iter_complete(dffull, K = K, 
                      outcome_name = outcome_name, ymodel = ymodel,
                      outcome_mintime = 0, 
                      intervention_name = intervention_name,
                      intervention = rep(0, K),
                      cov_mintimes = cov_mintimes,
                      histvars = histvars, 
                      histvals = histvals,
                      seed = myseed)
t3 <- Sys.time()
rst4 <- 
  gform_iter_match(dffull, K = K, J = 5,
                   outcome_name = outcome_name, ymodel = ymodel,
                   outcome_mintime = 0, 
                   intervention_name = intervention_name,
                   intervention = rep(0, K),
                   histvars = histvars, 
                   histvals = histvals,
                   seed = myseed)
t4 <- Sys.time()
# rst_pack <- gformula(obs_data = dffull, id = 'id',
#                      time_name = 't0', covnames = c(covnames, intervention_name),
#                      outcome_name = 'Y', outcome_type = 'survival',
#                      censor_name = censor_name, censor_model = censor_model,
#                      covtypes = covtypes,
#                      covparams = list(covmodels = covmodels),
#                      basecovs = base_covnames,
#                      ymodel = ymodel,
#                      histories = c(lagged), histvars = list(histvars),
#                      # intvars = list(intervention_name, intervention_name),
#                      # interventions = list(list(c(static, rep(0, K))),
#                      #                      list(c(static, rep(1, K)))),
#                      # int_descript = c('Never treat', 'Always treat'),
#                      intvars = list(intervention_name),
#                      interventions = list(list(c(static, rep(0, K)))),
#                      int_descript = c('Always treat'),
#                      model_fits = TRUE,
#                      seed = myseed)
# t5 <- Sys.time()
# rst5 <- list(); rst5$risks <- rst_pack$result$`g-form risk`[rst_pack$result$Interv.==1]
# rst5$risks


df_notreat <- rbind(data.frame(Time = 1:K, Risk = rst1$risks, Method = 'noniter_complete'), 
                    data.frame(Time = 1:K, Risk = rst2$risks, Method = 'noniter_match'),
                    data.frame(Time = 1:K, Risk = rst3$risks, Method = 'iter_complete'),
                    data.frame(Time = 1:K, Risk = rst4$risks, Method = 'iter_match'))
# data.frame(Time = 1:K, Risk = rst5$risks, Method = 'noniter_complete_package'))
time_notreat <- rbind(data.frame(Time = t1-t0, Method = 'noniter_complete'), 
                      data.frame(Time = t2-t1, Method = 'noniter_match'),
                      data.frame(Time = t3-t2, Method = 'iter_complete'),
                      data.frame(Time = t4-t3, Method = 'iter_match'))
# data.frame(Time = t5-t4, Method = 'noniter_complete_package'))


final_rst <- list(df_treat = df_treat,
                  df_notreat = df_notreat,
                  time_treat = time_treat,
                  time_notreat = time_notreat,
                  prev =  mean(dffull$Y, na.rm=TRUE))
