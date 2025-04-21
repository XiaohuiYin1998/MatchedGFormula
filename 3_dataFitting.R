source('matchedGFormula.R')

time_name = 't0'
id_name = 'id'
base_covnames <- c('Gender', 'Age', 'Race', 'Marital')

ymodel <- as.formula(sprintf('Y~%s', paste0(c(base_covnames, covnames, intervention_name), collapse = '+')))
censor_model <- as.formula(sprintf('C~%s', paste0(c(covnames, intervention_name), collapse = '+')))

covtypes <- rep('binary', length(covnames) + 1)
### the last model for treatment doesn't include previous A?
covmodels <- lapply(seq(covnames), function(j){
  as.formula(sprintf('%s~%s', c(covnames, intervention_name)[j], 
                     paste0(paste0('lag1_', c(covnames, intervention_name)[-j]), collapse = '+')))
})
### include previous treatment in the model for treatment
covmodels[['A']] = as.formula(sprintf('%s~%s', 
                                      intervention_name,
                                      paste0(paste0('lag1_',c(intervention_name, covnames)), collapse = '+')))
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
rst2<- 
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


df_treat <- rbind(data.frame(Time = 1:K, Risk = rst1$risks, Method = 'noniter_complete'), 
                    data.frame(Time = 1:K, Risk = rst2$risks, Method = 'noniter_match'),
                    data.frame(Time = 1:K, Risk = rst3$risks, Method = 'iter_complete'),
                    data.frame(Time = 1:K, Risk = rst4$risks, Method = 'iter_match'))
time_treat <- rbind(data.frame(Time = t1-t0, Method = 'noniter_complete'), 
                      data.frame(Time = t2-t1, Method = 'noniter_match'),
                      data.frame(Time = t3-t2, Method = 'iter_complete'),
                      data.frame(Time = t4-t3, Method = 'iter_match'))


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
  gform_noniter_match(dffull, K = K, J = 10,
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
                      id_name = 'id',
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
  gform_iter_match(dffull, K = K, Js = c(5,5,10,10,20,20),
                   outcome_name = outcome_name, ymodel = ymodel,
                   outcome_mintime = 0, 
                   base_covnames = base_covnames,
                   intervention_name = intervention_name,
                   intervention = rep(0, K),
                   histvars = histvars, 
                   histvals = histvals,
                   seed = myseed)
t4 <- Sys.time()


df_notreat <- rbind(data.frame(Time = 1:K, Risk = rst1$risks, Method = 'noniter_complete'), 
                    data.frame(Time = 1:K, Risk = rst2$risks, Method = 'noniter_match'),
                    data.frame(Time = 1:K, Risk = rst3$risks, Method = 'iter_complete'),
                    data.frame(Time = 1:K, Risk = rst4$risks, Method = 'iter_match'))
time_notreat <- rbind(data.frame(Time = t1-t0, Method = 'noniter_complete'), 
                      data.frame(Time = t2-t1, Method = 'noniter_match'),
                      data.frame(Time = t3-t2, Method = 'iter_complete'),
                      data.frame(Time = t4-t3, Method = 'iter_match'))

final_rst <- list(df_treat = df_treat,
                  df_notreat = df_notreat,
                  time_treat = time_treat,
                  time_notreat = time_notreat,
                  prev =  mean(dffull$Y, na.rm=TRUE))
saveRDS(final_rst, file = sprintf('../Results/repli_%d.RDS', simuID))



