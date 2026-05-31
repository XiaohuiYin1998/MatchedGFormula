library(fastDummies)
source('regularizedMatchedGformula.R')

n     <- as.integer(Sys.getenv('GFORM_N',     unset = '30000'))
K     <- as.integer(Sys.getenv('GFORM_K',     unset = '6'))
setup <- as.integer(Sys.getenv('GFORM_SETUP', unset = '5'))

source('0_regularizedParameterSetup.R')

set.seed(2024)
myseeds <- sample(1:1e5, size = 500)

simuID <- Sys.getenv("SLURM_ARRAY_TASK_ID")
simuID <- as.integer(simuID)
if (is.na(simuID)) simuID <- 1
myseed <- myseeds[simuID]

data_dir <- 'data'
out_dir  <- 'results_regularized'
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
load(sprintf('%s/Data_setup%d_replicate%d.rda', data_dir, setup, simuID))

time_name     <- 't0'
id_name       <- 'id'
base_covnames <- c('Gender', 'Age', 'Race', 'Marital')

dffull_dummy <- dummy_cols(dffull,
                           select_columns = base_covnames,
                           remove_selected_columns = TRUE,
                           remove_first_dummy = TRUE)
base_covnames <- setdiff(colnames(dffull_dummy), colnames(dffull))

ymodel       <- as.formula(sprintf('Y~%s',
                                   paste0(c(base_covnames, covnames, intervention_name),
                                          collapse = '+')))
censor_model <- as.formula(sprintf('C~%s',
                                   paste0(c(covnames, intervention_name),
                                          collapse = '+')))

covtypes <- rep('binary', length(covnames) + 1)
covmodels <- lapply(seq(covnames), function(j) {
  as.formula(sprintf('%s~%s', c(covnames, intervention_name)[j],
                     paste0(paste0('lag1_', c(covnames, intervention_name)[-j]),
                            collapse = '+')))
})
covmodels[['A']] <- as.formula(sprintf('%s~%s',
                                       intervention_name,
                                       paste0(paste0('lag1_',
                                                     c(intervention_name, covnames)),
                                              collapse = '+')))

cov_mintimes <- rep(1, length(covnames) + 1)
histvars     <- c(covnames, intervention_name)
histvals     <- 1

run_all_methods <- function(intervention_vec) {
  t0 <- Sys.time()
  rst1 <- gform_noniter_complete(dffull_dummy, K = K,
                                 time_name = 't0', id_name = 'id',
                                 outcome_name = outcome_name, ymodel = ymodel,
                                 outcome_mintime = 0,
                                 censor_name = censor_name, censor_model = censor_model,
                                 censor_mintime = 0,
                                 intervention_name = intervention_name,
                                 intervention = intervention_vec,
                                 covnames = covnames, covtypes = covtypes,
                                 covmodels = covmodels,
                                 base_covnames = base_covnames,
                                 cov_mintimes = cov_mintimes,
                                 histvars = histvars, histvals = histvals,
                                 seed = myseed)
  t1 <- Sys.time()
  rst2 <- gform_noniter_match(dffull_dummy, K = K, J = 10,
                              outcome_name = outcome_name, ymodel = ymodel,
                              outcome_mintime = 0,
                              censor_name = censor_name, censor_model = censor_model,
                              censor_mintime = 0,
                              intervention_name = intervention_name,
                              intervention = intervention_vec,
                              covnames = covnames, covtypes = covtypes,
                              covmodels = covmodels,
                              base_covnames = base_covnames,
                              cov_mintimes = cov_mintimes,
                              histvars = histvars, histvals = histvals,
                              seed = myseed)
  t2 <- Sys.time()
  rst3 <- gform_iter_complete(dffull_dummy, K = K, id_name = 'id',
                              outcome_name = outcome_name, ymodel = ymodel,
                              outcome_mintime = 0,
                              intervention_name = intervention_name,
                              intervention = intervention_vec,
                              cov_mintimes = cov_mintimes,
                              histvars = histvars, histvals = histvals,
                              seed = myseed)
  t3 <- Sys.time()
  rst4 <- gform_iter_match(dffull_dummy, K = K, Js = c(5, 5, 10, 10, 20, 20),
                           outcome_name = outcome_name, ymodel = ymodel,
                           outcome_mintime = 0,
                           base_covnames = base_covnames,
                           intervention_name = intervention_name,
                           intervention = intervention_vec,
                           histvars = histvars, histvals = histvals,
                           seed = myseed)
  t4 <- Sys.time()

  list(
    risks = rbind(data.frame(Time = 1:K, Risk = rst1$risks, Method = 'noniter_complete'),
                  data.frame(Time = 1:K, Risk = rst2$risks, Method = 'noniter_match'),
                  data.frame(Time = 1:K, Risk = rst3$risks, Method = 'iter_complete'),
                  data.frame(Time = 1:K, Risk = rst4$risks, Method = 'iter_match')),
    times = rbind(data.frame(Time = t1 - t0, Method = 'noniter_complete'),
                  data.frame(Time = t2 - t1, Method = 'noniter_match'),
                  data.frame(Time = t3 - t2, Method = 'iter_complete'),
                  data.frame(Time = t4 - t3, Method = 'iter_match'))
  )
}

treat   <- run_all_methods(rep(1, K))
notreat <- run_all_methods(rep(0, K))

final_rst <- list(df_treat     = treat$risks,
                  df_notreat   = notreat$risks,
                  time_treat   = treat$times,
                  time_notreat = notreat$times,
                  prev         = mean(dffull$Y, na.rm = TRUE))

saveRDS(final_rst,
        file = sprintf('%s/regularized_repli_%d.RDS', out_dir, simuID))
