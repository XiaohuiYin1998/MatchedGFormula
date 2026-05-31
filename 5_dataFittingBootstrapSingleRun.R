library(dplyr)
library(tidyr)
source('updatedMatchedGformula.R')

n     <- as.integer(Sys.getenv('GFORM_N',     unset = '30000'))
K     <- as.integer(Sys.getenv('GFORM_K',     unset = '6'))
setup <- as.integer(Sys.getenv('GFORM_SETUP', unset = '5'))
J     <- as.integer(Sys.getenv('GFORM_J',     unset = '5'))

source('0_regularizedParameterSetup.R')

simuID_start <- 0

set.seed(2024)
myseeds <- sample(1:1e5, size = 500)

simuID <- Sys.getenv("SLURM_ARRAY_TASK_ID")
simuID <- as.integer(simuID)
if (is.na(simuID)) simuID <- 1
simuID <- simuID + simuID_start

# Each replicate spans 800 array tasks: 4 estimator combos x 2 treatments x 100 bootstraps
dataID <- floor((simuID - 1) / 800) + 1
argID  <- (simuID - 1) %% 800 + 1
argDF  <- rbind(data.frame(itergform = c(0, 0, 1, 1),
                           matched   = c(0, 1, 0, 1),
                           treatment = 1),
                data.frame(itergform = c(0, 0, 1, 1),
                           matched   = c(0, 1, 0, 1),
                           treatment = 0))
argDF_rep <- bind_rows(lapply(1:100, function(i) argDF %>% mutate(b = i)))
itergform <- argDF_rep$itergform[argID]
matched   <- argDF_rep$matched[argID]
treatment <- argDF_rep$treatment[argID]
b         <- argDF_rep$b[argID]

myseed <- myseeds[dataID]

data_dir <- 'data'
out_dir  <- sprintf('results_bootstrap_single_setup%d', setup)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

load(sprintf('%s/Data_setup%d_replicate%d.rda', data_dir, setup, dataID))
print(dffull %>% group_by(t0) %>% summarise(mean(Y, na.rm = TRUE)))

time_name     <- 't0'
id_name       <- 'id'
base_covnames <- c('Gender', 'Age', 'Race', 'Marital')

ymodel        <- as.formula(sprintf('Y~%s',
                                    paste0(c(base_covnames, covnames, intervention_name),
                                           collapse = '+')))
censor_model  <- as.formula(sprintf('C~%s',
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

seed_b <- myseed + b
set.seed(seed_b)
sampled_ids <- sample(unique(dffull$id), size = length(unique(dffull$id)),
                      replace = TRUE)
boot_id_map <- tibble(old_id = sampled_ids, new_id = seq_along(sampled_ids))
dfboot <- boot_id_map %>%
  left_join(dffull, by = c("old_id" = "id"), relationship = "many-to-many") %>%
  mutate(id = new_id) %>%
  select(-old_id, -new_id)
dfboot <- as.data.table(dfboot)

if (itergform == 0 && matched == 0) {
  rst <- gform_noniter_complete(dfboot, K = K, time_name = time_name, id_name = id_name,
                                outcome_name = outcome_name, ymodel = ymodel,
                                outcome_mintime = 0, censor_name = censor_name,
                                censor_model = censor_model, censor_mintime = 0,
                                intervention_name = intervention_name,
                                intervention = rep(treatment, K),
                                covnames = covnames, covtypes = covtypes,
                                covmodels = covmodels,
                                base_covnames = base_covnames,
                                cov_mintimes = cov_mintimes,
                                histvars = histvars, histvals = histvals,
                                seed = seed_b)
} else if (itergform == 0 && matched == 1) {
  rst <- gform_noniter_match(dfboot, K = K, J = J,
                             outcome_name = outcome_name, ymodel = ymodel,
                             outcome_mintime = 0, censor_name = censor_name,
                             censor_model = censor_model, censor_mintime = 0,
                             intervention_name = intervention_name,
                             intervention = rep(treatment, K),
                             covnames = covnames, covtypes = covtypes,
                             covmodels = covmodels,
                             base_covnames = base_covnames,
                             cov_mintimes = cov_mintimes,
                             histvars = histvars, histvals = histvals,
                             seed = seed_b)
} else if (itergform == 1 && matched == 0) {
  rst <- gform_iter_complete(dfboot, K = K, id_name = 'id',
                             outcome_name = outcome_name, ymodel = ymodel,
                             outcome_mintime = 0,
                             intervention_name = intervention_name,
                             intervention = rep(treatment, K),
                             cov_mintimes = cov_mintimes,
                             histvars = histvars, histvals = histvals,
                             seed = seed_b)
} else {
  rst <- gform_iter_match(dfboot, K = K, J = J,
                          outcome_name = outcome_name, ymodel = ymodel,
                          outcome_mintime = 0,
                          base_covnames = base_covnames,
                          intervention_name = intervention_name,
                          intervention = rep(treatment, K),
                          histvars = histvars, histvals = histvals,
                          seed = seed_b)
}

df_risks <- data.frame(Risk             = rst$risks,
                       Time             = 1:K,
                       iterative        = itergform,
                       matched          = matched,
                       treatment        = treatment,
                       dataID           = dataID,
                       b                = b,
                       J                = J,
                       Computation_Data = as.numeric(rst$time_data),
                       Computation_Fit  = as.numeric(rst$time_fit))

print(df_risks)
out_file <- sprintf('%s/n_%d_setup%d_J%d_iter%d_matched%d_treat%d_simu%d_boots%d.txt',
                    out_dir, n, setup, J, itergform, matched, treatment, dataID, b)
cat(out_file, "\n")
write.table(df_risks, file = out_file,
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
