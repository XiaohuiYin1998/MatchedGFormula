library(dplyr)
library(tidyr)
source('updatedMatchedGformula.R')
n <- 30000
K <- 6
setup <- 5
source('0_regularizedParameterSetup.R')

simuID_start <- 0

set.seed(2024)
myseeds <- sample(1:1e5, size = 500)
simuID <- Sys.getenv("SLURM_ARRAY_TASK_ID")
simuID <- as.integer(simuID)
if(is.na(simuID)){simuID <- 1}

simuID <- simuID + simuID_start

dataID <- floor((simuID -1)/ 800) + 1
argID <- (simuID - 1) %% 800 + 1
argDF <- rbind(data.frame(itergform = c(0,0,1,1), matched = c(0,1,0,1), treatment = 1),
               data.frame(itergform = c(0,0,1,1), matched = c(0,1,0,1), treatment = 0))
argDF_rep <- bind_rows(lapply(1:100, function(i) {
  argDF %>% mutate(b = i)}))
itergform <- argDF_rep$itergform[argID] # 1 for iterative g formula, 0 for non-iterative gformula
matched <- argDF_rep$matched[argID] # 1 for matched, 0 for complete
treatment <- argDF_rep$treatment[argID] # 1 for always treatment, 0 for never treatment
b <- argDF_rep$b[argID] # bootstrap 
myseed <- myseeds[dataID]

load(sprintf('../Data_Bootstrap/Data_setup%d_replicate%d.rda', setup, dataID))
print(dffull %>% group_by(t0) %>% summarise(mean(Y, na.rm = T)))

time_name = 't0'
id_name = 'id'
base_covnames <- c('Gender', 'Age', 'Race', 'Marital')

ymodel <- as.formula(sprintf('Y~%s', paste0(c(base_covnames, covnames, intervention_name), collapse = '+')))
censor_model <- as.formula(sprintf('C~%s', paste0(c(covnames, intervention_name), collapse = '+')))

covtypes <- rep('binary', length(covnames) + 1)
covmodels <- lapply(seq(covnames), function(j){
  as.formula(sprintf('%s~%s', c(covnames, intervention_name)[j], 
                     paste0(paste0('lag1_', c(covnames, intervention_name)[-j]), collapse = '+')))
})
covmodels[['A']] <- as.formula(sprintf('%s~%s', 
                                       intervention_name,
                                       paste0(paste0('lag1_',c(intervention_name, covnames)), collapse = '+')))

cov_mintimes <- rep(1, length(covnames) + 1)
histvars <- c(covnames, intervention_name)
histvals <- 1


if(itergform == 0 & matched == 0){
  
  seed_b <- myseed + b
  set.seed(seed_b)  
  
  sampled_ids <- sample(unique(dffull$id), size = length(unique(dffull$id)), replace = TRUE)
  
  boot_id_map <- tibble(
    old_id = sampled_ids,
    new_id = seq_along(sampled_ids)
  )
  dfboot <- boot_id_map %>%
    left_join(dffull, by = c("old_id" = "id"), relationship = "many-to-many") %>%
    mutate(id = new_id) %>%
    select(-old_id, -new_id)
  
  dfboot <- as.data.table(dfboot)
  
  rst <- gform_noniter_complete(dfboot, K = K, time_name = time_name, id_name = id_name,
                                outcome_name = outcome_name, ymodel = ymodel,
                                outcome_mintime = 0, censor_name = censor_name,
                                censor_model = censor_model, censor_mintime = 0,
                                intervention_name = intervention_name, intervention = rep(treatment, K),
                                covnames = covnames, covtypes = covtypes, covmodels = covmodels,
                                base_covnames = base_covnames, cov_mintimes = cov_mintimes,
                                histvars = histvars, histvals = histvals, seed = seed_b)
}
if(itergform == 0 & matched == 1){
  seed_b <- myseed + b
  set.seed(seed_b)  
  
  sampled_ids <- sample(unique(dffull$id), size = length(unique(dffull$id)), replace = TRUE)
  
  boot_id_map <- tibble(
    old_id = sampled_ids,
    new_id = seq_along(sampled_ids)
  )
  dfboot <- boot_id_map %>%
    left_join(dffull, by = c("old_id" = "id"), relationship = "many-to-many") %>%
    mutate(id = new_id) %>%
    select(-old_id, -new_id)
  
  dfboot <- as.data.table(dfboot)
  
  rst <- gform_noniter_match(dfboot, K = K, J = 5,
                             outcome_name = outcome_name, ymodel = ymodel,
                             outcome_mintime = 0, censor_name = censor_name, 
                             censor_model = censor_model, censor_mintime = 0, 
                             intervention_name = intervention_name,
                             intervention = rep(treatment, K),
                             covnames = covnames, covtypes = covtypes, covmodels = covmodels, 
                             base_covnames = base_covnames, cov_mintimes = cov_mintimes,
                             histvars = histvars, histvals = histvals,
                             seed = seed_b)
  
}
if(itergform == 1 & matched == 0){
  
  seed_b <- myseed + b
  set.seed(seed_b)  
  
  sampled_ids <- sample(unique(dffull$id), size = length(unique(dffull$id)), replace = TRUE)
  
  boot_id_map <- tibble(
    old_id = sampled_ids,
    new_id = seq_along(sampled_ids)
  )
  dfboot <- boot_id_map %>%
    left_join(dffull, by = c("old_id" = "id"), relationship = "many-to-many") %>%
    mutate(id = new_id) %>%
    select(-old_id, -new_id)
  
  dfboot <- as.data.table(dfboot)
  
  rst <- gform_iter_complete(dfboot, K = K, 
                             id_name = 'id',
                             outcome_name = outcome_name, ymodel = ymodel,
                             outcome_mintime = 0, 
                             intervention_name = intervention_name,
                             intervention = rep(treatment, K),
                             cov_mintimes = cov_mintimes,
                             histvars = histvars, 
                             histvals = histvals,
                             seed = seed_b)
}
if(itergform == 1 & matched == 1){
  
  
  seed_b <- myseed + b
  set.seed(seed_b)  
  
  sampled_ids <- sample(unique(dffull$id), size = length(unique(dffull$id)), replace = TRUE)
  
  boot_id_map <- tibble(
    old_id = sampled_ids,
    new_id = seq_along(sampled_ids)
  )
  dfboot <- boot_id_map %>%
    left_join(dffull, by = c("old_id" = "id"), relationship = "many-to-many") %>%
    mutate(id = new_id) %>%
    select(-old_id, -new_id)
  
  
  dfboot <- as.data.table(dfboot)
  
  rst <- gform_iter_match(dfboot, K = K, J = 5,
                          outcome_name = outcome_name, ymodel = ymodel,
                          outcome_mintime = 0, 
                          base_covnames = base_covnames,
                          intervention_name = intervention_name,
                          intervention = rep(treatment, K),
                          histvars = histvars, 
                          histvals = histvals,
                          seed = seed_b)
}

df_risks <- data.frame(Risk = rst$risks, 
                       Time = 1:K,
                       iterative = itergform,
                       matched = matched, 
                       treatment = treatment,
                       dataID = dataID,
                       b = b,
                       Computation_Data = as.numeric(rst$time_data),
                       Computation_Fit = as.numeric(rst$time_fit))
# write.table(df_risks, 
#             file = sprintf('../Results_Boots_Single/Results_n_%d_setup%d.txt', 
#                            n, setup), 
#             append = TRUE, row.names = FALSE, col.names = FALSE, 
#             quote = FALSE, sep = "\t")

print(df_risks)
cat(sprintf('../Results_Boots_Single_Setup5/n_%d_setup%d_iter%d_matched%d_treat%d_simu%d_boots%d.txt', 
            n, setup, itergform, matched, treatment, dataID, b))
write.table(df_risks, 
            file = sprintf('../Results_Boots_Single_Setup5/n_%d_setup%d_iter%d_matched%d_treat%d_simu%d_boots%d.txt', 
                           n, setup, itergform, matched, treatment, dataID, b), 
            row.names = FALSE, col.names = TRUE, 
            quote = FALSE, sep = "\t")

# saveRDS(rst, 
#         file = sprintf('../Results_Boots_Single/n_%d_setup%d_iter%d_matched%d_treat%d_simu%d_boots%d.RDS', 
#                        n, setup, itergform, matched, treatment, dataID, b))



