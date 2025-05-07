library(dplyr)
library(tidyr)
source('updatedMatchedGformula.R')

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

nboot <- 100  # Set number of bootstrap replicates
boot_results <- vector("list", nboot)
pb <- txtProgressBar(min = 1, max = nboot + 1, style = 3) 

if(itergform == 0 & matched == 0){
  time_start <- Sys.time()
  for (b in 1:nboot) {
    
    setTxtProgressBar(pb, b) 
    
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
    
    
    rst <- gform_noniter_complete(dfboot, K = K, time_name = time_name, id_name = id_name,
                                  outcome_name = outcome_name, ymodel = ymodel,
                                  outcome_mintime = 0, censor_name = censor_name,
                                  censor_model = censor_model, censor_mintime = 0,
                                  intervention_name = intervention_name, intervention = rep(1, K),
                                  covnames = covnames, covtypes = covtypes, covmodels = covmodels,
                                  base_covnames = base_covnames, cov_mintimes = cov_mintimes,
                                  histvars = histvars, histvals = histvals, seed = seed_b)
    
    boot_results[[b]] <- rst
    
    
    setTxtProgressBar(pb, b); 
    cat(sprintf("Iterative %d; Matched; %d; Treatment %d; Replicate %3d finished at %s\n", 
                itergform, matched, treatment, b, Sys.time()))
  }
  time_end <- Sys.time()
}
if(itergform == 0 & matched == 1){
  time_start <- Sys.time()
  for (b in 1:nboot) {
    
    setTxtProgressBar(pb, b) 
    
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
    
    rst <- gform_noniter_match(dfboot, K = K, J = 10,
                               outcome_name = outcome_name, ymodel = ymodel,
                               outcome_mintime = 0, censor_name = censor_name, 
                               censor_model = censor_model, censor_mintime = 0, 
                               intervention_name = intervention_name,
                               intervention = rep(1, K),
                               covnames = covnames, covtypes = covtypes, covmodels = covmodels, 
                               base_covnames = base_covnames, cov_mintimes = cov_mintimes,
                               histvars = histvars, histvals = histvals,
                               seed = seed_b)
    
    boot_results[[b]] <- rst
    
    
    setTxtProgressBar(pb, b); 
    cat(sprintf("Iterative %d; Matched; %d; Treatment %d; Replicate %3d finished at %s\n", 
                itergform, matched, treatment, b, Sys.time()))
  }
  time_end <- Sys.time()
}
if(itergform == 1 & matched == 0){
  time_start <- Sys.time()
  for (b in 1:nboot) {
    
    setTxtProgressBar(pb, b) 
    
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
    
    rst <- gform_iter_complete(dfboot, K = K, 
                               id_name = 'id',
                               outcome_name = outcome_name, ymodel = ymodel,
                               outcome_mintime = 0, 
                               intervention_name = intervention_name,
                               intervention = rep(1, K),
                               cov_mintimes = cov_mintimes,
                               histvars = histvars, 
                               histvals = histvals,
                               seed = seed_b)
    
    boot_results[[b]] <- rst
    
    
    setTxtProgressBar(pb, b); 
    cat(sprintf("Iterative %d; Matched; %d; Treatment %d; Replicate %3d finished at %s\n", 
                itergform, matched, treatment, b, Sys.time()))
  }
  time_end <- Sys.time()
}
if(itergform == 1 & matched == 1){
  time_start <- Sys.time()
  for (b in 1:nboot) {
    
    setTxtProgressBar(pb, b) 
    
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
    
    rst <- gform_iter_match(dfboot, K = K, J = 10,
                            outcome_name = outcome_name, ymodel = ymodel,
                            outcome_mintime = 0, 
                            base_covnames = base_covnames,
                            intervention_name = intervention_name,
                            intervention = rep(1, K),
                            histvars = histvars, 
                            histvals = histvals,
                            seed = seed_b)
    
    boot_results[[b]] <- rst
    
    
    setTxtProgressBar(pb, b); 
    cat(sprintf("Iterative %d; Matched; %d; Treatment %d; Replicate %3d finished at %s\n", 
                itergform, matched, treatment, b, Sys.time()))
  }
  time_end <- Sys.time()
}


saveRDS(list(boot_results = boot_results, 
             event_prev =  dffull %>% group_by(t0) %>% summarise(avg = mean(Y, na.rm = TRUE)) %>% pull(avg),
             censor_prev = dffull %>% group_by(t0) %>% summarise(avg = mean(C, na.rm = TRUE)) %>% pull(avg),
             time = time_end - time_start), 
        file = sprintf('../Results_Updated/n_%d_setup%d_iter%d_matched%d_treat%d_simu%d.RDS', 
                       n, setup, itergform, matched, treatment, dataID))