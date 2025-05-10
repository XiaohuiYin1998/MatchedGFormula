library(dplyr)
set.seed(2024)
myseeds <- sample(1:1e5, size = 500)
simuID <- Sys.getenv("SLURM_ARRAY_TASK_ID")
simuID <- as.integer(simuID)
if(is.na(simuID)){simuID <- 1}

dataID <- simuID
# argID <- (simuID - 1) %% 8 + 1
# argDF <- rbind(data.frame(itergform = c(0,1,0,1), matched = c(0,1,0,1), treatment = 1),
#                 data.frame(itergform = c(0,1,0,1), matched = c(0,1,0,1), treatment = 0))
# itergform <- argDF$itergform[argID] # 1 for iterative g formula, 0 for non-iterative gformula
# matched <- argDF$matched[argID] # 1 for matched, 0 for complete
# treatment <- argDF$treatment[argID] # 1 for always treatment, 0 for never treatment



myseed <- myseeds[dataID]
set.seed(myseed)

min_count_Y <- 1
while(min_count_Y < 5){
  cat('generating data...\n')
  n <- 30000;
  df <- lapply(as.list(1:n), FUN=function(ind){
    va_simu_data(ind, K=K,
                 covmodels_star_coef_A = covmodels_star_coef_A,
                 covmodels_star_coef = covmodels_star_coef,
                 censor_model_star_coef = censor_model_star_coef,
                 outcome_model_star_coef = outcome_model_star_coef)
  })
  dffull <- do.call('rbind', df)
  colnames(dffull)[grep('X', colnames(dffull))] <- covnames
  colnames(dffull)[colnames(dffull) == 'A'] <- intervention_name
  
  dffull <- as.data.table(dffull)
  min_count_Y <- min(table(dffull$t0,dffull$Y))
}
print(dffull %>% group_by(t0) %>% summarise(mean(Y, na.rm = T)))
save(dffull, file = sprintf('../Data_Bootstrap/Data_setup%d_replicate%d.rda', setup, dataID))
