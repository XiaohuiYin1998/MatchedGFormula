set.seed(2024)
myseeds <- sample(1:1e5, size = 500)
simuID <- Sys.getenv("SLURM_ARRAY_TASK_ID")
simuID <- as.integer(simuID)
if(is.na(simuID)){simuID <- 23}
myseed <- myseeds[simuID]
set.seed(myseed)

min_count_Y <- 1
while(min_count_Y < 5){
  n <- 10000;
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



