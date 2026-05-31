library(dplyr)

set.seed(2024)
myseeds <- sample(1:1e5, size = 500)

simuID <- Sys.getenv("SLURM_ARRAY_TASK_ID")
simuID <- as.integer(simuID)
if (is.na(simuID)) simuID <- 1

dataID  <- simuID
myseed  <- myseeds[dataID]
set.seed(myseed)

out_dir <- 'data'
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

n <- 30000
min_count_Y <- 1
while (min_count_Y < 5) {
  cat('generating data...\n')
  df <- lapply(as.list(1:n), FUN = function(ind) {
    va_simu_data(ind, K = K,
                 covmodels_star_coef_A   = covmodels_star_coef_A,
                 covmodels_star_coef     = covmodels_star_coef,
                 censor_model_star_coef  = censor_model_star_coef,
                 outcome_model_star_coef = outcome_model_star_coef)
  })
  dffull <- do.call('rbind', df)
  colnames(dffull)[grep('X', colnames(dffull))]      <- covnames
  colnames(dffull)[colnames(dffull) == 'A']          <- intervention_name

  dffull <- as.data.table(dffull)
  min_count_Y <- min(table(dffull$t0, dffull$Y))
}

print(dffull %>% group_by(t0) %>% summarise(mean(Y, na.rm = TRUE)))
save(dffull, file = sprintf('%s/Data_setup%d_replicate%d.rda', out_dir, setup, dataID))
