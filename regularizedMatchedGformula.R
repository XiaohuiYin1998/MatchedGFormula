library(glmnet)
mylagged <- function(pool, histvars, histvals, time_name, t, id_name){
  for (i in histvals){
    if (t < i){
      lapply(histvars, FUN = function(histvar){
        classtmp <- class(pool[pool[[time_name]] == t][[histvar]])
        myclass <- paste('as.', classtmp, sep = '')
        if (is.factor(pool[pool[[time_name]] == t][[histvar]])){
          reflevel <- levels(pool[pool[[time_name]] == t][[histvar]])[1]
          pool[pool[[time_name]] == t, (paste('lag', i, '_', histvar, sep = '')) := reflevel]
        } else {
          pool[pool[[time_name]] == t, (paste('lag', i, '_', histvar, sep = '')) := get(myclass)(0)]
        }
      })
    } else {
      current_ids <- unique(pool[pool[[time_name]] == t][[id_name]])
      lapply(histvars, FUN = function(histvar){
        tlag <- t - i
        classtmp <- class(pool[pool[[time_name]] == t][[histvar]])
        myclass <- paste('as.', classtmp, sep = '')
        pool[pool[[time_name]] == t, (paste('lag', i, '_', histvar, sep = '')) :=
               get(myclass)(pool[pool[[time_name]] == tlag & get(id_name) %in% current_ids][[histvar]])]
      })
    }
  }
}

mycumavg <- function(pool, histvars, time_name, t, id_name){
  denom <- t + 1
  if (t == 0){
    # At first time point, set all cumulative averages equal to the actual value of the
    # variable
    lapply(histvars, FUN = function(histvar){
      pool[get(time_name) == t, (paste('cumavg_', histvar, sep = '')) :=
             as.double(pool[get(time_name) == t][[histvar]])]
    })
  } else {
    # At subsequent time points, create new column containing calculated cumulative
    # average until that time point
    current_ids <- unique(pool[get(time_name) == t][[id_name]])
    colnam <- colnames(pool)
    if (!(paste('cumavg_', '_', histvars[1], sep = '') %in% colnam)){
      # The cumulative average variable was not created yet
      # Therefore, cannot use recursive formula
      id_factor <- is.factor(pool[[id_name]])
      if (id_factor){
        lapply(histvars, FUN = function(histvar){
          pool[get(time_name) == t, (paste('cumavg_', histvar, sep = '')) :=
                 as.double(tapply(pool[get(id_name) %in% current_ids &
                                         get(time_name) <= t][[histvar]],
                                  droplevels(pool[get(id_name) %in% current_ids &
                                                    get(time_name) <= t][[id_name]]),
                                  FUN = mean))]
        })
      } else {
        lapply(histvars, FUN = function(histvar){
          pool[get(time_name) == t, (paste('cumavg_', histvar, sep = '')) :=
                 as.double(tapply(pool[get(id_name) %in% current_ids &
                                         get(time_name) <= t][[histvar]],
                                  pool[get(id_name) %in% current_ids &
                                         get(time_name) <= t][[id_name]], FUN = mean))]
        })
      }
    } else {
      # The cumulative average variable was already created
      # Therefore, can use recursive formula
      for (histvar in histvars){
        pool[get(time_name) == t, (paste('cumavg_', histvar, sep = '')) :=
               as.double(
                 (pool[get(id_name) %in% current_ids & get(time_name) == (t-1)][[paste('cumavg_', histvar, sep = '')]] * (denom - 1)
                  + pool[get(id_name) %in% current_ids & get(time_name) == t][[histvar]]) / denom
               )
        ]
      }
    }
  }
}

get_predictor_response_name <- function(fomu){
  terms_obj <- terms(fomu)
  response_name <- all.vars(fomu)[attr(terms_obj, 'response')]
  predictor_name <- attr(terms_obj, 'term.labels')
  list(response_name = response_name, predictor_name = predictor_name)
}

fit_regularized_model <- function(pool, response_name, predictor_name, 
                                  var_type, weight_name = NULL, seed, nfold = 5){
  
  y <- as.numeric(pool[[response_name]])
  x <- as.matrix(pool[, ..predictor_name])
  non_na_idx <- which(!is.na(y))
  
  if(!is.null(weight_name)){
    ws <- as.numeric(pool[[weight_name]])
    non_na_idx <- which(!is.na(y) & !is.na(ws))
  }else{ ws = NULL }
  
  set.seed(seed)
  foldid <- sample(rep(1:nfold, length.out = length(y)))
  
  if (var_type == 'binary') {
    if(sum(y!=0 | y!=1, na.rm = TRUE) > 0){
      mod <- glmnet::cv.glmnet(x = x[non_na_idx, , drop = FALSE], 
                               y = cbind(1-y, y)[non_na_idx,], 
                               weights = ws[non_na_idx],
                               family = 'binomial', 
                               foldid = foldid[non_na_idx])
    }else{
      mod <- glmnet::cv.glmnet(x = x[non_na_idx, , drop = FALSE], 
                               y = y[non_na_idx], 
                               weights = ws[non_na_idx],
                               family = 'binomial', 
                               foldid = foldid[non_na_idx])
    }
    
  }else if (var_type == 'normal') {
    mod <- glmnet::cv.glmnet(x = x[non_na_idx, , drop = FALSE], 
                             y = y[non_na_idx], 
                             weights = ws[non_na_idx],
                             family = 'gaussian', 
                             foldid = foldid[non_na_idx])
  }else if (var_type == 'categorical') {
    mod <- glmnet::cv.glmnet(x = x[non_na_idx, , drop = FALSE], 
                             y = y[non_na_idx], 
                             weights = ws,
                             family = 'multinomial', 
                             foldid = foldid[non_na_idx])
  }
  return(mod)
}

predict_regularized_model <- function(model, newx, predictor_name, 
                                      type = 'response', s = 'lambda.min') {
  newx <- as.matrix(newx[, ..predictor_name])
  predict(model, newx = newx, type = type, s = s)
}


gform_noniter_complete <- function(obs_data, K, 
                                   time_name = 't0', id_name = 'id', 
                                   outcome_name, ymodel, outcome_mintime = 0, 
                                   censor_name, censor_model, censor_mintime = 0, 
                                   intervention_name, 
                                   intervention = rep(0, K),
                                   covnames, covtypes, covmodels, cov_mintimes,
                                   base_covnames = NULL,
                                   histvars = NULL, histvals = NULL,
                                   seed = 2024, 
                                   npool = NULL){
  obs_data = copy(obs_data)
  obs_data <- as.data.table(obs_data)
  
  ## Create lagged data
  if(!is.null(histvals)){
    for(t_hist in 0:K){
      mylagged(obs_data, histvars = histvars, histvals = histvals, time_name = time_name, t = t_hist, id_name = id_name)
    }
  }
  
  ## Fit model 
  # fitY <- fit_model(formu = ymodel, pool = obs_data, var_type = 'binary')
  fitY <- fit_regularized_model(pool = obs_data, 
                                response_name = get_predictor_response_name(ymodel)$response_name, 
                                predictor_name = get_predictor_response_name(ymodel)$predictor_name, 
                                weight_name = NULL,
                                seed = 2025,
                                var_type = 'binary')
  
  # fitC <- fit_model(formu = censor_model, pool = obs_data, var_type = 'binary')
  # fitC <- fit_regularized_model(pool = obs_data, 
  #                               response_name = get_predictor_response_name(censor_model)$response_name, 
  #                               predictor_name = get_predictor_response_name(censor_model)$predictor_name, 
  #                               weight_name = NULL,
  #                               seed = 2025,
  #                               var_type = 'binary')
  fitC <- NULL
  fitcov <- lapply(seq(covnames), function(j){
    
    # fit_model(formu = covmodels[[j]], 
    #           pool = obs_data[obs_data[[time_name]]>=cov_mintimes[j]], 
    #           var_type = covtypes[j])
    fit_regularized_model(pool = obs_data, 
                          response_name = get_predictor_response_name(covmodels[[j]])$response_name, 
                          predictor_name = get_predictor_response_name(covmodels[[j]])$predictor_name, 
                          weight_name = NULL,
                          seed = 2025,
                          var_type = 'binary')
  })
  
  names(fitcov) <- covnames
  
  mod_list <- list('outcome' = fitY, 'censor' = fitC, 'cov' = fitcov)
  
  ## Sample pseudo pool and pseudo trajectory
  
  if(!is.null(npool)){
    ids <- as.data.table(sample(unique(obs_data[[id_name]]), 
                                npool, replace = TRUE))
    ids[, `:=`('sid', 1:npool)]
    colnames(ids) <- c(id_name, 'sid')
    obs_data <- obs_data[J(ids), allow.cartesian = TRUE, on = id_name]
    obs_data[, `:=`(id_name, sid)]
    obs_data <- obs_data[, -'sid']
  }
  
  t <- 0
  ids_unique <- unique(obs_data[[id_name]])
  data_len <- length(ids_unique)
  pool <- obs_data[obs_data[[time_name]] <= t, .SD, .SDcols = c(base_covnames, covnames, time_name)]
  set(pool, j = id_name, value = rep(ids_unique, each = 1))
  setcolorder(pool, c(id_name, time_name, base_covnames, covnames))
  set(pool, j = 'eligible_pt', value = TRUE)
  
  newdf <- pool[pool[[time_name]] == t]
  newdf[, (intervention_name) := intervention[1]]
  
  intervened <- rep(0, times = nrow(newdf))
  intervened <- intervened + 
    (abs(pool[[intervention_name]] - newdf[[intervention_name]]) > 1e-06)
  intervened <- ifelse(newdf$eligible_pt, intervened >= 1, NA)
  set(newdf, j = 'intervened', value = intervened)
  
  if (ncol(newdf) > ncol(pool)) {
    pool <- rbind(pool[pool[[time_name]] < t], newdf, 
                  fill = TRUE)
    pool <- pool[order(get(id_name), get(time_name))]
  }else {
    pool[pool[[time_name]] == t] <- newdf
  }
  
  if(!is.null(histvals)){
    for(t_hist in 0:K){
      mylagged(pool, histvars = histvars, histvals = histvals, time_name = time_name, t = t_hist, id_name = id_name)
    }
  }
  newdf <- pool[pool[[time_name]] == t]
  # set(newdf, j = 'Py', value = stats::predict(fitY, type = 'response', newdata = newdf))
  set(newdf, j = 'Py', value = 
        predict_regularized_model(model = fitY, newx = newdf, 
                                  predictor_name = get_predictor_response_name(ymodel)$predictor_name))
  set(newdf, j = 'Y', value = stats::rbinom(data_len, 1, newdf$Py))
  set(newdf, j = 'D', value = 0)
  set(newdf, j = 'prodp1', value = newdf$Py)
  set(newdf[newdf$D == 1], j = 'Y', value = NA)
  set(newdf, j = 'prodp0', value = 1 - newdf$Py)
  set(newdf, j = 'poprisk', value = newdf$prodp1)
  pool <- rbind(pool[pool[[time_name]] < t], newdf, 
                fill = TRUE)
  pool <- pool[order(get(id_name), get(time_name))]
  col_types <- sapply(pool, class)
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  for(t in 1:(K-1)){
    newdf <- pool[pool[[time_name]] == t - 1]
    set(newdf, j = time_name, value = rep(t, data_len))
    pool <- rbind(newdf, pool)
    
    if(!is.null(histvals)){
      for(t_hist in 0:K){
        mylagged(pool, histvars = histvars, histvals = histvals, time_name = time_name, t = t_hist, id_name = id_name)
      }
    }
    newdf <- pool[pool[[time_name]] == t]
    
    for (i in seq_along(covnames)) {
      cast <- get(paste0('as.', unname(col_types[covnames[i]])))
      if (covtypes[i] == 'binary') {
        set(newdf, j = covnames[i], 
            value = cast(stats::rbinom(data_len, 1, 
                                       predict_regularized_model(model = fitcov[[i]], newx = newdf, 
                                                                 predictor_name = get_predictor_response_name(covmodels[[i]])$predictor_name))
                         ))
      }else if (covtypes[i] == 'normal') {
        set(newdf, j = covnames[i], 
            value = cast(stats::rnorm(data_len, mean = stats::predict(fitcov[[i]], type = 'response', newdata = newdf), sd = fitcov[[i]]$rmse)))
      }else if (covtypes[i] == 'categorical') {
        set(newdf, j = covnames[i], 
            value = cast(stats::predict(fitcov[[i]], type = 'class', newdata = newdf)))
      }
      pool[pool[[time_name]] == t] <- newdf
      if (covnames[i] %in% histvars) {
        ind <- unlist(lapply(histvars, FUN = function(x) {
          covnames[i] %in% x
        }))
        mylagged(pool, histvars = rep(histvars, ind), histvals = 1, time_name = 't0', t = t, id_name = 'id')
        newdf <- pool[pool[[time_name]] == t]
      }
    }
    
    newdf <- pool[pool[[time_name]] == t]
    newdf[, (intervention_name) := intervention[t+1]]
    
    intervened <- rep(0, times = nrow(newdf))
    intervened <- intervened + 
      (abs(pool[[intervention_name]] - newdf[[intervention_name]]) > 1e-06)
    intervened <- ifelse(newdf$eligible_pt, intervened >= 1, NA)
    set(newdf, j = 'intervened', value = intervened)
    pool[pool[[time_name]] == t] <- newdf
    
    set(newdf, j = 'D', value = 0)
    set(newdf, j = 'Py', value = 
          predict_regularized_model(model = fitY, newx = newdf, 
                                    predictor_name = get_predictor_response_name(ymodel)$predictor_name))
    set(newdf, j = 'Y', value = stats::rbinom(data_len, 
                                              1, newdf$Py))
    newdf[newdf$D == 1, `:=`('Y', NA)]
    set(newdf, j = 'prodp1', value = newdf$Py * 
          pool[pool[[time_name]] == t - 1, ]$prodp0)
    set(newdf, j = 'prodp0', value = (1 - newdf$Py) * 
          pool[pool[[time_name]] == t - 1, ]$prodp0)
    set(newdf, j = 'poprisk', value = pool[pool[[time_name]] == 
                                             t - 1, ]$poprisk + newdf$prodp1)
    
    pool[pool[[time_name]] == t] <- newdf
    pool <- pool[pool[[time_name]] >= 0]
    pool[, `:=`('survival', 1 - pool$poprisk)]
  }
  risks <- tapply(pool$poprisk, pool[[time_name]], FUN = mean)
  return(list(risks = risks, mod_list = mod_list))
}

gform_noniter_match <- function(obs_data, K, J = 5, 
                                time_name = 't0', id_name = 'id', 
                                outcome_name, ymodel, outcome_mintime = 0, 
                                censor_name, censor_model, censor_mintime = 0, 
                                intervention_name, 
                                intervention = rep(0, K),
                                covnames, covtypes, covmodels, cov_mintimes,
                                base_covnames = NULL,
                                histvars = NULL, histvals = NULL,
                                seed = 2024, 
                                npool = NULL){
  obs_data <- copy(obs_data)
  set.seed(seed)
  colnames(obs_data)[colnames(obs_data) == id_name] <- "id"
  for(t in (K:1)-1){
    
    id_cases <- obs_data[obs_data[[outcome_name]] == 1 & obs_data[[time_name]] == t][["id"]]
    id_controls <- obs_data[obs_data[[outcome_name]] == 0 & obs_data[[time_name]] == t][["id"]]
    n_cases <- length(id_cases)
    n_controls <- length(id_controls)
    
    id_matched_controls <- sample(x = id_controls, size = J * n_cases, replace = TRUE)
    
    weight_cases <- n_cases / (n_cases + n_controls)
    weight_controls <- n_controls / (n_cases + n_controls) / J
    
    # If there are controls selected multiple time,
    # the weights are summed up.
    dt_control_weights <- data.table(id=id_matched_controls, Weights = weight_controls)
    dt_control_weights <- dt_control_weights[, Count := .N, by = id]
    dt_control_weights <- unique(dt_control_weights)
    
    obs_data[obs_data[[time_name]] == t, Weights := as.numeric(NA)]
    obs_data[obs_data[[outcome_name]] == 1 & obs_data[[time_name]] == t, 
             `:=`(Weights = weight_cases,
                  Count = 1)]
    obs_data[obs_data[[outcome_name]] == 0 & obs_data[[time_name]] == t, 
             `:=`(Weights = dt_control_weights$Weights[match(id, dt_control_weights[['id']])], 
                  Count = dt_control_weights$Count[match(id, dt_control_weights[['id']])])]
  }
  obs_data_match <- obs_data[, if (any(!is.na(Count) & Count != 0)) .SD, by = id]
  colnames(obs_data)[colnames(obs_data) == "id"] <- id_name
  colnames(obs_data_match)[colnames(obs_data_match) == "id"] <- id_name
  obs_data_match <- as.data.table(obs_data_match)
  ## Create lagged data
  if(!is.null(histvals)){
    for(t_hist in 0:K){
      mylagged(obs_data_match, histvars = histvars, histvals = histvals, time_name = time_name, t = t_hist, id_name = id_name)
    }
  }
  obs_data_match <- obs_data_match[rep(1:.N, ifelse(is.na(Count), 1, Count))]
  
  
  ## Fit model 
  # fitY <- fit_model(formu = ymodel, pool = obs_data_match, var_type = 'binary', weight_name = 'Weights')
  fitY <- fit_regularized_model(pool = obs_data_match, 
                                response_name = get_predictor_response_name(ymodel)$response_name, 
                                predictor_name = get_predictor_response_name(ymodel)$predictor_name, 
                                weight_name = 'Weights',
                                seed = 2025,
                                var_type = 'binary')
  
  # fitC <- fit_model(formu = censor_model, pool = obs_data_match, var_type = 'binary', weight_name = 'Weights')
  # fitC <- fit_regularized_model(pool = obs_data_match,
  #                               response_name = get_predictor_response_name(censor_model)$response_name,
  #                               predictor_name = get_predictor_response_name(censor_model)$predictor_name,
  #                               weight_name = 'Weights',
  #                               seed = 2025,
                                # var_type = 'binary')
  fitC <- NULL
  
  fitcov <- lapply(seq(covnames), function(j){
    # fit_model(formu = covmodels[[j]], 
    #           pool = obs_data_match[obs_data_match[[time_name]]>=cov_mintimes[j] ], 
    #           var_type = covtypes[[j]], weight_name = 'Weights')
    fit_regularized_model(pool = obs_data_match, 
                          response_name = get_predictor_response_name(covmodels[[j]])$response_name, 
                          predictor_name = get_predictor_response_name(covmodels[[j]])$predictor_name, 
                          weight_name = 'Weights',
                          seed = 2025,
                          var_type = covtypes[j])
    
  })
  names(fitcov) <- covnames
  
  mod_list <- list('outcome' = fitY, 'censor' = fitC, 'cov' = fitcov)
  
  ## Sample pseudo pool and pseudo trajectory
  
  if(!is.null(npool)){
    ids <- as.data.table(sample(unique(obs_data[[id_name]]), 
                                npool, replace = TRUE))
    ids[, `:=`("sid", 1:npool)]
    colnames(ids) <- c(id_name, "sid")
    obs_data <- obs_data[J(ids), allow.cartesian = TRUE, on = id_name]
    obs_data[, `:=`("id", sid)]
    obs_data <- obs_data[, -'sid']
  }
  
  t <- 0
  ids_unique <- unique(obs_data[[id_name]])
  data_len <- length(ids_unique)
  pool <- obs_data[obs_data[[time_name]] <= t, .SD, .SDcols = c(base_covnames, covnames, time_name)]
  set(pool, j = id_name, value = rep(ids_unique, each = 1))
  setcolorder(pool, c(id_name, time_name, base_covnames, covnames))
  set(pool, j = "eligible_pt", value = TRUE)
  
  newdf <- pool[pool[[time_name]] == t]
  newdf[, (intervention_name) := intervention[1]]
  
  intervened <- rep(0, times = nrow(newdf))
  intervened <- intervened + 
    (abs(pool[[intervention_name]] - newdf[[intervention_name]]) > 1e-06)
  intervened <- ifelse(newdf$eligible_pt, intervened >= 1, NA)
  set(newdf, j = "intervened", value = intervened)
  
  if (ncol(newdf) > ncol(pool)) {
    pool <- rbind(pool[pool[[time_name]] < t], newdf, 
                  fill = TRUE)
    pool <- pool[order(get(id_name), get(time_name))]
  }else {
    pool[pool[[time_name]] == t] <- newdf
  }
  
  if(!is.null(histvals)){
    for(t_hist in 0:K){
      mylagged(pool, histvars = histvars, histvals = histvals, time_name = time_name, t = t_hist, id_name = id_name)
    }
  }
  newdf <- pool[pool[[time_name]] == t]
  set(newdf, j = "Py", 
      value = 
        predict_regularized_model(model = fitY, newx = newdf, 
                                  predictor_name = get_predictor_response_name(ymodel)$predictor_name))
      # value = stats::predict(fitY, type = "response", newdata = newdf))
  set(newdf, j = "Y", value = stats::rbinom(data_len, 1, newdf$Py))
  set(newdf, j = "D", value = 0)
  set(newdf, j = "prodp1", value = newdf$Py)
  set(newdf[newdf$D == 1], j = "Y", value = NA)
  set(newdf, j = "prodp0", value = 1 - newdf$Py)
  set(newdf, j = "poprisk", value = newdf$prodp1)
  pool <- rbind(pool[pool[[time_name]] < t], newdf, 
                fill = TRUE)
  pool <- pool[order(get(id_name), get(time_name))]
  col_types <- sapply(pool, class)
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  for(t in 1:(K-1)){
    newdf <- pool[pool[[time_name]] == t - 1]
    set(newdf, j = time_name, value = rep(t, data_len))
    pool <- rbind(newdf, pool)
    
    if(!is.null(histvals)){
      for(t_hist in 0:K){
        mylagged(pool, histvars = histvars, histvals = histvals, time_name = time_name, t = t_hist, id_name = id_name)
      }
    }
    newdf <- pool[pool[[time_name]] == t]
    
    for (i in seq_along(covnames)) {
      cast <- get(paste0("as.", unname(col_types[covnames[i]])))
      if (covtypes[i] == "binary") {
        set(newdf, j = covnames[i], 
            value = cast(stats::rbinom(data_len, 1, 
                                       predict_regularized_model(model = fitcov[[i]], newx = newdf, 
                                                                 predictor_name = get_predictor_response_name(covmodels[[i]])$predictor_name))))
                                       # stats::predict(fitcov[[i]], type = "response", newdata = newdf))))
      }else if (covtypes[i] == "normal") {
        set(newdf, j = covnames[i], 
            value = cast(stats::rnorm(data_len, mean = stats::predict(fitcov[[i]], type = "response", newdata = newdf), sd = fitcov[[i]]$rmse)))
      }else if (covtypes[i] == "categorical") {
        set(newdf, j = covnames[i], 
            value = cast(stats::predict(fitcov[[i]], type = "class", newdata = newdf)))
      }
      pool[pool[[time_name]] == t] <- newdf
      if (covnames[i] %in% histvars) {
        ind <- unlist(lapply(histvars, FUN = function(x) {
          covnames[i] %in% x
        }))
        mylagged(pool, histvars = rep(histvars, ind), histvals = 1, time_name = 't0', t = t, id_name = 'id')
        newdf <- pool[pool[[time_name]] == t]
      }
    }
    
    newdf <- pool[pool[[time_name]] == t]
    newdf[, (intervention_name) := intervention[t+1]]
    
    intervened <- rep(0, times = nrow(newdf))
    intervened <- intervened + 
      (abs(pool[[intervention_name]] - newdf[[intervention_name]]) > 1e-06)
    intervened <- ifelse(newdf$eligible_pt, intervened >= 1, NA)
    set(newdf, j = "intervened", value = intervened)
    pool[pool[[time_name]] == t] <- newdf
    
    set(newdf, j = "D", value = 0)
    set(newdf, j = "Py", 
        value = 
          predict_regularized_model(model = fitY, newx = newdf, 
                                    predictor_name = get_predictor_response_name(ymodel)$predictor_name))
        # value = stats::predict(fitY, 
        #                                         type = "response", newdata = newdf))
    set(newdf, j = "Y", value = stats::rbinom(data_len, 
                                              1, newdf$Py))
    newdf[newdf$D == 1, `:=`("Y", NA)]
    set(newdf, j = "prodp1", value = newdf$Py * 
          pool[pool[[time_name]] == t - 1, ]$prodp0)
    set(newdf, j = "prodp0", value = (1 - newdf$Py) * 
          pool[pool[[time_name]] == t - 1, ]$prodp0)
    set(newdf, j = "poprisk", value = pool[pool[[time_name]] == 
                                             t - 1, ]$poprisk + newdf$prodp1)
    
    pool[pool[[time_name]] == t] <- newdf
    colnames(pool)[colnames(pool) == time_name] <- "t0"
    colnames(pool)[colnames(pool) == id_name] <- "id"
    setorder(pool, id, t0)
    colnames(pool)[colnames(pool) == "t0"] <- time_name
    colnames(pool)[colnames(pool) == "id"] <- id_name
    pool <- pool[pool[[time_name]] >= 0]
    pool[, `:=`("survival", 1 - pool$poprisk)]
  }
  risks <- tapply(pool$poprisk, pool[[time_name]], FUN = mean)
  return(list(risks = risks, mod_list = mod_list))
}
gform_iter_complete_inner <- function(obs_data, K, 
                                      time_name = 't0', id_name = 'id', 
                                      outcome_name, ymodel, outcome_mintime = 0, 
                                      intervention_name, 
                                      intervention = rep(0, K),
                                      histvars = NULL, histvals = NULL,
                                      seed = 2024, ...){
  obs_data = copy(obs_data)
  obs_data <- obs_data[obs_data[[time_name]]<=K-1]
  
  ymodel <- as.formula(paste0('YPred ~', as.character(ymodel)[3]))
  mod_list_complete <- list()
  
  obs_data[, YPred := as.numeric(NA)]
  
  t <- K - 1
  obs_data[obs_data[[time_name]] == t, YPred := as.numeric(Y)]
  yfit <- fit_regularized_model(pool = obs_data[obs_data[[time_name]] == t & obs_data[[censor_name]] == 0], 
                                response_name = get_predictor_response_name(ymodel)$response_name, 
                                predictor_name = get_predictor_response_name(ymodel)$predictor_name, 
                                weight_name = NULL,
                                seed = 2025,
                                var_type = 'binary')
  # yfit <- glm(ymodel, 
  #             data = obs_data[obs_data[[time_name]] == t & obs_data[[censor_name]] == 0], 
  #             family = binomial())
  tmpdf <- obs_data[obs_data[[time_name]] == t]
  tmpdf[, (intervention_name) := intervention[t+1]]
  tmpdf[tmpdf[[time_name]] == t, YPred := 
          predict_regularized_model(model = yfit, newx = tmpdf, 
                                    predictor_name = get_predictor_response_name(ymodel)$predictor_name)]
  
  
  obs_data[obs_data[[time_name]] == t - 1,
           YPred := tmpdf$YPred[match(get(id_name), tmpdf[[id_name]])]]
  obs_data[obs_data[[time_name]] == t - 1 & obs_data[[outcome_name]] == 1, YPred := 1]
  mod_list_complete[[t+1]] <- yfit
  
  # return(tmpdf)
  
  if(K == 1){
    return(list(risk = mean(tmpdf$YPred), mod_list = mod_list_complete))
  }
  ### better not use t
  for(t in (K-2):0){
    # yfit <- glm(ymodel, 
    #             data = obs_data[obs_data[[time_name]] == t & obs_data[[censor_name]] == 0], 
    #             family = binomial())
    yfit <- fit_regularized_model(pool = obs_data[obs_data[[time_name]] == t & obs_data[[censor_name]] == 0], 
                                  response_name = get_predictor_response_name(ymodel)$response_name, 
                                  predictor_name = get_predictor_response_name(ymodel)$predictor_name, 
                                  weight_name = NULL,
                                  seed = 2025,
                                  var_type = 'binary')
    # yfit <- glm(ymodel, 
    #             data = obs_data[obs_data[[time_name]] == t & obs_data[[censor_name]] == 0], 
    #             family = quasibinomial())
    
    tmpdf <- obs_data[obs_data[[time_name]] == t]
    tmpdf[, (intervention_name) := intervention[t+1]]
    tmpdf[tmpdf[[time_name]] == t, YPred := 
            predict_regularized_model(model = yfit, newx = tmpdf, 
                                      predictor_name = get_predictor_response_name(ymodel)$predictor_name)]
            # predict(yfit, tmpdf, type = 'response')]
    
    obs_data[obs_data[[time_name]] == t - 1,
             YPred := tmpdf$YPred[match(get(id_name), tmpdf[[id_name]])]]
    obs_data[obs_data[[time_name]] == t - 1 & obs_data[[outcome_name]] == 1, YPred := 1]
    mod_list_complete[[t+1]] <- yfit
  }
  
  return(list(risk = mean(tmpdf$YPred), mod_list = mod_list_complete))
}

gform_iter_complete <- function(obs_data, K, 
                                time_name = 't0', id_name = 'id', 
                                outcome_name, ymodel, outcome_mintime = 0, 
                                intervention_name, 
                                intervention = rep(0, K),
                                histvars = NULL, histvals = NULL,
                                seed = 2024, ...){
  obs_data <- copy(obs_data)
  rst <- lapply(1:K, function(k){
    gform_iter_complete_inner(obs_data, K = k, 
                              time_name = time_name, id_name = id_name,
                              outcome_name, ymodel, outcome_mintime = 0, 
                              intervention_name = intervention_name, 
                              intervention = intervention[1:k],
                              histvars = histvars, histvals = histvals,
                              seed = 2024)
  })
  
  risks <- sapply(1:K, function(k)rst[[k]]$risk)
  mod_list <- lapply(1:K, function(k)rst[[k]]$mod_list)
  return(list(risks = risks, mod_list = mod_list))
  
}

gform_iter_match_inner <- function(obs_data, K, 
                                   time_name = 't0', id_name = 'id', 
                                   outcome_name, ymodel, outcome_mintime = 0, 
                                   intervention_name, 
                                   intervention = rep(0, K),
                                   histvars = NULL, histvals = NULL,
                                   seed = 2024, ...){
  obs_data <- copy(obs_data)
  obs_data <- as.data.table(obs_data)
  ymodel <- as.formula(paste0('YPred ~', as.character(ymodel)[3]))
  
  mod_list_match <- list()
  t <- K-1
  obs_data[, YPred := as.numeric(NA)]
  tmpdf <- obs_data[obs_data[[time_name]] == t & obs_data[[censor_name]] == 0]
  tmpdf[, YPred := as.numeric(Y)]
  # yfit <- glm(ymodel, weights = Weights, 
  #             data = tmpdf, 
  #             family = quasibinomial())
  yfit <- fit_regularized_model(pool = tmpdf, 
                        response_name = get_predictor_response_name(ymodel)$response_name, 
                        predictor_name = get_predictor_response_name(ymodel)$predictor_name, 
                        weight_name = 'Weights',
                        seed = 2025,
                        var_type = 'binary')
  
  tmpdf <- obs_data[obs_data[[time_name]] == t]
  tmpdf[, (intervention_name) := intervention[t+1]]
  tmpdf[tmpdf[[time_name]] == t, 
        YPred := predict_regularized_model(model = yfit, newx = tmpdf, 
                                  predictor_name = get_predictor_response_name(ymodel)$predictor_name)]
        # YPred := predict_new_level(yfit, tmpdf, type = 'response')]
  
  obs_data[obs_data[[time_name]] == t - 1,
           YPred := tmpdf$YPred[match(get(id_name), tmpdf[[id_name]])]]
  obs_data[obs_data[[time_name]] == t - 1 & obs_data[[outcome_name]] == 1, YPred := 1]
  mod_list_match[[t+1]] <- yfit
  
  if(K == 1){
    return(list(risk = sum(tmpdf$YPred * tmpdf$Weights, na.rm = TRUE) / 
                  sum(tmpdf$Weights, na.rm = TRUE), 
                mod_list = mod_list_match))
  }
  
  for(t in (K-2):0){
    tmpdf <- obs_data[obs_data[[time_name]] == t & obs_data[[censor_name]] == 0]
    # yfit <- glm(ymodel, weights = Weights, 
    #             data = tmpdf, 
    #             family = quasibinomial())
    yfit <- fit_regularized_model(pool = tmpdf, 
                                  response_name = get_predictor_response_name(ymodel)$response_name, 
                                  predictor_name = get_predictor_response_name(ymodel)$predictor_name, 
                                  weight_name = 'Weights',
                                  seed = 2025,
                                  var_type = 'binary')
    
    tmpdf <- obs_data[obs_data[[time_name]] == t]
    tmpdf[, (intervention_name) := intervention[t+1]]
    tmpdf[tmpdf[[time_name]] == t, 
          YPred := 
            predict_regularized_model(model = yfit, newx = tmpdf, 
                                      predictor_name = get_predictor_response_name(ymodel)$predictor_name)]
          # YPred := predict_new_level(yfit, tmpdf, type = 'response')]
    
    obs_data[obs_data[[time_name]] == t - 1,
             YPred := tmpdf$YPred[match(get(id_name), tmpdf[[id_name]])]]
    obs_data[obs_data[[time_name]] == t - 1 & obs_data[[outcome_name]] == 1, YPred := 1]
    mod_list_match[[t+1]] <- yfit
  }
  
  return(list(risk = sum(tmpdf$YPred * tmpdf$Weights, na.rm = TRUE) / 
                sum(tmpdf$Weights, na.rm = TRUE), 
              mod_list = mod_list_match))
  
}

gform_iter_match <- function(obs_data, K, J, Js = NULL, 
                             time_name = 't0', id_name = 'id', 
                             outcome_name, ymodel, outcome_mintime = 0, 
                             base_covnames = NULL,
                             intervention_name, 
                             intervention = rep(0, K),
                             histvars = NULL, histvals = NULL,
                             seed = 2024, ...){
  obs_data <- copy(obs_data)
  obs_data <- as.data.table(obs_data)
  ## matching 
  if(!is.null(seed)){
    set.seed(seed)
  }
  for(t in (K:1)-1){
    
    if(!is.null(Js)){
      J <- Js[t+1]
    }
    
    id_cases <- obs_data[obs_data[[outcome_name]] == 1 & obs_data[[time_name]] == t][[id_name]]
    id_controls <- obs_data[obs_data[[outcome_name]] == 0 & obs_data[[time_name]] == t][[id_name]]
    # id_controls <- obs_data[(obs_data[[outcome_name]] == 0 | obs_data[[censor_name]] == 1) & obs_data[[time_name]] == t][[id_name]]
    
    n_cases <- length(id_cases)
    n_controls <- length(id_controls)
    
    if(is.null(base_covnames)){
      id_matched_controls <- sample(x = id_controls, size = J * n_cases, replace = TRUE)
    } else{
      num_levels <- 1
      while(num_levels == 1){
        id_matched_controls <- sample(x = id_controls, size = J * n_cases, replace = TRUE)
        num_levels <- sapply(base_covnames, function(v){
          length(unique(obs_data[obs_data[[id_name]] %in% c(id_cases, id_matched_controls)][[v]]))
        })
        num_levels <- min(num_levels)
      }
    }
    
    weight_cases <- n_cases / (n_cases + n_controls)
    weight_controls <- n_controls / (n_cases + n_controls) / J
    
    # If there are controls selected multiple time,
    # the weights are summed up.
    dt_control_weights <- data.table(id=id_matched_controls, Weights = weight_controls)
    dt_control_weights <- dt_control_weights[, Count := .N, by = id]
    dt_control_weights <- unique(dt_control_weights)
    
    obs_data[obs_data[[time_name]] == t, Weights := as.numeric(NA)]
    obs_data[obs_data[[outcome_name]] == 1 & obs_data[[time_name]] == t, 
             `:=`(Weights = weight_cases,
                  Count = 1)]
    obs_data[(obs_data[[outcome_name]] == 0 | 
                obs_data[[censor_name]] == 1) & obs_data[[time_name]] == t, 
             `:=`(Weights = dt_control_weights$Weights[match(get(id_name), dt_control_weights[['id']])], 
                  Count = dt_control_weights$Count[match(get(id_name), dt_control_weights[['id']])])]
  }
  match_data <- obs_data[, if (any(!is.na(Count) & Count != 0)) .SD, by = get(id_name)]
  match_data <- match_data[rep(1:.N, ifelse(is.na(Count), 1, Count))]
  
  rst <- lapply(1:K, function(k){
    gform_iter_match_inner(obs_data = match_data, K = k, 
                           time_name = time_name, id_name = id_name,
                           outcome_name, ymodel, outcome_mintime = 0, 
                           intervention_name = intervention_name, 
                           intervention = intervention[1:k],
                           histvars = histvars, histvals = histvals,
                           seed = seed)
  })
  tmpdf_time0 <- obs_data[obs_data[[time_name]] == 0, ]
  
  ### changes
  risks <- sapply(1:K, function(k)rst[[k]]$risk)
  # risks <- sapply(1:K, function(k){
  #   mean(predict(rst[[k]]$mod_list[[1]], tmpdf_time0, type = 'response'))
  # })
  mod_list <- lapply(1:K, function(k)rst[[k]]$mod_list)
  return(list(risks = risks, mod_list = mod_list))
}

