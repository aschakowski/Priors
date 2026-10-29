################################################################################
#
# TITLE: Run Model comparison spot selection
#
# AUTHORS: Schakowski, A. (schakowski@mpib-berlin.mpg.de)
#
# DESCRIPTION: This script takes prepared stan data and runs a series of 
# spot selection models.
#
# LAST UPDATED (DD/MM/YYYY): 10/10/2026
#
################################################################################
run_model_comparison_spot_selection <- function(){
  
  ################################################################################
  # helper function
  
  # fit models
  model_fit <- function(model_location, save_location, data_list, warmup, iterations, chains, fixed_param, stan_seed){
    # load model
    stan_model = cmdstan_model(model_location, 
                               cpp_options = list(stan_threads = TRUE))
    
    # fit models
    stan_model_fit <- stan_model$sample(
      data = data_list,
      chains = chains,
      parallel_chains = chains,
      refresh = 50,
      iter_sampling = iterations,
      iter_warmup = warmup,
      threads_per_chain = round(80/chains),
      fixed_param = fixed_param,
      seed = stan_seed
    )
    
    # save
    stan_model_fit$save_object(file = save_location)
    
    # empty return
    return()
    
  }
  
  # extract log lik
  compute_waic <- function(fit_location, save_location){
    
    # load model fit
    model_fit = readRDS(file = fit_location)
    
    # compute loo
    # models too large for loo, use WAIC instead
    #loo <- model_fit$loo(cores = 40)
    
    ll = model_fit$draws(variables = "log_lik")
    waic = waic(ll)
    
    # save loo
    saveRDS(waic, file = save_location)
    
    # keep workspace clear
    gc()
    
    # empty return
    return()
    
  }
  
  # fit models
  model_locations = c("stan/spot_selection_models/a_random.stan",
                      "stan/spot_selection_models/b_locality.stan",
                      "stan/spot_selection_models/c_personal.stan",
                      "stan/spot_selection_models/d_social.stan",
                      "stan/spot_selection_models/e_prior.stan",
                      "stan/spot_selection_models/f_cond_personal.stan",
                      "stan/spot_selection_models/g_cond_social.stan",
                      "stan/spot_selection_models/h_confidence.stan",
                      "stan/spot_selection_models/i_temporal.stan")
  
  save_locations = c("data/processed_data/model_fit/spot_selection_model/a_random_fit.rds",
                     "data/processed_data/model_fit/spot_selection_model/b_locality_fit.rds",
                     "data/processed_data/model_fit/spot_selection_model/c_personal_fit.rds",
                     "data/processed_data/model_fit/spot_selection_model/d_social_fit.rds",
                     "data/processed_data/model_fit/spot_selection_model/e_prior_fit.rds",
                     "data/processed_data/model_fit/spot_selection_model/f_cond_personal_fit.rds",
                     "data/processed_data/model_fit/spot_selection_model/g_cond_social_fit.rds",
                     "data/processed_data/model_fit/spot_selection_model/h_confidence_fit.rds",
                     "data/processed_data/model_fit/spot_selection_model/i_temporal_fit.rds")
  
  for (i in 1:length(model_locations)){
    
    if (i == 1){
      
      stan_data = readRDS(file = "data/processed_data/stan_data_spot_selection.rds")
      
      # fit model
      model_fit(model_location = model_locations[i],
                save_location = save_locations[i],
                data_list = stan_data,
                warmup = 1000,
                iterations = 2000,
                chains = 4,
                fixed_param = T, 
                stan_seed = stan_seed)
      
    } else {
      
      stan_data = readRDS(file = "data/processed_data/stan_data_spot_selection.rds")
      
      # fit model
      model_fit(model_location = model_locations[i],
                save_location = save_locations[i],
                data_list = stan_data,
                warmup = 2000,
                iterations = 2000,
                chains = 4,
                fixed_param = F,
                stan_seed = stan_seed)
      
    }
    
    gc()
    gc()
    gc()
  }
  
  # run model comparison
  fit_locations = save_locations
  save_locations = c("data/processed_data/model_fit/spot_selection_model/comparison/a_random_fit_waic.rds",
                     "data/processed_data/model_fit/spot_selection_model/comparison/b_locality_fit_waic.rds",
                     "data/processed_data/model_fit/spot_selection_model/comparison/c_personal_fit_waic.rds",
                     "data/processed_data/model_fit/spot_selection_model/comparison/d_social_fit_waic.rds",
                     "data/processed_data/model_fit/spot_selection_model/comparison/e_prior_fit_waic.rds",
                     "data/processed_data/model_fit/spot_selection_model/comparison/f_cond_personal_fit_waic.rds",
                     "data/processed_data/model_fit/spot_selection_model/comparison/g_cond_social_fit_waic.rds",
                     "data/processed_data/model_fit/spot_selection_model/comparison/h_confidence_fit_waic.rds",
                     "data/processed_data/model_fit/spot_selection_model/comparison/i_temporal_fit_waic.rds")
  
  
  
  for (i in 1:length(fit_locations)){
    if (i != 8){
      compute_waic(fit_location = fit_locations[i],
                   save_location = save_locations[i])
    }
  }
  
  # conduct model comparison
  loo_list = list(readRDS(save_locations[1]), 
                  readRDS(save_locations[2]), 
                  readRDS(save_locations[3]), 
                  readRDS(save_locations[4]),
                  readRDS(save_locations[5]),
                  readRDS(save_locations[6]),
                  readRDS(save_locations[7]),
                  readRDS(save_locations[9]))
  comparison = loo_compare(loo_list)
  
  # save results
  saveRDS(comparison, file = "data/processed_data/model_fit/spot_selection_model/comparison/comparison.rds")
}

################################################################################
# END
################################################################################
