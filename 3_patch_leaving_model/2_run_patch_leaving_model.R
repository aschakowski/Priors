################################################################################
#
# Title: Run and compare patch leaving models
#
# Authors: Schakowski, A.
#
# Last updated: 27/07/2024 (DD/MM/YYYY)
#
################################################################################
run_model_comparison_patch_leaving <- function(){
  
  ################################################################################
  # helper function
  
  # fit models
  model_fit <- function(model_location, save_location, data_list, warmup, iterations, chains, stan_seed){
    # load model
    stan_model = cmdstan_model(model_location, 
                               cpp_options = list(stan_threads = TRUE))
    
    # fit models
    stan_model_fit <- stan_model$sample(
      data = data_list,
      chains = chains,
      parallel_chains = chains,
      refresh = 5,
      iter_sampling = iterations,
      iter_warmup = warmup,
      threads_per_chain = round(100/chains),
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
    
    ll = model_fit$draws(variables = "log_lik")
    waic = waic(ll)
    
    # save loo
    saveRDS(waic, file = save_location)
    
    # keep workspace clear
    gc()
    
    # empty return
    return()
    
  }
  
  ################################################################################
  # load stan data
  stan_data = readRDS(file = "data/processed_data/stan_data_patch_leaving.rds")
  
  # fit models
  model_locations = c("stan/patch_leaving_models/a_patch_leaving_model_baseline.stan",
                      "stan/patch_leaving_models/b_patch_leaving_model_time_wo_catch.stan",
                      "stan/patch_leaving_models/c_patch_leaving_model_patch_discovery.stan",
                      "stan/patch_leaving_models/d_patch_leaving_model_personal_spatial_features.stan",
                      "stan/patch_leaving_models/e_patch_leaving_model_social_spatial_features.stan",
                      "stan/patch_leaving_models/f_patch_leaving_model_prior_spatial_features.stan",
                      "stan/patch_leaving_models/g_patch_leaving_model_confidence.stan")
  save_locations = c("data/processed_data/model_fit/patch_leaving_models/a_patch_leaving_model_baseline_fit.rds",
                     "data/processed_data/model_fit/patch_leaving_models/b_patch_leaving_model_time_wo_catch_fit.rds",
                     "data/processed_data/model_fit/patch_leaving_models/c_patch_leaving_model_patch_discovery_fit.rds",
                     "data/processed_data/model_fit/patch_leaving_models/d_patch_leaving_model_personal_spatial_features_fit.rds",
                     "data/processed_data/model_fit/patch_leaving_models/e_patch_leaving_model_social_spatial_fit.rds",
                     "data/processed_data/model_fit/patch_leaving_models/f_patch_leaving_model_prior_spatial_fit.rds",
                     "data/processed_data/model_fit/patch_leaving_models/g_patch_leaving_model_confidence_fit.rds")
  
  for (i in 1:length(model_locations)){
    model_fit(model_location = model_locations[i],
              save_location = save_locations[i],
              data_list = stan_data,
              warmup = 1000,
              iterations = 2000,
              chains = 4, 
              stan_seed = stan_seed)
  }
  
  # run model comparison
  fit_locations = save_locations
  save_locations = c("data/processed_data/model_fit/patch_leaving_models/comparison/a_loo_baseline.rds",
                     "data/processed_data/model_fit/patch_leaving_models/comparison/b_loo_time_wo_catch.rds",
                     "data/processed_data/model_fit/patch_leaving_models/comparison/c_loo_patch_discovery.rds",
                     "data/processed_data/model_fit/patch_leaving_models/comparison/d_loo_personal.rds",
                     "data/processed_data/model_fit/patch_leaving_models/comparison/e_loo_social.rds",
                     "data/processed_data/model_fit/patch_leaving_models/comparison/f_loo_prior.rds")
  
  for (i in 1:length(fit_locations)){
    compute_waic(fit_location = fit_locations[i],
                 save_location = save_locations[i])
  }
  
  # conduct model comparison
  loo_list = list(readRDS(save_locations[1]), 
                  readRDS(save_locations[2]), 
                  readRDS(save_locations[3]), 
                  readRDS(save_locations[4]), 
                  readRDS(save_locations[5]),
                  readRDS(save_locations[6]))
  
  comparison = loo_compare(loo_list)
  
  # save results
  saveRDS(comparison, file = "data/processed_data/model_fit/patch_leaving_models/comparison/comparison.rds")
  
}
################################################################################
# END
################################################################################