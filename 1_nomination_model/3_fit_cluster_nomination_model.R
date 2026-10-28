################################################################################
#
# TITLE: Fit area nomination model
# 
# AUTHORS: Schakowski, A. (schakowski@mpib-berlin.mpg.de)
#
# DESCRIPTION: This script fits a series of different model versions predicting
# which areas were nominated based on their spatial features.
#
# LAST UPDATED (DD/MM/YYYY): 06/10/2026
#
################################################################################

fit_nomination_model <- function(stan_seed){

# helper functions
  
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
    refresh = 5,
    iter_sampling = iterations,
    iter_warmup = warmup,
    threads_per_chain = round(100/chains),
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
  
  # compute loo (data too large to do loo, do
  # waic instead)
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

# load cluster data
clusters = readRDS(file = "data/output/simulated_clusters_eco.rds")
clusters$day = as.numeric(clusters$day)

# load catch data
catch_data = read_csv("data/catch_data.csv")

# load survey data
survey_data = read_csv("data/survey_data.csv")

# join survey and catch data
catch_data = left_join(catch_data, survey_data)

# join id information to cluster data
clusters = left_join(clusters, catch_data)
clusters$nomination = ifelse(clusters$condition == "simulated", 0, 1)

# choice identifiers
clusters$trip = as.numeric(as.factor(paste(clusters$day, clusters$camera_id)))
clusters$id = as.numeric(as.factor(clusters$participant_id))
clusters$lake = as.numeric(as.factor(clusters$day))

# arrange in correct order
clusters = clusters %>% arrange(unique_id)

# chosen location always last. 
#View(clusters %>% filter(unique_id == 1))

# save model data
saveRDS(clusters, file = "data/output/model_data.rds")
model_data = readRDS(file = "data/output/model_data.rds")

################################################################################
# prepare stan data
T = length(unique(model_data$unique_id))
N_clusters = max(model_data$iter, na.rm = T)+1
trip = model_data$trip
N_trips = length(unique(model_data$trip))
id = model_data$id
N_ids = length(unique(model_data$id))
lake = model_data$lake
N_lakes = length(unique(model_data$lake))
confidence = pull(model_data %>% group_by(unique_id) %>% slice(1) %>% dplyr::select(confidence) %>% st_drop_geometry())
roughness = matrix(NA, nrow = N_clusters, ncol = T)
displacement = matrix(NA, nrow = N_clusters, ncol = T)
depth = matrix(NA, nrow = N_clusters, ncol = T)
boundary = matrix(NA, nrow = N_clusters, ncol = T)
shore = matrix(NA, nrow = N_clusters, ncol = T)

for (i in 1:T){
  roughness[,i] = pull(model_data[model_data$unique_id == i, "roughness_scaled"] %>% st_drop_geometry())
  displacement[,i] = pull(model_data[model_data$unique_id == i, "displacement_scaled"] %>% st_drop_geometry())
  depth[,i] = pull(model_data[model_data$unique_id == i, "depth_scaled"] %>% st_drop_geometry())
  boundary[,i] = pull(model_data[model_data$unique_id == i, "boundary_scaled"] %>% st_drop_geometry())
  shore[,i] = pull(model_data[model_data$unique_id == i, "shore_scaled"] %>% st_drop_geometry())
}

choices = 1:T
grainsize = 1

################################################################################
# model locations
model_locations = c("stan/nomination_model/a_nomination_model_random.stan",
                    "stan/nomination_model/b_nomination_model_displacement.stan",
                    "stan/nomination_model/c_nomination_model_boundary.stan",
                    "stan/nomination_model/d_nomination_model_depth.stan",
                    "stan/nomination_model/e_nomination_model_shore.stan",
                    "stan/nomination_model/f_nomination_model_roughness.stan",
                    "stan/nomination_model/g_nomination_model_confidence.stan")
save_locations = c("data/processed_data/model_fit/nomination_models/a_random_fit.rds",
                   "data/processed_data/model_fit/nomination_models/b_displacement_fit.rds",
                   "data/processed_data/model_fit/nomination_models/c_boundary_fit.rds",
                   "data/processed_data/model_fit/nomination_models/d_depth_fit.rds",
                   "data/processed_data/model_fit/nomination_models/e_shore_fit.rds",
                   "data/processed_data/model_fit/nomination_models/f_roughness_fit.rds",
                   "data/processed_data/model_fit/nomination_models/g_confidence_fit.rds")

# load data
stan_data = list(T = T,
                 N_clusters = N_clusters,
                 trip = trip[seq(from = 1, to = T*N_clusters, by = N_clusters)],
                 N_trips = N_trips,
                 id = id[seq(from = 1, to = T*N_clusters, by = N_clusters)],
                 N_ids = N_ids,
                 lake = lake[seq(from = 1, to = T*N_clusters, by = N_clusters)],
                 N_lakes = N_lakes,
                 roughness = roughness,
                 depth = depth,
                 displacement = displacement,
                 boundary = boundary,
                 shore = shore,
                 choices = choices,
                 confidence = confidence,
                 grainsize = grainsize)

# fit model
for (i in 1:length(model_locations)){
  start_time = Sys.time()
  if (i == 1){
    
    model_fit(model_location = model_locations[i],
              save_location = save_locations[i],
              data_list = stan_data,
              warmup = 500,
              iterations = 1000,
              chains = 4,
              fixed_param = T,
              stan_seed = stan_seed)
  } else {
    
    model_fit(model_location = model_locations[i],
              save_location = save_locations[i],
              data_list = stan_data,
              warmup = 500,
              iterations = 1000,
              chains = 4,
              fixed_param = F,
              stan_seed = stan_seed)
    
  }
  end_time = Sys.time()
  time_diff = end_time-start_time
  print(paste("Model", i, "finished in", time_diff))
}

# Model comparison
# run model comparison
fit_locations = save_locations
save_locations = c("data/processed_data/model_fit/nomination_models/comparison/a_random.rds",
                   "data/processed_data/model_fit/nomination_models/comparison/b_displacement.rds",
                   "data/processed_data/model_fit/nomination_models/comparison/c_boundary.rds",
                   "data/processed_data/model_fit/nomination_models/comparison/d_depth.rds",
                   "data/processed_data/model_fit/nomination_models/comparison/e_shore.rds",
                   "data/processed_data/model_fit/nomination_models/comparison/f_roughness.rds",
                   "data/processed_data/model_fit/nomination_models/comparison/g_confidence.rds")

for (i in 1:length(fit_locations)){
  compute_waic(fit_location = fit_locations[i],
               save_location = save_locations[i])
}

# conduct model comparison with confidence model
loo_list = list(readRDS(save_locations[1]), 
                readRDS(save_locations[2]), 
                readRDS(save_locations[3]), 
                readRDS(save_locations[4]),
                readRDS(save_locations[5]),
                readRDS(save_locations[6]),
                readRDS(save_locations[7]))
comparison = loo_compare(loo_list)

# save results
saveRDS(comparison, file = "data/processed_data/model_fit/nomination_models/comparison/comparison.rds")

}

################################################################################
# END
################################################################################