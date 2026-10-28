################################################################################
#
# Title: Run success models
#
# Description: Run analyses and visualize results
#
# Authors: Schakowski, A.
#
# Last updated: 28/09/2026 (DD/MM/YYYY)
#
################################################################################
run_success_models <- function(){

# load competition results
catch_data <- read_csv("data/catch_data.csv")

# Load successful/unsuccessful spots
spot_selection_data = fread(file = "data/spot_selection_data_cleaned.csv")
spot_selection_data = spot_selection_data %>% group_by(day, year, camera_id) %>% summarize(perc_success = mean(Return),
                                                                                           n_spots = n())
# load model key
model_key_spot_selection = readRDS("data/processed_data/identifiers_spot_selection_model_key.rds")

# add success info
model_key_spot_selection = model_key_spot_selection %>% group_by(day, year, camera_id, id, trip, lake) %>% 
  slice(1) %>% 
  dplyr::select(day, year, camera_id, id, trip,lake)
model_key_spot_selection$day = as.numeric(model_key_spot_selection$day)
model_key_spot_selection$year = as.numeric(model_key_spot_selection$year)
catch_results_trip = left_join(model_key_spot_selection, spot_selection_data)
catch_results_trip = left_join(catch_results_trip, catch_data)

################################################################################
# 1: Spot selection strategies and success
################################################################################
# load model
spot_selection_model_fit = readRDS("data/processed_data/model_fit/spot_selection_model/g_full_model_fit.rds")

# load parameters
main_effects_spot_selection = spot_selection_model_fit %>% 
  spread_draws(weights[parameter]) %>% 
  group_by(parameter) %>% 
  summarize(weights = mean(weights))
lake_effects_spot_selection = spot_selection_model_fit %>% 
  spread_draws(v_lake[lake, parameter]) %>% 
  group_by(lake, parameter) %>% 
  summarize(v_lake = mean(v_lake))
trip_effects_spot_selection = spot_selection_model_fit %>% 
  spread_draws(v_trip[trip, parameter]) %>% 
  group_by(trip, parameter) %>% 
  summarize(v_trip = mean(v_trip))
id_effects_spot_selection = spot_selection_model_fit %>% 
  spread_draws(v_id[id, parameter]) %>% 
  group_by(id, parameter) %>% 
  summarize(v_id = mean(v_id))

# join to one dataset
spot_selection_parameters = left_join(trip_effects_spot_selection, model_key_spot_selection)
spot_selection_parameters = left_join(spot_selection_parameters, id_effects_spot_selection)
spot_selection_parameters = left_join(spot_selection_parameters, lake_effects_spot_selection)
spot_selection_parameters = left_join(spot_selection_parameters, main_effects_spot_selection)

# add catch results
catch_results_trip = left_join(spot_selection_parameters, catch_results_trip) %>% 
  group_by(trip, parameter) %>% 
  mutate(estimate = v_trip + v_id + v_lake + weights) %>% 
  dplyr::select(-c(v_id, v_lake, v_trip, weights)) %>% 
  pivot_wider(names_from = parameter, values_from = estimate, names_prefix = "parameter_")

# lognormal regression model

# priors
if (!file.exists("data/processed_data/spot_selection_success_prediction_catch.rds")){
priors_success_model1 <- c(prior(normal(0,1), class = b),
                           prior(exponential(1), class = sd))
model_ind_success = brms::brm(data = catch_results_trip, 
                              formula = brmsformula(catch ~ 1 + parameter_3 + parameter_4 + parameter_5 + parameter_6 + 
                                                      (1 + parameter_3 + parameter_4 + parameter_5 + parameter_6 | lake) + 
                                                      (1 + parameter_3 + parameter_4 + parameter_5 + parameter_6 | id)
                                                    , sigma ~ 1,
                                                    hu ~ 1 + parameter_3 + parameter_4 + parameter_5 + parameter_6 + 
                                                      (1 + parameter_3 + parameter_4 + parameter_5 + parameter_6 | lake) + 
                                                      (1 + parameter_3 + parameter_4 + parameter_5 + parameter_6 | id)),
                              family = hurdle_lognormal(),
                              chains = 4,
                              cores = 4,
                              prior = priors_success_model1)
}

if(!file.exists("data/processed_data/spot_selection_success_prediction_perc.rds")){
model_ind_success_perc =brms::brm(data = catch_results_trip, 
                                    formula = brmsformula(perc_success ~ 1 + parameter_3 + parameter_4 + parameter_5 + parameter_6 + 
                                                            (1 + parameter_3 + parameter_4 + parameter_5 + parameter_6 | lake) + 
                                                            (1 + parameter_3 + parameter_4 + parameter_5 + parameter_6 | id)
                                                          , phi ~ 1,# The precision
                                                          zi ~ 1 + parameter_3 + parameter_4 + parameter_5 + parameter_6 + 
                                                            (1 + parameter_3 + parameter_4 + parameter_5 + parameter_6 | lake) + 
                                                            (1 + parameter_3 + parameter_4 + parameter_5 + parameter_6 | id)),# The one-inflated part
                                    family = zero_inflated_beta(),
                                    chains = 4,
                                    cores = 4)
}

################################################################################
# 2: Patch leaving strategies and success
################################################################################

# load model
patch_leaving_model_fit = readRDS("data/processed_data/model_fit/patch_leaving_models/f_patch_leaving_model_prior_spatial_fit.rds")
identifiers_patch_leaving = readRDS("data/processed_data/identifiers_patch_leaving.rds")
model_key_patch_leaving = identifiers_patch_leaving %>% group_by(camera_id, day, year, id = model_id, trip = model_trips, lake = model_lakes) %>% slice(1) %>% dplyr::select(-c(spot_id, bin))
catch_results_trips_patch_leaving = left_join(model_key_patch_leaving, catch_results_trip %>% ungroup() %>% dplyr::select(day, year, camera_id, perc_success))

# load parameters
main_effects_patch_leaving = patch_leaving_model_fit %>% 
  spread_draws(beta[parameter]) %>% 
  group_by(parameter) %>% 
  summarize(beta = mean(beta))
lake_effects_patch_leaving = patch_leaving_model_fit %>% 
  spread_draws(v_lake[lake, parameter]) %>% 
  group_by(lake, parameter) %>% 
  summarize(v_lake = mean(v_lake))
trip_effects_patch_leaving = patch_leaving_model_fit %>% 
  spread_draws(v_trip[trip, parameter]) %>% 
  group_by(trip, parameter) %>% 
  summarize(v_trip = mean(v_trip))
id_effects_patch_leaving = patch_leaving_model_fit %>% 
  spread_draws(v_id[id, parameter]) %>% 
  group_by(id, parameter) %>% 
  summarize(v_id = mean(v_id))

# join to one dataset
patch_leaving_parameters = left_join(trip_effects_patch_leaving, catch_results_trips_patch_leaving)
patch_leaving_parameters = left_join(patch_leaving_parameters, id_effects_patch_leaving)
patch_leaving_parameters = left_join(patch_leaving_parameters, lake_effects_patch_leaving)
patch_leaving_parameters = left_join(patch_leaving_parameters, main_effects_patch_leaving)

# add catch results
catch_results_trips_patch_leaving = left_join(patch_leaving_parameters, catch_results_trips_patch_leaving) %>% 
  group_by(trip, parameter) %>% 
  mutate(estimate = v_trip + v_id + v_lake + beta) %>% 
  dplyr::select(-c(v_id, v_lake, v_trip, beta)) %>% 
  pivot_wider(names_from = parameter, values_from = estimate, names_prefix = "parameter_")

# run models 
# lognormal regression model

# priors
priors_success_model1 <- c(prior(normal(0,1), class = b),
                           prior(exponential(1), class = sd))
if (!file.exists("data/processed_data/patch_leaving_success_prediction_catch.rds")){
model_ind_success_patch_leaving = brms::brm(data = catch_results_trips_patch_leaving, 
                                            formula = brmsformula(catch ~ 1 + parameter_1 + parameter_2 + parameter_3 + parameter_4 + parameter_5 + parameter_6 + parameter_7 +
                                                                    (1 + parameter_1 + parameter_2 + parameter_3 + parameter_4 + parameter_5 + parameter_6 + parameter_7 | lake) + 
                                                                    (1 + parameter_1 + parameter_2 + parameter_3 + parameter_4 + parameter_5 + parameter_6 + parameter_7 | id)
                                                                  , sigma ~ 1,
                                                                  hu ~ 1 + parameter_1 + parameter_2 + parameter_3 + parameter_4 + parameter_5 + parameter_6 + parameter_7 + 
                                                                    (1 + parameter_1 + parameter_2 + parameter_3 + parameter_4 + parameter_5 + parameter_6 + parameter_7 | lake) + 
                                                                    (1 + parameter_1 + parameter_2 + parameter_3 + parameter_4 + parameter_5 + parameter_6 + parameter_7 | id)),
                                            family = hurdle_lognormal(),
                                            chains = 4,
                                            cores = 4,
                                            prior = priors_success_model1)
}

priors_success_model2 <- c(prior(normal(0,1), class = b),
                           prior(exponential(1), class = sd))
if (!file.exists("data/processed_data/patch_leaving_success_prediction_perc.rds")){
model_ind_success_perc_patch_leaving =brms::brm(data = catch_results_trip, 
                                                formula = brmsformula(perc_success ~ 1 + parameter_1 + parameter_2 + parameter_3 + parameter_4 + parameter_5 + parameter_6 + parameter_7 + 
                                                                        (1 + parameter_1 + parameter_2 + parameter_3 + parameter_4 + parameter_5 + parameter_6 + parameter_7 | lake) + 
                                                                        (1 + parameter_1 + parameter_2 + parameter_3 + parameter_4 + parameter_5 + parameter_6 + parameter_7 | id)
                                                                      , phi ~ 1,# The precision
                                                                      zi ~ 1 + parameter_1 + parameter_2 + parameter_3 + parameter_4 + parameter_5 + parameter_6 + parameter_7 + 
                                                                        (1 + parameter_1 + parameter_2 + parameter_3 + parameter_4 + parameter_5 + parameter_6 + parameter_7 | lake) + 
                                                                        (1 + parameter_1 + parameter_2 + parameter_3 + parameter_4 + parameter_5 + parameter_6 + parameter_7 | id)),# The one-inflated part
                                                family = zero_inflated_beta(),
                                                chains = 4,
                                                cores = 4,
                                                prior = priors_success_model2)
}
# save models
saveRDS(model_ind_success, file = "data/processed_data/spot_selection_success_prediction_catch.rds")
saveRDS(model_ind_success_perc, file = "data/processed_data/spot_selection_success_prediction_perc.rds")
saveRDS(model_ind_success_patch_leaving, file = "data/processed_data/patch_leaving_success_prediction_catch.rds")
saveRDS(model_ind_success_perc_patch_leaving, file = "data/processed_data/patch_leaving_success_prediction_perc.rds")

# summaries for tables
summary(model_ind_success)
summary(model_ind_success_perc)
summary(model_ind_success_patch_leaving)
summary(model_ind_success_perc_patch_leaving)

}
################################################################################
# END
################################################################################
