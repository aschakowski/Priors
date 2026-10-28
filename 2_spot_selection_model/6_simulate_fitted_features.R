################################################################################
#
# TITLE: Simulate fitted features
#
# AUTHORS: Schakowski, A. (schakowski@mpib-berlin.mpg.de)
#
# DESCRIPTION: This script simulates spatial features based on the estimated 
# bandwidth parameters.
#
# LAST UPDATED (DD/MM/YYYY): 10/10/2026
#
################################################################################
simulate_fitted_features<-function(){
  
# load model
spot_selection_fit = readRDS(file = "data/processed_data/model_fit/spot_selection_model/i_temporal_fit.rds")

# load data that includes features for first spots
stan_data = readRDS(file = "data/processed_data/stan_data_spot_selection_first_spots.rds")

# avg bandwidth parameters for each lake and each feature
main_lambdas = spot_selection_fit %>% 
  spread_draws(lambdas[feature]) %>% 
  group_by(feature) %>% summarize(log_lambda = mean(lambdas))
lake_lambdas = spot_selection_fit %>% 
  spread_draws(v_lake[lake,feature]) %>% 
  filter(feature %in% c(13, 14)) %>% 
  mutate(feature = feature - 12) %>% 
  group_by(feature, lake) %>% 
  summarize(log_lambda_offset = mean(v_lake))

# join to one dataframe
lambdas = left_join(main_lambdas, lake_lambdas)
lambdas$lake_lambda = exp(lambdas$log_lambda_offset + lambdas$log_lambda)

# helper
square <- function(x){x^2}

# loop through all choices and compute feature values
options(warn=2)
fitted_features = list()
for (t in 1:stan_data$T){#stan_data$T){
  
  # preallocate results
  success_feature = vector()
  loss_feature = vector()
  social_feature = vector()
  knowledge_feature = vector()
  lake = vector()
  id = vector()
  trip = vector()
  choice = vector()
  spot_id = vector()
  chosen = vector()
  confidence = vector()
  
  # success feature
  if (stan_data$Distances_successes[t,1,1] < 0){
    success_feature = rep(0, stan_data$N_spots)
  } else {
    if (stan_data$success_columns[t] == 1){
      success_feature = exp( - (stan_data$Distances_successes[t,, 1:stan_data$success_columns[t]] ^ 2) /  (2*square(lambdas$lake_lambda[lambdas$lake == stan_data$lake[t] & lambdas$feature == 2])))
    } else {
      success_feature = (exp( - (stan_data$Distances_successes[t,, 1:stan_data$success_columns[t]] ^ 2) /  (2*square(lambdas$lake_lambda[lambdas$lake == stan_data$lake[t] & lambdas$feature == 2]))) %*% rep(1.0, stan_data$success_columns[t]));  
    }
  }
  
  # loss feature
  if (stan_data$Distances_losses[t,1,1] < 0){
    loss_feature = rep(0, stan_data$N_spots);
  } else {
    if (stan_data$loss_columns[t] == 1){
      loss_feature = exp( - (stan_data$Distances_losses[t,, 1:stan_data$loss_columns[t]] ^ 2) /  (2*square(lambdas$lake_lambda[lambdas$lake == stan_data$lake[t] & lambdas$feature == 2])));  
    } else {
      loss_feature = (exp( - (stan_data$Distances_losses[t,, 1:stan_data$loss_columns[t]] ^ 2) /  (2*square(lambdas$lake_lambda[lambdas$lake == stan_data$lake[t] & lambdas$feature == 2]))) %*% rep(1.0, stan_data$loss_columns[t])); 
    }
  }
  
  # knowledge feature
  expectation_feature = stan_data$Distances_ecological_knowledge[,t]/1000;
  
  # social feature
  social_feature = rowSums(exp( - (stan_data$Distances_social[t,,] ^ 2) / (2*square(lambdas$lake_lambda[lambdas$lake == stan_data$lake[t] & lambdas$feature == 1]))) * (1.0 / (1 + exp(-20 * (stan_data$Distances_social[t,,] - 5)))));
  
  # locality feature
  locality_feature = stan_data$Locality[, t]/1000;
  
  # identifiers
  lake = stan_data$lake[t]
  id = stan_data$id[t]
  trip = stan_data$trip[t]
  choice = stan_data$choices[t]
  chosen = c(1, rep(0, stan_data$N_spots-1))
  spot_id = stan_data$spot_id[t]
  
  # time since start
  tss = stan_data$Time_since_start[t]
  
  # confidence
  confidence = stan_data$confidence[,t]
  
  # make tidy dataframe
  fitted_features_df = data.frame(success_feature = success_feature,
                                  loss_feature = loss_feature,
                                  social_feature = social_feature,
                                  locality_feature = locality_feature,
                                  knowledge_feature = expectation_feature,
                                  time_since_start = tss,
                                  lake = lake,
                                  id = id,
                                  trip = trip,
                                  choice = choice,
                                  spot_id = spot_id,
                                  chosen = chosen,
                                  confidence = confidence)
  
  # save to list
  fitted_features[[t]]<-fitted_features_df
  
  # progress
  if (t %% 100 == 0){
    print(paste(t, "/", stan_data$T))
  }
}

# make df from list
fitted_features = do.call(rbind, fitted_features)

# load ids to match back to dataset
identifiers = readRDS(file = "data/processed_data/identifiers_spot_selection_model_key.rds")

spot_selection_data = read.csv("data/spot_selection_data_cleaned.csv")

identifiers$day = as.numeric(identifiers$day)
identifiers$year = as.numeric(identifiers$year)
identifiers$spot_id = as.numeric(identifiers$spot_id)
identifiers = left_join(identifiers %>% dplyr::select(-time_since_start), spot_selection_data %>% dplyr::select(day, year, camera_id, Return, spot_id))

fitted_features = left_join(fitted_features, identifiers)
saveRDS(fitted_features, file = "data/processed_data/fitted_features.rds")
}

################################################################################
# END
################################################################################
