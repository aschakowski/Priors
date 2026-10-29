################################################################################
#
# TITLE: Wrapper script for: Competitive ice-fishers rely more on dynamic cues 
# than prior beliefs
# 
# AUTHORS: Schakowski, A. (schakowski@mpib-berlin.mpg.de)
#
# DESCRIPTION: This script loads helper scripts in succession to reproduce the
# analyses reported in the paper.
#
# LAST UPDATED (DD/MM/YYYY): 06/10/2026
#
################################################################################
# set working directory
# on linux server
if (Sys.info()[1] == "Linux"){
  setwd("/mnt/home/schakowski/MPIB-SRT/1136-icefish/private/repo_upload_priors")
} else {
  # on windows server
  setwd("P:/IceFishing/repo_upload_priors")
}

# check
getwd()

# load library and functions
source("functions/library.R")

################################################################################
# set seed
################################################################################
# global seed
set.seed(20260610)

# local seed for stan models
stan_seed = 20260610

################################################################################
# Step 0: Extract data from maps
################################################################################

# run data preparation script; prior areas
source("0_data_preparation/1_extract_ecological_knowledge.R")
extract_ecological_knowledge_data()

# create figure 1
source("0_data_preparation/2_create_figure_1.R")
create_figure_1()

################################################################################
# Step 1: Area nomination model: What features drive nomination?
################################################################################

# simulate alternative nominations
source("1_nomination_model/1_simulate_alternative_nominations.R")
simulate_alternative_nominations()

# extract ecological features of nominated areas
source("1_nomination_model/2_extract_ecological_properties.R")
extract_ecological_properties()

# fit cluster nomination model
source("1_nomination_model/3_model_comparison_area_nomination.R")
fit_nomination_model()

# model feature validity
source("1_nomination_model/4_success_modeling.R")
run_visualize_success_model()

# compute social consensus
source("1_nomination_model/5_consensus_simulation.R")
estimate_consensus()

# make visualizations
source("1_nomination_model/6_visualize_model_results.R")
visualize_nomination_model()

################################################################################
# Step 2: Load movement data and run spot selection model
################################################################################

# simulate alternative locations spot selection data
source("2_spot_selection_model/1_simulate_alternative_spots.R")
simulate_alternative_spots()

# extract ecological properties spot selection data
source("2_spot_selection_model/2_compute_features_spot_selection.R")
extract_features()

# run spot selection models
source("2_spot_selection_model/3_prepare_stan_data_spot_selection.R")
prepare_stan_data()

# run spot selection models
source("2_spot_selection_model/4_run_model_comparison.R")
run_models_spot_selection()

# visualize results
source("2_spot_selection_model/5_visualize_model_results.R")
visualize_model_results()

################################################################################
# Step 3: Load patch leaving data and run patch leaving model
################################################################################

# prepare stan data
source("3_patch_leaving_model/1_prepare_stan_data_patch_leaving.R")
prepare_stan_data_patch_leaving()

# run patch leaving models
source("3_patch_leaving_model/2_run_patch_leaving_models.R")
run_models_patch_leaving()

# visualize patch leaving model
source("3_patch_leaving_model/3_visualize_patch_leaving_model.R")
visualize_model_results_patch_leaving()

################################################################################
# Step 4: Association between foraging success and ecological knowledge
################################################################################

# run success prediction models
source("4_success_prediction/1_run_success_models.R")
run_success_models()

################################################################################
# Step 5: Descriptive statistics in text
################################################################################
# sample size
catch_data = fread("data/catch_data.csv")
nrow(catch_data)
catch_data %>% 
  group_by(day) %>% 
  summarize(n = n())

survey_data = fread("data/survey_data.csv")
left_join(catch_data, survey_data) %>% 
  summarize(median((years_icefishing), na.rm = T),
            sd(years_icefishing, na.rm = T))
length(unique(catch_data$participant_id))
catch_data %>% 
  group_by(participant_id) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  summarize(mean(n),sd(n))

# table s1
left_join(catch_data, survey_data) %>% 
  group_by(day) %>% 
  summarize(`Mean catch (g)` = mean(catch),
            `SD catch (g)` = sd(catch),
            `Median catch (g)` = median(catch),
            `Missing (N) 1` = sum(is.na(catch)),
            
            `Mean age (y)` = mean(2023-year_of_birth, na.rm = T),
            `SD age (y)` = sd(2023-year_of_birth, na.rm = T),
            `Median age (y)` = median(2022-year_of_birth, na.rm = T),
            
            `Missing (N) 2` = sum(is.na(year_of_birth)),
            
            `Female (%)` = mean(sex == "f", na.rm = T),
            `Missing (N) 3` = sum(is.na(sex)),
            
            `Mean experience (y)` = mean(as.numeric(years_icefishing)[!is.na(years_icefishing)], na.rm = T),
            `SD experience (y)` = sd(as.numeric(years_icefishing)[!is.na(years_icefishing)], na.rm = T),
            `Median experience (y)` = median(as.numeric(years_icefishing)[!is.na(years_icefishing)], na.rm = T),
            `Missing (N) 4` = sum(is.na(years_icefishing)))

# table s2
catch_data %>% 
  group_by(participant_id) %>% 
  summarize(n = n()) %>% 
  group_by(n) %>% 
  summarize(n = n())

# area size
areas = readRDS("data/study_lakes.rds")
size = vector()
for (i in 1:5){
  size[i] = (st_area(areas[[i]]) / (1000*1000))
}
median(size)
range(size)

# descriptive data on foraging trips/spots
spot_selection_data = fread("data/spot_selection_data_cleaned.csv")
nrow(spot_selection_data)
length(unique(spot_selection_data$track_id))
length(unique(spot_selection_data$participant_id))
spot_selection_data %>% 
  group_by(day, camera_id) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  summarize(mean(n), sd(n))
mean(spot_selection_data$angling_time)
mean(spot_selection_data$total_catch>0)

# nominated areas
ecological_knowledge_polygons = readRDS("data/ecological_knowledge_polygons.rds")
nrow(ecological_knowledge_polygons)
ecological_knowledge_polygons %>% 
  mutate(area = st_area(x)) %>%
  summarize(area = median(area) / (1000*1000))
ecological_knowledge_polygons %>% 
  group_by(confidence) %>% 
  summarize(n = n() / length(unique(paste(ecological_knowledge_polygons$camera_id, ecological_knowledge_polygons$day))))
ecological_knowledge_polygons %>% 
  summarize(n = n() / length(unique(paste(ecological_knowledge_polygons$camera_id, ecological_knowledge_polygons$day))))

################################################################################
# Step 6: Model estimates for tables
################################################################################
nomination_model = readRDS("data/processed_data/model_fit/nomination_models/f_roughness_fit.rds")
print(nomination_model$summary(posterior::default_summary_measures()[1],quantiles = ~quantile(.,probs = c(0.025, 0.975)), variables = c("weights", "sigma_lake", "sigma_id", "sigma_trip"), "rhat"), n = 50)
stan_data = readRDS("data/processed_data/stan_data_spot_selection.rds")
stan_data$N_ids
stan_data$N_lakes
stan_data$N_trips

print(spot_selection_model_fit, max_rows= 20)
print(spot_selection_model_fit$summary(posterior::default_summary_measures()[1],quantiles = ~quantile(.,probs = c(0.025, 0.975)), variables = c("weights","lambdas", "sigma_id", "sigma_trip", "sigma_lake")), n = 50)

nomination_model = readRDS("data/processed_data/model_fit/nomination_models/g_confidence_fit.rds")
print(nomination_model$summary(posterior::default_summary_measures()[1],quantiles = ~quantile(.,probs = c(0.025, 0.975)), variables = c("weights"), "rhat"), n = 50)

# explored area models
area_model1 = readRDS("data/processed_data/area_model_1.rds")
area_model2 = readRDS("data/processed_data/area_model_2.rds")
area_model3 = readRDS("data/processed_data/area_model_3.rds")

# validation model
validation_model = readRDS("data/processed_data/validation_model_results.rds")

# spot selection model
spot_selection_fit = readRDS("data/processed_data/model_fit/spot_selection_model/i_temporal_fit.rds")
spot_selection_fit$code()
# 1: locality
# 2: locality * Time
# 3: success
# 4: loss
# 5: social
# 6: prior
# 7: success * Time
# 8: loss * Time
# 9: social * Time
# 10: prior * Time
# 11: prior * Success
# 12: prios * Loss
# 13: prior * social
print(spot_selection_fit$summary(posterior::default_summary_measures()[1],quantiles = ~quantile(.,probs = c(0.025, 0.975)), variables = c("weights","lambdas", "sigma_lake", "sigma_id", "sigma_trip"), "rhat"), n = 100)

# confidence
spot_selection_fit_conf = readRDS("data/processed_data/model_fit/spot_selection_model/h_confidence_fit.rds")
spot_selection_fit_conf$code()
print(spot_selection_fit_conf$summary(posterior::default_summary_measures()[1],quantiles = ~quantile(.,probs = c(0.025, 0.975)), variables = c("weights","lambdas", "sigma_lake", "sigma_id", "sigma_trip"), "rhat"), n = 100)

# spot leaving
patch_leaving_model = readRDS("data/processed_data/model_fit/patch_leaving_models/f_patch_leaving_model_prior_spatial_fit.rds")
patch_leaving_model$code()
print(patch_leaving_model$summary(posterior::default_summary_measures()[1],quantiles = ~quantile(.,probs = c(0.025, 0.975)), variables = c("beta","sigma_lake", "sigma_id", "sigma_trip"), "rhat"), n = 190)

patch_leaving_model_conf = readRDS("data/processed_data/model_fit/patch_leaving_models/g_patch_leaving_model_confidence_fit.rds")
patch_leaving_model_conf$code()
print(patch_leaving_model_conf$summary(posterior::default_summary_measures()[1],quantiles = ~quantile(.,probs = c(0.025, 0.975)), variables = c("beta","sigma_lake", "sigma_id", "sigma_trip"), "rhat"), n = 190)

# success predictions
success_model_ind = readRDS("data/processed_data/success_model_ind.rds")
success_model_coll = readRDS("data/processed_data/success_model_coll.rds")

# success models predictions
success_spot_selection = readRDS("data/processed_data/spot_selection_success_prediction_catch.rds")
success_spot_selection_perc = readRDS("data/processed_data/spot_selection_success_prediction_perc.rds")

success_patch_leaving = readRDS("data/processed_data/patch_leaving_success_prediction_catch.rds")
success_patch_leaving_perc = readRDS("data/processed_data/patch_leaving_success_prediction_perc.rds")

nomination_model_conf = readRDS("data/processed_data/model_fit/nomination_model/e_confidence_fit.rds")
print(nomination_model_conf$summary(posterior::default_summary_measures()[1],quantiles = ~quantile(.,probs = c(0.025, 0.975)), variables = c("v_confidence")), n = 50)
nomination_model_conf$code()

# save models
model_ind_success = readRDS("data/processed_data/spot_selection_success_prediction_catch.rds")
model_ind_success_perc = readRDS("data/processed_data/spot_selection_success_prediction_perc.rds")
model_ind_success_patch_leaving = readRDS("data/processed_data/patch_leaving_success_prediction_catch.rds")
model_ind_success_perc_patch_leaving = readRDS("data/processed_data/patch_leaving_success_prediction_perc.rds")

# summaries for tables
summary(model_ind_success)
summary(model_ind_success_perc)
summary(model_ind_success_patch_leaving)
summary(model_ind_success_perc_patch_leaving)


################################################################################
# END
################################################################################
