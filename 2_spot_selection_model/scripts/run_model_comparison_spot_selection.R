################################################################################
#
# Title: Run and compare spot selection models
#
# Authors: Schakowski, A.
#
# Last updated: 27/07/2024 (DD/MM/YYYY)
#
################################################################################
run_model_comparison_spot_selection <- function(){
  
################################################################################
# helper function

# fit models
model_fit <- function(model_location, save_location, data_list, warmup, iterations, chains, fixed_param){
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
    threads_per_chain = round(32/chains),
    fixed_param = fixed_param
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

################################################################################
# load stan data

# fit models
model_locations = c("stan/spot_selection_models/a_random.stan",
                    "stan/spot_selection_models/b_no_personal.stan",
                    "stan/spot_selection_models/c_no_social.stan",
                    "stan/spot_selection_models/d_no_knowledge.stan",
                    "stan/spot_selection_models/e_no_cond_personal.stan",
                    "stan/spot_selection_models/f_no_cond_social.stan",
                    "stan/spot_selection_models/g_full_model.stan",
                    "stan/spot_selection_models/e_spot_selection_interaction.stan",
                    "stan/spot_selection_models/spot_selection_confidence.stan",
                    "stan/spot_selection_models/spot_selection_knowledge.stan",
                    "stan/spot_selection_models/spot_selection_confidence_v2.stan")

save_locations = c("data/processed_data/model_fit/a_random_fit.rds",
                   "data/processed_data/model_fit/b_no_personal_fit.rds",
                   "data/processed_data/model_fit/c_no_social_fit.rds",
                   "data/processed_data/model_fit/d_no_knowledge_fit.rds",
                   "data/processed_data/model_fit/e_no_cond_personal_fit.rds",
                   "data/processed_data/model_fit/f_no_cond_social_fit.rds",
                   "data/processed_data/model_fit/g_full_model_fit.rds",
                   "data/processed_data/model_fit/e_spot_selection_model_fit.rds",
                   "data/processed_data/model_fit/h_spot_selection_confidence_fit.rds",
                   "data/processed_data/model_fit/i_spot_selection_knowledge_fit.rds",
                   "data/processed_data/model_fit/j_spot_selection_confidence_v2_fit.rds")
for (i in 1:length(model_locations)){

  if (i == 1){
    
    stan_data = readRDS(file = "data/processed_data/stan_data_spot_selection.rds")
    
    # standardized time since start from 0 to 1
    stan_data$Time_since_start = as.numeric(stan_data$Time_since_start)
    #stan_data$Time_since_start = stan_data$Time_since_start / max(stan_data$Time_since_start)
    
    # any 0s in dataset? If so, replace with small value
    stan_data$Distances_social[stan_data$Distances_social == 0] <- .01
    stan_data$Distances_roughness[stan_data$Distances_roughness == 0] <- .01
    stan_data$Distances_success[stan_data$Distances_success == 0] <- .01
    stan_data$Distances_loss[stan_data$Distances_loss == 0] <- .01
    stan_data$Locality[stan_data$Locality == 0] <- .01
    min(stan_data$Time_since_start)
    max(stan_data$Time_since_start)
    sum(is.na(stan_data$Time_since_start))
    sum(is.na(stan_data$Distances_ecological_knowledge))
    
    # any missings in dataset?
    stan_data$Distances_social[is.na(stan_data$Distances_social)]
    stan_data$Distances_roughness[is.na(stan_data$Distances_roughness)]
    stan_data$Distances_success[is.na(stan_data$Distances_success)]
    stan_data$Distances_loss[is.na(stan_data$Distances_loss)]
    stan_data$Locality[is.na(stan_data$Locality)]
    
    #stan_data = stan_data[-22]
    stan_data$spot_ids = as.numeric(stan_data$spot_ids)
    
    
    model_fit(model_location = model_locations[i],
              save_location = save_locations[i],
              data_list = stan_data,
              warmup = 2000,
              iterations = 2000,
              chains = 4,
              fixed_param = T)
  } else if (i == 10){
    
    stan_data = readRDS(file = "data/processed_data/stan_data_spot_selection_prior_exposure.rds")
    
    # standardized time since start from 0 to 1
    stan_data$Time_since_start = as.numeric(stan_data$Time_since_start)
    #stan_data$Time_since_start = stan_data$Time_since_start / max(stan_data$Time_since_start)
    
    # any 0s in dataset? If so, replace with small value
    stan_data$Distances_social[stan_data$Distances_social == 0] <- .01
    stan_data$Distances_roughness[stan_data$Distances_roughness == 0] <- .01
    stan_data$Distances_success[stan_data$Distances_success == 0] <- .01
    stan_data$Distances_loss[stan_data$Distances_loss == 0] <- .01
    stan_data$Locality[stan_data$Locality == 0] <- .01
    min(stan_data$Time_since_start)
    max(stan_data$Time_since_start)
    sum(is.na(stan_data$Time_since_start))
    sum(is.na(stan_data$Distances_ecological_knowledge))
    
    # any missings in dataset?
    stan_data$Distances_social[is.na(stan_data$Distances_social)]
    stan_data$Distances_roughness[is.na(stan_data$Distances_roughness)]
    stan_data$Distances_success[is.na(stan_data$Distances_success)]
    stan_data$Distances_loss[is.na(stan_data$Distances_loss)]
    stan_data$Locality[is.na(stan_data$Locality)]
    
    #stan_data = stan_data[-22]
    stan_data$spot_ids = as.numeric(stan_data$spot_ids)
    
    model_fit(model_location = model_locations[i],
              save_location = save_locations[i],
              data_list = stan_data,
              warmup = 50,
              iterations = 50,
              chains = 4,
              fixed_param = F)
  } else {
    
    stan_data = readRDS(file = "data/processed_data/stan_data_spot_selection.rds")
    
    # standardized time since start from 0 to 1
    stan_data$Time_since_start = as.numeric(stan_data$Time_since_start)
    #stan_data$Time_since_start = stan_data$Time_since_start / max(stan_data$Time_since_start)
    
    # any 0s in dataset? If so, replace with small value
    stan_data$Distances_social[stan_data$Distances_social == 0] <- .01
    stan_data$Distances_roughness[stan_data$Distances_roughness == 0] <- .01
    stan_data$Distances_success[stan_data$Distances_success == 0] <- .01
    stan_data$Distances_loss[stan_data$Distances_loss == 0] <- .01
    stan_data$Locality[stan_data$Locality == 0] <- .01
    min(stan_data$Time_since_start)
    max(stan_data$Time_since_start)
    sum(is.na(stan_data$Time_since_start))
    sum(is.na(stan_data$Distances_ecological_knowledge))
    
    # any missings in dataset?
    stan_data$Distances_social[is.na(stan_data$Distances_social)]
    stan_data$Distances_roughness[is.na(stan_data$Distances_roughness)]
    stan_data$Distances_success[is.na(stan_data$Distances_success)]
    stan_data$Distances_loss[is.na(stan_data$Distances_loss)]
    stan_data$Locality[is.na(stan_data$Locality)]
    
    #stan_data = stan_data[-22]
    stan_data$spot_ids = as.numeric(stan_data$spot_ids)
    
    
    model_fit(model_location = model_locations[i],
              save_location = save_locations[i],
              data_list = stan_data,
              warmup = 50,
              iterations = 50,
              chains = 4,
              fixed_param = F)
    
  }
  
}
gc()

model = readRDS(save_locations[11])
print(model, max_rows = 20)
x=model$summary(variables = "weights", ~quantile(.x, probs = c(0.025, 0.975)))
model$code()
plot(x$`2.5%`[3:6], ylim = c(-1.5, 1.5))
points(x$`97.5%`[3:6])
abline(a = 0, b = 0)
exp(2.7)

dist = seq(from = 0, to = 250, by = 1)
k = exp(-dist^2 / (2 * 90^2))
k2 = exp(-dist^2 / (2 * 14^2))
plot(k~dist)
points(k2~dist)
ecological_knowledge_polygons = readRDS("data/ecological_knowledge_polygons.rds")
table(ecological_knowledge_polygons$label)

model %>% 
  spread_draws(v_lake[lake, parameter]) %>% 
  filter(parameter >= 3) %>% 
  ggplot(aes(x = parameter, y = v_lake, group = lake)) + 
  stat_pointinterval(position = position_dodge(.5))


model %>% 
  spread_draws(v_trip[trip, parameter]) %>% 
  ggplot(aes(x = trip, y = v_trip)) + 
  stat_pointinterval() + 
  facet_wrap(~parameter, scales = "free")

ests = model %>% 
  spread_draws(v_trip[trip, parameter]) %>% 
  group_by(trip, parameter) %>% 
  summarize(est = mean(v_trip)) %>% 
  pivot_wider(names_from = parameter, names_prefix = "parameter_", values_from = est)
ests %>% 
  ggplot(aes(x = parameter_2, y = parameter_4)) + 
  geom_point() + 
  geom_smooth(method = "lm")

# get standardized estimates

#1 social
#2 expectation
#3 success
#4 loss
#5 locality
#6 locality x Time
#d = seq(from = 0, to = 250, by = .5)
#plot(exp(-d^2 / (2 * exp(3.8))^2)~d)

model %>% 
  spread_draws(v_trip[trip, parameter]) %>%
  ggplot(aes(x = trip, y = v_trip)) + 
  stat_pointinterval() + 
  facet_wrap(~parameter, scales = "free")

# run model comparison
fit_locations = save_locations
save_locations = c("data/processed_data/model_fit/comparison/a_random_fit_waic.rds",
                   "data/processed_data/model_fit/comparison/b_no_personal_fit_waic.rds",
                   "data/processed_data/model_fit/comparison/c_no_social_fit_waic.rds",
                   "data/processed_data/model_fit/comparison/d_no_knowledge_fit_waic.rds",
                   "data/processed_data/model_fit/comparison/e_no_cond_personal_fit_waic.rds",
                   "data/processed_data/model_fit/comparison/f_no_cond_social_fit_waic.rds",
                   "data/processed_data/model_fit/comparison/g_full_model_fit_waic.rds",
                   "data/processed_data/model_fit/comparison/e_spot_selection_model_fit_waic.rds")

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
                readRDS(save_locations[6]),
                readRDS(save_locations[7]),
                readRDS(save_locations[8]))
comparison = loo_compare(loo_list)
plot(comparison)
# save results
saveRDS(comparison, file = "data/processed_data/model_fit/comparison/comparison.rds")
print(readRDS(fit_locations[5]), max_rows = 20)
model$code()
}


################################################################################
# END
################################################################################
comparison = readRDS(comparison, file = "data/processed_data/model_fit/comparison/comparison.rds")
comparison = as.data.frame(comparison)
comparison$name = c("+Conditional\nSocial",
                    "+Conditional\nPersonal",
                    "+Knowledge",
                    "+Social",
                    "+Personal",
                    "+Locality",
                    "Random")
comparison$model_id = 1:7
comparison_spot_selection = comparison %>%
  ggplot(aes(x = elpd_diff, y = fct_reorder(name, model_id))) + 
  geom_col(color = "black", width = .6) +
  geom_errorbarh(aes(xmin = elpd_diff - 4*se_diff, xmax = elpd_diff + 4*se_diff), height = .3, size = .8) + 
  ggtitle("Model Comparison Spot-Selection") + 
  xlab(expression(Delta~"ELPD (rel. to best)")) + 
  ylab("Model") + 
  #labs(tag = "a")+
  geom_vline(xintercept = 0)+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7, angle = 65, hjust = .5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
comparison_spot_selection
ggsave(comparison_spot_selection, file = "comparison_spot_selection.png", 
       width = 8, height = 8, bg = "white", dpi = 1200, units = "cm")

################################################################################
#