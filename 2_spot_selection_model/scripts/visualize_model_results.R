################################################################################
#
# Title: Visualize model results
#
# Description: Run analyses and visualize results
#
# Authors: Schakowski, A.
#
# Last updated: 02/04/2024 (DD/MM/YYYY)
#
################################################################################
visualize_spot_selection_model <- function(){

################################################################################
# load data
################################################################################
# load spatial features
fitted_features = readRDS(file = "data/fitted_features.rds")

spot_selection_model_fit = readRDS("data/processed_data/model_fit/g_full_model_fit.rds")
spot_selection_model_fit %>% 
  spread_draws(v_trip[trip, parameter]) %>% 
  ggplot(aes(x = trip, y = v_trip)) + 
  stat_pointinterval() + 
  facet_wrap(~parameter, scales = "free")
spot_selection_model_fit$code()
# spot data containing all angling spots
spot_selection_data = fread(file = "data/spot_selection_data.csv")
spot_selection_data = spot_selection_data %>% filter(year == 2023) %>% filter(!exclude)

# add unique ID as track identifier
spot_selection_data$unique_ID = paste(spot_selection_data$day, spot_selection_data$year, spot_selection_data$camera_id, sep = "")

# load ecological knowledge
ecological_knowledge = readRDS(file = "data/ecological_knowledge_polygons.rds")
ecological_knowledge = ecological_knowledge %>% filter(day %in% c(1, 2, 3, 4, 5))
################################################################################
# Visualize feature weights
################################################################################
lake_effects = spot_selection_model_fit %>% 
  spread_draws(weights[feature], v_lake[lake, feature]) %>% 
  mutate(lake_effect = weights + v_lake)
lake_effects$feature = lake_effects$feature %>% recode("Locality", "Locality\nxTime", "Successful","Unsuccessful", "Social","Prior", "Prior x\nSuccessful", "Prior x\nUnsuccessful", "Prior x\nSocial")

# compute standardized predictors
scale_features = fitted_features %>% 
  unnest() %>% 
  group_by(trip, lake) %>% 
  mutate(time_since_start = as.numeric(time_since_start) / max(as.numeric(time_since_start))) %>% 
  summarize(sd_social = sd(social_feature),
            sd_success = sd(success_feature),
            sd_loss = sd(loss_feature),
            sd_locality = sd(locality_feature),
            sd_locality_time = sd(locality_feature) * sd(time_since_start),
            sd_prior = sd(knowledge_feature),
            sd_prior_success = sd(knowledge_feature)*sd(success_feature),
            sd_prior_loss = sd(knowledge_feature)*sd(loss_feature),
            sd_prior_social = sd(knowledge_feature)*sd(social_feature)) %>% 
  ungroup() %>% 
  group_by(lake) %>% 
  summarize(sd_social = mean(sd_social),
            sd_prior = mean(sd_prior),
            sd_success = mean(sd_success),
            sd_loss = mean(sd_loss),
            sd_locality = mean(sd_locality),
            sd_locality_time = mean(sd_locality_time),
            sd_prior_success = mean(sd_prior_success),
            sd_prior_loss = mean(sd_prior_loss),
            sd_prior_social = mean(sd_prior_social)) %>% 
  unnest()

# compute standardized predictors
lake_effects = left_join(lake_effects, scale_features)

lake_effects = lake_effects %>% 
  mutate(standardized_weight = ifelse(feature == "Social", lake_effect*sd_social, 
                                      ifelse(feature == "Prior", lake_effect*sd_prior, 
                                             ifelse(feature == "Successful", lake_effect*sd_success, 
                                                    ifelse(feature == "Unsuccessful", lake_effect*sd_loss,
                                                           ifelse(feature == "Prior x\nUnsuccessful", lake_effect*sd_prior_loss,
                                                                  ifelse(feature == "Prior x\nSuccessful", lake_effect*sd_prior_success,
                                                                         ifelse(feature == "Prior x\nSocial", lake_effect*sd_prior_social, 
                                                                                ifelse(feature == "Distance", lake_effect*sd_locality, 0)))))))),
         standardized_main_effect = ifelse(feature == "Social", (weights)*mean(scale_features$sd_social), 
                                           ifelse(feature == "Prior", (weights)*mean(scale_features$sd_prior), 
                                                  ifelse(feature == "Successful", (weights)*mean(scale_features$sd_success), 
                                                         ifelse(feature == "Unsuccessful", weights*mean(scale_features$sd_loss),
                                                                ifelse(feature == "Prior x\nUnsuccessful", weights * mean(scale_features$sd_prior_loss),
                                                                       ifelse(feature == "Prior x\nSuccessful", weights * mean(scale_features$sd_prior_success),
                                                                              ifelse(feature == "Prior x\nSocial", weights*mean(scale_features$sd_prior_social), 
                                                                                     ifelse(feature == "Distance", weights*mean(scale_features$sd_locality), 0)))))))))
# visualize estimates
model_results = lake_effects %>% 
  filter(feature %in% c("Successful", "Unsuccessful", "Prior", "Social")) %>% 
  mutate(name = factor(feature, levels=c("Prior", "Social", "Unsuccessful", "Successful"))) %>%
  ggplot(aes(x = standardized_weight, y = name, group = as.factor(lake))) + 
  geom_vline(xintercept = 0, linetype = "dashed", size = .4) + 
  # .95
  annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful"], .025), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful"], .975), ymin = 3.65, ymax = 4.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful"], .17), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful"], .83), ymin = 3.65, ymax = 4.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful"], .485), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful"], .515), ymin = 3.65, ymax = 4.35,
           alpha = 1,fill = "black")+
  
  annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social"], .025), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social"], .975), ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "#BAE4B3")+
  annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social"], .17), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social"], .83),ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "#31A354")+
  annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social"], .49), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social"], .51),ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "black")+
  
  annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior"], .025), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior"], .975), ymin = 0.65, ymax = 1.35,
           alpha = 1,fill = "#BAE4B3")+
  annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior"], .17), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior"], .83), ymin = 0.65, ymax = 1.35,
           alpha = 1,fill = "#31A354")+
  annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior"], .49), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior"], .51),ymin = 0.65, ymax = 1.35,
           alpha = 1,fill = "black")+
  
  annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful"], .025), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful"], .975), ymin = 2.65, ymax = 3.35,
           alpha = 1,fill = "#BAE4B3")+
  annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful"], .17), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful"], .83), ymin = 2.65, ymax = 3.35,
           alpha = 1,fill = "#31A354")+ 
  annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful"], .49), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful"], .51), ymin = 2.65, ymax = 3.35,
           alpha = 1,fill = "black")+
  ylab("Feature")+
  stat_pointinterval(position = position_dodge(width = .8), color = "black", alpha = .10, size = 0) +
  
  #stat_pointinterval(color = "black", alpha = .15) +
  labs(tag = "D", color = "value")+
  ggtitle("Feature Weights")+
  xlab("Posterior Est. (95%CrI)\nSpot Selection ~ Feature")+
  scale_x_continuous(breaks = seq(from = -3, to = 6, by = 1))+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5, angle = 55, hjust = .5),#, hjust = .5
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
model_results
################################################################################
# Visualize interaction effects
################################################################################
model_results_interaction = lake_effects %>% 
  filter(feature %in% c("Prior x\nUnsuccessful", "Prior x\nSuccessful", "Prior x\nSocial")) %>% 
  mutate(name = factor(feature, levels=c("Prior x\nUnsuccessful", "Prior x\nSuccessful", "Prior x\nSocial"))) %>%
  ggplot(aes(x = standardized_weight, y = name, group = as.factor(lake))) + 
  geom_vline(xintercept = 0, linetype = "dashed", size = .4) + 
  # .95
  annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nUnsuccessful"], .025), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nUnsuccessful"], .975), ymin = .65, ymax = 1.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nUnsuccessful"], .17), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nUnsuccessful"], .83), ymin = .65, ymax = 1.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nUnsuccessful"], .5)-.0005, xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nUnsuccessful"], .5)+.0005, ymin = .65, ymax = 1.35,
           alpha = 1,fill = "black")+
  
  annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSuccessful"], .025), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSuccessful"], .975), ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "#BAE4B3")+
  annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSuccessful"], .17), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSuccessful"], .83),  ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "#31A354")+
  annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSuccessful"], .5)-.0005, xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSuccessful"], .5)+.0005,  ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "black")+
  stat_pointinterval(position = position_dodge(width = .8), color = "black", alpha = .10, size = 0)  +
  
  annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSocial"], .025), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSocial"], .975), ymin = 2.65, ymax = 3.35,
           alpha = 1,fill = "#BAE4B3")+
  annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSocial"], .17), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSocial"], .83),  ymin = 2.65, ymax = 3.35,
           alpha = 1,fill = "#31A354")+
  annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSocial"], .5)-.0005, xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSocial"], .5)+.0005,  ymin = 2.65, ymax = 3.35,
           alpha = 1,fill = "black")+
  stat_pointinterval(position = position_dodge(width = .8), color = "black", alpha = .10, size = 0)  +
  
  
  ylab("Feature")+
  labs(tag = "E", color = "value")+
  ggtitle("Conditional Weights")+
  xlab("Posterior Est. (95%CrI)\nSpot Selection ~ Feature")+
  scale_x_continuous(breaks = seq(from = -2, to = 2, by = .2))+
  scale_y_discrete(labels = c("Unsuccessful x\nPrior", "Successful x\nPrior", "Social x\nPrior"))+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5, angle = 55, hjust = .5),#, hjust = .5
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
model_results_interaction

################################################################################
# visualize model comparison
comparison = readRDS(file = "data/processed_data/model_fit/comparison/comparison.rds")
comparison = as.data.frame(comparison)
comparison$name = c("+Conditional\nSocial",
                    "+Conditional\nPersonal",
                    "+Prior",
                    "+Social",
                    "+Personal",
                    "+Locality",
                    "Random")
comparison$model_id = 1:7
comparison_spot_selection = comparison %>%
  ggplot(aes(x = elpd_diff, y = fct_reorder(name, model_id))) + 
  geom_col(color = "black", width = .6) +
  geom_errorbarh(aes(xmin = elpd_diff - 4*se_diff, xmax = elpd_diff + 4*se_diff), height = .3, size = .8) + 
  ggtitle("Model Comparison") + 
  xlab(expression(Delta~"ELPD (rel. to best)")) + 
  ylab("Model") + 
  labs(tag = "F")+
  geom_vline(xintercept = 0)+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 5, angle = 55, hjust = .5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
comparison_spot_selection

################################################################################
# make panel
# panel
panel_model_results = gridExtra::arrangeGrob(
  grobs = list(model_results,
               model_results_interaction,
               comparison_spot_selection),
  widths = c(1, 1, 1),
  layout_matrix = rbind(c(1, 2, 3))
)

ggsave(panel_model_results, file = "6_plots/panel_spatial_results.png", 
       width = 18, height = 8, bg = "white", dpi = 600, units = "cm")
ggsave(panel_model_results, file = "6_plots/panel_spatial_results.svg", 
       width = 18, height = 8, bg = "white", dpi = 600, units = "cm", device = svglite::svglite)

################################################################################
################################################################################
conf_model = readRDS("data/processed_data/model_fit/j_spot_selection_confidence_v2_fit.rds")
conf_model$code()
plot_s1_b = conf_model %>% 
  spread_draws(weights[parameter]) %>% 
  filter(parameter %in% c(6, 7, 8)) %>% 
  mutate(weights = weights) %>% 
  mutate(parameter = ifelse(parameter == 6, "Sure", ifelse(parameter == 7, "Very Likely", "Likely"))) %>% 
  mutate(parameter = factor(parameter, levels = rev(c("Sure", "Very Likely", "Likely")))) %>% 
  ggplot(aes(y = parameter, x = weights, color = parameter)) + 
  stat_pointinterval(position = position_dodge(width = .8))+
  scale_color_brewer(type = "seq", palette = 3)+
  #labs(tag = "A")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  ggtitle("Confidence in Prior Beliefs")+
  xlab("Posterior Est. (95%CrI)\nSpot Selection ~ Prior | Confidence")+
  scale_x_continuous(breaks = seq(from = -6, to = 6, by = 1))+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5, angle = 65, hjust = .5),#, hjust = .5
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "bottom",
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 5, margin=margin(0,0,0,0)),
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
        text = element_text(family = "sans"),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key.height = unit(.1, "cm"),
        plot.margin=unit(c(.1, .2, .1, .2), "cm"))
plot_s1_b
ggsave(plot_s1_b, file = "6_plots/plot_confidence_spatial.png", 
       width = 8, height = 8, bg = "white", dpi = 300, units = "cm")


################################################################################
# visualize what model is doing: distribution of trip level fw prior
trip_estimates = spot_selection_model_fit %>% 
  spread_draws(v_trip[trip, parameter]) %>% 
  group_by(trip, parameter) %>% 
  summarize(est = mean(v_trip)) %>% 
  pivot_wider(names_from = parameter, names_prefix = "parameter_", values_from = "est")

Feature_matrix_trip_selection = readRDS("data/processed_data/feature_matrix_trip_selection.rds")
identifiers_trip_selection = readRDS("data/processed_data/identifiers_trip_selection.rds")
identifiers_spot_selection = readRDS("data/processed_data/identifiers_spot_selection.rds")
ecological_knowledge_distances_trip_selection = readRDS(file = "data/processed_data/knowledge_distance_matrix_trip_selection.rds")

res = matrix(nrow = 251, ncol = 1)
for (j in 1){
  eco = vector()
  which_eco = vector()
  max_eco = vector()
  social = vector()
  #max_social = vector()
  for (i in 1:length(Feature_matrix_trip_selection)){
    
    eco[i] = min(ecological_knowledge_distances_trip_selection[[i]][j,])
    which_eco[i] = which.min(ecological_knowledge_distances_trip_selection[[i]][j,])
    max_eco[i] = length(ecological_knowledge_distances_trip_selection[[i]][j,])
    
  }
  
  identifiers_trip_selection$eco_knowledge = eco
  identifiers_trip_selection$which_eco_knowledge = which_eco
  identifiers_trip_selection$max_eco_knowledge = max_eco
  
  # descriptives: how many trips use priors?
  
  c=1
  for (epsilon in 100){
    
    use = (identifiers_trip_selection %>% 
             arrange(as.numeric(unique_trip_id)) %>% 
             mutate(use = ifelse(eco_knowledge <= epsilon, 1, 0)) %>% 
             group_by(as.numeric(unique_trip_id)) %>% 
             summarize(use = mean(use)))$use
    
    res[, j] = use
    c=c+1
  }
  #points(res~seq(from = 0, to = 500, by = 5))
}

# compare users vs parameter estimates
trip_estimates$prior_use = res[,1]

regression_model1 = brm(parameter_3 ~ parameter_6,
                        data = trip_estimates, 
                        family = "gaussian")
regression_model2 = brm(parameter_4 ~ parameter_6,
                        data = trip_estimates, 
                        family = "gaussian")

identifiers_spot_selection$day = as.integer(identifiers_spot_selection$day)
identifiers_spot_selection$year = as.integer(identifiers_spot_selection$year)
spot_selection_data = left_join(spot_selection_data, identifiers_spot_selection %>% group_by(day, camera_id, year) %>% slice(1) %>% dplyr::select(trip_id, day, year, ID = camera_id, participant_id))
spot_selection_data = left_join(spot_selection_data, trip_estimates %>% dplyr::select(trip_id = trip, prior_use = parameter_6, loss_attraction = parameter_4, social_attraction = parameter_5))

areas_data = spot_selection_data %>% 
  group_by(day, year, trip_id, participant_id) %>% 
  mutate(n_spots = n()) %>% 
  filter(n_spots > 2) %>% 
  summarize(area = (cbind(x[chull(x,y)], y[chull(x, y)])),
            prior_use = mean(prior_use),
            loss_attraction = mean(loss_attraction),
            social_attraction = mean(social_attraction),
            n_spots = n_spots[1]) %>% 
  group_by(day, year, trip_id, participant_id) %>% 
  mutate(n = n()) %>% 
  filter(n >= 4) %>% 
  summarize(area = Polygon(area)@area,
            prior_use = mean(prior_use),
            loss_attraction = mean(loss_attraction),
            social_attraction = mean(social_attraction),
            n = mean(n_spots)) %>% 
  mutate(area = area / (1000 * 1000))
areas_data$lake = paste(areas_data$day, areas_data$year)

#distance_data = spot_selection_data %>% 
#  group_by(day, year, trip_id) %>% 
#  summarize(distance = sum(step, na.rm = T))

avg_ret = spot_selection_data %>% 
  group_by(day, year, trip_id) %>% 
  summarize(return = mean(Return))
areas_data = left_join(areas_data, avg_ret)
model_area = brms::brm(data = areas_data, 
                            formula = brmsformula(area ~ 1 + prior_use + loss_attraction + return + n +
                                                    (1 + prior_use + loss_attraction + return + n | lake) + 
                                                    1 + prior_use + loss_attraction + return + n +
                                                    (1 + prior_use + loss_attraction + return + n | participant_id)
                                                  , sigma ~ 1),
                            family = lognormal(),
                            chains = 4,
                            cores = 4)
model_area_1 = brms::brm(data = areas_data, 
                       formula = brmsformula(area ~ 1 + prior_use +
                                               (1 + prior_use | lake), sigma ~ 1),
                       family = lognormal(),
                       chains = 4,
                       cores = 4)

model_area_2 = brms::brm(data = areas_data, 
                       formula = brmsformula(area ~ 1 + loss_attraction + 
                                               (1 + loss_attraction | lake), sigma ~ 1),
                       family = lognormal(),
                       chains = 4,
                       cores = 4)

# Plot posterior effects
posterior_predictions = posterior_epred(model_area_1, newdata = data.frame(prior_use = seq(from = -1.5, to = 1.5, by = .05)), re_formula = NA)
posterior_predictions = as.data.frame(posterior_predictions)

posterior_predictions$draw = 1:4000
plot_posterior = posterior_predictions %>% 
  pivot_longer(1:(ncol(posterior_predictions)-1)) %>% 
  group_by(draw) %>% 
  mutate(prior_use = seq(from = -1.5, to = 1.5, by = .05)) %>% 
  ggplot(aes(x = -prior_use, y = value)) + 
  stat_lineribbon(.width = .95, fill = "darkgreen", alpha = .3) + 
  xlab("Prior Beliefs Use") + 
  labs(tag = "A") +
  geom_point(data = areas_data, aes(x = -prior_use, y = area), alpha = .4, size = 1, color = "black")+
  ylab("Explored Area [sqkm]") +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
plot_posterior

##################
posterior_predictions2 = posterior_epred(model_area_2, newdata = data.frame(loss_attraction = seq(from = -1.5, to = 1.5, by = .05)), re_formula = NA)
posterior_predictions2 = as.data.frame(posterior_predictions2)

posterior_predictions2$draw = 1:4000
plot_posterior2 = posterior_predictions2 %>% 
  pivot_longer(1:(ncol(posterior_predictions2)-1)) %>% 
  group_by(draw) %>% 
  mutate(loss_attraction = seq(from = -1.5, to = 1.5, by = .05)) %>% 
  ggplot(aes(x = -loss_attraction, y = value)) + 
  stat_lineribbon(.width = .95, fill = "orange2", alpha = .3) + 
  ylab("Explored Area [sqkm]") + 
  labs(tag = "B")+
  geom_point(data = areas_data, aes(x = -loss_attraction, y = area), alpha = .4, size = 1, color = "black")+
  #stat_summary_bin(data = areas_data, aes(x = -loss_attraction, y = area), alpha = .4, fun.args = list(mult = 2), fun.data = mean_se)+
  xlab("Avoidance of Unsuccessful Spots") +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
plot_posterior2

panel_area = gridExtra::arrangeGrob(
  grobs = list(plot_posterior,
               plot_posterior2),
  widths = c(1, 1),
  layout_matrix = rbind(c(1,2))
)

ggsave(panel_area, file = "6_plots/panel_area.png", 
       width = 12, height = 5, bg = "white", dpi = 600, units = "cm")
ggsave(panel_area, file = "6_plots/panel_area.svg", 
       width = 12, height = 5, bg = "white", dpi = 600, units = "cm", device = svglite::svglite)


################################################################################

parameter_regression_1 = brm(parameter_3 ~ parameter_6,
                        data = trip_estimates, 
                        family = "gaussian")
parameter_regression_2 = brm(parameter_4 ~ parameter_6,
                             data = trip_estimates, 
                             family = "gaussian")

posterior_predictions_1 = posterior_epred(parameter_regression_1, newdata = data.frame(parameter_6 = seq(from = -1.5, to = 1.5, by = .05)), re_formula =NA)
posterior_predictions_1 = as.data.frame(posterior_predictions_1)
posterior_predictions_1$draw = 1:4000
posterior_predictions_1 = posterior_predictions_1 %>% 
  pivot_longer(1:(ncol(posterior_predictions_1)-1)) %>% 
  group_by(draw) %>% 
  mutate(prior_attraction = seq(from = -1.5, to = 1.5, by = .05))

posterior_predictions_2 = posterior_epred(parameter_regression_2, newdata = data.frame(parameter_6 = seq(from = -1.5, to = 1.5, by = .05)))
posterior_predictions_2 = as.data.frame(posterior_predictions_2)
posterior_predictions_2$draw = 1:4000
posterior_predictions_2 = posterior_predictions_2 %>% 
  pivot_longer(1:(ncol(posterior_predictions_2)-1)) %>% 
  group_by(draw) %>% 
  mutate(loss_avoidance = seq(from = -1.5, to = 1.5, by = .05))

trip_estimate_plot1 = trip_estimates %>%
  ggplot(aes(y = -parameter_3, x = -parameter_6)) + 
  geom_point(size = .7, alpha = .6, color = "lightgrey") + 
  stat_lineribbon(data = posterior_predictions_1, aes(x = prior_attraction, y = value), .width= .95, fill = "darkgreen", alpha = .3) + 
  coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  xlab("Prior Belief Use") + 
  labs(tag = "A")+
  geom_point(aes(y = -trip_estimates$parameter_3[trip_estimates$trip == 103], x = -trip_estimates$parameter_6[trip_estimates$trip == 103]),  alpha = .6, color = "red", size = 1.5) + 
  geom_point(aes(y = -trip_estimates$parameter_3[trip_estimates$trip == 115], x = -trip_estimates$parameter_6[trip_estimates$trip == 115]),  alpha = .6, color = "red", size = 1.5) + 
  annotate("text", y = -trip_estimates$parameter_3[trip_estimates$trip == 103], x = -trip_estimates$parameter_6[trip_estimates$trip == 103] - .25, label = "C")+
  #annotate("text", y = -trip_estimates$parameter_4[trip_estimates$trip == 108], x = -trip_estimates$parameter_6[trip_estimates$trip == 108], label = "(ii)")+
  annotate("text", y = -trip_estimates$parameter_3[trip_estimates$trip == 115], x = -trip_estimates$parameter_6[trip_estimates$trip == 115] - .25, label = "D")+
  ylab("Attraction to Successful Spots") + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
trip_estimate_plot1

trip_estimate_plot2 = trip_estimates %>%
  ggplot(aes(y = -parameter_4, x = -parameter_6)) + 
  geom_point(size = .7, alpha = .6, color = "lightgrey") + 
  stat_lineribbon(data = posterior_predictions_2, aes(x = loss_avoidance, y = value), .width= .95, fill = "darkgreen", alpha = .3) + 
  #coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  xlab("Prior Belief Use") + 
  labs(tag = "B")+
  geom_point(aes(y = -trip_estimates$parameter_4[trip_estimates$trip == 103], x = -trip_estimates$parameter_6[trip_estimates$trip == 103]),  alpha = .6, color = "red", size = 1.5) + 
  geom_point(aes(y = -trip_estimates$parameter_4[trip_estimates$trip == 115], x = -trip_estimates$parameter_6[trip_estimates$trip == 115]),  alpha = .6, color = "red", size = 1.5) + 
  annotate("text", y = -trip_estimates$parameter_4[trip_estimates$trip == 103], x = -trip_estimates$parameter_6[trip_estimates$trip == 103] - .25, label = "C")+
  #annotate("text", y = -trip_estimates$parameter_4[trip_estimates$trip == 108], x = -trip_estimates$parameter_6[trip_estimates$trip == 108], label = "(ii)")+
  annotate("text", y = -trip_estimates$parameter_4[trip_estimates$trip == 115], x = -trip_estimates$parameter_6[trip_estimates$trip == 115] - .25, label = "D")+
  ylab("Avoidance of Unsuccessful Spots") + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
trip_estimate_plot2

# Step 2: Load competition areas
competition_areas = readRDS("data/study_lakes.rds")
for (i in 1:10){
  st_crs(competition_areas[[i]]) <- "epsg:3067"
}

depth_profiles = readRDS("data/depth_profiles.rds")

# Visualize trips for illustration
trips = spot_selection_data %>% dplyr::select(trip_id) #filter(day == 4) %>% 

trip_estimates %>% filter(trip %in% pull(trips)) %>% filter(parameter_6 > .7 & parameter_4 < -.7)
trip_estimates %>% filter(trip %in% pull(trips)) %>% filter(abs(parameter_6) < .1 & abs(parameter_4) < .1)
trip_estimates %>% filter(-parameter_6 > .9 & -parameter_4 < -.9)
trips2 = trip_estimates %>% filter(-parameter_6 > .9 & -parameter_4 < -.9) %>% dplyr::select(trip)
spot_selection_data %>% filter(trip_id %in% pull(trips2)) %>% group_by(trip_id) %>% summarize(r = mean(Return),
                                                                                             day = day[1])

trips1 = trip_estimates %>% filter(-parameter_6 < -.9 & -parameter_4 > .9) %>% dplyr::select(trip)
spot_selection_data %>% filter(trip_id %in% pull(trips1)) %>% group_by(trip_id) %>% summarize(r = mean(Return),
                                                                                             day = day[1])
ggplot() + 
  geom_point(data = spot_selection_data %>% 
               filter(day == 5), aes(x = x, y = y)) + 
  geom_sf(data = depth_profiles[[8]])

profile_transformed = st_transform(depth_profiles[[6]], "epsg:3067")

spot_selection_data %>% 
  filter(day == 4) %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_point()

ex1=ggplot() + 
  geom_sf(data = competition_areas[[3]], alpha = .5) + 
  geom_sf(data = st_crop(profile_transformed, st_buffer(competition_areas[[3]], 100)), color = "lightgrey", alpha = .8) + 
  geom_sf(data = ecological_knowledge %>% filter(camera_id == 21 & day == 3), fill = "lightblue") + 
  geom_path(data = spot_selection_data %>% filter(trip_id == 111), aes(x = x, y = y), color = "darkgrey") + 
  geom_point(data = spot_selection_data %>% filter(trip_id == 111) %>% filter(Return == 0), aes(x = x, y = y), color = "lightgreen", size = 1, alpha = .5) + 
  geom_point(data = spot_selection_data %>% filter(trip_id == 111) %>% filter(Return == 1), aes(x = x, y = y), color = "darkgreen", size = 1) + 
  geom_point(data = spot_selection_data %>% filter(trip_id == 111) %>% filter(spot_id == 1), aes(x = x, y = y), color = "black", shape = 4, size = 3) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(tag = "C") +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
ex1
ex2=ggplot() + 
  geom_sf(data = competition_areas[[3]], alpha = .5) + 
  geom_sf(data = st_crop(profile_transformed, st_buffer(competition_areas[[3]], 100)), color = "lightgrey", alpha = .8) + 
  geom_sf(data = ecological_knowledge %>% filter(camera_id == 19 & day == 3), fill = "lightblue") + 
  geom_path(data = spot_selection_data %>% filter(trip_id == 108), aes(x = x, y = y), color = "darkgrey") + 
  geom_point(data = spot_selection_data %>% filter(trip_id == 108) %>% filter(Return == 0), aes(x = x, y = y), color = "lightgreen", size = 1, alpha = .5) + 
  geom_point(data = spot_selection_data %>% filter(trip_id == 108) %>% filter(Return == 1), aes(x = x, y = y), color = "darkgreen", size = 1) + 
  geom_point(data = spot_selection_data %>% filter(trip_id == 108) %>% filter(spot_id == 1), aes(x = x, y = y), color = "black", shape = 4, size = 3) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(tag = "D") +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
ex2
spot_selection_data %>% filter(trip_id == 111)
ex3=ggplot() + 
  geom_sf(data = competition_areas[[3]], alpha = .5) + 
  geom_sf(data = st_crop(profile_transformed, st_buffer(competition_areas[[3]], 100)), color = "lightgrey", alpha = .8) + 
  geom_sf(data = ecological_knowledge %>% filter(camera_id == 26 & day == 3), fill = "lightblue") + 
  geom_path(data = spot_selection_data %>% filter(trip_id == 115), aes(x = x, y = y), color = "darkgrey") + 
  geom_point(data = spot_selection_data %>% filter(trip_id == 115) %>% filter(Return == 0), aes(x = x, y = y), color = "lightgreen", size = 1, alpha = .5) + 
  geom_point(data = spot_selection_data %>% filter(trip_id == 115) %>% filter(Return == 1), aes(x = x, y = y), color = "darkgreen", size = 1) + 
  geom_point(data = spot_selection_data %>% filter(trip_id == 115) %>% filter(spot_id == 1), aes(x = x, y = y), color = "black", shape = 4, size = 3) + 
  #geom_segment(aes(x=(spot_selection_data %>% filter(trip_id == 115) %>% filter(spot_id == 1))$x, y=(spot_selection_data %>% filter(trip_id == 115) %>% filter(spot_id == 1))$y, xend=652605.5+380, yend=6917742), linetype = "dashed", size = .2)+
  #geom_segment(aes(x=(spot_selection_data %>% filter(trip_id == 115) %>% filter(spot_id == 1))$x, y=(spot_selection_data %>% filter(trip_id == 115) %>% filter(spot_id == 1))$y, xend=652605.5+800, yend=6917742-500), linetype = "dashed", size = .2)+
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(tag = "D") +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
ex3
extent(st_crop(profile_transformed, st_buffer(competition_areas[[3]], 100)))
ex4 = ggplot() + 
  geom_sf(data = competition_areas[[3]], alpha = .5) + 
  geom_sf(data = st_crop(profile_transformed, st_buffer(competition_areas[[3]], 100)), color = "lightgrey", alpha = .8) + 
  geom_sf(data = ecological_knowledge %>% filter(camera_id == 26 & day == 3), fill = "lightblue") + 
  coord_sf(xlim = c(min((spot_selection_data %>% filter(trip_id == 115))$x - 5), max((spot_selection_data %>% filter(trip_id == 115))$x) + 5),
           ylim = c(min((spot_selection_data %>% filter(trip_id == 115))$y - 5), max((spot_selection_data %>% filter(trip_id == 115))$y) + 5))+
  geom_path(data = spot_selection_data %>% filter(trip_id == 115), aes(x = x, y = y), color = "darkgrey") + 
  geom_point(data = spot_selection_data %>% filter(trip_id == 115) %>% filter(Return == 0), aes(x = x, y = y), color = "lightgreen", size = 1, alpha = 1) + 
  geom_point(data = spot_selection_data %>% filter(trip_id == 115) %>% filter(Return == 1), aes(x = x, y = y), color = "darkgreen", size = 1) + 
  geom_point(data = spot_selection_data %>% filter(trip_id == 115) %>% filter(spot_id == 1), aes(x = x, y = y), color = "black", shape = 4, size = 3) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(tag = "") +
  theme_nothing()
ex4
library(cowplot)
plot.with.inset <-
  ggdraw() +
  draw_plot(ex3) +
  draw_plot(ex4, x = 0.65, y = .7, width = .25, height = .25)
ggsave(plot.with.inset, file = "6_plots/panel_illustration_test.png", 
       width = 12, height = 12, bg = "white", dpi = 300, units = "cm")

# make panel
# panel
panel_illustrations = gridExtra::arrangeGrob(
  grobs = list(trip_estimate_plot1,
               trip_estimate_plot2,
               ex1,
               #ex2,
               ex3),
  widths = c(1, 1),
  layout_matrix = rbind(c(1, 2),
                        c(3, 4))
)

ggsave(panel_illustrations, file = "6_plots/panel_illustration.png", 
       width = 12, height = 12, bg = "white", dpi = 600, units = "cm")
ggsave(panel_illustrations, file = "6_plots/panel_illustration.svg", 
       width = 12, height = 12, bg = "white", dpi = 600, units = "cm", device = svglite::svglite)
ggsave(ex4, file = "6_plots/panel_illustration_inset.png", 
       width = 2.5, height = 2.5, bg = "white", dpi = 600, units = "cm")

# posterior correlation
spot_selection_model_fit %>% 
  spread_draws(Rho_trip[parameter_1, parameter_2]) %>% 
  filter(parameter_2 == 3 & parameter_1 == 6) %>% 
  summarize(quantile(Rho_trip, probs = c(0.025, .975)),
            mean = mean(Rho_trip))
spot_selection_model_fit %>% 
  spread_draws(Rho_trip[parameter_1, parameter_2]) %>% 
  filter(parameter_2 == 4 & parameter_1 == 6) %>% 
  summarize(quantile(Rho_trip, probs = c(0.025, .975)),
            mean = mean(Rho_trip))

################################################################################
# predict prior use from trip estimates
validation_model = brms::brm(bf(prior_use ~ parameter_6,
                                phi ~ parameter_6,
                                zoi ~ parameter_6,
                                coi ~ parameter_6),
                             family = zero_one_inflated_beta(),
                             data = trip_estimates, 
                             cores = 4,
                             chains = 4)

posterior_predictions_zoi = posterior_epred(validation_model, newdata = data.frame(parameter_6 = seq(from = -1.5, to = 1.5, by = .01)))
posterior_predictions_zoi = as.data.frame(posterior_predictions_zoi)
posterior_predictions_zoi$draw = 1:4000
plot_posterior_zoi = posterior_predictions_zoi %>% 
  pivot_longer(1:(ncol(posterior_predictions_zoi)-1)) %>% 
  group_by(draw) %>% 
  mutate(prior_attraction = seq(from = -1.5, to = 1.5, by = .01)) %>% 
  ggplot(aes(x = -prior_attraction, y = value)) + 
  stat_lineribbon(.width = .95, fill = "darkolivegreen", alpha = .3) + 
  ylab("Proportion of spots within prior area") + 
  geom_point(data = trip_estimates, aes(x = -parameter_6, y = prior_use), alpha = .4, size = 1, color = "black")+
  #stat_summary_bin(data = areas_data, aes(x = -loss_attraction, y = area), alpha = .4, fun.args = list(mult = 2), fun.data = mean_se)+
  xlab("Model estimate: Prior Belief Use") +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
plot_posterior_zoi
ggsave(plot_posterior_zoi, file = "6_plots/model_validation.png", 
       width = 8, height = 8, bg = "white", dpi = 600, units = "cm")
################################################################################
# Are prior users more successful?
catch_data <- read_delim("data/catch_data.csv",delim = ";", escape_double = FALSE, trim_ws = TRUE)
catch_data = catch_data %>% 
  filter(year == 2023 & competition == "individual")
spot_selection_data = left_join(spot_selection_data, catch_data)

spot_selection_data %>%
  group_by(day, ID) %>% 
  slice(1) %>% 
  ggplot(aes(x = prior_use, y = catch)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~day)

model_dat = spot_selection_data %>%
  group_by(day, ID) %>% 
  slice(1)


success_model_priors = brms::brm(data = model_dat, 
          formula = brmsformula(catch ~ 1 + prior_use +
                                  (1 + prior_use | day) + 
                                  (1 + prior_use | participant_id), 
                                sigma ~ 1,
                                hu ~ 1 + prior_use +
                                  (1 + prior_use | day) + 
                                  (1 + prior_use | participant_id)),
          family = hurdle_lognormal(),
          chains = 4,
          cores = 4)

################################################################################
# descriptive measure of prior use

# Compute descriptive measure of prior use and compare to model estimates

# load alternative trips
Feature_matrix = readRDS("data/processed_data/feature_matrix_trip_selection.rds")

# compute distance to prior area
knowledge_distance_matrix = readRDS(file = "data/processed_data/knowledge_distance_matrix_trip_selection.rds")

# load catch data for catch success
catch_data = read_delim("data/catch_data.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
catch_data = catch_data %>% filter(year == 2023 & competition == "individual")

# compute prior distance for each location
prior_use = list()
for (i in 1:length(knowledge_distance_matrix)){
  prior_use[[i]] = apply(knowledge_distance_matrix[[i]], 1, min)
}
rm(c)
prior_use = do.call(c, prior_use)

Feature_matrix = Feature_matrix %>% dplyr::bind_rows()
Feature_matrix$prior_distance = prior_use

compute_descriptive_measure <- function(Feature_matrix, eps){
  # compute binary measure
  Feature_matrix$prior_use = ifelse(Feature_matrix$prior_distance <= eps, 1, 0)
  
  # compute trip level estimates
  prior_use_trip = Feature_matrix %>% 
    group_by(day, camera_id, unique_trip_id, iteration) %>% 
    summarize(prior_use = mean(prior_use))
  
  # add trip-level estimates
  model = readRDS("data/processed_data/model_fit/g_full_model_fit.rds")
  
  # join trip level estimates to trip level estimates of prior use
  trip_estimates = 
    model %>% spread_draws(v_trip[unique_trip_id, parameter]) %>% 
    filter(parameter %in% c(6)) %>% 
    group_by(unique_trip_id, parameter) %>% 
    summarize(v_trip = mean(v_trip))
  
  trip_ests = left_join(trip_estimates, prior_use_trip %>% filter(iteration == 0))
  
  #trip_ests[[9]] %>% 
  #  ggplot(aes(x = prior_use, y = -v_trip)) + 
  #  geom_point() + 
  #  geom_smooth() +
  #  facet_wrap(~parameter)
  return(trip_ests)
}

trip_ests=compute_descriptive_measure(Feature_matrix, 100)
trip_ests$prior_attraction = -trip_ests$v_trip
trip_ests = left_join(trip_ests, catch_data)
trip_ests %>% 
  ggplot(aes(x = prior_use, y = catch)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~day)

mod = brms::brm(data = trip_ests,
                formula = bf(prior_use ~ prior_attraction,
                             phi ~ prior_attraction,
                             zoi ~ prior_attraction,
                             coi ~ prior_attraction),
                family = zero_one_inflated_beta(),
                iter = 2000,
                chains = 4,
                cores = 4)
posterior_predictions = posterior_epred(mod, newdata = data.frame(prior_attraction = seq(from = -1.5, to = 1.5, by = .1)))
posterior_predictions = as.data.frame(posterior_predictions)
posterior_predictions$draw = 1:4000
plot_posterior = posterior_predictions %>% 
  pivot_longer(1:length(seq(from = -1.5, to = 1.5, by = .1))) %>% 
  group_by(draw) %>% 
  mutate(v_trip = seq(from = -1.5, to = 1.5, by = .1))

# make plot
model_validation = ggplot() + 
  geom_point(data = trip_ests, aes(x = -v_trip, y = prior_use), size = .4, alpha = .5) + 
  stat_lineribbon(data = plot_posterior, aes(x = v_trip, y = value), .width = c(.95), fill = "darkolivegreen", alpha = .4) +
  xlab("Model Estimate: Attraction to Prior") + 
  ylab("Fraction of locations within prior area") + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),#, hjust = .5
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "bottom",
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 5, margin=margin(0,0,0,0)),
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
        text = element_text(family = "sans"),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key.height = unit(.1, "cm"),
        plot.margin=unit(c(.1, .2, .1, .2), "cm"))
model_validation

ggsave(model_validation, file = "6_plots/model_validation.png", 
       width = 6, height = 6, bg = "white", dpi = 300, units = "cm")


# make distribution plot: prior use
Feature_matrix$prior_use = ifelse(Feature_matrix$prior_distance <= 100, 1, 0)
descr_plot = Feature_matrix %>% 
  group_by(unique_trip_id, iteration) %>% 
  summarize(prior_use = mean(prior_use)) %>% 
  ungroup() %>% 
  mutate(quantile_prior_use = cut(prior_use, breaks = seq(from = 0, to = 1, by = .1), include.lowest = T, right = T)) %>% 
  group_by(quantile_prior_use, iteration) %>% 
  summarize(n = n()) 

colors = c("Random" = "darkolivegreen","Data" = "black")

hists = ggplot() + 
  stat_interval(data = descr_plot %>% filter(iteration > 0), aes(x = as.numeric(quantile_prior_use)/10 - .05, y = n, color = "Random"), .width = c(.5,.95), alpha = .7) + 
  #stat_summary(data = descr_plot %>% filter(iteration > 0), aes(x = as.numeric(quantile_prior_use), y = n), geom = "line", color = "green3") + 
  geom_point(data = descr_plot %>% filter(iteration == 0), aes(x = as.numeric(quantile_prior_use)/10 - .05, y = n, color = "Data"), shape = 15) + 
  #stat_summary(data = descr_plot %>% filter(iteration == 0), aes(x = as.numeric(quantile_prior_use), y = n), geom = "line", color = "black") +
  labs(tag = "C") +
  annotate("text", x = .8, y = 175, label = expression(epsilon == 100), size = 7/.pt)+
  xlab("Proportion of spots within prior area")+
  ylab("N trips")+
  scale_x_continuous(breaks = seq(from = 0, to = 1, by = .1))+
  scale_color_manual(values = colors)+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),#, hjust = .5
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 5, margin=margin(0,0,0,0)),
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
        text = element_text(family = "sans"),
        legend.box.margin = margin(-.2, -.2, -.2, -.2),
        legend.key.height = unit(.01, "cm"),
        plot.margin=unit(c(.1, .2, .1, .2), "cm"))
hists
# prior depending on navigational uncertainty
compute_descriptive_measure_eps <- function(Feature_matrix, eps){
  # compute binary measure
  Feature_matrix$prior_use = ifelse(Feature_matrix$prior_distance <= eps, 1, 0)
  
  # compute trip level estimates
  prior_use_trip = Feature_matrix %>% 
    group_by(day, camera_id, unique_trip_id, iteration) %>% 
    summarize(prior_use = mean(prior_use))
  prior_use_trip$eps = eps
  #trip_ests[[9]] %>% 
  #  ggplot(aes(x = prior_use, y = -v_trip)) + 
  #  geom_point() + 
  #  geom_smooth() +
  #  facet_wrap(~parameter)
  return(prior_use_trip)
}

res = list()
c=1
for (eps in seq(from = 0, to = 500, by = 10)){
  res[[c]] = compute_descriptive_measure_eps(Feature_matrix, eps)
  c = c + 1
}
res = res %>% dplyr::bind_rows()

lines = ggplot() + 
  #stat_summary(data = res %>% filter(iteration == 0), aes(x = eps, y = prior_use)) + 
  stat_lineribbon(data = res %>% filter(iteration > 0) %>% group_by(iteration, eps) %>% summarize(prior_use = mean(prior_use)), aes(x = eps, y = prior_use, color = "Random"), .width = c(.95) , alpha = .5, fill = "darkolivegreen") + 
  geom_path(data = res %>% filter(iteration == 0) %>% group_by(iteration, eps) %>% summarize(prior_use = mean(prior_use)), aes(x = eps, y = prior_use, color = "Data")) +
  #scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+
  coord_cartesian(ylim = c(0, 1))+
  geom_vline(xintercept = 100, linetype ="dashed")+
  labs(tag = "B")+
  xlab(expression("Navigational error"~epsilon))+
  ylab("Avg. proportion of spots in prior area")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),#, hjust = .5
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 5, margin=margin(0,0,0,0)),
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
        text = element_text(family = "sans"),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key.height = unit(.1, "cm"),
        plot.margin=unit(c(.1, .2, .1, .2), "cm"))

lines2 = ggplot() + 
  #stat_summary(data = res %>% filter(iteration == 0), aes(x = eps, y = prior_use)) + 
  stat_lineribbon(data = res %>% filter(iteration > 0) %>% group_by(iteration, eps) %>% summarize(prior_use = mean(prior_use>0)), aes(x = eps, y = prior_use, color = "Random"), .width = c(.95) , alpha = .5, fill = "darkolivegreen") + 
  geom_path(data = res %>% filter(iteration == 0) %>% group_by(iteration, eps) %>% summarize(prior_use = mean(prior_use>0)), aes(x = eps, y = prior_use, color = "Data")) +
  #scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+
  labs(tag = "A")+
  xlab(expression("Navigational error"~epsilon))+
  ylab("Proportion of trips using prior")+
  coord_cartesian(ylim = c(0, 1))+
  geom_vline(xintercept = 100, linetype ="dashed")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),#, hjust = .5
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.7,.2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 5, margin=margin(0,0,0,0)),
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
        text = element_text(family = "sans"),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key.height = unit(.1, "cm"),
        plot.margin=unit(c(.1, .2, .1, .2), "cm"))

# make panel
panel_descr_measure = gridExtra::arrangeGrob(
  grobs = list(lines2,
               lines,
               hists),
  widths = c(1, 1, 1),
  layout_matrix = rbind(c(1, 2, 3))
)

ggsave(panel_descr_measure, file = "6_plots/panel_descr_measure.png", 
       width = 18, height = 6, bg = "white", dpi = 600, units = "cm")
ggsave(panel_descr_measure, file = "6_plots/panel_descr_measure.svg", 
       width = 18, height = 6, bg = "white", dpi = 600, units = "cm", device = svglite::svglite)

################################################################################ 
# make full panel figure spatial search

panel_model_results = gridExtra::arrangeGrob(
  grobs = list(lines2,
               lines,
               hists,
               model_results,
               model_results_interaction,
               comparison_spot_selection),
  widths = c(1, 1, 1),
  heights = c(1, 1.1),
  layout_matrix = rbind(c(1, 2, 3),
                        c(4, 5, 6))
)

ggsave(panel_model_results, file = "6_plots/panel_spatial_results.png", 
       width = 18, height = 12, bg = "white", dpi = 600, units = "cm")
ggsave(panel_model_results, file = "6_plots/panel_spatial_results.svg", 
       width = 18, height = 12, bg = "white", dpi = 600, units = "cm", device = svglite::svglite)

################################################################################
print(spot_selection_model_fit, max_rows= 20)
print(spot_selection_model_fit$summary(posterior::default_summary_measures()[1],quantiles = ~quantile(.,probs = c(0.025, 0.975)), variables = c("weights","lambdas", "sigma_id", "sigma_trip", "sigma_lake")), n = 50)

conf_model = readRDS("data/processed_data/model_fit/j_spot_selection_confidence_v2_fit.rds")
print(conf_model, max_rows= 20)
print(conf_model$summary(posterior::default_summary_measures()[1],quantiles = ~quantile(.,probs = c(0.025, 0.975)), variables = c("weights","lambdas", "sigma_id", "sigma_trip", "sigma_lake")), n = 50)

patch_leaving_model = readRDS("data/processed_data/model_fit/patch_leaving_models/e_patch_leaving_model_spatial_features_fit.rds")
print(patch_leaving_model$summary(posterior::default_summary_measures()[1],quantiles = ~quantile(.,probs = c(0.025, 0.975)), variables = c("beta", "sigma_id", "sigma_trip", "sigma_lake")), n = 50)

nomination_model = readRDS("data/processed_data/model_fit/nomination_model/d_full_fit.rds")
print(nomination_model$summary(posterior::default_summary_measures()[1],quantiles = ~quantile(.,probs = c(0.025, 0.975)), variables = c("weights", "sigma_id", "sigma_trip", "sigma_lake")), n = 50)

nomination_model_conf = readRDS("data/processed_data/model_fit/nomination_model/e_confidence_fit.rds")
print(nomination_model_conf$summary(posterior::default_summary_measures()[1],quantiles = ~quantile(.,probs = c(0.025, 0.975)), variables = c("v_confidence")), n = 50)
nomination_model_conf$code()
################################################################################
test = left_join(trip_ests,fitted_features %>% group_by(day, year, camera_id) %>% summarize(loss_feature = mean(loss_feature),
                                                                                     social_feature = mean(social_feature),
                                                                                     knowledge_feature = mean(knowledge_feature),
                                                                                     n = n()))
test %>% 
  ggplot(aes(x = loss_feature, y = knowledge_feature)) + 
  geom_point() + 
  geom_smooth()
################################################################################
# correlate trip level estimates and success
trip_estimates = spot_selection_model_fit %>% 
  spread_draws(v_trip[trip_id, parameter]) %>% 
  filter(parameter == 6) %>% 
  group_by(trip_id) %>% 
  summarize(trip_estimate = mean(v_trip))

trip_estimates_patch_leaving = patch_leaving_model %>% 
  spread_draws(v_trip[trip_id, parameter]) %>% 
  filter(parameter == 6) %>% 
  group_by(trip_id) %>% 
  summarize(trip_estimate = mean(v_trip))

################################################################################
# patch choice
ids = identifiers_spot_selection %>% group_by(trip_id) %>% slice(1)
ids = left_join(ids, spot_selection_data %>% group_by(camera_id = ID, day, year) %>% summarize(return = mean(total_catch > 0)))
ids = left_join(ids, trip_estimates)
ids = left_join(ids, catch_data %>% filter(year == 2023 & competition == "individual") %>% dplyr::select(day, year, camera_id, catch))

# lognormal regression model
model_ind_success = brms::brm(data = ids, 
                       formula = brmsformula(catch ~ 1 + trip_estimate + 
                                               (1 + trip_estimate | lake_id) + 
                                               (1 + trip_estimate | participant_id)
                                             , sigma ~ 1,
                                             hu ~ 1 + trip_estimate + 
                                               (1 + trip_estimate | lake_id) + 
                                               (1 + trip_estimate | participant_id)),
                       family = hurdle_lognormal(),
                       chains = 4,
                       cores = 4)

model_ind_success_perc = brms::brm(data = ids, 
                              formula = brmsformula(return ~ 1 + trip_estimate + 
                                                      (1 + trip_estimate | lake_id) + 
                                                      (1 + trip_estimate | participant_id), 
                                                    phi ~ 1,
                                                    zi ~ 1 + trip_estimate + 
                                                      (1 + trip_estimate | lake_id) + 
                                                      (1 + trip_estimate | participant_id)),
                              family = zero_inflated_beta(),
                              chains = 4,
                              cores = 4)

################################################################################
# patch leaving
identifiers_patch_leaving = readRDS(file = "data/processed_data/identifiers_patch_leaving.rds")
ids = identifiers_patch_leaving %>% group_by(trip_id = model_trips) %>% slice(1)
ids = left_join(ids, spot_selection_data %>% group_by(camera_id = ID, day, year) %>% summarize(return = mean(total_catch > 0)))
ids = left_join(ids, trip_estimates_patch_leaving)
ids = left_join(ids, catch_data %>% filter(year == 2023 & competition == "individual") %>% dplyr::select(day, year, camera_id, catch))

# lognormal regression model
model_ind_success_pl = brms::brm(data = ids, 
                              formula = brmsformula(catch ~ 1 + trip_estimate + 
                                                      (1 + trip_estimate | day) + 
                                                      (1 + trip_estimate | unique_id)
                                                    , sigma ~ 1,
                                                    hu ~ 1 + trip_estimate + 
                                                      (1 + trip_estimate | day) + 
                                                      (1 + trip_estimate | unique_id)),
                              family = hurdle_lognormal(),
                              chains = 4,
                              cores = 4)

model_ind_success_perc_pl = brms::brm(data = ids, 
                                   formula = brmsformula(return ~ 1 + trip_estimate + 
                                                           (1 + trip_estimate | day) + 
                                                           (1 + trip_estimate | unique_id), 
                                                         phi ~ 1,
                                                         zi ~ 1 + trip_estimate + 
                                                           (1 + trip_estimate | day) + 
                                                           (1 + trip_estimate | unique_id)),
                                   family = zero_inflated_beta(),
                                   chains = 4,
                                   cores = 4)

