################################################################################
#
# TITLE: Foraging success and nominated areas
# 
# AUTHORS: Schakowski, A. (schakowski@mpib-berlin.mpg.de)
#
# DESCRIPTION: This script analyses whether nominated areas are more or less
# successful than locations outside those areas.
#
# LAST UPDATED (DD/MM/YYYY): 06/10/2026
#
################################################################################
run_visualize_success_model <- function(stan_seed){
  
# Step 1: Load Nominated Areas
ecological_knowledge_polygons = readRDS("data/ecological_knowledge_polygons_cleaned.rds")

# Load foraging locations
spot_selection_data = fread(file = "data/spot_selection_data_cleaned.csv")

# compute distance of each foraging location to nominated areas
n_contained_collective_conf1 = n_contained_collective_conf2 = n_contained_collective_conf3 = n_contained_collective = contained_individual = confidence_closest_individual = vector()
for (i in 1:nrow(spot_selection_data)){
  
  # subset foraging location
  spot = spot_selection_data[i,]
  
  # subset polygons of correct day
  polygons_collective = ecological_knowledge_polygons %>% filter(day == spot$day & camera_id != spot$camera_id)
  polygons_individual = ecological_knowledge_polygons %>% filter(day == spot$day & camera_id == spot$camera_id)
  
  # compute distance to polygon
  coords = st_as_sf(spot[, c("x", "y")], coords = c("x","y"))
  st_crs(coords) <- crs(polygons_collective)
  
  # compute distance
  dist_collective = st_distance(coords, polygons_collective)
  dist_collective1 = st_distance(coords, polygons_collective %>% filter(confidence == 1))
  dist_collective2 = st_distance(coords, polygons_collective %>% filter(confidence == 2))
  dist_collective3 = st_distance(coords, polygons_collective %>% filter(confidence == 3))
  
  # compute distance to own priors
  dist_individual = st_distance(coords, polygons_individual)
  
  # make binary distinction based on epsilon = 100
  dist_collective = sum(ifelse(dist_collective %>% drop_units() < 100,1,0))
  dist_collective1 = sum(ifelse(dist_collective1 %>% drop_units() < 100,1,0))
  dist_collective2 = sum(ifelse(dist_collective2 %>% drop_units() < 100,1,0))
  dist_collective3 = sum(ifelse(dist_collective3 %>% drop_units() < 100,1,0))
  dist_individual = ifelse(min(as.vector(dist_individual))<100, 1, 0)
  
  # get confidence values
  confidence_individual = polygons_individual$confidence[which.min(as.vector(dist_individual))]
  
  # save results
  n_contained_collective[i] = dist_collective
  n_contained_collective_conf1[i] = dist_collective1
  n_contained_collective_conf2[i] = dist_collective2
  n_contained_collective_conf3[i] = dist_collective3
  contained_individual[i] = dist_individual
  confidence_closest_individual[i] = confidence_individual
  
  # progress bar
  if (i %% 10 == 0){
    print(paste(i, "/", nrow(spot_selection_data))) 
  } else if (i == nrow(spot_selection_data)){
    print("Completed.")
  }
  
}

spot_selection_data$n_contained_collective =  n_contained_collective
spot_selection_data$n_contained_collective_conf1 = n_contained_collective_conf1
spot_selection_data$n_contained_collective_conf2 = n_contained_collective_conf2
spot_selection_data$n_contained_collective_conf3 = n_contained_collective_conf3
spot_selection_data$contained_individual = contained_individual
spot_selection_data$confidence_closest_individual = confidence_closest_individual

write_csv(spot_selection_data, file = "data/spot_selection_data_contained.csv")
spot_selection_data = fread("data/spot_selection_data_contained.csv")

# set priors for models
bprior <- c(prior_string("normal(0,1)", class = "b"),
            prior_string("exponential(1)", class = "sd"))

# predict success from distance to priors
success_model_ind = brms::brm(data = spot_selection_data,
                          formula = bf(total_catch ~ contained_individual + 
                                         (1+contained_individual|day) + 
                                         (1+contained_individual|participant_id) + 
                                         (1+contained_individual|track_id),
                                       zi ~ contained_individual + 
                                         (1+contained_individual|day) + 
                                         (1+contained_individual|participant_id) + 
                                         (1+contained_individual|track_id)),
                          family = zero_inflated_poisson(),
                          cores = 4, chains = 4,
                          seed = stan_seed)
saveRDS(success_model_ind, file = "data/processed_data/success_model_ind.rds")

success_model_coll = brms::brm(data = spot_selection_data,
                              formula = bf(total_catch ~ n_contained_collective + 
                                             (1+n_contained_collective|day) + 
                                             (1+n_contained_collective|participant_id) + 
                                             (1+n_contained_collective|track_id),
                                           zi ~ n_contained_collective + 
                                             (1+n_contained_collective|day) + 
                                             (1+n_contained_collective|participant_id) + 
                                             (1+n_contained_collective|track_id)),
                              family = zero_inflated_poisson(),
                              cores = 4, chains = 4,
                              seed = stan_seed)

saveRDS(success_model_coll, file = "data/processed_data/success_model_coll.rds")

# Confidence model
spot_selection_data$confidence_closest_individual = as.factor(spot_selection_data$confidence_closest_individual)
success_model_ind_confidence = brms::brm(data = spot_selection_data,
                              formula = bf(total_catch ~ contained_individual*confidence_closest_individual + 
                                             (1+contained_individual*confidence_closest_individual|day) + 
                                             (1+contained_individual*confidence_closest_individual|participant_id) + 
                                             (1+contained_individual*confidence_closest_individual|track_id),
                                           zi ~ contained_individual*confidence_closest_individual + 
                                             (1+contained_individual*confidence_closest_individual|day) + 
                                             (1+contained_individual*confidence_closest_individual|participant_id) + 
                                             (1+contained_individual*confidence_closest_individual|track_id)),
                              family = zero_inflated_poisson(),
                              cores = 4, chains = 4,
                              seed = stan_seed)
saveRDS(success_model_ind_confidence, file = "data/processed_data/success_model_ind_conf.rds")

success_model_collective_confidence = brms::brm(data = spot_selection_data,
                          formula = bf(total_catch ~ n_contained_collective_conf1 + n_contained_collective_conf2 + n_contained_collective_conf3 + 
                                         (1+n_contained_collective_conf1 + n_contained_collective_conf2 + n_contained_collective_conf3|participant_id) + 
                                         (1+n_contained_collective_conf1 + n_contained_collective_conf2 + n_contained_collective_conf3|track_id),
                                       zi ~ n_contained_collective_conf1 + n_contained_collective_conf2 + n_contained_collective_conf3 + 
                                         (1+n_contained_collective_conf1 + n_contained_collective_conf2 + n_contained_collective_conf3|participant_id) + 
                                         (1+n_contained_collective_conf1 + n_contained_collective_conf2 + n_contained_collective_conf3|track_id)),
                          family = zero_inflated_poisson(),
                          cores = 4, chains = 4,
                          seed = stan_seed)
saveRDS(success_model_collective_confidence, file = "data/processed_data/success_model_coll_conf.rds")

# make visualization
ind_parameters = success_model_ind %>% 
  spread_draws(b_zi_contained_individual, r_day__zi[day, name],
               b_contained_individual, r_day[day, name]) %>% 
  mutate(`N fish | P(Catch>0)` = b_contained_individual + r_day,
         `P(Catch>0)` = b_zi_contained_individual + r_day__zi) %>% 
  filter(name != "Intercept") %>% 
  pivot_longer(cols = c(`N fish | P(Catch>0)`, `P(Catch>0)`), names_to = "parameter", values_to = "lake_effect")

plot_success_ind = ind_parameters %>% 
  ggplot(aes(x = lake_effect, y = parameter, group = day)) + 
  coord_cartesian(xlim = c(-2,2)) + 
  geom_vline(xintercept = 0, linetype = "dashed", size = .4) + 
  # .95
  annotate("rect", xmin = quantile(ind_parameters$b_zi_contained_individual, .025, na.rm = T), xmax =  quantile(ind_parameters$b_zi_contained_individual, .975, na.rm = T), ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = quantile(ind_parameters$b_zi_contained_individual, .17, na.rm = T), xmax =  quantile(ind_parameters$b_zi_contained_individual, .83, na.rm = T), ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = quantile(ind_parameters$b_zi_contained_individual, .5, na.rm = T)-0.005, xmax =  quantile(ind_parameters$b_zi_contained_individual, .5, na.rm = T) + 0.005, ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "black")+
  
  annotate("rect", xmin = quantile(ind_parameters$b_contained_individual, .025, na.rm = T), xmax =  quantile(ind_parameters$b_contained_individual, .975, na.rm = T), ymin = 0.65, ymax = 1.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = quantile(ind_parameters$b_contained_individual, .17, na.rm = T), xmax =  quantile(ind_parameters$b_contained_individual, .83, na.rm = T), ymin = 0.65, ymax = 1.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = quantile(ind_parameters$b_contained_individual, .5, na.rm = T)-0.005, xmax =  quantile(ind_parameters$b_contained_individual, .5, na.rm = T) + 0.005, ymin = 0.65, ymax = 1.35,
           alpha = 1,fill = "black")+
  stat_pointinterval(position = position_dodge(width = .8), color = "black", alpha = .10, size = 0)  +
  ylab("Parameter")+
  labs(tag = "A", color = "value")+
  xlab("Posterior Est. (95%CrI)\nParameter ~ Individual Prior Belief")+
  scale_x_continuous(breaks = seq(from = -2, to = 2, by = 1))+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5, angle = 65, hjust = .5),#, hjust = .5
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
plot_success_ind

# collective
coll_parameters = success_model_coll %>% 
  spread_draws(b_zi_n_contained_collective, r_day__zi[day, name],
               b_n_contained_collective, r_day[day, name]) %>% 
  mutate(`N fish | P(Catch>0)` = b_n_contained_collective + r_day,
         `P(Catch>0)` = b_zi_n_contained_collective + r_day__zi) %>% 
  filter(name != "Intercept") %>% 
  pivot_longer(cols = c(`N fish | P(Catch>0)`, `P(Catch>0)`), names_to = "parameter", values_to = "lake_effect")

plot_success_coll = coll_parameters %>% 
  ggplot(aes(x = lake_effect, y = parameter, group = day)) + 
  #coord_cartesian(xlim = c(-2,2)) + 
  geom_vline(xintercept = 0, linetype = "dashed", size = .4) + 
  # .95
  annotate("rect", xmin = quantile(coll_parameters$b_zi_n_contained_collective, .025, na.rm = T), xmax =  quantile(coll_parameters$b_zi_n_contained_collective, .975, na.rm = T), ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = quantile(coll_parameters$b_zi_n_contained_collective, .17, na.rm = T), xmax =  quantile(coll_parameters$b_zi_n_contained_collective, .83, na.rm = T), ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = quantile(coll_parameters$b_zi_n_contained_collective, .5, na.rm = T)-0.0005, xmax =  quantile(coll_parameters$b_zi_n_contained_collective, .5, na.rm = T) + 0.0005, ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "black")+
  
  annotate("rect", xmin = quantile(coll_parameters$b_n_contained_collective, .025, na.rm = T), xmax =  quantile(coll_parameters$b_n_contained_collective, .975, na.rm = T), ymin = 0.65, ymax = 1.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = quantile(coll_parameters$b_n_contained_collective, .17, na.rm = T), xmax =  quantile(coll_parameters$b_n_contained_collective, .83, na.rm = T), ymin = 0.65, ymax = 1.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = quantile(coll_parameters$b_n_contained_collective, .5, na.rm = T)-0.0005, xmax =  quantile(coll_parameters$b_n_contained_collective, .5, na.rm = T) + 0.0005, ymin = 0.65, ymax = 1.35,
           alpha = 1,fill = "black")+
  stat_pointinterval(position = position_dodge(width = .8), color = "black", alpha = .10, size = 0)  +
  ylab("Parameter")+
  labs(tag = "B", color = "value")+
  xlab("Posterior Est. (95%CrI)\nParameter ~ Collective Prior Belief")+
  #scale_x_continuous(breaks = seq(from = -2, to = 2, by = 1))+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5, angle = 65, hjust = .5),#, hjust = .5
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
plot_success_coll

# save plot
panel_success_model = gridExtra::arrangeGrob(
  grobs = list(plot_success_ind,
               plot_success_coll),
  widths = c(1, 1),
  layout_matrix = rbind(c(1, 2))
)

ggsave(panel_success_model, file = "6_plots/panel_success_model.png", 
       width = 12, height = 6, bg = "white", dpi = 600, units = "cm")
ggsave(panel_success_model, file = "6_plots/panel_success_model.svg", 
       width = 12, height = 6, bg = "white", dpi = 600, units = "cm")
ggsave(panel_success_model, file = "6_plots/panel_success_model.pdf", 
       width = 12, height = 6, bg = "white", dpi = 600, units = "cm")

# split by confidence


################################################################################
}
