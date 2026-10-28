################################################################################
#
# TITLE: Create figures nomination model
# 
# AUTHORS: Schakowski, A. (schakowski@mpib-berlin.mpg.de)
#
# DESCRIPTION: Creates main text and supplemental figures of the results
# of the area nomination analyses
#
# LAST UPDATED (DD/MM/YYYY): 06/10/2026
#
################################################################################

# Figure 1: Visualize nomination model results
visualize_nomination_model <- function(){
  
# Load nominated areas and depth profiles
alternative_clusters = readRDS(file = "data/output/simulated_clusters.rds")
competition_areas = readRDS(file = "data/study_lakes.rds")
depth_profiles = readRDS("data/depth_profiles.rds")

# load model fit
model = readRDS(file = "data/processed_data/model_fit/nomination_models/f_roughness_fit.rds")

# take example lake
example_lake = depth_profiles[[3]]
example_area = competition_areas[[3]]

# sample points for visualization
sampled_points = st_sample(example_area, 20000, type = "regular")
example_lake <- st_transform(example_lake, st_crs(sampled_points))

# compute ecological features for sampled points
# depth
depth = example_lake$depth[apply(st_distance(sampled_points, example_lake), 1, which.min)]

# distance to closest depth curve
rough = apply(st_distance(sampled_points, example_lake) %>% drop_units(), 1, min)

# distance to starting point
gps_data <- fread("data/gps_data.csv")
starting_locations = gps_data %>% 
  group_by(day, year, camera_id) %>% 
  slice(1) %>% 
  group_by(day, year) %>% 
  filter(year == 2023) %>% 
  summarize(starting_x = mean(position_long),
            starting_y = mean(position_lat))
starting_coords = LongLatToUTM(starting_locations$starting_x, starting_locations$starting_y, zone = 35)  
starting_locations$starting_x = starting_coords[, 2]
starting_locations$starting_y = starting_coords[, 3]
starting_locations = st_as_sf(starting_locations, coords = c("starting_x", "starting_y"))
st_crs(starting_locations) = "epsg:3067"
displacement = apply(st_distance(sampled_points, starting_locations), 1, min)

# distance to shoreline (depth curve with depth 0)
shore = st_distance(x = sampled_points, y = example_lake[example_lake$depth == 0,])
shore = (apply(shore, 1, FUN = min))

# distance to boundary (boundary of competition area)
# extract competition boundary
boundary = example_area
boundary_transformed = st_transform(boundary, "epsg:3067")
boundary_transformed = boundary_transformed %>% st_cast("MULTIPOINT")
boundary_dist = st_distance(x = sampled_points, y = boundary_transformed)
boundary_dist = (apply(boundary_dist, 1, FUN = min))

# add to df
sampled_points = st_sf(sampled_points)
sampled_points$depth = depth
sampled_points$roughness = rough
sampled_points$displacement = displacement
sampled_points$boundary = boundary_dist
sampled_points$shore = shore

# get example individual and lake
pars = model %>% 
  spread_draws(weights[parameter]) %>% 
  group_by(parameter) %>% 
  summarize(weights = mean(weights))
#1=displacement, 2=boundary, 3=depth, 4 = shore, 5 =roughness

# visualize
sampled_points = st_intersection(example_area, sampled_points)
coords = st_coordinates(sampled_points)
sampled_points = sampled_points %>% st_drop_geometry()
sampled_points$x = coords[,1]
sampled_points$y = coords[,2]

softmax <- function(par){
  n.par <- length(par)
  par1 <- sort(par, decreasing = TRUE)
  Lk <- par1[1]
  for (k in 1:(n.par-1)) {
    Lk <- max(par1[k+1], Lk) + log1p(exp(-abs(par1[k+1] - Lk))) 
  }
  val <- exp(par - Lk)
  return(val)
}

sampled_points$prediction = softmax(pars$weights[1] * scale(sampled_points$displacement) + 
          pars$weights[2] * scale(sampled_points$boundary) + 
          pars$weights[3] * scale(sampled_points$depth) + 
          pars$weights[4] * scale(sampled_points$shore) + 
          pars$weights[5] * scale(sampled_points$roughness))

colors = c("Data" = "red", "Simulated" = "lightgrey")
plot_a1_displacement=ggplot() + 
  geom_point(data = sampled_points, aes(x = x, y = y, color = displacement), alpha = .6, size = .75) +
  geom_sf(data = st_intersection(example_area,example_lake)) + 
  geom_sf(data = alternative_clusters %>% filter(day == 3 & condition == "data" & camera_id == 10), aes(fill = "Data"), alpha = .8) +
  geom_sf(data = alternative_clusters %>% filter(day == 3 & condition == "simulated" & camera_id == 10), aes(fill = "Simulated"), alpha = .6, linetype = "dotted") +
  scale_fill_manual(values = colors, guide = F)+
  scale_color_viridis_c() +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(size = 5, margin=margin(0,0,0,0)),
        plot.tag = element_text(size = 7, face = "bold", margin=margin(0,0,0,0)),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        plot.tag.position = "topleft",
        plot.margin=unit(c(0, 0, 0, 0), "cm"))+
  labs(tag = "A", color = "Displacement (m) ", fill = "Condition") + 
  ggtitle("Displacement")
plot_a1_displacement
plot_a2_boundary=ggplot() + 
  geom_point(data = sampled_points, aes(x = x, y = y, color = -boundary), alpha = .6, size = .75) +
  geom_sf(data = st_intersection(example_area,example_lake)) + 
  geom_sf(data = alternative_clusters %>% filter(day == 3 & condition == "data" & camera_id == 10), aes(fill = "Data"), alpha = .8) +
  geom_sf(data = alternative_clusters %>% filter(day == 3 & condition == "simulated" & camera_id == 10), aes(fill = "Simulated"), alpha = .6, linetype = "dotted") +
  scale_fill_manual(values = colors, guide = F)+
  scale_color_viridis_c() +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(size = 5, margin=margin(0,0,0,0)),
        plot.tag = element_text(size = 7, face = "bold", margin=margin(0,0,0,0)),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        plot.tag.position = "topleft",
        plot.margin=unit(c(0, 0, 0, 0), "cm"))+
  labs(tag = "B", color = "Boundary (m)", fill = "Condition") + 
  ggtitle("Boundary")
plot_a3_depth=ggplot() + 
  geom_point(data = sampled_points, aes(x = x, y = y, color = depth), alpha = .6, size = .75) +
  geom_sf(data = st_intersection(example_area,example_lake)) + 
  geom_sf(data = alternative_clusters %>% filter(day == 3 & condition == "data" & camera_id == 10), aes(fill = "Data"), alpha = .8) +
  geom_sf(data = alternative_clusters %>% filter(day == 3 & condition == "simulated" & camera_id == 10), aes(fill = "Simulated"), alpha = .6, linetype = "dotted") +
  scale_fill_manual(values = colors, guide = F)+
  scale_color_viridis_c() +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(size = 5, margin=margin(0,0,0,0)),
        plot.tag = element_text(size = 7, face = "bold", margin=margin(0,0,0,0)),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        plot.tag.position = "topleft",
        plot.margin=unit(c(0, 0, 0, 0), "cm"))+
  labs(tag = "C", color = "Depth Curve (m)", fill = "Condition") + 
  ggtitle("Depth")
plot_a4_shore=ggplot() + 
  geom_point(data = sampled_points, aes(x = x, y = y, color = -shore), alpha = .6, size = .75) +
  geom_sf(data = st_intersection(example_area,example_lake)) + 
  geom_sf(data = alternative_clusters %>% filter(day == 3 & condition == "data" & camera_id == 10), aes(fill = "Data"), alpha = .8) +
  geom_sf(data = alternative_clusters %>% filter(day == 3 & condition == "simulated" & camera_id == 10), aes(fill = "Simulated"), alpha = .6, linetype = "dotted") +
  scale_fill_manual(values = colors)+
  scale_color_viridis_c() +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(size = 5, margin=margin(0,0,0,0)),
        plot.tag = element_text(size = 7, face = "bold", margin=margin(0,0,0,0)),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        plot.tag.position = "topleft",
        plot.margin=unit(c(0, 0, 0, 0), "cm"))+
  labs(tag = "D", color = "Probability", fill = "Condition") + 
  ggtitle("Shore")
plot_a5_roughness=ggplot() + 
  geom_point(data = sampled_points, aes(x = x, y = y, color = -roughness), alpha = .6, size = .75) +
  geom_sf(data = st_intersection(example_area,example_lake)) + 
  geom_sf(data = alternative_clusters %>% filter(day == 3 & condition == "data" & camera_id == 10), aes(fill = "Data"), alpha = .8) +
  geom_sf(data = alternative_clusters %>% filter(day == 3 & condition == "simulated" & camera_id == 10), aes(fill = "Simulated"), alpha = .6, linetype = "dotted") +
  scale_fill_manual(values = colors)+
  scale_color_viridis_c() +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(size = 5, margin=margin(0,0,0,0)),
        plot.tag = element_text(size = 7, face = "bold", margin=margin(0,0,0,0)),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        plot.tag.position = "topleft",
        plot.margin=unit(c(0, 0, 0, 0), "cm"))+
  labs(tag = "E", color = "Probability", fill = "Condition") + 
  ggtitle("Roughness")
plot_a6_prediction=ggplot() + 
  geom_point(data = sampled_points, aes(x = x, y = y, color = prediction), alpha = .6, size = .75) +
  geom_sf(data = st_intersection(example_area,example_lake)) + 
  geom_sf(data = alternative_clusters %>% filter(day == 3 & condition == "data" & camera_id == 10), aes(fill = "Data"), alpha = .8) +
  geom_sf(data = alternative_clusters %>% filter(day == 3 & condition == "simulated" & camera_id == 10), aes(fill = "Simulated"), alpha = .6, linetype = "dotted") +
  scale_fill_manual(values = colors)+
  scale_color_viridis_c() +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(size = 5, margin=margin(0,0,0,0)),
        plot.tag = element_text(size = 7, face = "bold", margin=margin(0,0,0,0)),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        plot.tag.position = "topleft",
        plot.margin=unit(c(0, 0, 0, 0), "cm"))+
  labs(tag = "F", color = "Probability", fill = "Condition") + 
  ggtitle("Prediction")
plot_a6_prediction

# Parameter estimates
lake_effects = model %>% 
  spread_draws(v_lake[lake, parameter]) %>% 
  #mutate(parameter = ifelse(parameter == 1, "Roughness", ifelse(parameter == 2, "Displacement", "Depth")))
  mutate(parameter = ifelse(parameter == 1, "Displacement", ifelse(parameter == 2, "Boundary", ifelse(parameter == 3, "Depth", ifelse(parameter == 4, "Shore", "Roughness")))))
main_effects = model %>% 
  spread_draws(weights[parameter]) %>% 
  #mutate(parameter = ifelse(parameter == 1, "Roughness", ifelse(parameter == 2, "Displacement", "Depth")))
  mutate(parameter = ifelse(parameter == 1, "Displacement", ifelse(parameter == 2, "Boundary", ifelse(parameter == 3, "Depth", ifelse(parameter == 4, "Shore", "Roughness")))))
full_effects = left_join(lake_effects, main_effects)
full_effects$lake_effect = full_effects$weights + full_effects$v_lake

# visualize model estimates
plot_b=full_effects %>% 
  mutate(parameter = factor(parameter, levels = rev(c("Displacement", "Boundary","Depth","Shore", "Roughness")))) %>% 
  mutate(lake_effect = ifelse(parameter != "Displacement", -lake_effect, lake_effect)) %>% 
  ggplot(aes(x = lake_effect, y = parameter, group = as.factor(lake))) + 
  geom_vline(xintercept = 0, linetype = "dashed", size = .4) + 
  annotate("rect", xmin = quantile(main_effects$weights[main_effects$parameter == "Displacement"], .025), xmax =  quantile(main_effects$weights[main_effects$parameter == "Displacement"], .975), ymin = 4.65, ymax = 5.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = quantile(main_effects$weights[main_effects$parameter == "Displacement"], .17), xmax =  quantile(main_effects$weights[main_effects$parameter == "Displacement"], .83), ymin = 4.65, ymax = 5.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = quantile(main_effects$weights[main_effects$parameter == "Displacement"], .5)-0.005, xmax =  quantile(main_effects$weights[main_effects$parameter == "Displacement"], .50) + 0.005, ymin = 4.65, ymax = 5.35,
           alpha = 1,fill = "black") + 
  
  annotate("rect", xmin = quantile(-main_effects$weights[main_effects$parameter == "Boundary"], .025), xmax =  quantile(-main_effects$weights[main_effects$parameter == "Boundary"], .975), ymin = 3.65, ymax = 4.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = quantile(-main_effects$weights[main_effects$parameter == "Boundary"], .17), xmax =  quantile(-main_effects$weights[main_effects$parameter == "Boundary"], .83), ymin = 3.65, ymax = 4.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = quantile(-main_effects$weights[main_effects$parameter == "Boundary"], .5)-0.005, xmax =  quantile(-main_effects$weights[main_effects$parameter == "Boundary"], .5)+0.005, ymin = 3.65, ymax = 4.35,
           alpha = 1,fill = "black") + 
  
  annotate("rect", xmin = quantile(main_effects$weights[main_effects$parameter == "Depth"], .025), xmax =  quantile(main_effects$weights[main_effects$parameter == "Depth"], .975), ymin = 2.65, ymax = 3.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = quantile(main_effects$weights[main_effects$parameter == "Depth"], .17), xmax =  quantile(main_effects$weights[main_effects$parameter == "Depth"], .83), ymin = 2.65, ymax = 3.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = quantile(main_effects$weights[main_effects$parameter == "Depth"], .5)-0.005, xmax =  quantile(main_effects$weights[main_effects$parameter == "Depth"], .5)+0.005, ymin = 2.65, ymax = 3.35,
           alpha = 1,fill = "black") + 
  
  annotate("rect", xmin = quantile(-main_effects$weights[main_effects$parameter == "Shore"], .025), xmax =  quantile(-main_effects$weights[main_effects$parameter == "Shore"], .975), ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = quantile(-main_effects$weights[main_effects$parameter == "Shore"], .17), xmax =  quantile(-main_effects$weights[main_effects$parameter == "Shore"], .83), ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = quantile(-main_effects$weights[main_effects$parameter == "Shore"], .5)-0.005, xmax =  quantile(-main_effects$weights[main_effects$parameter == "Shore"], .5)+0.005, ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "black") + 
  
  annotate("rect", xmin = quantile(-main_effects$weights[main_effects$parameter == "Roughness"], .025), xmax =  quantile(-main_effects$weights[main_effects$parameter == "Roughness"], .975), ymin = .65, ymax = 1.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = quantile(-main_effects$weights[main_effects$parameter == "Roughness"], .17), xmax =  quantile(-main_effects$weights[main_effects$parameter == "Roughness"], .83), ymin = .65, ymax = 1.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = quantile(-main_effects$weights[main_effects$parameter == "Roughness"], .5)-0.005, xmax =  quantile(-main_effects$weights[main_effects$parameter == "Roughness"], .5)+0.005, ymin = .65, ymax = 1.35,
           alpha = 1,fill = "black") + 
  
  stat_pointinterval(position = position_dodge(width = .8), color = "black", alpha = .10, size = 0)+
  labs(tag = "G", color = "value")+
  ggtitle("Feature Weights")+
  xlab("Posterior Est. (95%CrI)\nArea Selection ~ Feature")+
  scale_x_continuous(breaks = seq(from = -6, to = 6, by = 1))+
  #coord_cartesian(ylim = c(1, 6))+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5, angle = 65, hjust = .5),#, hjust = .5
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 5, margin=margin(0,0,0,0)),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"),
        plot.margin=unit(c(.1, .2, .1, .2), "cm"))
plot_b 

# Model comparison
comparison = readRDS(file = "data/processed_data/model_fit/nomination_models/comparison/comparison.rds")
comparison = as.data.frame(comparison)
comparison$name = c("+Confidence",
                    "+Roughness",
                    "+Shore",
                    "+Depth",
                    "+Boundary",
                    "+Displacement",
                    "Random")
comparison$model_id = 1:7
plot_c = comparison %>%
  ggplot(aes(x = elpd_diff, y = fct_reorder(name, model_id))) + 
  geom_col(color = "black", width = .6) +
  geom_errorbarh(aes(xmin = elpd_diff - 4*se_diff, xmax = elpd_diff + 4*se_diff), height = .3, size = .8) + 
  ggtitle("Model Comparison") + 
  xlab(expression(Delta~"ELPD (rel. to best)")) + 
  ylab("Model") + 
  #labs(tag = "F")+
  scale_x_continuous(labels = scales::comma_format(big.mark = ','))+
  geom_vline(xintercept = 0)+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5, angle = 65, hjust = .5),#, hjust = .5
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 5, margin=margin(0,0,0,0)),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"),
        plot.margin=unit(c(.1, .2, .1, .2), "cm"))
plot_c

# Load Consensus Estimates
consensus_df = readRDS("data/output/consensus.rds")
consensus_df = consensus_df %>% 
  mutate(model = ifelse(model_id == "random", "Random", 
                        ifelse(model_id == "data", "Data", 
                               ifelse(model_id == "+roughness", "+Roughness", 
                                      ifelse(model_id == "+depth", "+Depth", 
                                             ifelse(model_id == "+displacement", "+Displacement", 
                                                    ifelse(model_id == "+boundary", "+Boundary",
                                                           ifelse(model_id == "+shore", "+Shore", NA))))))))
ns = alternative_clusters %>%
  st_drop_geometry() %>% 
  group_by(day) %>% 
  filter(iter == 0) %>% 
  summarize(n = n())
consensus_df$day = as.character(consensus_df$day)
consensus_df = left_join(consensus_df, ns)
plot_sim = consensus_df %>% 
  group_by(model, iteration, day) %>% 
  summarize(consensus = mean(consensus)) %>% #n[1]
  mutate(model = factor(model, levels = rev(c("Data","Random", "+Displacement", "+Boundary", "+Depth", "+Shore", "+Roughness"))))
plot_sim_main = consensus_df %>% 
  group_by(model, iteration) %>% 
  summarize(consensus = mean(consensus)) %>% #/mean(n)
  mutate(model = factor(model, levels = rev(c("Data","Random", "+Displacement", "+Boundary", "+Depth", "+Shore", "+Roughness"))))
plot_data = consensus_df %>% 
  filter(iteration == 0) %>% 
  group_by(model) %>% 
  summarize(mean = mean(consensus),
            upper = mean(consensus) + 2 * sd(consensus)/sqrt(n()),
            lower = mean(consensus) - 2 * sd(consensus)/sqrt(n()))

plot_data_day = consensus_df %>% 
  filter(iteration == 0) %>% 
  group_by(model, day) %>% 
  summarize(mean = mean(consensus),
            upper = mean(consensus) + 2 * sd(consensus)/sqrt(n()),
            lower = mean(consensus) - 2 * sd(consensus)/sqrt(n()))


plot_d = plot_sim %>% filter(model != "Data") %>% 
  ggplot(aes(y = model, x = consensus, group = day)) + 
  
  stat_pointinterval(.width = .95, position = position_dodge(.5), color = "white")+
  
  annotate("rect", xmin = plot_data$lower, xmax = plot_data$upper, ymin = -Inf, ymax = Inf,
           alpha = 1,fill = "lightgrey")+
  #annotate("segment", x = plot_data$lower, y = 4, xend = plot_data$upper, yend = 4,
  #         arrow = arrow(type = "closed", length = unit(0.02, "npc"), ends = "both"))+
  annotate(geom="text", x= (plot_data$lower + plot_data$upper) / 2, y=1.5, label="Data (95%CI)",
           color="black", size = 7/.pt)+
  
  annotate("rect", xmin = quantile(plot_sim_main$consensus[plot_sim_main$model == "Random"], .025), xmax =  quantile(plot_sim_main$consensus[plot_sim_main$model == "Random"], .975), ymin = 5.65, ymax = 6.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = quantile(plot_sim_main$consensus[plot_sim_main$model == "Random"], .17), xmax =  quantile(plot_sim_main$consensus[plot_sim_main$model == "Random"], .83), ymin = 5.65, ymax = 6.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = quantile(plot_sim_main$consensus[plot_sim_main$model == "Random"], .5)-.00005, xmax =  quantile(plot_sim_main$consensus[plot_sim_main$model == "Random"], .5)+.00005, ymin = 5.65, ymax = 6.35,
           alpha = 1,fill = "black") + 
  
  annotate("rect", xmin = quantile(plot_sim_main$consensus[plot_sim_main$model == "+Displacement"], .025), xmax =  quantile(plot_sim_main$consensus[plot_sim_main$model == "+Displacement"], .975), ymin = 4.65, ymax = 5.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = quantile(plot_sim_main$consensus[plot_sim_main$model == "+Displacement"], .17), xmax =  quantile(plot_sim_main$consensus[plot_sim_main$model == "+Displacement"], .83), ymin = 4.65, ymax = 5.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = quantile(plot_sim_main$consensus[plot_sim_main$model == "+Displacement"], .5)-.00005, xmax =  quantile(plot_sim_main$consensus[plot_sim_main$model == "+Displacement"], .5)+.00005, ymin = 4.65, ymax = 5.35,
           alpha = 1,fill = "black") + 
  
  annotate("rect", xmin = quantile(plot_sim_main$consensus[plot_sim_main$model == "+Boundary"], .025), xmax =  quantile(plot_sim_main$consensus[plot_sim_main$model == "+Boundary"], .975), ymin = 3.65, ymax = 4.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = quantile(plot_sim_main$consensus[plot_sim_main$model == "+Boundary"], .17), xmax =  quantile(plot_sim_main$consensus[plot_sim_main$model == "+Boundary"], .83), ymin = 3.65, ymax = 4.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = quantile(plot_sim_main$consensus[plot_sim_main$model == "+Boundary"], .5)-.00005, xmax =  quantile(plot_sim_main$consensus[plot_sim_main$model == "+Boundary"], .5)+.00005, ymin = 3.65, ymax = 4.35,
           alpha = 1,fill = "black") + 
  
  annotate("rect", xmin = quantile(plot_sim_main$consensus[plot_sim_main$model == "+Depth"], .025), xmax =  quantile(plot_sim_main$consensus[plot_sim_main$model == "+Depth"], .975), ymin = 2.65, ymax = 3.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = quantile(plot_sim_main$consensus[plot_sim_main$model == "+Depth"], .17), xmax =  quantile(plot_sim_main$consensus[plot_sim_main$model == "+Depth"], .83), ymin = 2.65, ymax = 3.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = quantile(plot_sim_main$consensus[plot_sim_main$model == "+Depth"], .5)-.00005, xmax =  quantile(plot_sim_main$consensus[plot_sim_main$model == "+Depth"], .5)+.00005, ymin = 2.65, ymax = 3.35,
           alpha = 1,fill = "black") + 
  
  annotate("rect", xmin = quantile(plot_sim_main$consensus[plot_sim_main$model == "+Shore"], .025), xmax =  quantile(plot_sim_main$consensus[plot_sim_main$model == "+Shore"], .975), ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = quantile(plot_sim_main$consensus[plot_sim_main$model == "+Shore"], .17), xmax =  quantile(plot_sim_main$consensus[plot_sim_main$model == "+Shore"], .83), ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = quantile(plot_sim_main$consensus[plot_sim_main$model == "+Shore"], .5)-.00005, xmax =  quantile(plot_sim_main$consensus[plot_sim_main$model == "+Shore"], .5)+.00005, ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "black") + 
  
  annotate("rect", xmin = quantile(plot_sim_main$consensus[plot_sim_main$model == "+Roughness"], .025), xmax =  quantile(plot_sim_main$consensus[plot_sim_main$model == "+Roughness"], .975), ymin = 0.65, ymax = 1.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = quantile(plot_sim_main$consensus[plot_sim_main$model == "+Roughness"], .17), xmax =  quantile(plot_sim_main$consensus[plot_sim_main$model == "+Roughness"], .83), ymin = 0.65, ymax = 1.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = quantile(plot_sim_main$consensus[plot_sim_main$model == "+Roughness"], .5)-.00005, xmax =  quantile(plot_sim_main$consensus[plot_sim_main$model == "+Roughness"], .5)+.00005, ymin = 0.65, ymax = 1.35,
           alpha = 1,fill = "black") + 
  
  
  stat_pointinterval(position = position_dodge(width = .8), color = "black", alpha = .10, size = 0, width = .95)+
  geom_errorbar(data = plot_data_day, aes(y = 1, x = mean, xmin = lower, xmax = upper), width = .4, position = position_dodge(.8))+
  
  ggtitle("Shared Beliefs")+
  coord_cartesian(xlim = c(0, 0.08))+
  xlab("Consensus of prior beliefs\n(N intersections/sqm)")+
  ylab("Model")+
  labs(tag = "H")+
  scale_x_continuous(breaks = seq(from = 0, to = 0.08, by = .02))+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5, angle = 65, hjust = .5),#, hjust = .5
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 5, margin=margin(0,0,0,0)),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"),
        plot.margin=unit(c(.1, .2, .1, .2), "cm"))
plot_d

################################################################################
# arrange panel
panel_figure_1 = gridExtra::arrangeGrob(
  grobs = list(plot_a1_displacement,
               plot_a2_boundary,
               plot_a3_depth,
               plot_a4_shore,
               plot_a5_roughness,
               plot_a6_prediction,
               plot_b,
               plot_d),
  heights = c(1, 1.5),
  widths = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
  layout_matrix = rbind(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
                        c(NA, 7, 7, 7, 7, 7, 7, 7, NA,NA, 8,8, 8, 8, 8, 8, 8, NA))
  #layout_matrix = rbind(c(1, 1, 2, 2, 5, 5),
  #                      c(3, 3, 4, 4, 6, 6))
)

ggsave(panel_figure_1, file = "6_plots/panel_figure_1.png", 
       width = 18, height = 10, bg = "white", dpi = 600, units = "cm")
ggsave(panel_figure_1, file = "6_plots/panel_figure_1.pdf", 
       width = 18, height = 10, bg = "white", dpi = 600, units = "cm")
ggsave(panel_figure_1, file = "6_plots/panel_figure_1.svg", 
       width = 18, height = 10, bg = "white", dpi = 600, units = "cm", device = svglite::svglite)

################################################################################
# Supplemental figures
# Confidence levels 
full_model = readRDS("data/processed_data/model_fit/nomination_models/g_confidence_fit.rds")
plot_s1_a = full_model %>% 
  spread_draws(weights[parameter, condition]) %>% 
  mutate(condition = ifelse(condition == 1, "Likely", ifelse(condition == 2, "Very Likely", "Sure"))) %>% 
  mutate(condition = factor(condition, levels = c("Likely", "Very Likely", "Sure"))) %>% 
  mutate(parameter = ifelse(parameter == 1, "Displacement", ifelse(parameter == 2, "Boundary", ifelse(parameter == 3, "Depth", ifelse(parameter == 4, "Shore", "Roughness"))))) %>% 
  mutate(parameter = factor(parameter, levels = rev(c("Displacement", "Boundary", "Depth", "Shore", "Roughness")))) %>% 
  ggplot(aes(y = parameter, x = -weights, group = condition, color = condition)) + 
  geom_vline(xintercept = 0, linetype = "dotted")+
  stat_pointinterval(position = position_dodge(width = .8))+
  scale_color_brewer(type = "seq", palette = 3)+
  #labs(tag = "A")+
  ggtitle("Confidence in Prior Beliefs")+
  xlab("Posterior Est. (95%CrI)\nArea Selection ~ Feature")+
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
plot_s1_a

ggsave(plot_s1_a, file = "6_plots/confidence_model_area_nomination.png", 
       width = 6, height = 6, bg = "white", dpi = 300, units = "cm")
ggsave(plot_s1_a, file = "6_plots/confidence_model_area_nomination.pdf", 
       width = 6, height = 6, bg = "white", dpi = 300, units = "cm")
ggsave(plot_s1_a, file = "6_plots/confidence_model_area_nomination.svg", 
       width = 6, height = 6, bg = "white", dpi = 300, units = "cm", device = svglite::svglite)

# model comparison
plot_c
ggsave(plot_c, file = "6_plots/model_comparison_area_nomination.png", 
       width = 6, height = 6, bg = "white", dpi = 300, units = "cm")
ggsave(plot_c, file = "6_plots/model_comparison_area_nomination.pdf", 
       width = 6, height = 6, bg = "white", dpi = 300, units = "cm")
ggsave(plot_c, file = "6_plots/model_comparison_area_nomination.svg", 
       width = 6, height = 6, bg = "white", dpi = 300, units = "cm", device = svglite::svglite)

################################################################################
# descriptive histograms
ecological_knowledge_polygons = readRDS(file = "data/ecological_knowledge_polygons_cleaned.rds")
hist_n = ecological_knowledge_polygons %>% 
  filter(day %in% 1:5) %>% 
  group_by(camera_id, day) %>% 
  summarize(n = n()) %>% 
  ggplot(aes(x = n)) + 
  geom_histogram() + 
  facet_wrap(~day, nrow = 1) + 
  scale_x_continuous(breaks = seq(from = 0, to = 10, by = 1))+
  scale_y_continuous(breaks = seq(from = 0, to = 25, by = 5))+
  labs(x = "N nominated areas", y = "N trips")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.ticks = element_blank(),
        title = element_text(size = 5, margin=margin(0,0,0,0)),
        plot.tag = element_text(size = 7, face = "bold", margin=margin(0,0,0,0)),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        plot.tag.position = "topleft",
        legend.key.size = unit(2, "lines"),
        legend.text = element_text(size = 3),
        legend.title = element_text(size = 3, margin=margin(0,0,0,0)),
        legend.background = element_rect(colour = "black", fill = "white"),
        plot.margin=unit(c(.1, .2, .1, .2), "cm"))
ggsave(hist_n, file = "6_plots/hist_nominated_areas.png", 
       width = 18, height = 5, bg = "white", dpi = 300, units = "cm")
ggsave(hist_n, file = "6_plots/hist_nominated_areas.svg", 
       width = 18, height = 5, bg = "white", dpi = 300, units = "cm", device = svglite::svglite)
ggsave(hist_n, file = "6_plots/hist_nominated_areas.pdf", 
       width = 18, height = 5, bg = "white", dpi = 300, units = "cm")

}
################################################################################
# END
################################################################################