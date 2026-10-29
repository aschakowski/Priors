################################################################################
#
# TITLE: Create figure 1
# 
# AUTHORS: Schakowski, A. (schakowski@mpib-berlin.mpg.de)
#
# DESCRIPTION: This script loads data and plots figure 1
#
# LAST UPDATED (DD/MM/YYYY): 06/10/2026
#
################################################################################
create_figure_1 <- function(){
# Panel A: screenshot of competition scene.
scene <- readPNG("data/figure1/competition_scene.png")
map <- readPNG("data/figure1/bathymetric_map.png")

# plot with picture as layer
#ggplot(mapping = aes(1:10, 1:10)) +
#  annotation_raster(scene, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
#  geom_point()

# make facet from each picture
scene_im <- ggplot(mapping = aes(x = 1:max(dim(scene)), y = 1:max(dim(scene)))) +
  annotation_raster(scene, xmin = 1, xmax = dim(scene)[2], ymin = 1, ymax = dim(scene)[1]) +
  geom_point(alpha = 0) + 
  coord_fixed(ratio = 1, xlim = c(1, dim(scene)[2]),ylim = c(1, dim(scene)[1]), expand = 0) + 
  labs(tag = "A") + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        axis.ticks = element_blank(),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))

scene_im

# make facet from each picture
map_im <- ggplot(mapping = aes(x = 1:max(dim(map)), y = 1:max(dim(map)))) +
  annotation_raster(map, xmin = 1, xmax = dim(map)[2], ymin = 1, ymax = dim(map)[1]) +
  geom_point(alpha = 0) + 
  coord_fixed(ratio = 1, xlim = c(1, dim(map)[2]),ylim = c(1, dim(map)[1]), expand = 0) + 
  labs(tag = "B") + 
  xlab("Longitude") +
  ylab("Latitude") +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
map_im

# illustrate spot selection model
competition_areas = readRDS("data/study_lakes.rds")
depth_profiles = readRDS("data/depth_profiles.rds")

# Visualize trips for illustration
profile_transformed = st_transform(depth_profiles[[3]], "epsg:3067")

#ecological_knowledge_polygons %>% 
#  filter(day == 3) %>% 
#  filter(camera_id %in% c(39:44)) %>% 
#  ggplot(aes(fill = (camera_id))) + 
#  geom_sf(alpha = .3)

# select example trip
spot_selection_data = fread(file = "data/spot_selection_data_cleaned.csv")
trip_tmp = 287
spot_selection_data %>% 
  filter(day == 3) %>% 
  filter(camera_id == 42)
trip = spot_selection_data %>% filter(track_id == trip_tmp)
spot_selection_data %>% filter(day == 3) %>% 
  reframe(unique(track_id))
Feature_matrix = readRDS("data/processed_data/feature_matrix_spot_selection.rds")
Feature_matrix = Feature_matrix %>% dplyr::bind_rows()
Feature_matrix_tmp = Feature_matrix %>% filter(day == trip$day[1] & camera_id == trip$camera_id[1])
spots_ill = ggplot() + 
  geom_sf(data = competition_areas[[3]], alpha = .5) + 
  geom_sf(data = st_crop(profile_transformed, st_buffer(competition_areas[[3]], trip_tmp)), color = "lightgrey", alpha = .8) + 
  geom_sf(data = ecological_knowledge_polygons %>% filter(camera_id == trip$camera_id[1] & day == trip$day[1]), fill = "lightblue") + 
  geom_path(data = spot_selection_data %>% filter(track_id == trip_tmp), aes(x = x, y = y), color = "darkgrey") +
  geom_point(data = Feature_matrix_tmp, aes(x = x, y = y), size = 1, color = "darkgrey", alpha = .4) + 
  geom_point(data = spot_selection_data %>% filter(track_id == trip_tmp), aes(x = x, y = y), color = "black") +
  geom_point(data = spot_selection_data %>% filter(track_id == trip_tmp & spot_id == 2), aes(x = x, y = y), color = "red", size = 2) +
  geom_point(data = spot_selection_data %>% filter(track_id == trip_tmp & spot_id == 58), aes(x = x, y = y), color = "red", size = 2) +
  annotate("text", x = pull(spot_selection_data[spot_selection_data$track_id == trip_tmp & spot_selection_data$spot_id == 2, "x"])+100, 
           y = pull(spot_selection_data[spot_selection_data$track_id == trip_tmp & spot_selection_data$spot_id == 2, "y"]) + 70, label = "D")+
  annotate("text", x = pull(spot_selection_data[spot_selection_data$track_id == trip_tmp & spot_selection_data$spot_id == 58, "x"])+100, 
           y = pull(spot_selection_data[spot_selection_data$track_id == trip_tmp & spot_selection_data$spot_id == 58, "y"]) + 70, label = "E")+
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
        axis.ticks = element_blank(),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
spots_ill

## visualize patch leaving
# angling plot
patch_leaving_data = fread("data/patch_leaving_data.csv")
patch_leaving_data = patch_leaving_data %>% filter(!exclude & year == 2023)

#View(patch_leaving_data %>% 
#  filter(day == trip$day[1] & camera_id == trip$ID[1]) %>% 
#  group_by(sequence) %>% 
#    summarize(sum(reward),
#              max(bin*10/60)))

angling_sub = patch_leaving_data %>% 
  filter(day == trip$day[1] & camera_id == trip$camera_id[1]) %>% 
  filter(sequence == 2*2) %>% 
  mutate(cumulative_r = cumsum(reward)) %>% 
  mutate(bin = bin - 1) %>% 
  mutate(bin = 10 * bin / 60)
max(angling_sub$bin)


#View(patch_leaving_data %>% 
#  filter(day == trip$day[1] & camera_id == trip$ID[1]) %>% 
#  group_by(sequence) %>% 
#  summarize(t = max(bin)*10/60))

angling_sub_2 = patch_leaving_data %>% 
  filter(day == trip$day[1] & camera_id == trip$camera_id[1]) %>% 
  filter(sequence == 29*2) %>% 
  mutate(cumulative_r = cumsum(reward)) %>% 
  mutate(bin = bin - 1) %>% 
  mutate(bin = 10 * bin / 60)

angling_1 = ggplot() + 
  geom_point(data = angling_sub %>% filter(reward == 1), aes(x = bin, y = cumulative_r)) + 
  geom_path(data = angling_sub, aes(x = bin, y = cumulative_r)) + 
  geom_vline(data = angling_sub %>% filter(reward == 1), aes(xintercept = bin), linetype = "dotted", size = .6) + 
  geom_vline(data = angling_sub %>% filter(bin == 0 | bin == max(bin)), aes(xintercept = bin), size = .6) + 
  ylab("Cumulative Catch (n)") + 
  xlab("Angling Time (min)") + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.8, .2),
        legend.text = element_text(size = 5),
        plot.tag = element_text(size = 7, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white")) + 
  coord_fixed(ratio = .75, xlim = c(0, max(angling_sub$bin)))+
  scale_x_continuous(breaks = seq(from = 0, to = 10, by = 1))+
  scale_y_continuous(breaks = seq(from = 0, to = 10, by = 1))+
  labs(tag = "D")+
  geom_rect(aes(xmin = angling_sub$bin[angling_sub$reward==1][4], xmax = max(angling_sub$bin), ymin = -Inf, ymax = Inf), fill = "lightgrey", size = .6, alpha = .3)+
  geom_rect(aes(xmin = 0-.2, xmax = 0+.2, ymin = 5-.2*1/.75, ymax = 5+.2*1/.75), color = "black",fill = "white", size = .6)+
  geom_rect(aes(xmin = max(angling_sub$bin)-.2, xmax = max(angling_sub$bin)+.2, ymin = 5-.2*1/.75, ymax = 5+.2*1/.75), color = "black",fill = "white", size = .6)+
  annotate("text", x = 0, y=5, label = "S", size = 7/.pt) +
  add_fishape(family = "Percidae",option = "Perca_flavescens",
              xmin = angling_sub$bin[angling_sub$reward==1][1]-.4, xmax = angling_sub$bin[angling_sub$reward==1][1]+.4,
              ymin = 4, fill = "cadetblue",alpha = .9)+
  add_fishape(family = "Percidae",option = "Perca_flavescens",
              xmin = angling_sub$bin[angling_sub$reward==1][2]-.4, xmax = angling_sub$bin[angling_sub$reward==1][2]+.4,
              ymin = 4, fill = "cadetblue",alpha = .9)+
  add_fishape(family = "Percidae",option = "Perca_flavescens",
              xmin = angling_sub$bin[angling_sub$reward==1][3]-.4, xmax = angling_sub$bin[angling_sub$reward==1][3]+.4,
              ymin = 4, fill = "cadetblue",alpha = .9)+
  add_fishape(family = "Percidae",option = "Perca_flavescens",
              xmin = angling_sub$bin[angling_sub$reward==1][4]-.4, xmax = angling_sub$bin[angling_sub$reward==1][4]+.4,
              ymin = 4, fill = "cadetblue",alpha = .9)+
  annotate("text", x = max(angling_sub$bin), y=5, label = "L", size = 7/.pt)

angling_2 = ggplot() + 
  geom_point(data = angling_sub_2 %>% filter(reward == 1), aes(x = bin, y = cumulative_r)) + 
  geom_path(data = angling_sub_2, aes(x = bin, y = cumulative_r)) + 
  geom_vline(data = angling_sub_2 %>% filter(reward == 1), aes(xintercept = bin), linetype = "dotted", size = .6) + 
  geom_vline(data = angling_sub_2 %>% filter(bin == 0 | bin == max(bin)), aes(xintercept = bin), size = .6) + 
  ylab("Cumulative Catch (n)") + 
  xlab("Angling Time (min)") + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.8, .2),
        legend.text = element_text(size = 5),
        plot.tag = element_text(size = 7, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white")) + 
  coord_fixed(ratio = .75, xlim = c(0, max(angling_sub$bin)))+
  scale_x_continuous(breaks = seq(from = 0, to = 10, by = 1))+
  scale_y_continuous(breaks = seq(from = 0, to = 10, by = 1))+
  labs(tag = "E")+
  geom_rect(aes(xmin = angling_sub_2$bin[angling_sub_2$reward==1][4], xmax = max(angling_sub_2$bin), ymin = -Inf, ymax = Inf), fill = "lightgrey", size = .6, alpha = .3)+
  geom_rect(aes(xmin = 0-.2, xmax = 0+.2, ymin = 5-.2*1/.75, ymax = 5+.2*1/.75), color = "black",fill = "white", size = .6)+
  geom_rect(aes(xmin = max(angling_sub_2$bin)-.2, xmax = max(angling_sub_2$bin)+.2, ymin = 5-.2*1/.75, ymax = 5+.2*1/.75), color = "black",fill = "white", size = .6)+
  annotate("text", x = 0, y=5, label = "S", size = 7/.pt) +
  add_fishape(family = "Percidae",option = "Perca_flavescens",
              xmin = angling_sub_2$bin[angling_sub_2$reward==1][1]-.4, xmax = angling_sub_2$bin[angling_sub_2$reward==1][1]+.4,
              ymin = 4, fill = "cadetblue",alpha = .9)+
  add_fishape(family = "Percidae",option = "Perca_flavescens",
              xmin = angling_sub_2$bin[angling_sub_2$reward==1][2]-.4, xmax = angling_sub_2$bin[angling_sub_2$reward==1][2]+.4,
              ymin = 4, fill = "cadetblue",alpha = .9)+
  add_fishape(family = "Percidae",option = "Perca_flavescens",
              xmin = angling_sub_2$bin[angling_sub_2$reward==1][3]-.4, xmax = angling_sub_2$bin[angling_sub_2$reward==1][3]+.4,
              ymin = 4, fill = "cadetblue",alpha = .9)+
  add_fishape(family = "Percidae",option = "Perca_flavescens",
              xmin = angling_sub_2$bin[angling_sub_2$reward==1][4]-.4, xmax = angling_sub_2$bin[angling_sub_2$reward==1][4]+.4,
              ymin = 4, fill = "cadetblue",alpha = .9)+
  annotate("text", x = max(angling_sub_2$bin), y=5, label = "L", size = 7/.pt)

# make panel
# panel
panel_illustrations_rq = gridExtra::arrangeGrob(
  grobs = list(scene_im,
               map_im,
               spots_ill,
               angling_1,
               angling_2),
  widths = c(1.1, 1.7, 1.6),
  layout_matrix = rbind(c(1, 3, 4),
                        c(2, 3, 5))
)

ggsave(panel_illustrations_rq, file = "5_plots/panel_illustration_rq.png", 
       width = 16, height = 8, bg = "white", dpi = 600, units = "cm")
ggsave(panel_illustrations_rq, file = "5_plots/panel_illustration_rq.pdf", 
       width = 16, height = 8, bg = "white", dpi = 600, units = "cm")
ggsave(panel_illustrations_rq, file = "5_plots/panel_illustration_rq.svg", 
       width = 16, height = 8, bg = "white", dpi = 600, units = "cm", device = svglite::svglite)

# supplementary figure: map extraction visualization
left <- readPNG("data/figure1/visualization_filled_out_map.png")
middle <- readPNG("data/figure1/visualization_intensified_colors.png")
right <- readPNG("data/figure1/digitized_map.png")

# make facet from each picture
left_plot <- ggplot(mapping = aes(x = 1:max(dim(left)), y = 1:max(dim(left)))) +
  annotation_raster(left, xmin = 1, xmax = dim(left)[2], ymin = 1, ymax = dim(left)[1]) +
  geom_point(alpha = 0) + 
  coord_fixed(ratio = 1, xlim = c(1, dim(left)[2]),ylim = c(1, dim(left)[1]), expand = 0) + 
  labs(tag = "A") + 
  #ggtitle("Competition Scene")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        axis.ticks = element_blank(),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
middle_plot <- ggplot(mapping = aes(x = 1:max(dim(middle)), y = 1:max(dim(middle)))) +
  annotation_raster(middle, xmin = 1, xmax = dim(middle)[2], ymin = 1, ymax = dim(middle)[1]) +
  geom_point(alpha = 0) + 
  coord_fixed(ratio = 1, xlim = c(1, dim(middle)[2]),ylim = c(1, dim(middle)[1]), expand = 0) + 
  labs(tag = "B") + 
  #ggtitle("Competition Scene")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        axis.ticks = element_blank(),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
right_plot <- ggplot(mapping = aes(x = 1:max(dim(right)), y = 1:max(dim(right)))) +
  annotation_raster(right, xmin = 1, xmax = dim(right)[2], ymin = 1, ymax = dim(right)[1]) +
  geom_point(alpha = 0) + 
  coord_fixed(ratio = 1, xlim = c(1, dim(right)[2]),ylim = c(1, dim(right)[1]), expand = 0) + 
  labs(tag = "C") + 
  #ggtitle("Competition Scene")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        axis.ticks = element_blank(),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))

panel_illustrations_process = gridExtra::arrangeGrob(
  grobs = list(left_plot, middle_plot, right_plot),
  widths = c(1, 1, 1),
  layout_matrix = rbind(c(1, 2, 3))
)

ggsave(panel_illustrations_process, file = "5_plots/panel_illustrations_process.png", 
       width = 18, height = 6, bg = "white", dpi = 300, units = "cm")
ggsave(panel_illustrations_process, file = "5_plots/panel_illustrations_process.pdf", 
       width = 18, height = 6, bg = "white", dpi = 300, units = "cm")
ggsave(panel_illustrations_process, file = "5_plots/panel_illustrations_process.svg", 
       width = 18, height = 6, bg = "white", dpi = 300, units = "cm", device = svglite::svglite)

}
################################################################################
# END
################################################################################