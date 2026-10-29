################################################################################
#
# TITLE: Simulate alternative foraging locations
#
# AUTHORS: Schakowski, A. (schakowski@mpib-berlin.mpg.de)
#
# DESCRIPTION: This script simulate n random locations for each observed foraging
# spot for comparison.
#
# LAST UPDATED (DD/MM/YYYY): 10/10/2026
#
################################################################################

simulate_alternative_spots <- function(){
  
  ################################################################################
  # load data
  ################################################################################
  
  # data contains information for all foraging locations
  spot_selection_data = fread(file = "data/spot_selection_data_cleaned.csv")
  
  # Import competition areas
  competition_area = readRDS("data/study_lakes.rds")
  
  # polygon lakes for faster execution time
  point_lakes = list()
  for (i in 1:5){
    point_lakes[[i]] = competition_area[[i]] %>% 
      st_buffer(10) %>% 
      st_cast("POINT") %>% 
      as.data.frame()
    point_lakes[[i]] = point_lakes[[i]] %>% 
      mutate(x = unlist(map(point_lakes[[i]]$polygons,1)),
             y = unlist(map(point_lakes[[i]]$polygons,2))) %>% 
      dplyr::select(x,y)
  }
  
  # Import GPS Data
  gps_data <- fread("data/gps_data.csv")
  
  # arrange data chronologically
  spot_selection_data = spot_selection_data %>% 
    arrange(year, day, camera_id, spot_id)
  
  # how many alternative spots to simulate
  N_spots = 50
  
  # source helper functions
  sourceCpp("functions/fast_dist.cpp")
  sourceCpp("functions/ptinpoly.cpp")
  
  # check if any coordinates are outside competition area
  check = vector()
  c = 1
  competition_area_new=list()
  for (i in 1:5){
    
    area = point_lakes[[i]]
    
    coords = spot_selection_data %>% filter(day == i) %>% dplyr::select(x, y)
    
    if (nrow(coords) == 0){
      
      competition_area_new[[i]] = area
      
    } else {
      
      # convex hull of angling spots
      hull = coords[coords %>% chull(), ]
      colnames(hull) <- c("x", "y")
      hull2 = rbind(area, hull)
      hull3 = hull2[hull2 %>% chull(),c("x", "y")]
      
      for (j in 1:nrow(coords)){
        check[c] = pnpoly(cbind(coords$x[j], coords$y[j]), as.matrix(hull3)) 
        c = c+1
      }
      print(sum(!check))
      
      competition_area_new[[i]] <- hull3
      
    }
    
  }
  
  # compute starting point of each competition
  starting_point = gps_data %>% 
    group_by(day, year, camera_id) %>% 
    arrange(timestamp) %>% 
    slice(1) %>% 
    group_by(day, year) %>%
    drop_na(year) %>% 
    summarize(position_long = mean(position_long),
              position_lat = mean(position_lat))
  
  # convert to xy
  starting_points_xy = LongLatToUTM(x = starting_point$position_long,
                                    y = starting_point$position_lat,
                                    zone = 35)
  
  starting_point$x = starting_points_xy[, 3]
  starting_point$y = starting_points_xy[, 3]
  
  # how many spots
  T = nrow(spot_selection_data)
  
  # preallocate features
  Feature_matrix = list()
  
  # start simulation
  # takes ~ 10min
  # initialize step length distribution
  steps = spot_selection_data$step[spot_selection_data$step > -99 & spot_selection_data$angle > -99]
  steps = ifelse(steps == 0, .1, steps)
  angles = spot_selection_data$angle[spot_selection_data$step > -99 & spot_selection_data$angle > -99]
  
  # plot for supplement
  vis = data.frame(steps = steps,
                   angles = angles)
  
  a=ggplot(vis, aes(x = steps)) +
    geom_density(color = "black", fill = "darkgrey") + 
    coord_cartesian(xlim = c(0, 300)) + 
    ylab("Density")+
    labs(tag = "A", color = "value")+
    xlab("Relocation Distance (m)")+
    scale_x_continuous(breaks = seq(from = 0, to = 300, by = 100))+
    scale_y_continuous(breaks = seq(from = 0, to = .03, by = .01))+
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
  b=ggplot(vis, aes(x = angles)) +
    geom_density(color = "black", fill = "darkgrey")+
    #geom_histogram(bins = 50, color = "black", fill = "darkgrey") + 
    #coord_cartesian(xlim = c(-3.5, 3.5)) + 
    ylab("Density")+
    labs(tag = "B", color = "value")+
    xlab("Turning Angle (rad.)")+
    scale_x_continuous(breaks = seq(from = -4, to = 4, by = 2))+
    scale_y_continuous(breaks = seq(from = 0, to = .4, by = .1))+
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
  
  panel_distributions = gridExtra::arrangeGrob(
    grobs = list(a, b),
    widths = c(1, 1),
    layout_matrix = matrix(1:2, nrow = 1, ncol = 2)
  )
  
  ggsave(panel_distributions, file = "5_plots/panel_distributions.png", 
         width = 9, height = 5, bg = "white", dpi = 1200, units = "cm")
  ggsave(panel_distributions, file = "5_plots/panel_distributions.svg", 
         width = 9, height = 5, bg = "white", dpi = 1200, units = "cm")
  ggsave(panel_distributions, file = "5_plots/panel_distributions.pdf", 
         width = 9, height = 5, bg = "white", dpi = 1200, units = "cm")
  
  # start sim
  for (t in 1:T){
    
    # get spot data
    chosen_spot = spot_selection_data[t,]
    
    # subset all choices of trip for return vector
    choices_id = spot_selection_data[spot_selection_data$track_id == chosen_spot$track_id, ]
    first_choice = min(choices_id$unique_spot_id)
    
    # lake surface as polygon
    lake = competition_area[[chosen_spot$day]]
    
    # lake surface as points
    lake_points = competition_area_new[[chosen_spot$day]]
    
    # simulate coordinates of alternative spots
    if (chosen_spot$spot_id == 1){
      
      # simulate randomly from lake
      simulated_points = as.data.frame(st_sample(lake, N_spots))
      simulated_points = simulated_points %>% mutate(x = unlist(map(simulated_points$geometry,1)),
                                                     y = unlist(map(simulated_points$geometry,2)))
      
      simulated_x = simulated_points$x
      simulated_y = simulated_points$y
      
    } else {
      
      # compute heading angle:
      if(chosen_spot$spot_id == 2){
        
        heading = heading_angle(x = as.matrix((starting_point %>% filter(day == chosen_spot$day)))[, c("x","y")], 
                                y = as.matrix(spot_selection_data[t-1,c("x", "y")]))
        
      } else {
        
        heading = heading_angle(x = as.matrix(spot_selection_data[t-2, c("x", "y")]), 
                                y = as.matrix(spot_selection_data[t-1,c("x", "y")]))
        
      }
      
      # if simulated spots are not in competition area
      radius = sample(steps, N_spots, replace = T)
      
      # simulate turning angle from empirical turning angle distribution
      turning_angle = sample(angles, N_spots, replace = F)
      
      # add to coordinates of last spot
      simulated_x = pull(spot_selection_data[t-1, "x"]) + radius * cos(heading + turning_angle)
      simulated_y = pull(spot_selection_data[t-1, "y"]) + radius * sin(heading + turning_angle)
      
      # check if points are on lake surface, if not, resample  
      for (pt in 1:length(simulated_x)){
        counter = 1
        while(!pnpoly(c(simulated_x[pt], simulated_y[pt]), as.matrix(lake_points)) & counter < 100){
          a = sample(angles, 1)
          r = sample(steps, 1)
          simulated_x[pt] = pull(spot_selection_data[t-1, "x"] + r * cos(heading + a))
          simulated_y[pt] = pull(spot_selection_data[t-1, "y"] + r * sin(heading + a))
          counter = counter + 1
        }
      }
      
    }
    
    if (chosen_spot$spot_id == 1){
      # create features object
      features = data.frame(x = simulated_x,
                            y = simulated_y,
                            spot_id = spot_selection_data[t, "spot_id"],
                            unique_spot_id = spot_selection_data[t, "unique_spot_id"],
                            track_id = spot_selection_data[t, "track_id"],
                            participant_id = spot_selection_data[t, "participant_id"],
                            lake_id = spot_selection_data[t, "lake_id"],
                            time_start = spot_selection_data[t, "time_start"],
                            time_end = 1,
                            day = spot_selection_data[t, "day"], 
                            year = spot_selection_data[t, "year"], 
                            camera_id = spot_selection_data[t, "camera_id"],
                            chosen = 0,
                            timestamp_start = spot_selection_data[t, "timestamp_start"])
    } else {
      # create features object
      features = data.frame(x = simulated_x,
                            y = simulated_y,
                            spot_id = spot_selection_data[t, "spot_id"],
                            unique_spot_id = spot_selection_data[t, "unique_spot_id"],
                            track_id = spot_selection_data[t, "track_id"],
                            participant_id = spot_selection_data[t, "participant_id"],
                            lake_id = spot_selection_data[t, "lake_id"],
                            time_start = spot_selection_data[t, "time_start"],
                            time_end = spot_selection_data[t-1, "time_end"],
                            day = spot_selection_data[t, "day"], 
                            year = spot_selection_data[t, "year"], 
                            camera_id = spot_selection_data[t, "camera_id"],
                            chosen = 0,
                            timestamp_start = spot_selection_data[t, "timestamp_start"])
    }
    
    # chosen spot in first row
    features[N_spots + 1, ] <- as.data.frame(chosen_spot) %>% 
      dplyr::select(x, y, spot_id, unique_spot_id, track_id, participant_id, lake_id, time_start, time_end, day, year, camera_id, chosen, timestamp_start)
    features[N_spots + 1, "time_end"] <- features$time_end[1]
    features <- features[(N_spots+1):1,]
    
    # add to feature matrix
    Feature_matrix[[t]] <- features
    
    # progress
    if (t %% 100 == 0){
      print(paste(t,"/", T, sep = ""))
    } 
    
    if (t == T){
      print("Simulation Completed.")
    }
    
  }
  
  # save results
  saveRDS(Feature_matrix, file = "data/processed_data/feature_matrix_spot_selection.rds")
  
  # create identifier list
  identifiers = data.frame(day = vector(),
                           year = vector(),
                           spot_id = vector(),
                           unique_spot_id = vector(),
                           track_id = vector(),
                           participant_id = vector(),
                           lake_id = vector(),
                           camera_id = vector(),
                           t = vector(),
                           time_since_start = vector())
  for (t in 1:length(Feature_matrix)){
    day = Feature_matrix[[t]]$day[1]
    year = Feature_matrix[[t]]$year[1]
    spot_id = Feature_matrix[[t]]$spot_id[1]
    unique_spot_id = Feature_matrix[[t]]$unique_spot_id[1]
    track_id = Feature_matrix[[t]]$track_id[1]
    participant_id = Feature_matrix[[t]]$participant_id[1]
    lake_id = Feature_matrix[[t]]$lake_id[1]
    camera_id = Feature_matrix[[t]]$camera_id[1]
    time_since_start = Feature_matrix[[t]]$time_start[1]
    identifiers[t,] <- c(day, year, spot_id, unique_spot_id, track_id, participant_id, lake_id, camera_id, t, time_since_start)
  }
  
  # load catch data
  catch_data = read_csv("data/catch_data.csv")
  
  # rename participant id to avoid doubling
  names(catch_data)[names(catch_data) == 'participant_id'] <- 'original_participant_id'
  
  # load survey data
  survey_data = read_csv("data/survey_data.csv")
  names(survey_data)[names(survey_data) == 'participant_id'] <- 'original_participant_id'
  
  # join survey and catch data
  catch_data = left_join(catch_data, survey_data)
  
  # convert for joining
  catch_data$day = as.character(catch_data$day)
  catch_data$year = as.character(catch_data$year)
  
  # join id information to cluster data
  identifiers = left_join(identifiers, catch_data)
  
  # join avg resource return to data
  identifiers = left_join(identifiers, spot_selection_data %>% 
    mutate(day = as.character(day)) %>% 
    group_by(day, camera_id) %>% 
    summarize(average_success = mean(Return)))
  
  saveRDS(identifiers, "data/processed_data/identifiers_spot_selection.rds")
  
}