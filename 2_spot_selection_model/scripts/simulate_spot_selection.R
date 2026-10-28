################################################################################
#
# Title: Simulate spot selection
#
# Description: Simulates trajectories for each id based on posterior estimates
#
# Authors: Schakowski, A.
#
# Last updated: 02/04/2024 (DD/MM/YYYY)
#
################################################################################
simulate_spot_selection <- function(spot_selection_data, gps_data){

################################################################################
# compute gaussian kernel for kde2d
kern <- function(d, h, option, ncol){
  d = as.matrix(d[, 1:ncol])
  if(option == "social"){
    d = exp(-d^2 / (2 * h^2)) * (1.0 / (1 + exp(-20 * (d - 5))))
  } else {
    d = exp(-d^2 / (2 * h^2))
  }
  
  d = d %*% rep(1, ncol)
  
  return(d)
}

# softmax function to compute features
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

# function to compute heading angles
heading_angle <- function(x, y){
  v <- c(y[1]-x[1],y[2]-x[2])
  angle <- atan2(v[2],v[1])
  
  while(angle<=(-pi))
    angle <- angle + 2*pi
  while(angle>pi)
    angle <- angle - 2*pi
  
  return(angle)
}

# library
library(tidybayes)
library(Rcpp)

# source helper functions
sourceCpp("functions/fast_dist.cpp")
sourceCpp("functions/ptinpoly.cpp")

################################################################################
# load posterior estimates model
#spot_selection_fit = readRDS(file = "data/processed_data/model_fit/e_no_cond_personal_fit.rds")
spot_selection_fit = readRDS("data/processed_data/model_fit/g_full_model_fit.rds")

# data
stan_data = readRDS(file = "data/processed_data/stan_data_spot_selection.rds")
stan_data$Time_since_start = as.numeric(stan_data$Time_since_start)
stan_data$Time_since_start = stan_data$Time_since_start / max(stan_data$Time_since_start)

# compute individual level parameters
library(dplyr)
print(spot_selection_fit, max_rows = 20)
main_effects = spot_selection_fit %>% 
  spread_draws(weights[feature]) %>% 
  group_by(feature) %>% 
  summarize(weights = mean(weights))
trip_effects = spot_selection_fit %>% 
  spread_draws(v_trip[trip, feature]) %>% 
  group_by(trip, feature) %>% 
  summarize(weights = mean(v_trip))
lake_effects = spot_selection_fit %>% 
  spread_draws(v_lake[lake, feature]) %>% 
  group_by(lake, feature) %>% 
  filter(feature <=9) %>% 
  summarize(weights = mean(v_lake))
participant_effects = spot_selection_fit %>% 
  spread_draws(v_id[id, feature]) %>% 
  group_by(id, feature) %>% 
  summarize(weights = mean(v_id))

# compute kernel parameters
bandwidths = spot_selection_fit %>% 
  spread_draws(lambdas[feature]) %>% 
  group_by(feature) %>% 
  filter(feature %in% c(1, 2)) %>% 
  summarize(h = mean(lambdas))
bandwidths_offset = spot_selection_fit %>% 
  spread_draws(v_lake[lake, feature]) %>% 
  filter(feature > 9 & feature < 12) %>% 
  group_by(lake, feature) %>% 
  summarize(offset = mean(v_lake)) %>% 
  mutate(feature = feature - 9)
bandwidths = left_join(bandwidths, bandwidths_offset)
bandwidths = bandwidths %>% 
  mutate(bandwidth = exp(h + offset))

# load ids
identifiers = readRDS(file = "data/processed_data/identifiers_spot_selection.rds")

# add model ids
identifiers$lake = as.numeric(as.factor(identifiers$lake))
identifiers$trip = as.numeric(as.factor(identifiers$trip_id))
identifiers$id = as.numeric(as.factor(identifiers$participant_id))

# convert time
identifiers$time_since_start = as.numeric(identifiers$time_since_start) / max(as.numeric(identifiers$time_since_start))

# recreate trip index and participant id and lake id as used for model fitting
social_distances = readRDS(file = "data/processed_data/social_distance_matrix.rds")
social_distances = social_distances[(identifiers$spot_id!= 1)]

# subset modeled choices
sub_ids = identifiers[identifiers$spot_id != 1, ][1:length(social_distances),]

# create key matrix
sub_ids$trip_index = as.numeric(as.factor(sub_ids$trip_id))
sub_ids$lake_index = as.numeric(as.factor(sub_ids$lake_id))
sub_ids$id_index = as.numeric(as.factor(sub_ids$participant_id))
key_matrix = sub_ids %>% group_by(trip_index, lake_index, id_index) %>% slice(1) %>% dplyr::select(trip_index, lake_index, id_index, day, year, camera_id)

# add to model data
spot_selection_data$camera_id = spot_selection_data$ID
key_matrix$day = as.integer(key_matrix$day)
key_matrix$year = as.integer(key_matrix$year)
spot_selection_data = left_join(spot_selection_data, key_matrix)

# remove missing trips
spot_selection_data = spot_selection_data %>% filter(!is.na(trip_index))

# add time since start
identifiers$day = as.integer(identifiers$day)
identifiers$year = as.integer(identifiers$year)
identifiers$spot_id = as.integer(identifiers$spot_id)
spot_selection_data = left_join(spot_selection_data, identifiers %>% dplyr::select(day, year, spot_id, camera_id, time_since_start))

# rename
identifiers = spot_selection_data

################################################################################
# load lake data
# lake key to link day, year to unique lake index
lake_key = data.frame(day = c(1:5, 1:3, 5,4), 
                      year = rep(c(2023, 2022), each = 5),
                      lake_index = 1:10)
lake_key$year = as.numeric(lake_key$year)
lake_key$day = as.numeric(lake_key$day)
identifiers = left_join(identifiers, lake_key)

# Import Lake Data
lake_data = readRDS("data/depth_profiles.rds")
competition_area = readRDS("data/study_lakes.rds")

# get in same order
depth_profiles = list()
depth_profiles[[1]] <- lake_data[[7]]
depth_profiles[[2]] <- lake_data[[10]]
depth_profiles[[3]] <- lake_data[[6]]
depth_profiles[[4]] <- lake_data[[9]]
depth_profiles[[5]] <- lake_data[[8]]
depth_profiles[[6]] <- lake_data[[1]]
depth_profiles[[7]] <- lake_data[[2]]
depth_profiles[[8]] <- lake_data[[3]]
depth_profiles[[9]] <- lake_data[[5]]
depth_profiles[[10]] <- lake_data[[4]]

# polygon lakes for faster execution time
point_lakes = list()
for (i in 1:10){
  point_lakes[[i]] = competition_area[[i]] %>% 
    st_buffer(-10) %>% 
    st_cast("POINT") %>% 
    as.data.frame()
  point_lakes[[i]] = point_lakes[[i]] %>% 
    mutate(x = unlist(map(point_lakes[[i]]$polygons,1)),
           y = unlist(map(point_lakes[[i]]$polygons,2))) %>% 
    dplyr::select(x,y)
}

################################################################################
# compute competition area
check = vector()
c = 1
competition_area_new=list()
for (i in 1:5){
  
  # subset lake area
  area = point_lakes[[i]]
  
  # coordinates of spot choices
  coords = identifiers %>% filter(day == lake_key[i, "day"] & year == lake_key[i, "year"]) %>% dplyr::select(x, y) %>% filter(!is.na(x))
  
  # convex hull around spot choices
  hull = coords[coords %>% chull(),c("x", "y")]
  colnames(hull) <- c("x", "y")
  hull2 = rbind(area, hull)
  
  # hull encompassing lake area and spot choices
  hull3 = hull2[hull2 %>% chull(),c("x", "y")]
  
  # check whether any spots are outside competition area
  for (j in 1:nrow(coords)){
    check[c] = pnpoly(cbind(coords$x[j], coords$y[j]), as.matrix(hull3)) 
    c = c+1
  }
  
  # save hull
  competition_area_new[[i]] <- hull3
  
}

# 
depth_profile_points = list()
for (i in 1:5){
  tmp = st_line_sample(depth_profiles[[i]] %>% filter(depth>0), density = .1, type = "regular")
  depth_profile_points[[i]] = tmp[!st_is_empty(tmp)] %>% st_cast("POINT")
}

# load ecological knowledge
ecological_knowledge_polygons = readRDS("data/ecological_knowledge_polygons.rds")

################################################################################
# run simulation
################################################################################
# number of alternative spots (same as in model)
N_spots = 50

################################################################################
simulate_trajectory <- function(iter=1,
                                identifiers=identifiers,
                                competition_area_new=competition_area_new,
                                depth_profiles=depth_profiles,
                                gps_data=gps_data,
                                main_effects=main_effects,
                                trip_effects=trip_effects,
                                lake_effects=lake_effects,
                                participant_effects=participant_effects,
                                bandwidths=bandwidths,
                                N_spots=N_spots,
                                ecological_knowledge_polygons,
                                option = "brw"){
  
  # load c++ functions (needs to be here for future_map export)
  sourceCpp("functions/fast_dist.cpp")
  sourceCpp("functions/ptinpoly.cpp")
  
  # movement constraints
  steps = identifiers$step[identifiers$step > -99 & identifiers$angle>-99]
  steps = ifelse(steps == 0, .1, steps)
  
  turning_angles = identifiers$angle[identifiers$step > -99 & identifiers$angle>-99]
  
  c=1
  Result = list()
  
  for (trip_tmp in unique(identifiers$trip_index)){
    
    # subset spots for this participant
    spot_sub = identifiers %>% 
      filter(trip_index == trip_tmp)
    
    # select angling times
    times = spot_sub$timestamp_start
    
    # select successful spots
    return = spot_sub$Return
    
    
    # predict success from logistic model (to do)
    
    # time since start
    time_since_start = spot_sub$time_since_start
    
    # select coordinates of all spots
    coordinates = as.matrix(spot_sub[, c("x", "y")])
    
    # identifiers
    lake_tmp = spot_sub$lake_index[1]
    participant_tmp = spot_sub$id_index[1]
    day_tmp = spot_sub$day[1]
    year_tmp = spot_sub$year[1]
    camera_id_tmp = spot_sub$camera_id[1]
    
    # lake surface as points
    lake_points = competition_area_new[[spot_sub$lake_index[1]]]
    
    # depth profile
    depth_profile = depth_profiles[[spot_sub$lake_index[1]]]
    
    # social coordinates
    coordinates_social = gps_data %>% 
      filter(day == day_tmp &
               year == year_tmp &
               camera_id != camera_id_tmp &
               timestamp %in% times) %>% 
      drop_na(position_lat, position_long) %>% 
      ungroup() %>% 
      dplyr::select(x = position_lat,y = position_long, ID = camera_id, start = timestamp)
    
    ################################################################################
    # preallocate spot choices: start simulating from initial spot
    spot_choices = matrix(nrow = length(times),
                          ncol = 2)
    spot_choices[1,1:2] <- c(coordinates[1,1],coordinates[1,2])
    
    # preallocate feature vector
    social_features_choice = vector()
    success_features_choice = vector()
    roughness_features_choice = vector()
    locality_features_choice = vector()
    loss_features_choice = vector()
    knowledge_features_choice = vector()
    time_since_last_success_choice = vector()
    check = vector()
    chosen_angle = vector()
    chosen_radius = vector()
    
    # simulate alternative spots
    if (nrow(spot_sub)<2){
      
    } else {
      
     for (n in 2:nrow(spot_sub)){
     
      # compute heading angle
      if (n >= 3){
        heading = heading_angle(x = spot_choices[n-2,], y = spot_choices[n-1,])
      } else {
        heading = sample(turning_angles, 1, replace = T)
      }
      
      #plot(spot_choices[1:15,1], spot_choices[1:15,2])
      #lines(spot_choices[1:15,1], spot_choices[1:15,2])
      
      sampled_indices = sample(1:length(turning_angles), N_spots, replace = T)
      # simulate turning angle from empirical turning angle distribution
      angle = turning_angles[sampled_indices]
      radius = steps[sampled_indices]
      
      # add to coordinates of last spot
      simulated_x = spot_choices[n-1,1] + radius * cos(heading + angle)
      simulated_y = spot_choices[n-1,2] + radius * sin(heading + angle)
      
      # check if points are on lake surface, if not, resample  
      for (pt in 1:length(simulated_x)){
        counter = 1
        while(!pnpoly(c(simulated_x[pt], simulated_y[pt]), as.matrix(lake_points)) & counter < 200){
          
          sampled_index = sample(1:length(turning_angles), 1)
          angle[pt] = turning_angles[sampled_index]
          radius[pt] = steps[sampled_index]
          
          simulated_x[pt] =  spot_choices[n-1,1] + radius[pt] * cos(heading + angle[pt])
          simulated_y[pt] =  spot_choices[n-1,2] + radius[pt] * sin(heading + angle[pt])
          counter = counter + 1
        }
        
      }
      
      # alternative coordinates
      coord_tmp = cbind(simulated_x, simulated_y)
      
      # compute features for each alternative
      ####### Social Feature ########
      coordinates_social_sub = coordinates_social[coordinates_social$start == times[n],]
      coordinates_xy = LongLatToUTM(coordinates_social_sub$y, coordinates_social_sub$x, zone = 35)
      social_distances  = fastPdist2(as.matrix(coord_tmp), as.matrix(coordinates_xy[,2:3]))
      social_distances[is.na(social_distances)] <- sample(social_distances, 1)
      social_distances = t(apply(social_distances, 1, sort))[, 1:min(nrow(coordinates_xy), 15)]
      social_feature = kern(social_distances, h = pull(bandwidths[bandwidths$lake == lake_tmp & bandwidths$feature == 1, "bandwidth"]), option = "social", ncol = min(nrow(coordinates_xy), 15))
      
      ####### Roughness Feature ######
      #points = st_coordinates(depth_profile_points[[spot_sub$lake_index[1]]])
      #roughness_distances  = fastPdist2(as.matrix(coord_tmp), as.matrix(points[,1:2]))
      #roughness_distances[is.na(roughness_distances)] <- sample(roughness_distances, 1)
      #roughness_distances = t(apply(roughness_distances, 1, sort))[, 1]
      ##roughness_feature = kern(as.matrix(roughness_distances), h = pull(bandwidths[bandwidths$lake == lake_tmp & bandwidths$feature == 2, "bandwidth"]), option = "rough", ncol = 1)
      #roughness_distances = roughness_distances/1000
      ######### Success Feature #######
      if (is.null(dim(spot_choices[1:(n-1),]))){
        successful_spots = cbind(t(as.matrix(spot_choices[1:(n-1),])), return[1:(n-1)])
        successful_spots = successful_spots[successful_spots[, 3] == 1, 1:2]
      } else {
        successful_spots = cbind(spot_choices[1:(n-1),], return[1:(n-1)])
        successful_spots = successful_spots[successful_spots[, 3] == 1, 1:2]
      }
    
      time_since_last_success = ifelse(length(successful_spots)>0, n-max(which(return[1:(n-1)]==1)), NA)
      
      if(length(successful_spots) == 0){
        success_feature = rep(0, (N_spots))
      } else {
        
        if (is.null(dim(successful_spots))){
          result = fastPdist2(as.matrix(coord_tmp), t(as.matrix(successful_spots)))
        } else {
          result = fastPdist2(as.matrix(coord_tmp), as.matrix(successful_spots))
        }
        
        # sometimes missing, if distances very close
        result[is.na(result)]<-0
        
        # keep only 50 closest features for storage
        if (dim(result)[2] == 1){
          result = result
        } else {
          result = t(apply(result, 1, sort))[, 1:min(ncol(result), 15)] 
        }
        
        success_feature = kern(result, h = pull(bandwidths[bandwidths$lake == lake_tmp & bandwidths$feature == 2, "bandwidth"]), option = "success", ncol = min(ncol(result), 15))
        #success_feature = kern(result, h = 8, option = "success", ncol = min(ncol(result), 15))
        
      }
      
      ######## Loss Feature #######
      if (is.null(dim(spot_choices[1:(n-1),]))){
        unsuccessful_spots = cbind(t(as.matrix(spot_choices[1:(n-1),])), return[1:(n-1)])
        unsuccessful_spots = unsuccessful_spots[unsuccessful_spots[, 3] == 0, 1:2]
      } else {
        unsuccessful_spots = cbind(spot_choices[1:(n-1),], return[1:(n-1)])
        unsuccessful_spots = unsuccessful_spots[unsuccessful_spots[, 3] == 0, 1:2]
      }
      
      if(length(unsuccessful_spots) == 0){
        loss_feature = rep(0, N_spots)
      } else {
        
        if (is.null(dim(unsuccessful_spots))){
          result = fastPdist2(as.matrix(coord_tmp), t(as.matrix(unsuccessful_spots)))
        } else {
          result = fastPdist2(as.matrix(coord_tmp), as.matrix(unsuccessful_spots))
        }
        
        # sometimes missing, if distances very close
        result[is.na(result)]<-0
        
        # keep only 50 closest features for storage
        if (dim(result)[2] == 1){
          result = result
        } else {
          result = t(apply(result, 1, sort))[, 1:min(ncol(result), 15)] 
        }
        loss_feature = kern(result, h = pull(bandwidths[bandwidths$lake == lake_tmp & bandwidths$feature == 2, "bandwidth"]), option = "loss", ncol = min(ncol(result), 15))
        
      }
      
      # ecological knowledge feature
      polygons = ecological_knowledge_polygons %>% filter(day == day_tmp & camera_id == camera_id_tmp)
      
      # if no nomination, then everything set to 0
      if (nrow(polygons) == 0){
        knowledge_dist = rep(0, N_spots)
      } else {
        # convert fm coordinates to sf
        coords = st_as_sf(as.data.frame(coord_tmp[, c("simulated_x", "simulated_y")]), coords = c("simulated_x","simulated_y"))
        st_crs(coords) <- crs(polygons)
        
        # compute distance
        knowledge_dist = as.matrix(apply(st_distance(coords, polygons)/1000, 1, sort))
        if (dim(knowledge_dist)[1] != 50){
          knowledge_dist = t(knowledge_dist)
        }
        knowledge_dist = knowledge_dist[,1]
        #knowledge_dist = kern(as.matrix(knowledge_dist), h = pull(bandwidths[bandwidths$lake == lake_tmp & bandwidths$feature == 3, "bandwidth"]), option = "success", ncol = min(ncol(result), 1))
        
        # keep only minimum
        #knowledge_dist = apply(knowledge_dist, 1, min)
      }
      #knowledge_dist = knowledge_dist#/1000
      # locality
      locality_feature = radius/1000
      
      # matrix
      #f = cbind(social_feature,
      #          #roughness_feature,
      #          success_feature,
      #          loss_feature,
      #          locality_feature,
      #          locality_feature*time_since_start[n],
      #          knowledge_dist,
      #          knowledge_dist*success_feature,
      #          knowledge_dist*loss_feature,
      #          knowledge_dist*social_feature)
      
      f = cbind(locality_feature,
                locality_feature*time_since_start[n],
                success_feature,
                loss_feature,
                social_feature,
                knowledge_dist,
                knowledge_dist*success_feature,
                knowledge_dist*loss_feature,
                knowledge_dist*social_feature)
      
      
      w =  main_effects$weights + trip_effects$weights[trip_effects$trip == trip_tmp] + lake_effects$weights[lake_effects$lake==lake_tmp] + participant_effects$weights[participant_effects$id == participant_tmp]
      
      if (option == "high_prior"){
        
        #w[6] = w[6] - 4
        
      } else if (option == "brw"){
      
      } else if (option == "no_loss"){
        w[4] <- 0
      }
      
      # compute softmax
      s = softmax(f%*%w)
      
      # select spot acc to softmax
      if (sum(is.na(s))>0){
        ind_fw = sample(1:N_spots, 1, prob = rep(1/N_spots, N_spots))
      } else {
        ind_fw = sample(1:N_spots, 1, prob = s)
      }
      
      # save results
      spot_choices[n,] <- coord_tmp[ind_fw,]
      
      # features
      social_features_choice[n] = social_feature[ind_fw]
      success_features_choice[n] = success_feature[ind_fw]
      loss_features_choice[n] = loss_feature[ind_fw]
      #roughness_features_choice[n] = roughness_feature[ind_fw]
      locality_features_choice[n] = locality_feature[ind_fw]
      knowledge_features_choice[n] = knowledge_dist[ind_fw]
      time_since_last_success_choice[n] = time_since_last_success
      chosen_angle[n] = angle[ind_fw]
      chosen_radius[n] = radius[ind_fw]
      check[n] = sum(is.na(s))
      
    }
    
    # progress
    
    # save results tidy
    if (nrow(spot_choices)>2){
      res <- moveHMM::prepData(as.data.frame(spot_choices), type = 'UTM', coordNames = c('V1', "V2")) 
    } else {
      res = data.frame(x = spot_choices[,1], y = spot_choices[, 2], step = NA, angle = NA)
    }
    
    #hist(res$angle, breaks = 20)
    res$ID = camera_id_tmp
    res$id = participant_tmp
    res$day = day_tmp
    res$year = year_tmp
    res$spot_id = 1:n
    res$trip = trip_tmp
    res$prior_weight = w[6]
    res$loss_weight = w[4]
    res$social_feature = social_features_choice
    res$success_feature = success_features_choice
    res$loss_feature = loss_features_choice
    #res$roughness_feature = roughness_features_choice
    res$locality_feature = locality_features_choice
    res$knowledge_feature = knowledge_features_choice
    res$time_since_start = time_since_start
    res$Return = return
    res$time_since_last_success = time_since_last_success_choice
    res$iter = iter
    res$check = check
    res$chosen_angle = chosen_angle
    res$chosen_radius = chosen_radius
    Result[[c]] <- res
    c=c+1
    
  }

  }
  
  return(do.call(rbind, Result)) 
}
################################################################################
library(future)
library(furrr)
library(sf)
library(sp)
cl <- makeClusterPSOCK(15)
future::plan(future::cluster, workers = cl, gc = T)
options(future.globals.maxSize = 6000*1024^2)
model_sim = future_map(1:50, function(x) simulate_trajectory(x,
                                                            identifiers,
                                                            competition_area_new,
                                                            depth_profiles,
                                                            gps_data,
                                                            main_effects,
                                                            trip_effects,
                                                            lake_effects,
                                                            participant_effects,
                                                            bandwidths,
                                                            N_spots,
                                                            ecological_knowledge_polygons,
                                                            option = "brw"),
                       .progress = TRUE,
                       .options = furrr_options(seed = TRUE))
options(future.globals.maxSize = 500*1024^2)
predictions = do.call(rbind, model_sim)
parallel::stopCluster(cl)

# non loss walk
cl <- makeClusterPSOCK(15)
future::plan(cluster, workers = cl, gc = T)
options(future.globals.maxSize = 6000*1024^2)
nonloss_sim = future_map(1:50, function(x) simulate_trajectory(x,
                                                                 identifiers,
                                                                 competition_area_new,
                                                                 depth_profiles,
                                                                 gps_data,
                                                                 main_effects,
                                                                 trip_effects,
                                                                 lake_effects,
                                                                 participant_effects,
                                                                 bandwidths,
                                                                 N_spots,
                                                                 ecological_knowledge_polygons,
                                                                 option = "no_loss"),
                           .progress = TRUE,
                           .options = furrr_options(seed = TRUE))
options(future.globals.maxSize = 500*1024^2)
parallel::stopCluster(cl)
nonloss = do.call(rbind, nonloss_sim)

################################################################################
saveRDS(predictions, file = "data/processed_data/spot_selection_simulations_full.rds")
saveRDS(nonloss, file = "data/processed_data/spot_selection_simulations_no_loss.rds")
################################################################################
}
################################################################################
# END
################################################################################
# compute explored area
data_sub = predictions %>% 
  group_by(trip) %>% 
  slice(1) %>% 
  dplyr::select(ID, day, year, trip, prior_weight, loss_weight)
data_sub = left_join(spot_selection_data, data_sub)

areas_data = data_sub %>% 
  group_by(trip) %>% 
  mutate(n_spots = n()) %>% 
  filter(n_spots > 3) %>% 
  summarize(area = (cbind(x[chull(x,y)], y[chull(x, y)])),
            ratio = mean(Return > 0),
            prior_weight = mean(prior_weight),
            loss_weight = mean(loss_weight)) %>% 
  group_by(trip) %>% 
  mutate(n = n()) %>% 
  filter(n >= 4) %>% 
  summarize(area = Polygon(area)@area,
            ratio = mean(ratio),
            prior_weight = mean(prior_weight),
            loss_weight = mean(loss_weight)) %>% 
  mutate(area = area / (1000 * 1000)) %>% 
  mutate(prior_weight_binned = cut(-prior_weight, breaks = rev(-c(-4.8, -3.6, -3.1, -2.8, -2.5, -2.2, -1.9, -1.6, -1.2, -0.5, 0.4, Inf)),
                                   include.lowest = TRUE, right = TRUE))



test_data = 
  data_sub %>% 
  group_by(trip, day) %>% 
  mutate(n_spots = n()) %>% 
  filter(n_spots > 3) %>% 
  summarize(area = (cbind(x[chull(x,y)], y[chull(x, y)])),
            ratio = mean(Return > 0),
            prior_weight = mean(prior_weight),
            loss_weight = mean(loss_weight)) %>% 
  group_by(trip, day) %>% 
  mutate(n = n()) %>% 
  filter(n >= 4) %>% 
  summarize(area = Polygon(area)@area,
            ratio = mean(ratio),
            prior_weight = mean(prior_weight),
            loss_weight = mean(loss_weight)) %>% 
  mutate(area = area / (1000 * 1000))

mod = brms::brm(data = test_data,
                area ~ prior_weight + loss_weight + (prior_weight + loss_weight|day))

areas_data %>% 
  ggplot(aes(x = prior_weight_binned, y = area)) + 
  stat_summary()

areas_brw = predictions %>% 
  group_by(trip, iter) %>% 
  mutate(n_spots = n()) %>% 
  filter(n_spots > 3) %>% 
  summarize(area = (cbind(x[chull(x,y)], y[chull(x, y)])),
            ratio = mean(Return > 0),
            prior_weight = mean(prior_weight),
            loss_weight = mean(loss_weight)) %>% 
  group_by(trip, iter) %>% 
  mutate(n = n()) %>% 
  filter(n >= 4) %>% 
  summarize(area = Polygon(area)@area,
            ratio = mean(ratio),
            prior_weight = mean(prior_weight),
            loss_weight = mean(loss_weight)) %>% 
  mutate(area = area / (1000 * 1000)) %>% 
  mutate(prior_weight_binned = cut(-prior_weight, breaks = rev(-c(-4.8, -3.6, -3.1, -2.8, -2.5, -2.2, -1.9, -1.6, -1.2, -0.5, 0.4, Inf)),
                                   include.lowest = TRUE, right = TRUE)) %>% 
  group_by(iter, prior_weight_binned) %>% 
  summarize(area = mean(area))

areas_brw_loss_binned = predictions %>% 
  group_by(trip, iter) %>% 
  mutate(n_spots = n()) %>% 
  filter(n_spots > 3) %>% 
  summarize(area = (cbind(x[chull(x,y)], y[chull(x, y)])),
            ratio = mean(Return > 0),
            prior_weight = mean(prior_weight),
            loss_weight = mean(loss_weight)) %>% 
  group_by(trip, iter) %>% 
  mutate(n = n()) %>% 
  filter(n >= 4) %>% 
  summarize(area = Polygon(area)@area,
            ratio = mean(ratio),
            prior_weight = mean(prior_weight),
            loss_weight = mean(loss_weight)) %>% 
  mutate(area = area / (1000 * 1000)) %>% 
  mutate(prior_weight_binned = cut(-prior_weight, breaks = rev(-c(-4.8, -3.6, -3.1, -2.8, -2.5, -2.2, -1.9, -1.6, -1.2, -0.5, 0.4, Inf)),
                                   include.lowest = TRUE, right = TRUE)) %>% 
  group_by(iter, loss_weight_binned, prior_weight_binned) %>% 
  summarize(area = mean(area))


areas_brw %>% 
  ggplot(aes(x = prior_weight_binned, y = area)) + 
  stat_lineribbon(.width = .95)

areas_brw_loss_binned %>% 
  ggplot(aes(x = loss_weight_binned, y = prior_weight_binned, fill = area)) + 
  geom_tile()

areas_noloss = nonloss %>% 
  group_by(trip, iter) %>% 
  mutate(n_spots = n()) %>% 
  filter(n_spots > 3) %>% 
  summarize(area = (cbind(x[chull(x,y)], y[chull(x, y)])),
            ratio = mean(Return > 0),
            prior_weight = mean(prior_weight),
            loss_weight = mean(loss_weight)) %>% 
  group_by(trip, iter) %>% 
  mutate(n = n()) %>% 
  filter(n >= 4) %>% 
  summarize(area = Polygon(area)@area,
            ratio = mean(ratio),
            prior_weight = mean(prior_weight),
            loss_weight = mean(loss_weight)) %>% 
  mutate(area = area / (1000 * 1000)) %>% 
  mutate(prior_weight_binned = cut(-prior_weight, breaks = rev(-c(-4.8, -3.6, -3.1, -2.8, -2.5, -2.2, -1.9, -1.6, -1.2, -0.5, 0.4, Inf)),
                                   include.lowest = TRUE, right = TRUE)) %>% 
  group_by(iter, prior_weight_binned) %>% 
  summarize(area = mean(area))

areas_noloss %>% 
  ggplot(aes(x = prior_weight_binned, y = area)) + 
  stat_lineribbon(.width = .95)

# make plot
n_fun <- function(x){
  return(data.frame(y = mean(x), label = paste0(length(x),"\n\n")))
}

median_se = function(z) {
  data.frame(
    y = median(z),
    ymin = median(z) - 2*(1.253 * (sd(z)/sqrt(length(z)))),
    ymax = median(z) + 2*(1.253 * (sd(z)/sqrt(length(z))))
  )
}
colors <- c("Model Simulation" = "darkgreen", "Data" = "black", "No Loss Avoidance" = "cadetblue")
exploration_loss = ggplot() + 
  stat_summary(data = areas_data, aes(x = prior_weight_binned, y = (area), size = n), alpha = .2, show.legend = F, color = "black") + 
  stat_summary(data = areas_data, aes(x = prior_weight_binned, y = (area)), alpha = .6, geom = "line", color = "darkgrey") + 
  stat_summary(data = areas_data, aes(x = prior_weight_binned, y = (area)), fun.data = n_fun, geom = "text", size = 5/.pt) + 
  
  stat_lineribbon(data = areas_brw, aes(x = prior_weight_binned, y = (area), color = "Model Simulation"),.width = .95, alpha = .5, fill = "lightblue")+
  stat_lineribbon(data = areas_noloss, aes(x = prior_weight_binned, y = (area), color = "No Loss Avoidance"), .width = .95, alpha = .5, fill = "lightgreen")+
  scale_color_manual(values = colors) +
  scale_y_continuous(breaks = seq(from = 0, to = .3, by = .05))+
  ylab("Avg. explored area (sqkm)")+
  #ggtitle("Exploration")+
  #labs(tag = "c", color = "") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "lightgrey", size = 1)+
  xlab("Attraction to Priors")+
  scale_size(range=c(0.1, 2))+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x =  element_text(size = 7),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 6, angle = 55, vjust = .5),
        plot.title = element_text(size = 7),
        plot.tag = element_text(size = 7, face = "bold"),
        axis.text.y = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.65, .7),
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "white", fill = "white"))
ggsave(exploration_loss, file = "simulations.png", 
       width = 8, height = 8, bg = "white", dpi = 300, units = "cm")
