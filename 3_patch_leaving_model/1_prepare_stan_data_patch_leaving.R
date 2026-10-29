################################################################################
#
# TITLE: Prepare stan data for patch leaving models
# 
# AUTHORS: Schakowski, A. (schakowski@mpib-berlin.mpg.de)
#
# DESCRIPTION: This script prepares a stan data list to run the patch leaving
# models.
#
# LAST UPDATED (DD/MM/YYYY): 06/10/2026
#
################################################################################
prepare_stan_data_patch_leaving <- function(){
  
  
  # kernel function to simulate features
  kern <- function(d, h, option, ncol){
    if(option == "social"){
      d = exp(-d^2 / (2 * h^2)) * (1.0 / (1 + exp(-20 * (d - 5))))
    } else {
      d = exp(-d^2 / (2 * h^2))
    }
    
    if (is.null(dim(d))){
      d = sum(d)
    } else {
      d = d %*% rep(1, ncol(d))
    }
    
    return(d)
  }
  
  ################################################################################
  # load angling data
  patch_leaving_data <- fread("data/patch_leaving_data.csv")
  patch_leaving_data = patch_leaving_data %>% filter(!exclude) %>% filter(year == 2023)
  patch_leaving_data = patch_leaving_data %>% group_by(day, year, camera_id) %>% mutate(n_spots = length(unique(spot_id)))
  
  # load participant ids
  catch_data <- read_csv("data/catch_data.csv")
  patch_leaving_data = left_join(patch_leaving_data, catch_data %>% dplyr::select(day, year, camera_id, participant_id))
  
  # compute unique trip id
  patch_leaving_data$unique_trip_id = as.numeric(as.factor(paste(patch_leaving_data$day, patch_leaving_data$year, patch_leaving_data$camera_id)))
  
  # compute unique participant id
  patch_leaving_data$unique_participant_id = as.numeric(as.factor(paste(patch_leaving_data$participant_id)))
  
  # compute unique lake id
  patch_leaving_data$unique_lake_id = as.numeric(as.factor(paste(patch_leaving_data$day, patch_leaving_data$year)))
  
  # how many unique angling spots
  patch_leaving_data$unique_spot_id = as.numeric(as.factor(paste(patch_leaving_data$unique_trip_id, patch_leaving_data$spot_id)))
  
  # compute cumulative catch at spots
  patch_leaving_data = patch_leaving_data %>% 
    group_by(unique_spot_id) %>% 
    mutate(cumulative_catch = cumsum(reward))
  
  # compute time since event
  patch_leaving_data = patch_leaving_data %>% 
    arrange(day, year, camera_id, spot_id, bin)
  
  time_since_event = vector()
  time_since_event[1] = 0
  for (i in 2:nrow(patch_leaving_data)){
    
    if (patch_leaving_data$reward[i] == 1 || patch_leaving_data$unique_spot_id[i] != patch_leaving_data$unique_spot_id[i-1]){
      
      time_since_event[i] = 0
      
    } else {
      time_since_event[i] = 1 + time_since_event[i-1]
    }
    
  }
  
  patch_leaving_data$time_since_event = time_since_event
  
  # compute leaving state
  patch_leaving_data = 
    patch_leaving_data %>% 
    group_by(unique_spot_id) %>% 
    mutate(state = ifelse(bin == max(bin), 1, 0)) %>% 
    mutate(state = ifelse(last_spot == 1, 0, state))
  
  # compute angling time
  patch_leaving_data = 
    patch_leaving_data %>% 
    group_by(unique_spot_id) %>% 
    mutate(angling_time = 1:n())
  
  # compute catch
  patch_leaving_data$catch = patch_leaving_data$reward
  
  # for each trip compute start and end row
  patch_leaving_data = patch_leaving_data %>% 
    arrange(day, year, camera_id, spot_id, bin)
  patch_leaving_data$row_index = 1:nrow(patch_leaving_data)
  
  trip_start = patch_leaving_data %>% group_by(unique_trip_id) %>% summarize(t = min(row_index))
  trip_end = patch_leaving_data %>% group_by(unique_trip_id) %>% summarize(t = max(row_index))
  trip_start = trip_start %>% arrange(unique_trip_id)
  trip_end = trip_end %>% arrange(unique_trip_id)
  trip_start = as.data.frame(trip_start)
  trip_end = as.data.frame(trip_end)
  
  # add spatial features
  spatial_features = readRDS(file = "data/processed_data/fitted_features.rds")
  spatial_features = spatial_features %>% group_by(day, year, camera_id, spot_id) %>% slice(1) %>% dplyr::select(success_feature, loss_feature, social_feature, knowledge_feature, day, year, camera_id, spot_id, confidence)
  patch_leaving_data = left_join(patch_leaving_data, spatial_features)
  
  # compute contiuous social feature across whole angling time
  # load gps data
  gps_data <- fread("data/gps_data.csv")
  
  # fast distance computation
  sourceCpp("functions/fast_dist.cpp")
  
  # bandwidth estimates
  spot_selection_fit = readRDS(file = "data/processed_data/model_fit/spot_selection_model/i_temporal_fit.rds")
  bandwidths = spot_selection_fit %>% 
    spread_draws(lambdas[feature]) %>% 
    group_by(feature) %>% 
    summarize(h = mean(lambdas))
  bandwidths_offset = spot_selection_fit %>% 
    spread_draws(v_lake[lake, feature]) %>% 
    filter(feature %in% c(14, 15)) %>% 
    mutate(feature = feature - 13) %>% 
    group_by(lake, feature) %>% 
    summarize(offset = mean(v_lake))
  bandwidths = left_join(bandwidths, bandwidths_offset)
  bandwidths = bandwidths %>% 
    mutate(bandwidth = exp(h + offset))
  
  # project coordinates
  coords_xy = LongLatToUTM(gps_data$position_long, gps_data$position_lat, zone = 35)
  gps_data$position_x = coords_xy[, 2]
  gps_data$position_y = coords_xy[, 3]
  gps_data = gps_data %>% filter(year == 2023)
  gps_data = gps_data %>% filter(competition == "individual")
  
  # compute time
  patch_leaving_data$time_sec = patch_leaving_data$time_start + ( (patch_leaving_data$bin - 1) * 10)
  patch_leaving_data$timestamp_sec = patch_leaving_data$timestamp_start + ( (patch_leaving_data$bin - 1) * 10)
  
  # lake key
  lake_key = patch_leaving_data %>% group_by(day, year) %>% slice(1) %>% dplyr::select(day, year, lake_index = unique_lake_id)
  
  # test
  #View(patch_leaving_data[67650:67660,])
  
  # write as function
  compute_social_feature <- function(timestamp_sec,
                                     unique_lake_id,
                                     camera_id,
                                     lake_key,
                                     gps_data,
                                     bandwidths){
    
    # define cpp function for cluster export
    sourceCpp("functions/fast_dist.cpp")
    
    # get day and year
    tmp = lake_key[lake_key$lake_index == unique_lake_id, ]
    day_tmp = tmp$day
    year_tmp = tmp$year
    
    # tmp id
    id_tmp = camera_id
    
    target_coordinates = gps_data %>% 
      filter(camera_id == id_tmp & year == year_tmp & day == day_tmp & timestamp == timestamp_sec) %>% 
      dplyr::select(position_x, position_y)
    
    # sometimes at end of competition, no GPS available due to end of competition (happens for roughly 300 datapoints / 170000 = 0.2%)
    if (nrow(target_coordinates)==0){
      target_coordinates = gps_data %>% 
        filter(camera_id == id_tmp & year == year_tmp & day == day_tmp & timestamp == (timestamp_sec-(12*10))) %>% 
        dplyr::select(position_x, position_y)
      
    }
    
    # extract coordinates social
    social_coordinates = gps_data %>% 
      filter(camera_id != id_tmp & year == year_tmp & day == day_tmp & timestamp == timestamp_sec) %>% 
      dplyr::select(position_x, position_y)
    
    # take social feature 2 minutes prior
    if (nrow(social_coordinates)<15){
     social_coordinates = gps_data %>% 
        filter(camera_id != id_tmp & year == year_tmp & day == day_tmp & timestamp == (timestamp_sec-(12*10))) %>% 
        dplyr::select(position_x, position_y)
    }
    # compute distances
    distances = fastPdist2(as.matrix(target_coordinates), as.matrix(social_coordinates))
    
    # extract bandwidth 
    bandwidth_tmp = bandwidths$bandwidth[bandwidths$lake == unique_lake_id & bandwidths$feature == 1]
    
    # compute kde
    return(kern(d = sort(distances)[1:15], h = bandwidth_tmp, option = "social", ncol = 15))
    
  }
  
  # check if continous social feature already computed, if not recompute
  if(!file.exists("data/processed_data/social_feature_continuous.rds")){
    
    # takes roughly 10 minutes on 50 cores
    cl <- makeClusterPSOCK(50)
    future::plan(cluster, workers = cl, gc = T)
    options(future.globals.maxSize= 6000*1024^2)
    
    start_time <- Sys.time()
    social_feature_angling <- future_pmap(
      list(
        patch_leaving_data$timestamp_sec,
        patch_leaving_data$unique_lake_id,
        patch_leaving_data$camera_id
      ),
      function(x, y, z) {
        compute_social_feature(
          x, y, z,
          lake_key,
          gps_data,
          bandwidths
        )
      },
      .progress = TRUE
    )
    end_time = Sys.time()
    end_time - start_time
    
    # stop cluster
    parallel::stopCluster(cl)
    
    # save
    social_feature_angling = unlist(social_feature_angling)
    saveRDS(social_feature_angling, file = "data/processed_data/social_feature_continuous.rds")
  }
  
  # read continuous social feature
  social_feature_angling_vector = readRDS(file = "data/processed_data/social_feature_continuous.rds")
  
  # add to data
  patch_leaving_data$social_feature_cont = social_feature_angling_vector
  
  # lake and id 
  lakes = patch_leaving_data$unique_lake_id
  ids = patch_leaving_data$unique_participant_id
  
  # test
  #plot(patch_leaving_data$social_feature_cont[500:1500], patch_leaving_data$social_feature[500:1500])
  
  # compute time in competition
  time_start = as.numeric(patch_leaving_data$time_start)
  
  max_times = data.frame(lake = lakes, 
                         time = time_start) %>% 
    ungroup() %>% 
    dplyr::group_by(lake) %>% 
    dplyr::summarize(max_time = max(time),
                     min_time = min(time))
  times = data.frame(lake = lakes, 
                     time = time_start)
  times = left_join(times, max_times)
  time_since_start_normalized = (times$time-times$min_time) / (times$max_time-times$min_time)
  patch_leaving_data$time_since_start_normalized = time_since_start_normalized
  
  # make stan data
  stan_data <- list(T = nrow(patch_leaving_data),
                    state = patch_leaving_data$state,
                    catch = patch_leaving_data$catch,
                    time_since_event = patch_leaving_data$time_since_event,
                    cumulative_catch = patch_leaving_data$cumulative_catch,
                    angling_time = patch_leaving_data$angling_time,
                    success_feature = scale(patch_leaving_data$success_feature)[,1],
                    loss_feature = scale(patch_leaving_data$loss_feature)[,1],
                    knowledge_feature = scale(patch_leaving_data$knowledge_feature)[,1],
                    social_feature = scale(patch_leaving_data$social_feature_cont)[,1],
                    time_since_start_normalized = time_since_start_normalized,
                    grainsize = 1,
                    confidence = patch_leaving_data$confidence,
                    N_trips = length(unique(patch_leaving_data$unique_trip_id)),
                    trips = 1:length(unique(patch_leaving_data$unique_trip_id)),
                    trip_start = trip_start$t,
                    trip_end = trip_end$t,
                    lakes = lakes, 
                    ids = ids,
                    N_ids = length(unique(ids)),
                    N_lakes = length(unique(lakes)))
  saveRDS(stan_data, file = "data/processed_data/stan_data_patch_leaving.rds")
  
  # make identifiers as for spatial model
  # need day, year, camera id, time bin, ids
  identifiers = data.frame(camera_id = patch_leaving_data$camera_id,
                           day = patch_leaving_data$day,
                           spot_id = patch_leaving_data$spot_id,
                           year = patch_leaving_data$year,
                           bin = patch_leaving_data$bin,
                           model_id = ids,
                           model_trips = patch_leaving_data$unique_trip_id,
                           model_lakes = lakes)
  
  # add unique participant id
  identifiers = left_join(identifiers, catch_data %>% dplyr::select(camera_id, sex, day, catch, year_of_birth, year, unique_id = participant_id))
  
  # save
  saveRDS(identifiers, file = "data/processed_data/identifiers_patch_leaving.rds")
  
}
################################################################################
# END
################################################################################
