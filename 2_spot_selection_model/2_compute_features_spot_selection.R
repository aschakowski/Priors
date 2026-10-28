################################################################################
#
# TITLE: Compute features spot selection
#
# AUTHORS: Schakowski, A. (schakowski@mpib-berlin.mpg.de)
#
# DESCRIPTION: This script computes spatial features for observed and simulated
# foraging spots.
#
# LAST UPDATED (DD/MM/YYYY): 10/10/2026
#
################################################################################

extract_features <- function(){
  
  ############### HELPER FUNCTIONS ###############################################
  # Extract social distances
  extract_social <- function(sub, gps_data, time_option){
    
    # subset spot identifier
    day_tmp = sub[1, "day"]
    year_tmp = sub[1, "year"]
    id_tmp = sub[1, "camera_id"]
    coord_tmp = sub[, c("x", "y")]
    time_tmp = sub[1, "timestamp_start"]
    
    # subset coordinates of other participants
    coordinates = gps_data %>% 
      filter(day == day_tmp &
               year == year_tmp &
               camera_id != id_tmp &
               timestamp == time_tmp - time_option) %>% 
      drop_na(position_lat, position_long) %>% 
      ungroup() %>% 
      dplyr::select(x = position_lat,y = position_long, camera_id = camera_id)
    
    if (nrow(coordinates)== 0){
      coordinates = gps_data %>% 
        filter(day == day_tmp &
                 year == year_tmp &
                 camera_id != id_tmp &
                 timestamp == time_tmp - 3*60) %>% 
        drop_na(position_lat, position_long) %>% 
        ungroup() %>% 
        dplyr::select(x = position_lat,y = position_long, camera_id = camera_id)
      
    }
    
    coordinates_xy = LongLatToUTM(coordinates$y, coordinates$x, zone = 35)
    coordinates_xy$ID = coordinates$camera_id
    
    # compute distances to target
    tmp_distances = fastPdist2(as.matrix(coord_tmp), as.matrix(coordinates_xy[, 2:3]))
    
    return(tmp_distances)
    
  }
  
  # compute distance to successes and losses
  extract_spot_distance <- function(sub, spot_selection, option){ #option either 0 or 1
    
    # subset spot identifier
    day_tmp = sub[1, "day"]
    year_tmp = sub[1, "year"]
    id_tmp = sub[1, "camera_id"]
    coord_tmp = sub[, c("x", "y")]
    spot_id_tmp = sub[1, "spot_id"]
    
    # get coordinates of last spots
    if (option == "all"){
      spots = spot_selection %>% 
        ungroup() %>% 
        filter(day == day_tmp & year == year_tmp & camera_id == id_tmp & spot_id<spot_id_tmp) %>% 
        dplyr::select(x, y, spot_id)
    } else {
      spots = spot_selection %>% 
        ungroup() %>% 
        filter(day == day_tmp & year == year_tmp & camera_id == id_tmp & spot_id<spot_id_tmp & Return == option) %>% 
        dplyr::select(x, y, spot_id)
    }
    
    if(nrow(spots) == 0){
      result = rep(-99, nrow(sub))
    } else {
      result = fastPdist2(as.matrix(coord_tmp), as.matrix(spots[, c("x", "y")]))
      
      # sometimes missing, if distances very close
      result[is.na(result)]<-.00001
      
      # keep only 50 closest features for storage
      if (dim(result)[2] == 1){
        result = result
      } else {
        result = t(apply(result, 1, sort))[, 1:min(ncol(result), 50)] 
      }
      
    }
    
    return(result)
    
  }
  
  # extract distance to last spot
  extract_locality_feature <- function(sub, spot_selection_data, starting_locations){
    
    # define cpp function for cluster export
    #sourceCpp("functions/fast_dist.cpp")
    
    # subset spot identifier
    day_tmp = sub[1, "day"]
    year_tmp = sub[1, "year"]
    id_tmp = sub[1, "camera_id"]
    coord_tmp = sub[, c("x", "y")]
    spot_id_tmp = sub[1, "spot_id"]
    
    # get coordinates of last spots
    spots = spot_selection_data %>% 
      ungroup() %>% 
      filter(day == day_tmp & year == year_tmp & camera_id == id_tmp & spot_id==(spot_id_tmp-1)) %>% 
      dplyr::select(x, y, spot_id, Return)
    
    if(spot_id_tmp == 1){
      result = fastPdist2(as.matrix(coord_tmp), t(as.matrix(st_coordinates(starting_locations)[day_tmp,])))
    } else {
      result = fastPdist2(as.matrix(coord_tmp), as.matrix(spots[nrow(spots), c("x", "y")]))
    }
    
    # sometimes missing, if distances very close
    result[is.na(result)]<-0
    
    return(result)
    
  }
  
  # compute distance for prior nomination
  extract_prior_knowledge_feature <- function(sub, ecological_knowledge_polygons){
    
    # extract identifiers
    day_tmp = sub[1, "day"]
    year_tmp = sub[1, "year"]
    id_tmp = sub[1, "camera_id"]
    coord_tmp = sub[, c("x", "y")]
    spot_id_tmp = sub[1, "spot_id"]
    
    # extract polygons of ID
    polygons = ecological_knowledge_polygons %>% filter(day == day_tmp & camera_id == id_tmp)
    
    # filter by option
    polygons = polygons %>% filter(confidence %in% option)
    
    # if no nomination, then everything set to 0 (will not affect estimates)
    if (nrow(polygons) == 0){
      dist = rep(0, nrow(sub))
    } else {
      # convert fm coordinates to sf
      coords = st_as_sf(sub[, c("x", "y")], coords = c("x","y"))
      st_crs(coords) <- crs(polygons)
      
      # compute distance
      dist = st_distance(coords, polygons)
      dist = t(apply(dist, 1, sort))
      
      if (is.null(dim(dist))){
        dist = t(t(dist))
      }
      
      if (dim(dist)[1] == 1){
        dist = t(dist)
      }
      
    }
    
    # return
    return(dist)
    
  }
  
  
  # compute confidence
  extract_confidence <- function(sub, ecological_knowledge_polygons){
    
    # extract identifiers
    day_tmp = sub[1, "day"]
    id_tmp = sub[1, "camera_id"]
    coord_tmp = sub[, c("x", "y")]
    spot_id_tmp = sub[1, "spot_id"]
    
    # extract polygons of ID
    polygons = ecological_knowledge_polygons %>% filter(day == day_tmp & camera_id == id_tmp)
    
    # if no nomination, then confidence is missing (-99)
    if (nrow(polygons) == 0){
      dist = rep(0, nrow(sub))
      confidence = -99
    } else {
      # convert fm coordinates to sf
      coords = st_as_sf(sub[, c("x", "y")], coords = c("x","y"))
      st_crs(coords) <- crs(polygons)
      
      # compute distance
      dist = st_distance(coords, polygons)
      
      confidence = polygons$confidence[apply(dist, 1, which.min)]
      
    }
    
    # return
    return(confidence)
    
  }
  
  ################################################################################
  # load spot selection data
  spot_selection_data = fread(file = "data/spot_selection_data_cleaned.csv")
  
  # Load feature matrix
  Feature_matrix = readRDS(file = "data/processed_data/feature_matrix_spot_selection.rds")
  
  # load identifiers
  identifiers = readRDS("data/processed_data/identifiers_spot_selection.rds")
  
  # Import GPS Data
  gps_data <- fread("data/gps_data.csv")
  
  # arrange data chronologically
  spot_selection_data = spot_selection_data %>% 
    arrange(year, day, camera_id, spot_id)
  
  # define cpp function for cluster export
  sourceCpp("functions/fast_dist.cpp")
  
  ################################################################################
  # extract social distances
  ################################################################################
  # at time of joining
  social_distance_matrix = map(Feature_matrix, function(x) extract_social(x, gps_data, time_option = 0), .progress = TRUE)
  
  # save social distances
  saveRDS(social_distance_matrix, file = "data/processed_data/social_distance_matrix.rds")
  
  rm(social_distance_matrix)
  gc()
  gc()
  
  ################################################################################
  # Extract spot distances
  ################################################################################
  # extract feature
  success_distance_matrix = map(Feature_matrix, function(x) extract_spot_distance(x, spot_selection_data, option = 1), .progress = TRUE)
  loss_distance_matrix = map(Feature_matrix, function(x) extract_spot_distance(x, spot_selection_data, option = 0), .progress = TRUE)
  
  saveRDS(success_distance_matrix, file = "data/processed_data/success_distance_matrix.rds")
  saveRDS(loss_distance_matrix, file = "data/processed_data/loss_distance_matrix.rds")
  
  ################################################################################
  # extract locality feature
  ################################################################################
  # compute starting locations from gps data
  gps_data <- fread("data/gps_data.csv")
  
  starting_locations = gps_data %>% 
    group_by(day, camera_id) %>% 
    arrange(date) %>% 
    slice(1) %>% 
    group_by(day) %>% 
    summarize(starting_x = mean(position_long),
              starting_y = mean(position_lat))
  starting_coords = LongLatToUTM(starting_locations$starting_x, starting_locations$starting_y, zone = 35)  
  starting_locations$starting_x = starting_coords[, 2]
  starting_locations$starting_y = starting_coords[, 3]
  starting_locations = st_as_sf(starting_locations, coords = c("starting_x", "starting_y"))
  st_crs(starting_locations) = "epsg:3067"
  
  locality_matrix = list()
  for (i in 1:length(Feature_matrix)){
    locality_matrix[[i]] <- extract_locality_feature(Feature_matrix[[i]], spot_selection_data, starting_locations)
  }
  
  saveRDS(locality_matrix, file = "data/processed_data/locality_matrix.rds")
  
  rm(locality_matrix)
  gc()
  
  ################################################################################
  # extract ecological knowledge feature
  ################################################################################
  ecological_knowledge_polygons = readRDS("data/ecological_knowledge_polygons.rds")
  knowledge_distance_matrix = map(Feature_matrix, function(x) extract_prior_knowledge_feature(x, ecological_knowledge_polygons), .progress = TRUE)
  
  # save result
  saveRDS(knowledge_distance_matrix, file = "data/processed_data/knowledge_distance_matrix.rds")
  rm(knowledge_distance_matrix)
  gc()
  
  ##############################################################################
  # extract confidence matrix
  ##############################################################################
  confidence_matrix = map(Feature_matrix, function(x) extract_confidence(x, ecological_knowledge_polygons), .progress = TRUE)
  confidence_matrix = do.call(cbind, confidence_matrix)
  
  # save result
  saveRDS(confidence_matrix, file = "data/processed_data/confidence_matrix.rds")
  rm(confidence_matrix)
  gc()
  
  ##############################################################################
  # extract roughness feature
  ##############################################################################
  # Import Lake Data
  lake_data = readRDS("data/depth_profiles.rds")
  
  # depth profile points
  depth_profile_points = list()
  for (i in 1:5){
    tmp = st_line_sample(lake_data[[i]] %>% filter(depth>0), density = .1, type = "regular")
    depth_profile_points[[i]] = tmp[!st_is_empty(tmp)] %>% st_cast("POINT")
  }
  
  # extract feature
  roughness_distance_matrix = map(Feature_matrix, function(x) extract_roughness(x, depth_profile_points), .progress = TRUE)
  roughness_distance_matrix = do.call(cbind, roughness_distance_matrix)
  
  # save result
  saveRDS(roughness_distance_matrix, file = "data/processed_data/roughness_distance_matrix.rds")
  rm(roughness_distance_matrix)
  gc()
  
}

################################################################################
# END
################################################################################
