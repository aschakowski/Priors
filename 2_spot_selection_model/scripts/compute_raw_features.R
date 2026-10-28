################################################################################
#
# Title: Compute raw features
#
# Description: Computes raw feature distances based on feature matrix
#
# Authors: Schakowski, A.
#
# Last updated: 13/02/2024 (DD/MM/YYYY)
#
################################################################################
extract_features <- function(gps_data, spot_selection_data){
############### HELPER FUNCTIONS ###############################################
# Extract social distances
extract_social <- function(sub, gps_data, time_option){
  
  # subset spot identifier
  day_tmp = sub[1, "day"]
  year_tmp = sub[1, "year"]
  id_tmp = sub[1, "ID"]
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
    dplyr::select(x = position_lat,y = position_long, ID = camera_id)
  
  coordinates_xy = LongLatToUTM(coordinates$y, coordinates$x, zone = 35)
  coordinates_xy$ID = coordinates$ID
  
  # compute distances to target
  tmp_distances = fastPdist2(as.matrix(coord_tmp), as.matrix(coordinates_xy[, 2:3]))
  
  return(tmp_distances)
  
}

# extract roughness
extract_roughness <- function(sub, depth_profile_points, lake_key){
  
  sub = left_join(sub, lake_key, by = c("day", "year"))
  points = st_coordinates(depth_profile_points[[sub$lake_index[1]]])
  
  coord_tmp = sub[, c("x", "y")]
  depth_profile = depth_profiles[[sub$lake_index[1]]]
  
  # compute distances to target
  d = fastPdist2(as.matrix(coord_tmp), as.matrix(points[, 1:2]))
  
  # if very small sometimes na
  d[which(is.na(d))] <- .00001
  
  # keep only 50 closest features for storage
  d = t(apply(d, 1, sort))[, 1:50]
  return(d)
  
}

# compute distance to successes and losses
extract_spot_distance <- function(sub, spot_selection, option){ #option either 0 or 1
  
  # define cpp function for cluster export
  #sourceCpp("functions/fast_dist.cpp")
  
  # subset spot identifier
  day_tmp = sub[1, "day"]
  year_tmp = sub[1, "year"]
  id_tmp = sub[1, "ID"]
  coord_tmp = sub[, c("x", "y")]
  spot_id_tmp = sub[1, "spot_id"]
  
  # get coordinates of last spots
  
  if (option == "all"){
    spots = spot_selection %>% 
      ungroup() %>% 
      filter(day == day_tmp & year == year_tmp & ID == id_tmp & spot_id<spot_id_tmp) %>% 
      dplyr::select(x, y, spot_id)
  } else {
    spots = spot_selection %>% 
      ungroup() %>% 
      filter(day == day_tmp & year == year_tmp & ID == id_tmp & spot_id<spot_id_tmp & Return == option) %>% 
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
extract_locality_feature <- function(sub, spot_selection_data){
  
  # define cpp function for cluster export
  #sourceCpp("functions/fast_dist.cpp")
  
  # subset spot identifier
  day_tmp = sub[1, "day"]
  year_tmp = sub[1, "year"]
  id_tmp = sub[1, "ID"]
  coord_tmp = sub[, c("x", "y")]
  spot_id_tmp = sub[1, "spot_id"]
  
  # get coordinates of last spots
  spots = spot_selection_data %>% 
    ungroup() %>% 
    filter(day == day_tmp & year == year_tmp & ID == id_tmp & spot_id<spot_id_tmp) %>% 
    dplyr::select(x, y, spot_id, Return)
  
  if(nrow(spots) == 0){
    result = rep(-99, nrow(sub))
  } else {
    result = fastPdist2(as.matrix(coord_tmp), as.matrix(spots[nrow(spots), c("x", "y")]))
  }
  
  # sometimes missing, if distances very close
  result[is.na(result)]<-0
  
  return(result)
  
}

# compute distance for prior nomination
extract_prior_knowledge_feature <- function(sub, ecological_knowledge_polygons, option){
  
  # extract identifiers
  day_tmp = sub[1, "day"]
  year_tmp = sub[1, "year"]
  id_tmp = sub[1, "ID"]
  coord_tmp = sub[, c("x", "y")]
  spot_id_tmp = sub[1, "spot_id"]
  
  # extract polygons of ID
  polygons = ecological_knowledge_polygons %>% filter(day == day_tmp & camera_id == id_tmp)
  
  # filter by option
  polygons = polygons %>% filter(label %in% option)
  
  # if no nomination, then everything set to 0
  if (nrow(polygons) == 0){
    dist = rep(0, nrow(sub))
  } else {
    # convert fm coordinates to sf
    coords = st_as_sf(sub[, c("x", "y")], coords = c("x","y"))
    st_crs(coords) <- crs(polygons)
    
    # compute polygon centroids
    #centroids = list()
    #for (p in 1:nrow(polygons)){
    #  centroids[[p]] = st_centroid(polygons[p,])
    #}
    #centroids = do.call(rbind, centroids)
    
    # compute distance
    dist = st_distance(coords, polygons)
    dist = t(apply(dist, 1, sort))
    
    if (is.null(dim(dist))){
      dist = t(t(dist))
    }

  }
  
  # return
  return(dist)
  
}


# compute confidence
extract_confidence <- function(sub, ecological_knowledge_polygons){
  
  # extract identifiers
  day_tmp = sub[1, "day"]
  year_tmp = sub[1, "year"]
  id_tmp = sub[1, "ID"]
  coord_tmp = sub[, c("x", "y")]
  spot_id_tmp = sub[1, "spot_id"]
  
  # extract polygons of ID
  polygons = ecological_knowledge_polygons %>% filter(day == day_tmp & camera_id == id_tmp)
  
  # if no nomination, then everything set to 0
  if (nrow(polygons) == 0){
    dist = rep(0, nrow(sub))
    confidence = 1
  } else {
    # convert fm coordinates to sf
    coords = st_as_sf(sub[, c("x", "y")], coords = c("x","y"))
    st_crs(coords) <- crs(polygons)
    
    # compute polygon centroids
    #centroids = list()
    #for (p in 1:nrow(polygons)){
    #  centroids[[p]] = st_centroid(polygons[p,])
    #}
    #centroids = do.call(rbind, centroids)
    
    # compute distance
    dist = st_distance(coords, polygons)
    
    confidence = polygons$label[apply(dist, 1, which.min)]
    
  }
  
  
  # which is the closest nomination
  
  
  # return
  return(confidence)
  
}
################################################################################
# Load feature matrix
Feature_matrix = readRDS(file = "data/processed_data/feature_matrix_spot_selection.rds")

# load identifiers
identifiers = readRDS("data/processed_data/identifiers_spot_selection.rds")

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
# Extract Roughness
################################################################################
# Import Lake Data
lake_key = data.frame(day = c(1:5, 1:3, 5,4), 
                      year = rep(c(2023, 2022), each = 5),
                      lake_index = 1:10)

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

# depth profile points
depth_profile_points = list()
for (i in 1:10){
  tmp = st_line_sample(depth_profiles[[i]] %>% filter(depth>0), density = .1, type = "regular")
  depth_profile_points[[i]] = tmp[!st_is_empty(tmp)] %>% st_cast("POINT")
}

# extract feature
roughness_distance_matrix = map(Feature_matrix, function(x) extract_roughness(x, depth_profile_points, lake_key), .progress = TRUE)

# save result
saveRDS(roughness_distance_matrix, file = "data/processed_data/roughness_distance_matrix.rds")
rm(roughness_distance_matrix)
gc()

################################################################################
# Extract spot distances
################################################################################
# extract feature
success_distance_matrix = map(Feature_matrix, function(x) extract_spot_distance(x, spot_selection_data, option = 1), .progress = TRUE)
loss_distance_matrix = map(Feature_matrix, function(x) extract_spot_distance(x, spot_selection_data, option = 0), .progress = TRUE)
spot_distance_matrix = map(Feature_matrix, function(x) extract_spot_distance(x, spot_selection_data, option = "all"), .progress = TRUE)

saveRDS(success_distance_matrix, file = "data/processed_data/success_distance_matrix.rds")
saveRDS(loss_distance_matrix, file = "data/processed_data/loss_distance_matrix.rds")
saveRDS(spot_distance_matrix, file = "data/processed_data/spot_distance_matrix.rds")

################################################################################
# extract locality feature
################################################################################
locality_matrix = list()
for (i in 1:length(Feature_matrix)){
  locality_matrix[[i]] <- extract_locality_feature(Feature_matrix[[i]], spot_selection_data)
}

saveRDS(locality_matrix, file = "data/processed_data/locality_matrix.rds")

rm(locality_matrix)
gc()

################################################################################
# extract ecological knowledge feature
################################################################################
ecological_knowledge_polygons = readRDS("data/ecological_knowledge_polygons.rds")
knowledge_distance_matrix = map(Feature_matrix, function(x) extract_prior_knowledge_feature(x, ecological_knowledge_polygons, option = c(1, 2, 3)), .progress = TRUE)
nas = vector()
for (i in 1:length(knowledge_distance_matrix)){
  nas[i] = sum(!is.finite(knowledge_distance_matrix[[i]]))
}

# save result
saveRDS(knowledge_distance_matrix, file = "data/processed_data/knowledge_distance_matrix.rds")
rm(knowledge_distance_matrix)
gc()

# extract confidence matrix
confidence_matrix = map(Feature_matrix, function(x) extract_confidence(x, ecological_knowledge_polygons), .progress = TRUE)
confidence_matrix = do.call(cbind, confidence_matrix)

# save result
saveRDS(confidence_matrix, file = "data/processed_data/confidence_matrix.rds")
rm(confidence_matrix)
gc()


# confidence level 1
knowledge_distance_matrix_1 = map(Feature_matrix, function(x) extract_prior_knowledge_feature(x, ecological_knowledge_polygons, option = 1), .progress = TRUE)
nas = vector()
for (i in 1:length(knowledge_distance_matrix_1)){
  nas[i] = sum(!is.finite(knowledge_distance_matrix_1[[i]]))
}

# save result
saveRDS(knowledge_distance_matrix_1, file = "data/processed_data/knowledge_distance_matrix_1.rds")
rm(knowledge_distance_matrix_1)
gc()

# confidence level 2
knowledge_distance_matrix_2 = map(Feature_matrix, function(x) extract_prior_knowledge_feature(x, ecological_knowledge_polygons, option = 2), .progress = TRUE)
nas = vector()
for (i in 1:length(knowledge_distance_matrix_2)){
  nas[i] = sum(!is.finite(knowledge_distance_matrix_2[[i]]))
}

# save result
saveRDS(knowledge_distance_matrix_2, file = "data/processed_data/knowledge_distance_matrix_2.rds")
rm(knowledge_distance_matrix_2)
gc()

# confidence level 3
knowledge_distance_matrix_3 = map(Feature_matrix, function(x) extract_prior_knowledge_feature(x, ecological_knowledge_polygons, option = 3), .progress = TRUE)
nas = vector()
for (i in 1:length(knowledge_distance_matrix_3)){
  nas[i] = sum(!is.finite(knowledge_distance_matrix_3[[i]]))
}

# save result
saveRDS(knowledge_distance_matrix_3, file = "data/processed_data/knowledge_distance_matrix_3.rds")
rm(knowledge_distance_matrix_3)
gc()

}

################################################################################
# END
################################################################################
