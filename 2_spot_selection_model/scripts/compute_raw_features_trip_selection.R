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
    day_tmp = pull(sub[1, "day"])
    year_tmp = 2023
    id_tmp = pull(sub[1, "camera_id"])
    coord_tmp = sub[, c("x", "y")]
    time_tmp = pull(sub[1, "timestamp_start"])
    
    # subset coordinates of other participants
    coordinates = gps_data %>% 
      ungroup() %>% 
      dplyr::filter(day == day_tmp) %>% 
      filter(year == year_tmp) %>% 
      filter(camera_id != id_tmp) %>% 
      filter(timestamp == time_tmp - time_option) %>% 
      drop_na(position_lat, position_long) %>% 
      ungroup() %>% 
      dplyr::select(x = position_lat,y = position_long, ID = camera_id)
    
    coordinates_xy = LongLatToUTM(coordinates$y, coordinates$x, zone = 35)
    coordinates_xy$ID = coordinates$ID
    
    # compute distances to target
    tmp_distances = fastPdist2(as.matrix(coord_tmp), as.matrix(coordinates_xy[, 2:3]))
    
    return(tmp_distances)
    
  }
  
  # compute distance for prior nomination
  extract_prior_knowledge_feature <- function(sub, ecological_knowledge_polygons, sort = TRUE){
    
    # extract identifiers
    day_tmp = pull(sub[1, "day"])
    year_tmp = 2023
    id_tmp = pull(sub[1, "camera_id"])
    coord_tmp = sub[, c("x", "y")]
    spot_id_tmp = pull(sub[1, "spot_id"])
    
    # extract polygons of ID
    polygons = ecological_knowledge_polygons %>% dplyr::filter(day == day_tmp) %>% filter(camera_id == id_tmp)
    
    # if no nomination, then everything set to 0
    if (nrow(polygons) == 0){
      dist = rep(-99, nrow(sub))
    } else {
      # convert fm coordinates to sf
      coords = st_as_sf(sub[, c("x", "y")], coords = c("x","y"))
      st_crs(coords) <- crs(polygons)
      
      # compute distance
      dist = st_distance(coords, polygons) %>% drop_units()
      
      if (sort == TRUE){
        dist = t(apply(dist, 1, sort)) 
      }
      
    }
    
    # return
    return(dist)
    
  }
  ################################################################################
  # Load feature matrix
  Feature_matrix = readRDS(file = "data/processed_data/feature_matrix_trip_selection.rds")
  
  # define cpp function for cluster export
  sourceCpp("functions/fast_dist.cpp")
  
  ################################################################################
  # extract social distances
  ################################################################################
  # at time of joining
  social_distance_matrix = map(Feature_matrix, function(x) extract_social(x, gps_data, time_option = 0), .progress = TRUE)
  
  # save social distances
  saveRDS(social_distance_matrix, file = "data/processed_data/social_distance_matrix_trip_selection.rds")
  
  rm(social_distance_matrix)
  gc()
  gc()
  
  ################################################################################
  # extract ecological knowledge feature
  ################################################################################
  ecological_knowledge_polygons = readRDS("data/ecological_knowledge_polygons.rds")
  knowledge_distance_matrix = map(Feature_matrix, function(x) extract_prior_knowledge_feature(x, ecological_knowledge_polygons, sort = FALSE), .progress = TRUE)
  nas = vector()
  for (i in 1:length(knowledge_distance_matrix)){
    nas[i] = sum(!is.finite(knowledge_distance_matrix[[i]]))
  }
  
  # check dimensions
  dims = matrix(ncol = 2, nrow = length(knowledge_distance_matrix))
  for (i in 1:length(knowledge_distance_matrix)){
    
    if (is.null(dim(knowledge_distance_matrix[[i]]))){
      knowledge_distance_matrix[[i]] = as.matrix(knowledge_distance_matrix[[i]])
    } else if (dim(knowledge_distance_matrix[[i]])[2] == 51){
      knowledge_distance_matrix[[i]] = t(knowledge_distance_matrix[[i]])
    }
    dims[i,] <- dim(knowledge_distance_matrix[[i]])   
  }
  
  # save result
  saveRDS(knowledge_distance_matrix, file = "data/processed_data/knowledge_distance_matrix_trip_selection.rds")
  rm(knowledge_distance_matrix)
  gc()
  
}

################################################################################
# END
################################################################################
