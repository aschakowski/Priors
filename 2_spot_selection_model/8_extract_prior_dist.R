################################################################################
#
# TITLE: Extract prior distance
#
# AUTHORS: Schakowski, A. (schakowski@mpib-berlin.mpg.de)
#
# DESCRIPTION: This script computes prior distance for all actual and simulated
# trips to compute the descriptive measure of prior use
#
# LAST UPDATED (DD/MM/YYYY): 10/10/2026
#
################################################################################
extract_prior_dist <- function(){
  
  ############### HELPER FUNCTIONS ###############################################
  # compute distance for prior nomination
  extract_prior_knowledge_feature <- function(sub, ecological_knowledge_polygons){
    
    # extract identifiers
    day_tmp = sub[1, "day"]
    id_tmp = sub[1, "camera_id"]
    coord_tmp = sub[, c("x", "y")]
    spot_id_tmp = sub[1, "spot_id"]
    
    # extract polygons of ID
    polygons = ecological_knowledge_polygons %>% dplyr::filter(day == pull(day_tmp) & camera_id == pull(id_tmp))
    
    # if no nomination, then everything set to NA
    if (nrow(polygons) == 0){
      dist = rep(NA, nrow(sub))
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
  
  # compute confidence of prior nomination
  extract_confidence <- function(sub, ecological_knowledge_polygons){
    
    # extract identifiers
    day_tmp = sub[1, "day"]
    id_tmp = sub[1, "camera_id"]
    coord_tmp = sub[, c("x", "y")]
    spot_id_tmp = sub[1, "spot_id"]
    
    # extract polygons of ID
    polygons = ecological_knowledge_polygons %>% filter(day == pull(day_tmp) & camera_id == pull(id_tmp))
    
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
  Feature_matrix = readRDS(file = "data/processed_data/feature_matrix_trip_selection.rds")
  
  # define cpp function for cluster export
  sourceCpp("functions/fast_dist.cpp")
  
  ################################################################################
  # extract ecological knowledge feature
  ################################################################################
  ecological_knowledge_polygons = readRDS("data/ecological_knowledge_polygons.rds")
  knowledge_distance_matrix = map(Feature_matrix, function(x) extract_prior_knowledge_feature(x, ecological_knowledge_polygons), .progress = TRUE)
  
  # save result
  saveRDS(knowledge_distance_matrix, file = "data/processed_data/knowledge_distance_matrix_trip_selection.rds")
  rm(knowledge_distance_matrix)
  gc()
  
  ################################################################################
  # extract confidence feature
  ################################################################################
  confidence_matrix = map(Feature_matrix, function(x) extract_confidence(x, ecological_knowledge_polygons), .progress = TRUE)
  confidence_matrix = do.call(cbind, confidence_matrix)
  
  # save result
  saveRDS(confidence_matrix, file = "data/processed_data/confidence_matrix_trip_selection.rds")
  rm(confidence_matrix)
  gc()
  
}

################################################################################
# END
################################################################################
