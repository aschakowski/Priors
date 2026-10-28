################################################################################
#
# TITLE: Extract ecological properties
# 
# AUTHORS: Schakowski, A. (schakowski@mpib-berlin.mpg.de)
#
# DESCRIPTION: For each actual and randomly rotated area, this script computes
# a set of spatial features for comparison.
#
# LAST UPDATED (DD/MM/YYYY): 06/10/2026
#
################################################################################
extract_ecological_properties <- function(){
  
  # extract cluster properties from depth curves
  # depth profiles
  depth_profiles = readRDS("data/depth_profiles.rds")
  
  # load simulated clusters
  simulated_clusters = readRDS("data/output/simulated_clusters.rds")
  
  # check whether translation has worked
  #simulated_clusters %>% 
  #  filter(unique_id == 1) %>% 
  #  ggplot(aes(color = condition)) + 
  #  geom_sf()
  
  #simulated_clusters %>% 
  #  filter(unique_id == 25) %>% 
  #  ggplot(aes(color = condition)) + 
  #  geom_sf()
  
  # competition areas
  competition_areas = readRDS("data/study_lakes.rds")
  
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
  
  # for faster execution time: 
  dt_sim_clusters = as.data.table(simulated_clusters)
  
  # make function for parallization
  extract_ecology <- function(cluster, competition_areas, starting_locations, depth_profiles){
    
    library(dplyr)
    library(data.table)
    
    # day
    day_tmp = as.numeric(cluster$day[1])
    
    # extract depth profile and transform, if not correct crs
    profile = depth_profiles[[day_tmp]]
    profile_transformed = st_transform(profile, "epsg:3067")
    
    # extract competition boundary
    boundary = competition_areas[[day_tmp]]
    boundary_transformed = st_transform(boundary, "epsg:3067")
    boundary_transformed = boundary_transformed %>% st_cast("MULTIPOINT")
    
    # compute area
    area = drop_units(st_area(cluster$x))
    
    # sample points from polygon (1pt per 100m^2)
    sample_area = round(area/(100))
    if(sample_area < 5){
      sample_area = 5
    }
    
    samples = st_sample(cluster$x[1], size = sample_area, type = "regular")
    
    # plot cluster and profile as check
    #ggplot() + 
    #  geom_sf(data = competition_areas[[day_tmp]]) +
    # geom_sf(data = profile_transformed) + 
    # geom_sf(data = simulated_clusters %>% filter(unique_id == 1)) +
    # geom_sf(data = simulated_clusters[1,], color = "red") + 
    # geom_sf(data = samples, color = "blue", alpha = .2) + 
    #  geom_sf(data = starting_locations %>% filter(day ==1))
    
    # extract ecological properties of cluster coordinates from depth profile
    
    # depth: which is the closest curve
    profile_transformed = st_transform(profile_transformed, st_crs(samples))
    distance = st_distance(x = samples, y = profile_transformed)#rows are points, cols are curves
    distance_min = apply(distance, 1, FUN = which.min)#which col/curve contains minimum for each pt
    depth = mean(profile_transformed$depth[distance_min]) # depth as the average depth of the closest curve
    
    # distance to closest depth curve (= inv. roughness)
    roughness = st_distance(x = samples, y = profile_transformed)
    roughness = mean(apply(roughness, 1, FUN = min))
    
    # distance to starting point
    starting_locations = st_transform(starting_locations, st_crs(samples))
    displacement = st_distance(x = samples, y = starting_locations)
    displacement = mean(apply(displacement, 1, FUN = min))#avg distance to starting location across all pts
    
    # distance to shoreline (depth curve with depth 0)
    shore = st_distance(x = samples, y = profile_transformed[profile_transformed$depth == 0, ])
    shore = mean(apply(shore, 1, FUN = min))
    
    # distance to boundary (boundary of competition area)
    boundary_dist = st_distance(x = samples, y = boundary_transformed)
    boundary_dist = mean(apply(boundary_dist, 1, FUN = min))
    
    # add to original data frame
    cluster$depth = as.numeric(depth)
    cluster$roughness = as.numeric(roughness)
    cluster$shore = as.numeric(shore)
    cluster$displacement = as.numeric(displacement)
    cluster$boundary = as.numeric(boundary_dist)
    
    # add to return list
    return(cluster)
    
  }
  
  ################################################################################
  # Parellelize
  ################################################################################
  
  # split dataframe into list
  dt_sim_clusters$batch = 1:nrow(dt_sim_clusters)
  out <- split(dt_sim_clusters, dt_sim_clusters$batch) 
  
  # compute ecological variables
  cl <- makeClusterPSOCK(30)
  future::plan(future::cluster, workers = cl, gc = T)
  options(future.globals.maxSize = 500*1024^2)
  eco = future_map(out, function(x) extract_ecology(x, competition_areas, starting_locations, depth_profiles),
                         .progress = TRUE,
                         .options = furrr_options(seed = TRUE))
  options(future.globals.maxSize = 500*1024^2)
  clusters_full = dplyr::bind_rows(eco)
  parallel::stopCluster(cl)
  
  ################################################################################
  # Postprocessing
  ################################################################################
  # standardize predictors
  clusters_full = clusters_full %>% 
    group_by(day, camera_id) %>% 
    mutate(sd_depth = sd(depth),
           mean_depth = mean(depth),
           sd_roughness = sd(roughness),
           mean_roughness = mean(roughness),
           sd_shore = sd(shore),
           mean_shore = mean(shore),
           sd_boundary = sd(boundary),
           mean_boundary = mean(boundary),
           sd_displacement = sd(displacement),
           mean_displacement = mean(displacement)) %>% 
    mutate(depth_scaled = (depth - mean_depth) / sd_depth,
           roughness_scaled = (roughness - mean_roughness)/sd_roughness,
           shore_scaled = (shore - mean_shore) / sd_shore,
           displacement_scaled = (displacement - mean_displacement) / sd_displacement,
           boundary_scaled = (boundary - mean_boundary) / sd_boundary)
  saveRDS(clusters_full, file = "data/output/simulated_clusters_eco.rds")
  
}
################################################################################
# END
################################################################################
