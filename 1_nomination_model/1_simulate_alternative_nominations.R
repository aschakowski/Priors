################################################################################
#
# TITLE: Simulate alternative nominated areas 
#
# AUTHORS: Schakowski, A. (schakowski@mpib-berlin.mpg.de)
#
# DESCRIPTION: For each nominated area, this script creates a set of 50 random
# areas through random shifting and rotating.
#
# LAST UPDATED (DD/MM/YYYY): 06/10/2026
#
################################################################################
simulate_alternative_nominations <- function(){

  # helper rotate polygon from SpatialEco package
  rotate.polygon <- function(p, angle = 180, sp = FALSE,  
                             anchor = c("center", "lower.left", "upper.right")) {
    if(class(p)[1] != "sf") p <- sf::st_as_sf(p)
    p.coords <- sf::st_coordinates(p)[,1:2]
    if(anchor[1] == "center") {
      p.center <- suppressWarnings(sf::st_coordinates(sf::st_centroid(p)))
    } else if(anchor[1] == "lower.left") {
      p.center <- c(min(sf::st_coordinates(p)[,1]),
                    min(sf::st_coordinates(p)[,2]))
    } else if(anchor[1] == "upper.right") { 
      p.center <- c(max(sf::st_coordinates(p)[,1]),
                    max(sf::st_coordinates(p)[,2])) 
    }  
    rotate.coords <- function(xy, a, center) {
      co <- cos(-a * pi / 180)
      si <- sin(-a * pi / 180)
      adj <- matrix(rep(center, nrow(xy)), ncol=2, byrow=TRUE)
      xy <- xy-adj
      cbind(co * xy[,1] - si * xy[,2],si * xy[,1] + co * xy[,2]) + adj
    }
    p.rotate <- rotate.coords(p.coords, a = angle, center = c(p.center[1], 
                                                            p.center[2]))
    sf::st_geometry(p) <- sf::st_sfc(sf::st_polygon(list(p.rotate)))
    if(sp) p <- sf::as_Spatial(p)
    
    return( p )	
    
  }

  # Create random shift of nominated areas within competition area
  # Step 1: Load nominated areas polygons for each day
  areas = readRDS(file = "data/ecological_knowledge_polygons.rds")
  
  # assign missing values to confidence ratings
  areas = areas %>% 
    mutate(confidence = ifelse(camera_id == 37 & day == 2, NA, confidence),
           color = ifelse(camera_id == 37 & day == 2, NA, color)) %>% 
    mutate(confidence = ifelse(camera_id == 8 & day == 5, NA, confidence),
           color = ifelse(camera_id == 8 & day == 5, NA, color))
  
  # save to original file
  saveRDS(areas, file = "data/ecological_knowledge_polygons.rds")
  
  # filter out areas with missing confidence ratings
  areas = areas %>% 
    filter(!is.na(confidence))
  
  # filter out areas with too small areas
  area_size = st_area(areas) %>% drop_units() #in sqm
  areas = areas[area_size > 10,]
  
  # filter out areas that clearly lie outside the competition arena.
  # load competition areas
  competition_areas = readRDS("data/study_lakes.rds")
  
  intersection = vector()
  for (i in 1:nrow(areas)){
    
    # which day?
    day_tmp = areas[i, "day"] %>% st_drop_geometry()
    
    # competition arena
    arena = competition_areas[[as.numeric(day_tmp)]]
    tmp = (st_area(st_intersection(st_geometry(arena), areas[i,])) / st_area(areas[i,])) %>% drop_units() > .5
    
    if (length(tmp)==0){
      intersection[i] = FALSE
    } else {
      intersection[i] = TRUE
    }
    
  }
  
  # 3 areas filtered out
  table(intersection)
  areas = areas[intersection, ]
  
  # assign unique id to each area
  areas$unique_id = as.numeric(as.factor(paste(areas$camera_id, areas$day, areas$nomination_id)))
  
  # save cleaned dataset for later use
  saveRDS(areas, file = "data/ecological_knowledge_polygons_cleaned.rds")
  
  # number of iterations (how many alternative areas for each nominated area)
  n_iter = 50
  
  # preallocate results
  translated_clusters = list()
  c=1
  # loop through all areas
  # takes roughly 2hrs.
  for (j in 1:nrow(areas)){
    
    # subset area
    sub_clusters = areas[j,]
    
    # select correct competition area
    day_tmp = sub_clusters$day[1]
    unique_id_tmp = sub_clusters$unique_id[1]
    competition_area_tmp = competition_areas[[as.numeric(day_tmp)]]
    
    # confidence level
    confidence = sub_clusters$confidence #1=low, 3 = high confidence
    
    # define maximum vertical and horizontal shifts
    # to limit random simulation to sensible locations (i.e., shift in either
    # direction cannot be larger than arena extent.)
    extent = raster::extent(competition_area_tmp)
    
    # cluster centroid
    centroid = st_centroid(st_geometry(sub_clusters))
    
    # define maximum neg/pos vertical/horizontal shift
    # based on position of cluster centroid. 
    # shift can't be larger than distance of centroid to arena boundaries
    max_delta_x_left = st_coordinates(centroid)[,1] - extent@xmin
    max_delta_x_right =  extent@xmax - st_coordinates(centroid)[,1]
    max_delta_y_up = extent@ymax - st_coordinates(centroid)[,2]
    max_delta_y_down = st_coordinates(centroid)[,2] - extent@ymin
    
    # translate polygons
    for (i in 1:n_iter){
    
      # initialize test: cluster in competition area? (0=no, 1 = yes)
      test = 0
        
      # counter: how often tried to find suitable rotation
      counter = 0
        
      # how many clusters skipped
      skip = 0
        
      while(test == 0){
        
        # sample vertical and horizontal shift for area
        delta_x = runif(1, -max_delta_x_left, max_delta_x_right)
        delta_y = runif(1, -max_delta_y_down, max_delta_y_up)
        
        # sample rotation angle
        rotation_angle = runif(1, 0, 360)
        
        # translate cluster within extent: 1st rotate, 2nd add xy shift
        translated_cluster = st_geometry(rotate.polygon(st_geometry(sub_clusters), anchor = "center", angle = rotation_angle)) + c(delta_x, delta_y)
        
        # add original coordinate system
        st_crs(translated_cluster) <- st_crs(sub_clusters)
        
        # transform competition area, if in different coordinate system
        competition_area_tmp = st_transform(competition_area_tmp, st_crs(translated_cluster))
        
        # check if cluster is in bounding box
        #test = length(st_contains(st_geometry(competition_area_tmp), translated_cluster)[[1]])
        
        # if more than 90% of cluster are in bounding box, accept
        if(st_is_valid(translated_cluster)){
          test = (st_area(st_intersection(st_geometry(competition_area_tmp), translated_cluster)) / st_area(translated_cluster)) %>% drop_units() > .9
          
          if (length(test) == 0){
            test = 0
          }
        } else {
          test = length(st_contains(st_geometry(competition_area_tmp), translated_cluster)[[1]])
        }
        
        #ggplot() + 
        #  geom_sf(data = competition_area_tmp) + 
        #  geom_sf(data = translated_clusters %>% dplyr::bind_rows())
        
        # count up
        counter = counter + 1
        
        # if no solution found after 100 iterations, translated cluster is equal to original cluster
        # (conservative, this cluster will not affect model result bc of 0 variation across features)
        if (counter > 100){
          
          translated_cluster = sub_clusters$x          
          
          # test true
          test = 1
          
          # cluster count
          skip = 1
      
        }
        
       #print(i) 
       
      }
      
      # convert to dataframe
      translated_cluster = st_as_sf(translated_cluster)
      translated_cluster$iter = i
      translated_cluster$unique_id = unique_id_tmp
      translated_cluster$day = day_tmp
      translated_cluster$camera_id = sub_clusters$camera_id
      translated_cluster$successful_translation = 1-skip
      translated_cluster$confidence = confidence
      
      # add to list
      translated_clusters[[c]] <- translated_cluster
      c = c + 1
    
    }
    
    # progress
    print(paste(j, "/", nrow(areas)))
    
  }
  
  # merge to one df
  alternative_clusters = dplyr::bind_rows(translated_clusters)
  alternative_clusters$condition = "simulated"
  
  # add original clusters
  areas = areas %>% dplyr::select(unique_id, day, camera_id, confidence, x)
  areas$condition = "data"
  areas$iter = 0
  areas$successful_translation = 1
  
  # join to simulated
  alternative_clusters = bind_rows(alternative_clusters, areas)
  
  # save results
  saveRDS(alternative_clusters, file = "data/output/simulated_clusters.rds")
  
################################################################################
# END
################################################################################
}
