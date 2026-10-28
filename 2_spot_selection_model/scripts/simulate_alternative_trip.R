################################################################################
#
# Title: Simulate alternative foraging trips
#
# Authors: Schakowski, A.
#
# Last updated: 27/07/2024 (DD/MM/YYYY)
#
################################################################################
# helper rotate and shift foraging trip.
shift.rotate.trip <- function(trip, angle, shift.x, shift.y){
  
  # shift by fixed offset
  shifted_trip = trip
  
  # standardize at origin
  shifted_trip[, 1] = shifted_trip[,1]-shifted_trip[1,1]
  shifted_trip[, 2] = shifted_trip[,2]-shifted_trip[1,2]
  
  # now rotate
  rotated_trip = shifted_trip
  rotated_trip[, 1] = shifted_trip[,1]*cos(angle) - shifted_trip[,2]*sin(angle)
  rotated_trip[, 2] = shifted_trip[,1]*sin(angle) + shifted_trip[,2]*cos(angle)
  
  # backtranslate
  rotated_trip[,1] = rotated_trip[,1] + trip[1,1]
  rotated_trip[,2] = rotated_trip[,2] + trip[1,2]
  
  # shift
  rotated_trip[, 1] = rotated_trip[, 1] + shift.x
  rotated_trip[, 2] = rotated_trip[, 2] + shift.y
  
  # return shifted trip
  return(rotated_trip)
  
}

# source helper functions
sourceCpp("functions/fast_dist.cpp")
sourceCpp("functions/ptinpoly.cpp")

# lake key
lake_key = data.frame(day = c(1:5, 1:3, 5,4), 
                      year = rep(c(2023, 2022), each = 5),
                      lake_index = 1:10)

# add to data
spot_selection_data = fread("data/spot_selection_data_cleaned.csv")
spot_selection_data$chosen = 1
spot_selection_data$total_spot_id = 1:nrow(spot_selection_data)

# check if any coordinates are outside competition area
check = vector()
c = 1
competition_area_new=list()
for (i in 1:5){
  
  area = point_lakes[[i]]
  
  coords = spot_selection_data %>% filter(day == lake_key[i, "day"] & year == lake_key[i, "year"]) %>% dplyr::select(x, y)
  
  if (nrow(coords) == 0){
    
    competition_area_new[[i]] = area
    
  } else {
    # convex hull of angling spots
    hull = coords[coords[, 4:5] %>% chull(),c("x", "y")]
    colnames(hull) <- c("x", "y")
    hull2 = rbind(area, hull)
    hull3 = hull2[hull2 %>% chull(),c("x", "y")]
    
    for (j in 1:nrow(coords)){
      check[c] = pnpoly(cbind(coords$x[j], coords$y[j]), as.matrix(hull3)) 
      c = c+1
    }
    
    competition_area_new[[i]] <- hull3
    
  }
  
}

# compute starting point of each competition
starting_point = gps_data %>% 
  group_by(day, year, camera_id) %>% 
  filter(competition == "individual") %>% 
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

# create unique trip id
spot_selection_data$unique_trip_id = as.numeric(as.factor(paste(spot_selection_data$day, spot_selection_data$ID)))

# how many spots
N_trips = length(unique(spot_selection_data$unique_trip_id))

# preallocate features
Feature_matrix = list()

# number of iterations
N_iter = 50

# start simulation
# takes ~ 5min
# start sim
for (t in 1:N_trips){
  
  # get spot data
  trip = spot_selection_data %>% filter(unique_trip_id == t)
  trip_coords = trip[, c("x", "y")]
  
  # lake surface as polygon
  lake = competition_area[[trip$lake_index[1]]]
  
  # lake surface as points
  lake_points = competition_area_new[[trip$lake_index[1]]]
  
  # define max xy shift
  shift.x.max = extent(lake)[2]-extent(lake)[1]
  shift.y.max = extent(lake)[4]-extent(lake)[3]
  
  # simulate coordinates of alternative trip
  alt_trips = list()
  for (s in 1:N_iter){
    
    #plot(trip_coords, xlim = c(extent(lake)[1],extent(lake)[2]), ylim = c(extent(lake)[3],extent(lake)[4]))
    #lines(trip_coords)
    #lines(alt_trip)
    #points(lake_points)
    
    # check if trip on lake
    counter = 1
    check = FALSE
    while(!check & counter < 500){
      
      # sample angles and shifts
      angle = runif(1, -pi, pi)
      shift.x = runif(1, min = -shift.x.max, max = shift.x.max)
      shift.y = runif(1, min = -shift.y.max, max = shift.y.max)
      
      # shift rotate alternative trip
      alt_trip = shift.rotate.trip(as.matrix(trip_coords), angle = angle, shift.x = shift.x, shift.y = shift.y)
      
      # check if trip on lake
      sum = 0
      for (pt in 1:nrow(alt_trip)){
        sum = sum + !pnpoly(cbind(alt_trip[,1], alt_trip[,2])[pt,], as.matrix(lake_points))
      }
      
      # if all points on lake, then sum = 0
      check = ifelse(sum == 0, TRUE, FALSE)
      
      counter = counter + 1
      
    }
    
    if (counter >= 500){
      # sample randomly
      alt_trip = st_coordinates(st_sample(lake, size = nrow(alt_trip)))
      colnames(alt_trip) <- c("x", "y")
    }
    
    # convert to df
    alt_trip = as.data.frame(alt_trip)
    
    # add number of iteration
    alt_trip$iteration = s
    
    # add spot id
    alt_trip$spot_id = 1:nrow(alt_trip)
    
    # add unique trip id
    alt_trip$unique_trip_id = t
    
    # add time start
    alt_trip$timestamp_start = trip$timestamp_start
    
    # add day
    alt_trip$day = trip$day
    
    # add camera id
    alt_trip$camera_id = trip$ID
    
    # add Return
    alt_trip$Return = trip$Return
    
    # add chosen id
    alt_trip$chosen = 0
    
    # save trip
    alt_trips[[s]] = alt_trip
    
  }
  
  # convert to single df
  alt_trips = do.call(rbind, alt_trips)
  
  # add chosen trip
  chosen_trip = trip %>% ungroup() %>% dplyr::select(c(x, y, spot_id, unique_trip_id, timestamp_start, day, Return, chosen, camera_id = ID))
  
  chosen_trip$chosen = 1
  chosen_trip$iteration = 0
  
  # bind rows
  sim_trips = bind_rows(chosen_trip, alt_trips)
  
  # add to feature matrix
  Feature_matrix[[t]] <- sim_trips
  
  # progress
  print(paste(t,"/", N_trips, sep = ""))
  
}

# create identifiers
Feature_matrix = do.call(rbind, Feature_matrix)
Feature_matrix$unique_spot_id = as.numeric(as.factor(paste(Feature_matrix$unique_trip_id, Feature_matrix$spot_id)))
Feature_matrix = split(Feature_matrix, Feature_matrix$unique_spot_id)

# create identifiers
# create identifier list
identifiers = data.frame(day = vector(),
                         spot_id = vector(),
                         camera_id = vector(),
                         t = vector(),
                         unique_trip_id = vector(),
                         unique_spot_id = vector())
for (t in 1:length(Feature_matrix)){
  day = Feature_matrix[[t]]$day[1]
  spot_id = Feature_matrix[[t]]$spot_id[1]
  camera_id = Feature_matrix[[t]]$camera_id[1]
  unique_trip_id = Feature_matrix[[t]]$unique_trip_id[1]
  unique_spot_id = Feature_matrix[[t]]$unique_spot_id[1]
  identifiers[t,] <- c(day, spot_id, camera_id, t, unique_trip_id, unique_spot_id)
}

# load catch data
catch_data <- read_delim("data/catch_data.csv", 
                         delim = ";", escape_double = FALSE, trim_ws = TRUE)
catch_data = catch_data %>% filter(year == 2023) %>% filter(competition == "individual")

## compute unique participant id
identifiers$camera_id = as.character(identifiers$camera_id)
catch_data$day = as.character(catch_data$day)
identifiers = left_join(identifiers, catch_data %>% dplyr::select(camera_id, day, year, participant_id))

saveRDS(identifiers, "data/processed_data/identifiers_trip_selection.rds")

# save results
saveRDS(Feature_matrix, file = "data/processed_data/feature_matrix_trip_selection.rds")

################################################################################
Feature_matrix = readRDS("data/processed_data/feature_matrix_trip_selection.rds")
Feature_matrix = do.call(rbind, Feature_matrix)

Feature_matrix %>% 
  filter(day == 5) %>% 
  ggplot(aes(x = x, y = y, color = chosen, group = iteration)) + 
  geom_path() + 
  facet_wrap(~chosen)
