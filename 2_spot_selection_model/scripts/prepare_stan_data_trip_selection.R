################################################################################
#
# Title: Prepare Stan Data
#
# Description: Prepare stan data for model fit
#
# Authors: Schakowski, A.
#
# Last updated: 13/02/2024 (DD/MM/YYYY)
#
################################################################################
prepare_stan_data_spot_selection <- function(){
  
  ################################################################################
  # Helper Functions
  ################################################################################
  pad_matrix <- function(matrix, new_rows, new_columns, filler){
    
    if ( new_rows > nrow(matrix) & new_columns > ncol(matrix)){
      
      # add columns
      column_diff = new_columns - ncol(matrix)
      columns = replicate(n = column_diff, expr = rep(filler, nrow(matrix)))
      matrix = cbind(matrix, columns)
      
      # add rows
      row_diff = new_rows - nrow(matrix)
      rows = t(replicate(n = row_diff, expr = rep(filler, ncol(matrix))))
      matrix = rbind(matrix, rows)
      
    } else if (new_rows > nrow(matrix) & new_columns == ncol(matrix)){
      
      # add rows
      row_diff = new_rows - nrow(matrix)
      rows = t(replicate(n = row_diff, expr = rep(filler, ncol(matrix))))
      matrix = rbind(matrix, rows)
      
    } else if (new_rows == nrow(matrix) & new_columns > ncol(matrix) & new_rows > 1){
      
      # add columns
      column_diff = new_columns - ncol(matrix)
      columns = replicate(n = column_diff, expr = rep(filler, nrow(matrix)))
      matrix = cbind(matrix, columns)
      
    } else if (new_rows == 1){
      
      ## add columns
      column_diff = new_columns - ncol(matrix)
      columns = replicate(n = column_diff, expr = rep(filler, 1))
      matrix = c(matrix, columns)
      
    } else {
      # do nothing
      #print("New dimensions are old dimensions.")
    }
    
    return(matrix)
  }
  
  prepare_stan_data <- function(social_distances,
                                ecological_knowledge_distances,
                                grainsize, 
                                identifiers, 
                                ncolumns){
    
    # subset ids
    sub_ids = identifiers[1:length(social_distances),]
    
    # lake index
    sub_ids$lake_index = as.numeric(as.factor(sub_ids$day))
    lake = sub_ids$day
    N_lakes = length(unique(sub_ids$day))
    
    # trip index
    sub_ids$trip_index = sub_ids$unique_trip_id
    trip = sub_ids$unique_trip_id
    N_trips = length(unique(sub_ids$unique_trip_id))
    
    # participant index
    sub_ids$id_index = as.numeric(as.factor(sub_ids$participant_id))
    id = sub_ids$id_index
    N_ids = length(unique(sub_ids$id_index))
    
    # compute knowledge columns
    ecological_knowledge_columns = vector(mode = "numeric", length = length(ecological_knowledge_distances))
    for (i in 1:length(ecological_knowledge_distances)){
      
      if (is.null(dim(ecological_knowledge_distances[[i]]))){
        ecological_knowledge_columns[i] <- 1
      } else {
        ecological_knowledge_columns[i] <- ncol(ecological_knowledge_distances[[i]])
        ecological_knowledge_columns[i] <- min(ncolumns, ecological_knowledge_columns[i]) 
      }
    }
    
    # select closest ncol spots
    for (i in 1:length(social_distances)){
      
      # social
      social_distances[[i]][is.na(social_distances[[i]])] <- 0.0001
      social_distances[[i]] = t(apply(social_distances[[i]], 1, sort))
      social_distances[[i]] <- social_distances[[i]][, 1:ncolumns]
      
    }
    
    # pad matrix to max columns
    for (i in 1:length(social_distances)){
      
      social_distances[[i]] <- pad_matrix(social_distances[[i]],
                                          new_rows = nrow(social_distances[[i]]),
                                          new_columns = ncolumns,
                                          filler = -99)
      
      ecological_knowledge_distances[[i]] <- pad_matrix(as.matrix(ecological_knowledge_distances[[i]]),
                                                        new_rows = nrow(as.matrix(ecological_knowledge_distances[[i]])),
                                                        new_columns = ncolumns,
                                                        filler = -99)
      
    }
    
    # convert list of matrices to stan conform array
    social_feature <- abind(social_distances, along=3)
    social_feature <- aperm(social_feature, c(3,1,2)) #dim 1: time, dim2: alternative spots, dim3: distances to 55 players
    
    ecological_knowledge_feature <- abind(ecological_knowledge_distances, along=3)
    ecological_knowledge_feature <- aperm(ecological_knowledge_feature, c(3, 1, 2))
    
    # number of alternative spots
    N_spots = nrow(social_distances[[1]])
    
    # choice vector
    choices = 1:length(social_distances)
    
    # grainsize for within chain parallelization
    grainsize = grainsize
    
    # return padded distance matrix and columns
    return(list(Distances_social = social_feature,
                Distances_ecological_knowledge = ecological_knowledge_feature,
                ecological_knowledge_columns = ecological_knowledge_columns,
                T = length(social_distances), 
                N_spots = N_spots, 
                choices = choices, 
                grainsize = grainsize,
                columns = ncolumns,
                lake = lake,
                N_lakes = N_lakes,
                trip = trip,
                N_trips = N_trips, 
                id = id,
                N_ids = N_ids))
    
  }
  ################################################################################
  # load raw features
  social_distances = readRDS(file = "data/processed_data/social_distance_matrix_trip_selection.rds")
  ecological_knowledge_distances = readRDS(file = "data/processed_data/knowledge_distance_matrix_trip_selection.rds")
  
  # load identifiers
  identifiers = readRDS(file = "data/processed_data/identifiers_trip_selection.rds")
  
  # drop first choices
  #social_distances = social_distances[(identifiers$spot_id!= 1)]
  #roughness_distances = roughness_distances[(identifiers$spot_id!= 1)]
  #success_distances = success_distances[(identifiers$spot_id!= 1)]
  #loss_distances = loss_distances[(identifiers$spot_id!= 1)]
  #locality = locality[(identifiers$spot_id!= 1)]
  #time_since_start = time_since_start[identifiers$spot_id!=1]
  #ecological_knowledge_distances=ecological_knowledge_distances[identifiers$spot_id!=1]
  
  # compute feature matrix
  data_list = prepare_stan_data(social_distances = social_distances,
                                ecological_knowledge_distances = ecological_knowledge_distances,
                                grainsize = 1, 
                                identifiers,
                                ncolumns = 15)
  saveRDS(data_list, file = "data/processed_data/stan_data_trip_selection.rds")
  gc()
}


################################################################################
# END
################################################################################
ecological_knowledge_distances = readRDS(file = "data/processed_data/knowledge_distance_matrix_trip_selection.rds")
social_distances = readRDS(file = "data/processed_data/social_distance_matrix_trip_selection.rds")

dist = matrix(ncol = 7, nrow = length(ecological_knowledge_distances))
for (i in 1:length(ecological_knowledge_distances)){
  if (is.null(dim(ecological_knowledge_distances[[i]]))){
    dist[i,1] <- min(ecological_knowledge_distances[[i]][1])
    dist[i,2] <- mean(ecological_knowledge_distances[[i]][2:51])
    if (dist[i,1] == -99){
      dist[i,1] = dist[i,2] = NA
    }
  } else if (dim(ecological_knowledge_distances[[i]])[2]==1){
    dist[i,1] <- min(ecological_knowledge_distances[[i]][1,])
    dist[i,2] <- mean(ecological_knowledge_distances[[i]][2:51])
  } else {
    dist[i,1] <- min(ecological_knowledge_distances[[i]][1,])
    dist[i,2] <- mean(apply(ecological_knowledge_distances[[i]][2:51,], 1, min)) 
  }
  dist[i,3] = identifiers$spot_id[i]
  dist[i,4] = identifiers$day[i]
  dist[i,5] = identifiers$unique_trip_id[i]
  dist[i,6] = identifiers$camera_id[i]
  dist[i,7] = identifiers$participant_id[i]
}

dist = as.data.frame(dist)
colnames(dist) = c("data", "sim", "spot_id", "day", "trip_id", "camera_id", "participant_id")
dist %>% 
  group_by(trip_id, spot_id, camera_id, participant_id, day) %>% 
  summarize(data = mean(as.numeric(data)),
            sim = mean(as.numeric(sim))) %>% 
  mutate(diff = data - sim) %>% 
  group_by(trip_id, camera_id, participant_id, day) %>% 
  summarize(diff = mean(diff)) %>% 
  filter(abs(diff) < 10)

dist %>% 
  group_by(trip_id, spot_id) %>% 
  summarize(data = mean(as.numeric(data)),
            sim = mean(as.numeric(sim))) %>% 
  mutate(data = exp(-data^2 / (2 * exp(3.2)^2)),
         sim = exp(-sim^2 / 2*exp(3.2)^2)) %>% 
  pivot_longer(3:4) %>% 
  ggplot(aes(x = value, fill = name)) + 
  geom_histogram(position = position_identity(), alpha = .3) + 
  theme_minimal() 
mean(as.numeric(dist$data))
dist %>% 
  group_by(trip_id, spot_id, day) %>% 
  summarize(data = mean(exp(-as.numeric(data)^2 / (2*exp(3.2)^2)), na.rm =T),
            sim = mean(exp(-as.numeric(sim)^2 / (2*exp(3.2)^2)), na.rm = T)) %>% 
  pivot_longer(4:5) %>% 
  ggplot(aes(x = value, fill = name)) + 
  geom_histogram(position = position_identity(), alpha = .3) + 
  theme_minimal()  +
  facet_wrap(~day)

dist %>% 
  group_by(trip_id, spot_id, day) %>% 
  summarize(data = mean(as.numeric(data)),
            sim = mean(as.numeric(sim))) %>% 
  pivot_longer(4:5)

dist %>% 
  group_by(trip_id, day) %>% 
  mutate(spot_id = as.numeric(spot_id) / max(as.numeric(spot_id))) %>%
  group_by(trip_id, spot_id, day) %>% 
  summarize(data = mean(as.numeric(data)),
            sim = mean(as.numeric(sim))) %>% 
  pivot_longer(4:5) %>% 
  ggplot(aes(x = spot_id, y = value, color = name)) + 
  stat_summary_bin(
    fun = mean)+
  facet_wrap(~day)
  

dist = matrix(ncol = 5, nrow = length(social_distances))
for (i in 1:length(social_distances)){
  if (is.null(dim(social_distances[[i]]))){
    dist[i,1] <- social_distances[[i]][1]
    dist[i,2] <- mean(social_distances[[i]][2:51])
  } else {
    dist[i,1] <- social_distances[[i]][1,1]
    dist[i,2] <- median(social_distances[[i]][2:51,1]) 
  }
  dist[i,3] = identifiers$spot_id[i]
  dist[i,4] = identifiers$day[i]
  dist[i,5] = identifiers$unique_trip_id[i]
}

dist = as.data.frame(dist)
colnames(dist) = c("data", "sim", "spot_id", "day", "trip_id")

dist %>% 
  group_by(trip_id, spot_id) %>% 
  summarize(data = mean(as.numeric(data)),
            sim = mean(as.numeric(sim))) %>% 
  #mutate(data = exp(-data^2 / (2 * exp(3.2)^2)),
  #       sim = exp(-sim^2 / 2*exp(3.2)^2)) %>% 
  pivot_longer(3:4) %>% 
  ggplot(aes(x = value, fill = name)) + 
  geom_histogram(position = position_identity(), alpha = .3) + 
  theme_minimal() 



#  test
test = do.call(rbind, Feature_matrix)
plot_dat = test %>% 
  #filter(chosen == 1) %>% 
  filter(day == 1 & camera_id == 27) %>%
  arrange(spot_id)

ggplot() + 
  geom_path(data = plot_dat, aes(x = x, y = y, group = iteration, color = as.factor(chosen))) + 
  geom_sf(data = test2)

test2 = ecological_knowledge_polygons %>% 
  filter(day == 1 & camera_id == 27)
