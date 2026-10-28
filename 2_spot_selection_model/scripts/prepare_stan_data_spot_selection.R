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
                              roughness_distances, 
                              ecological_knowledge_distances,
                              ecological_knowledge_distances_1,
                              ecological_knowledge_distances_2,
                              ecological_knowledge_distances_3,
                              success_distances, 
                              loss_distances, 
                              spot_distances,
                              locality, 
                              #fished_before,
                              time_since_start,
                              grainsize, 
                              identifiers,
                              confidence_matrix,
                              ncolumns){
  
  # subset ids
  sub_ids = identifiers[identifiers$spot_id != 1, ][1:length(social_distances),]
  #sub_ids = identifiers[1:length(social_distances),]
  spot_ids = sub_ids$spot_id
  
  # lake index
  sub_ids$lake_index = as.numeric(as.factor(sub_ids$lake_id))
  lake = sub_ids$lake_index
  N_lakes = length(unique(sub_ids$lake_index))
  
  # trip index
  sub_ids$trip_index = as.numeric(as.factor(sub_ids$trip_id))
  trip = sub_ids$trip_index
  N_trips = length(unique(sub_ids$trip_index))
  
  # participant index
  sub_ids$id_index = as.numeric(as.factor(sub_ids$participant_id))
  id = sub_ids$id_index
  N_ids = length(unique(sub_ids$id_index))
  
  # compute success columns
  success_columns = vector(mode = "numeric", length = length(success_distances))
  for (i in 1:length(success_distances)){
    
    if (is.null(dim(success_distances[[i]]))){
      success_columns[i] <- 1
    } else {
      success_columns[i] <- ncol(success_distances[[i]])
      success_columns[i] <- min(ncolumns, success_columns[[i]]) 
    }
  }
  
  
  # compute loss columns
  loss_columns = vector(mode = "numeric", length = length(loss_distances))
  for (i in 1:length(loss_distances)){
    
    if (is.null(dim(loss_distances[[i]]))){
      loss_columns[i] <- 1
    } else {
      loss_columns[i] <- ncol(loss_distances[[i]])
      loss_columns[i] <- min(ncolumns, loss_columns[[i]]) 
    }
  }
  
  # compute spot columns
  spot_columns = vector(mode = "numeric", length = length(spot_distances))
  for (i in 1:length(spot_distances)){
    
    if (is.null(dim(spot_distances[[i]]))){
      spot_columns[i] <- 1
    } else {
      spot_columns[i] <- ncol(spot_distances[[i]])
      spot_columns[i] <- min(ncolumns, spot_columns[[i]]) 
    }
  }
  
  # compute knowledge columns
  ecological_knowledge_columns = vector(mode = "numeric", length = length(ecological_knowledge_distances))
  for (i in 1:length(ecological_knowledge_distances)){
    
    if (is.null(dim(ecological_knowledge_distances[[i]]))){
      ecological_knowledge_columns[i] <- 1
    } else {
      ecological_knowledge_columns[i] <- ncol(ecological_knowledge_distances[[i]])
      ecological_knowledge_columns[i] <- min(ncolumns, ecological_knowledge_columns[[i]]) 
    }
  }
  
  # select closest ncol spots
  for (i in 1:length(social_distances)){
    
    # social
    social_distances[[i]][is.na(social_distances[[i]])] <- 0.0001
    social_distances[[i]] = t(apply(social_distances[[i]], 1, sort))
    social_distances[[i]] <- social_distances[[i]][, 1:ncolumns]
    
    # roughness
    roughness_distances[[i]][is.na(roughness_distances[[i]])] <- 0.0001
    roughness_distances[[i]] = t(apply(roughness_distances[[i]], 1, sort))
    roughness_distances[[i]] <- as.matrix(roughness_distances[[i]])[, 1:ncolumns]
    
    
    # select last ncol spots
    if (success_columns[[i]] == 1){
      success_distances[[i]] = success_distances[[i]]
    } else {
      success_distances[[i]] = success_distances[[i]][, 1:success_columns[i]]
    }
    
    if (loss_columns[[i]] == 1){
      loss_distances[[i]] = loss_distances[[i]]
    } else {
      loss_distances[[i]] = loss_distances[[i]][, 1:loss_columns[i]] 
    }
    
    if (spot_columns[[i]] == 1){
      spot_distances[[i]] = spot_distances[[i]]
    } else {
      spot_distances[[i]] = spot_distances[[i]][, 1:spot_columns[i]] 
    }
    
  }
  
  # pad matrix to max columns
  for (i in 1:length(social_distances)){
    
    social_distances[[i]] <- pad_matrix(social_distances[[i]],
                                        new_rows = nrow(social_distances[[i]]),
                                        new_columns = ncolumns,
                                        filler = -99)
    
    success_distances[[i]] <- pad_matrix(as.matrix(success_distances[[i]]),
                                         new_rows = nrow(as.matrix(success_distances[[i]])),
                                         new_columns = ncolumns,
                                         filler = -99)
    
    spot_distances[[i]] <- pad_matrix(as.matrix(spot_distances[[i]]),
                                         new_rows = nrow(as.matrix(spot_distances[[i]])),
                                         new_columns = ncolumns,
                                         filler = -99)
    
    roughness_distances[[i]] <- pad_matrix(roughness_distances[[i]],
                                           new_rows = nrow(roughness_distances[[i]]),
                                           new_columns = ncolumns,
                                           filler = -99)
    
    loss_distances[[i]] <- pad_matrix(as.matrix(loss_distances[[i]]),
                                      new_rows = nrow(as.matrix(loss_distances[[i]])),
                                      new_columns = ncolumns,
                                      filler = -99)
    
    if (is.null(dim(ecological_knowledge_distances[[i]]))){
      ecological_knowledge_distances[[i]] <- pad_matrix(as.matrix(t(t(ecological_knowledge_distances[[i]]))),
                                                        new_rows = nrow(as.matrix(t(t(ecological_knowledge_distances[[i]])))),
                                                        new_columns = ncolumns,
                                                        filler = -99)
      
    } else if (nrow(ecological_knowledge_distances[[i]])<51){
      ecological_knowledge_distances[[i]] <- pad_matrix(as.matrix(t(ecological_knowledge_distances[[i]])),
                                                        new_rows = nrow(as.matrix(t(ecological_knowledge_distances[[i]]))),
                                                        new_columns = ncolumns,
                                                        filler = -99)
    } else {
      ecological_knowledge_distances[[i]] <- pad_matrix(as.matrix(ecological_knowledge_distances[[i]]),
                                                        new_rows = nrow(as.matrix(ecological_knowledge_distances[[i]])),
                                                        new_columns = ncolumns,
                                                        filler = -99)
      
    }
    
    # level 1
    if (is.null(dim(ecological_knowledge_distances_1[[i]]))){
      ecological_knowledge_distances_1[[i]] <- pad_matrix(as.matrix(t(t(ecological_knowledge_distances_1[[i]]))),
                                                        new_rows = nrow(as.matrix(t(t(ecological_knowledge_distances_1[[i]])))),
                                                        new_columns = ncolumns,
                                                        filler = -99)
      
    } else if (nrow(ecological_knowledge_distances_1[[i]])<51){
      ecological_knowledge_distances_1[[i]] <- pad_matrix(as.matrix(t(ecological_knowledge_distances_1[[i]])),
                                                        new_rows = nrow(as.matrix(t(ecological_knowledge_distances_1[[i]]))),
                                                        new_columns = ncolumns,
                                                        filler = -99)
    } else {
      ecological_knowledge_distances_1[[i]] <- pad_matrix(as.matrix(ecological_knowledge_distances_1[[i]]),
                                                        new_rows = nrow(as.matrix(ecological_knowledge_distances_1[[i]])),
                                                        new_columns = ncolumns,
                                                        filler = -99)
      
    }
    
    
    # level 2
    if (is.null(dim(ecological_knowledge_distances_2[[i]]))){
      ecological_knowledge_distances_2[[i]] <- pad_matrix(as.matrix(t(t(ecological_knowledge_distances_2[[i]]))),
                                                          new_rows = nrow(as.matrix(t(t(ecological_knowledge_distances_2[[i]])))),
                                                          new_columns = ncolumns,
                                                          filler = -99)
      
    } else if (nrow(ecological_knowledge_distances_2[[i]])<51){
      ecological_knowledge_distances_2[[i]] <- pad_matrix(as.matrix(t(ecological_knowledge_distances_2[[i]])),
                                                          new_rows = nrow(as.matrix(t(ecological_knowledge_distances_2[[i]]))),
                                                          new_columns = ncolumns,
                                                          filler = -99)
    } else {
      ecological_knowledge_distances_2[[i]] <- pad_matrix(as.matrix(ecological_knowledge_distances_2[[i]]),
                                                          new_rows = nrow(as.matrix(ecological_knowledge_distances_2[[i]])),
                                                          new_columns = ncolumns,
                                                          filler = -99)
      
    }
    
    # level 3
    if (is.null(dim(ecological_knowledge_distances_3[[i]]))){
      ecological_knowledge_distances_3[[i]] <- pad_matrix(as.matrix(t(t(ecological_knowledge_distances_3[[i]]))),
                                                          new_rows = nrow(as.matrix(t(t(ecological_knowledge_distances_3[[i]])))),
                                                          new_columns = ncolumns,
                                                          filler = -99)
      
    } else if (nrow(ecological_knowledge_distances_3[[i]])<51){
      ecological_knowledge_distances_3[[i]] <- pad_matrix(as.matrix(t(ecological_knowledge_distances_3[[i]])),
                                                          new_rows = nrow(as.matrix(t(ecological_knowledge_distances_3[[i]]))),
                                                          new_columns = ncolumns,
                                                          filler = -99)
    } else {
      ecological_knowledge_distances_3[[i]] <- pad_matrix(as.matrix(ecological_knowledge_distances_3[[i]]),
                                                          new_rows = nrow(as.matrix(ecological_knowledge_distances_3[[i]])),
                                                          new_columns = ncolumns,
                                                          filler = -99)
      
    }
}

  # convert list of matrices to stan conform array
  social_feature <- abind(social_distances, along=3)
  social_feature <- aperm(social_feature, c(3,1,2)) #dim 1: time, dim2: alternative spots, dim3: distances to 55 players
  
  success_feature <- abind(success_distances, along=3)
  success_feature <- aperm(success_feature, c(3,1,2)) #dim 1: time, dim2: alternative spots, dim3: distances to 55 players
  
  loss_feature <- abind(loss_distances, along=3)
  loss_feature <- aperm(loss_feature, c(3,1,2)) #dim 1: time, dim2: alternative spots, dim3: distances to 55 players
  
  #spot_feature <- abind(spot_distances, along=3)
  #spot_feature <- aperm(spot_feature, c(3,1,2)) #dim 1: time, dim2: alternative spots, dim3: distances to 55 players
  
  roughness_feature <- abind(roughness_distances, along=3)
  roughness_feature <- aperm(roughness_feature, c(3,1,2)) #dim 1: time, dim2: alternative spots, dim3: distances to 55 players
  
  ecological_knowledge_feature <- abind(ecological_knowledge_distances, along=3)
  ecological_knowledge_feature <- aperm(ecological_knowledge_feature, c(3, 1, 2))
  
  ecological_knowledge_feature_1 <- abind(ecological_knowledge_distances_1, along=3)
  ecological_knowledge_feature_1 <- aperm(ecological_knowledge_feature_1, c(3, 1, 2))
  
  ecological_knowledge_feature_2 <- abind(ecological_knowledge_distances_2, along=3)
  ecological_knowledge_feature_2 <- aperm(ecological_knowledge_feature_2, c(3, 1, 2))
  
  ecological_knowledge_feature_3 <- abind(ecological_knowledge_distances_3, along=3)
  ecological_knowledge_feature_3 <- aperm(ecological_knowledge_feature_3, c(3, 1, 2))
  
  # locality matrix: columns: all spot choices, rows: alternative spots
  locality_feature = do.call(cbind, locality)
  
  # number of alternative spots
  N_spots = nrow(social_distances[[1]])
  
  # choice vector
  choices = 1:length(social_distances)
  
  # grainsize for within chain parallelization
  grainsize = grainsize
  
  # confidence
  confidence = confidence_matrix
  
  # return padded distance matrix and columns
  return(list(Distances_social = social_feature,
              Distances_roughness = roughness_feature, 
              Distances_successes = success_feature,
              Distances_losses = loss_feature,
              #Distances_spots = spot_feature,
              Distances_ecological_knowledge = ecological_knowledge_feature,
              Distances_ecological_knowledge_1 = ecological_knowledge_feature_1,
              Distances_ecological_knowledge_2 = ecological_knowledge_feature_2,
              Distances_ecological_knowledge_3 = ecological_knowledge_feature_3,
              confidence = confidence,
              Locality = locality_feature,
              Time_since_start = as.numeric(time_since_start)/max(as.numeric(time_since_start)),
              success_columns = success_columns, 
              loss_columns = loss_columns,
              spot_columns = spot_columns,
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
              #prior_exposure = fished_before,
              N_ids = N_ids,
              spot_ids = spot_ids))
  
}
################################################################################
# spot ids
identifiers = readRDS(file = "data/processed_data/identifiers_spot_selection.rds")

# catch data
catch_data = read_delim("data/catch_data.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
catch_data = catch_data %>% filter(year == 2023 & competition == "individual") %>% dplyr::select(camera_id, day, year, participant_id, fished_before)
catch_data$fished_before = as.numeric(as.factor(catch_data$fished_before))
catch_data$day = as.character(catch_data$day)
catch_data$year = as.character(catch_data$year)
identifiers = left_join(identifiers, catch_data)

# load raw features
social_distances = readRDS(file = "data/processed_data/social_distance_matrix.rds")
roughness_distances = readRDS(file = "data/processed_data/roughness_distance_matrix.rds")
success_distances = readRDS(file = "data/processed_data/success_distance_matrix.rds")
loss_distances = readRDS(file = "data/processed_data/loss_distance_matrix.rds")
spot_distances = readRDS(file = "data/processed_data/spot_distance_matrix.rds")
locality = readRDS(file = "data/processed_data/locality_matrix.rds")
ecological_knowledge_distances = readRDS(file = "data/processed_data/knowledge_distance_matrix.rds")
ecological_knowledge_distances_1 = readRDS(file = "data/processed_data/knowledge_distance_matrix_1.rds")
ecological_knowledge_distances_2 = readRDS(file = "data/processed_data/knowledge_distance_matrix_2.rds")
ecological_knowledge_distances_3 = readRDS(file = "data/processed_data/knowledge_distance_matrix_3.rds")
confidence_matrix = readRDS(file = "data/processed_data/confidence_matrix.rds")[,identifiers$spot_id != 1]

# time since start as feature
time_since_start = identifiers$time_since_start[(identifiers$spot_id!= 1)]

# fished before
fished_before = identifiers$fished_before[(identifiers$spot_id!= 1)]

# drop first choices
social_distances = social_distances[(identifiers$spot_id!= 1)]
roughness_distances = roughness_distances[(identifiers$spot_id!= 1)]
success_distances = success_distances[(identifiers$spot_id!= 1)]
loss_distances = loss_distances[(identifiers$spot_id!= 1)]
locality = locality[(identifiers$spot_id!= 1)]
ecological_knowledge_distances=ecological_knowledge_distances[identifiers$spot_id!=1]
ecological_knowledge_distances_1=ecological_knowledge_distances_1[identifiers$spot_id!=1]
ecological_knowledge_distances_2=ecological_knowledge_distances_2[identifiers$spot_id!=1]
ecological_knowledge_distances_3=ecological_knowledge_distances_3[identifiers$spot_id!=1]

# compute feature matrix
data_list = prepare_stan_data(social_distances = social_distances,
                              roughness_distances = roughness_distances,
                              ecological_knowledge_distances = ecological_knowledge_distances,
                              ecological_knowledge_distances_1 = ecological_knowledge_distances_1,
                              ecological_knowledge_distances_2 = ecological_knowledge_distances_2,
                              ecological_knowledge_distances_3 = ecological_knowledge_distances_3,
                              success_distances = success_distances,
                              loss_distances = loss_distances,
                              spot_distances = spot_distances,
                              locality = locality, 
                              time_since_start = time_since_start,
                              #fished_before = fished_before,
                              grainsize = 1, identifiers = identifiers, confidence_matrix,ncolumns = 15)
saveRDS(data_list, file = "data/processed_data/stan_data_spot_selection.rds")
gc()

################################################################################
# spot ids
identifiers = readRDS(file = "data/processed_data/identifiers_spot_selection.rds")

# catch data
catch_data = read_delim("data/catch_data.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
catch_data = catch_data %>% filter(year == 2023 & competition == "individual") %>% dplyr::select(camera_id, day, year, participant_id, fished_before)
catch_data$fished_before = as.numeric(as.factor(catch_data$fished_before))
catch_data$day = as.character(catch_data$day)
catch_data$year = as.character(catch_data$year)
identifiers = left_join(identifiers, catch_data)

# load raw features
social_distances = readRDS(file = "data/processed_data/social_distance_matrix.rds")
roughness_distances = readRDS(file = "data/processed_data/roughness_distance_matrix.rds")
success_distances = readRDS(file = "data/processed_data/success_distance_matrix.rds")
loss_distances = readRDS(file = "data/processed_data/loss_distance_matrix.rds")
spot_distances = readRDS(file = "data/processed_data/spot_distance_matrix.rds")
locality = readRDS(file = "data/processed_data/locality_matrix.rds")
ecological_knowledge_distances = readRDS(file = "data/processed_data/knowledge_distance_matrix.rds")
ecological_knowledge_distances_1 = readRDS(file = "data/processed_data/knowledge_distance_matrix_1.rds")
ecological_knowledge_distances_2 = readRDS(file = "data/processed_data/knowledge_distance_matrix_2.rds")
ecological_knowledge_distances_3 = readRDS(file = "data/processed_data/knowledge_distance_matrix_3.rds")

# time since start as feature
time_since_start = identifiers$time_since_start[(identifiers$spot_id!= 1 & !is.na(identifiers$fished_before))]

# fished before
fished_before = identifiers$fished_before[(identifiers$spot_id!= 1 & !is.na(identifiers$fished_before))]

# drop first choices
social_distances = social_distances[(identifiers$spot_id!= 1 & !is.na(identifiers$fished_before))]
roughness_distances = roughness_distances[(identifiers$spot_id!= 1 & !is.na(identifiers$fished_before))]
success_distances = success_distances[(identifiers$spot_id!= 1 & !is.na(identifiers$fished_before))]
loss_distances = loss_distances[(identifiers$spot_id!= 1 & !is.na(identifiers$fished_before))]
locality = locality[(identifiers$spot_id!= 1 & !is.na(identifiers$fished_before))]
ecological_knowledge_distances=ecological_knowledge_distances[identifiers$spot_id!=1 & !is.na(identifiers$fished_before)]
ecological_knowledge_distances_1=ecological_knowledge_distances_1[identifiers$spot_id!=1 & !is.na(identifiers$fished_before)]
ecological_knowledge_distances_2=ecological_knowledge_distances_2[identifiers$spot_id!=1 & !is.na(identifiers$fished_before)]
ecological_knowledge_distances_3=ecological_knowledge_distances_3[identifiers$spot_id!=1 & !is.na(identifiers$fished_before)]

# compute feature matrix
data_list = prepare_stan_data(social_distances = social_distances,
                              roughness_distances = roughness_distances,
                              ecological_knowledge_distances = ecological_knowledge_distances,
                              ecological_knowledge_distances_1 = ecological_knowledge_distances_1,
                              ecological_knowledge_distances_2 = ecological_knowledge_distances_2,
                              ecological_knowledge_distances_3 = ecological_knowledge_distances_3,
                              success_distances = success_distances,
                              loss_distances = loss_distances,
                              spot_distances = spot_distances,
                              locality = locality, 
                              time_since_start = time_since_start,
                              fished_before = fished_before,
                              grainsize = 1, identifiers = identifiers, ncolumns = 15)
saveRDS(data_list, file = "data/processed_data/stan_data_spot_selection_prior_exposure.rds")
################################################################################


}
data_list$Locality[,1]

################################################################################
# END
################################################################################