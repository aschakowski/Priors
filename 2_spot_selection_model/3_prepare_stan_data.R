################################################################################
#
# TITLE: Prepare Stan Data
#
# AUTHORS: Schakowski, A. (schakowski@mpib-berlin.mpg.de)
#
# DESCRIPTION: This script prepares a stan data list for fitting the spot 
# selection models.
#
# LAST UPDATED (DD/MM/YYYY): 10/10/2026
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
                                success_distances, 
                                loss_distances, 
                                locality, 
                                grainsize, 
                                identifiers,
                                confidence_matrix,
                                ncolumns,
                                lake_characteristics,
                                option = "option"){
    
    # subset ids
    sub_ids = identifiers
    
    # filter first spots
    if (option != "first_spots"){
      sub_ids = sub_ids[sub_ids$spot_id > 1,]
    }
    
    # spot ids
    spot_ids = as.numeric(sub_ids$spot_id)
    
    # lake index
    lake = as.numeric(as.factor(sub_ids$day))
    N_lakes = length(unique(lake))
    
    # trip index
    trip = as.numeric(as.factor(sub_ids$track_id))
    N_trips = length(unique(trip))
    
    # participant index
    id = as.numeric(as.factor(sub_ids$participant_id))
    N_ids = length(unique(id))
    
    # make model ids key
    sub_ids$id = id
    sub_ids$trip = trip
    sub_ids$lake = lake
    
    # save model key
    saveRDS(sub_ids, file = "data/processed_data/identifiers_spot_selection_model_key.rds")
    
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
  
    # select closest ncol spots
    for (i in 1:length(social_distances)){
      
      # social
      social_distances[[i]][is.na(social_distances[[i]])] <- 0.0001
      social_distances[[i]] = t(apply(social_distances[[i]], 1, sort))
      social_distances[[i]] <- social_distances[[i]][, 1:ncolumns]
      
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
      
      loss_distances[[i]] <- pad_matrix(as.matrix(loss_distances[[i]]),
                                        new_rows = nrow(as.matrix(loss_distances[[i]])),
                                        new_columns = ncolumns,
                                        filler = -99)
    }
    
    # convert list of matrices to stan conform array
    social_feature <- abind(social_distances, along=3)
    social_feature <- aperm(social_feature, c(3,1,2)) #dim 1: time, dim2: alternative spots, dim3: distances to 55 players
    
    success_feature <- abind(success_distances, along=3)
    success_feature <- aperm(success_feature, c(3,1,2)) #dim 1: time, dim2: alternative spots, dim3: distances to 55 players
    
    loss_feature <- abind(loss_distances, along=3)
    loss_feature <- aperm(loss_feature, c(3,1,2)) #dim 1: time, dim2: alternative spots, dim3: distances to 55 players
    
    # locality matrix: columns: all spot choices, rows: alternative spots
    locality_feature = do.call(cbind, locality)
    
    # convert knowledge to matrix
    ecological_knowledge_feature = list()
    for (i in 1:length(ecological_knowledge_distances)){
      if (is.null(dim(ecological_knowledge_distances[[i]]))){
        ecological_knowledge_feature[[i]] = t(t(ecological_knowledge_distances[[i]]))
      } else {
        ecological_knowledge_feature[[i]] = ecological_knowledge_distances[[i]][,1]
      }
    }
    ecological_knowledge_feature = do.call(cbind, ecological_knowledge_feature)
    
    # number of alternative spots
    N_spots = nrow(social_distances[[1]])
    
    # choice vector
    choices = 1:length(social_distances)
    
    # grainsize for within chain parallelization
    grainsize = grainsize
    
    # confidence
    confidence = confidence_matrix
    
    # normalize time since start within each competition
    time_since_start = as.numeric(sub_ids$time_since_start)
    
    max_times = data.frame(lake = lake, 
               time = time_since_start) %>% 
      ungroup() %>% 
      dplyr::group_by(lake) %>% 
      dplyr::summarize(max_time = max(time),
                       min_time = min(time))
    times = data.frame(lake = lake, 
                       time = time_since_start)
    times = left_join(times, max_times)
    time_since_start_normalized = (times$time-times$min_time) / (times$max_time-times$min_time)
      
    # return padded distance matrix and columns
      return(list(Distances_social = social_feature,
                  Distances_successes = success_feature,
                  Distances_losses = loss_feature,
                  Distances_ecological_knowledge = ecological_knowledge_feature,
                  confidence = confidence,
                  Locality = locality_feature,
                  Time_since_start = time_since_start_normalized,
                  success_columns = success_columns, 
                  loss_columns = loss_columns,
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
                  N_ids = N_ids,
                  spot_id = spot_ids))
      
    
  }
  
  ##############################################################################
  # Run function
  ##############################################################################
  # spot ids
  identifiers = readRDS(file = "data/processed_data/identifiers_spot_selection.rds")
  
  # load raw features
  social_distances = readRDS(file = "data/processed_data/social_distance_matrix.rds")
  success_distances = readRDS(file = "data/processed_data/success_distance_matrix.rds")
  loss_distances = readRDS(file = "data/processed_data/loss_distance_matrix.rds")
  locality = readRDS(file = "data/processed_data/locality_matrix.rds")
  ecological_knowledge_distances = readRDS(file = "data/processed_data/knowledge_distance_matrix.rds")
  confidence_matrix = readRDS(file = "data/processed_data/confidence_matrix.rds")
  
  # filter out first spots
  social_distances = social_distances[identifiers$spot_id > 1]
  success_distances = success_distances[identifiers$spot_id > 1]
  loss_distances = loss_distances[identifiers$spot_id > 1]
  locality = locality[identifiers$spot_id > 1]
  ecological_knowledge_distances = ecological_knowledge_distances[identifiers$spot_id > 1]
  confidence_matrix = confidence_matrix[, identifiers$spot_id > 1]
  confidence_matrix[is.na(confidence_matrix)] <- -99
  
  # compute feature matrix
  data_list = prepare_stan_data(social_distances = social_distances,
                                ecological_knowledge_distances = ecological_knowledge_distances,
                                success_distances = success_distances,
                                loss_distances = loss_distances,
                                locality = locality, 
                                grainsize = 1, 
                                identifiers = identifiers, 
                                confidence_matrix = confidence_matrix,
                                ncolumns = 15)
  
  # any 0s in dataset? If so, replace with small value
  data_list$Distances_social[data_list$Distances_social == 0] <- .01
  data_list$Distances_success[data_list$Distances_success == 0] <- .01
  data_list$Distances_loss[data_list$Distances_loss == 0] <- .01
  data_list$Locality[data_list$Locality == 0] <- .01
  data_list$Distances_ecological_knowledge[data_list$Distances_ecological_knowledge == 0] <- 0.01
  data_list$confidence[is.na(data_list$confidence)]<- -99
  saveRDS(data_list, file = "data/processed_data/stan_data_spot_selection.rds")
  gc()
  
  ##############################################################################
  # prepare data for first spots (will be used later to calculate spatial features
  # for all spots including first)
  ##############################################################################
  # spot ids
  identifiers = readRDS(file = "data/processed_data/identifiers_spot_selection.rds")
  
  # load raw features
  social_distances = readRDS(file = "data/processed_data/social_distance_matrix.rds")
  success_distances = readRDS(file = "data/processed_data/success_distance_matrix.rds")
  loss_distances = readRDS(file = "data/processed_data/loss_distance_matrix.rds")
  locality = readRDS(file = "data/processed_data/locality_matrix.rds")
  ecological_knowledge_distances = readRDS(file = "data/processed_data/knowledge_distance_matrix.rds")
  confidence_matrix = readRDS(file = "data/processed_data/confidence_matrix.rds")
  
  # compute feature matrix
  data_list = prepare_stan_data(social_distances = social_distances,
                                ecological_knowledge_distances = ecological_knowledge_distances,
                                success_distances = success_distances,
                                loss_distances = loss_distances,
                                locality = locality, 
                                grainsize = 1, 
                                identifiers = identifiers, 
                                confidence_matrix = confidence_matrix,
                                ncolumns = 15,
                                option = "first_spots")
  
  # any 0s in dataset? If so, replace with small value
  data_list$Distances_social[data_list$Distances_social == 0] <- .01
  data_list$Distances_success[data_list$Distances_success == 0] <- .01
  data_list$Distances_loss[data_list$Distances_loss == 0] <- .01
  data_list$Locality[data_list$Locality == 0] <- .01
  data_list$Distances_ecological_knowledge[data_list$Distances_ecological_knowledge == 0] <- 0.01
  data_list$confidence[is.na(data_list$confidence)]<- -99
  saveRDS(data_list, file = "data/processed_data/stan_data_spot_selection_first_spots.rds")
  gc()
  
}
################################################################################
# END
################################################################################