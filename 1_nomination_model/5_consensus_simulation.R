################################################################################
#
# TITLE: Estimate social consensus of area nominations
# 
# AUTHORS: Schakowski, A. (schakowski@mpib-berlin.mpg.de)
#
# DESCRIPTION: This script simulates and estimates social consensus of area
# nominations as predicted by the respective models.
#
# LAST UPDATED (DD/MM/YYYY): 06/10/2026
#
################################################################################

estimate_consensus <- function(){
  
  # softmax function helper
  softmax <- function(par){
    n.par <- length(par)
    par1 <- sort(par, decreasing = TRUE)
    Lk <- par1[1]
    for (k in 1:(n.par-1)) {
      Lk <- max(par1[k+1], Lk) + log1p(exp(-abs(par1[k+1] - Lk))) 
    }
    val <- exp(par - Lk)
    return(val)
  }
  
  ################################################################################
  # selection of areas according to model parameters
  ################################################################################
  # loop through all areas
  select_areas <- function(iter, clusters, model, condition, model_id){
    
    # extract model parameters
    if (condition == "random"){
      # random choice
    } else {
      # choice according to model parameters
      main_effects = model %>% spread_draws(weights[feature]) %>% group_by(feature) %>% summarize(weights = mean(weights))
      lake_effects = model %>% spread_draws(v_lake[lake, feature]) %>% group_by(lake,feature) %>% summarize(weights = mean(v_lake))
      id_effects = model %>% spread_draws(v_id[id, feature]) %>% group_by(id, feature) %>% summarize(weights = mean(v_id))
      trip_effects = model %>% spread_draws(v_trip[trip, feature]) %>% group_by(trip, feature) %>% summarize(weights = mean(v_trip))
    }
    
    # loop through all clusters
    polygons = list()
    c = 1
    for (i in unique(clusters$unique_id)){
      
      # subset simulated clusters
      sub_clusters = clusters %>% filter(unique_id == i) %>% filter(condition == "simulated")
      
      # identifiers
      trip_tmp = sub_clusters$trip[1]
      lake_tmp = sub_clusters$lake[1]
      id_tmp = sub_clusters$id[1]
      
      # feature matrix
      f = sub_clusters[, c("displacement_scaled", "boundary_scaled", "depth_scaled", "shore_scaled", "roughness_scaled")] %>% st_drop_geometry()
      
      # compute feature weights
      if (condition == "random"){
        w = t(t(c(0, 0, 0, 0, 0)))
      } else {
        w = main_effects$weights + lake_effects[lake_effects$lake == lake_tmp,"weights"] + trip_effects[trip_effects$trip == trip_tmp,"weights"]  + id_effects[id_effects$id == id_tmp,"weights"]
      }
      
      # select polygon
      p = softmax(as.matrix(f)[, 1:nrow(w)] %*% as.matrix(w)) 
      
      # sample
      choice = sample(1:length(p), 1, prob = p)
      
      # select polygon
      polygons[[c]] = sub_clusters[choice,]
      
      # count
      c = c + 1
      
    }
    
    # to df
    polygons = dplyr::bind_rows(polygons)
    
    # compute consensus
    # preallocate results
    intersections_model = vector(mode = "numeric", length = max(clusters$unique_id))
    c = 1
    for (j in unique(polygons$unique_id)){
      
      # subset polygon
      target = polygons %>% filter(unique_id == j)
      target = st_as_sf(target)
      
      # subset reference set
      reference = polygons %>% filter(day == target$day & unique_id != target$unique_id)
      reference = st_as_sf(reference %>% dplyr::select(x))
      
      # sample points from polygon surface
      area = st_area(target)
      area = units::set_units(area, m^2)
      pts = st_sample(target, size = round(max(1, drop_units(area / 100))), type = "regular", replace = F)
      
      # how many intersections with other polygons for each point
      test = st_contains(reference, pts, sparse = T)
      
      # relative number of intersections
      intersections_model[c] = length(unlist(test)) / (area %>% drop_units() / 100)
      c = c + 1
      
    }
    
    # add identifiers to consensus estimates
    tidy_df = data.frame(unique_id = unique(polygons$unique_id),
                         consensus = intersections_model)
    
    identifiers = polygons %>% 
      group_by(unique_id) %>% 
      slice(1) %>% 
      dplyr::select(unique_id, 
                    day, 
                    camera_id,
                    confidence)
    
    tidy_df = left_join(tidy_df, identifiers)
    tidy_df$iteration = iter
    tidy_df$model_id = model_id
    
    return(tidy_df)
    
  }
  
  ################################################################################
  # load clusters
  clusters = readRDS(file = "data/output/model_data.rds")
  
  # loop through fitted models
  save_locations = c("data/processed_data/model_fit/nomination_models/a_random_fit.rds",
                     "data/processed_data/model_fit/nomination_models/b_displacement_fit.rds",
                     "data/processed_data/model_fit/nomination_models/c_boundary_fit.rds",
                     "data/processed_data/model_fit/nomination_models/d_depth_fit.rds",
                     "data/processed_data/model_fit/nomination_models/e_shore_fit.rds",
                     "data/processed_data/model_fit/nomination_models/f_roughness_fit.rds")
  
  consensus_list = list()
  options(future.globals.maxSize = 8000*1024^2)
  for (i in 1:length(save_locations)){
    
    # load model
    model = readRDS(save_locations[i])
    model_id = c("random", "+displacement", "+boundary", "+depth", "+shore", "+roughness")[i]
    
    if (i == 1){
      condition = "random"
    } else {
      condition = "not random"
    }
    
    # compute consensus
    cl <- makeClusterPSOCK(30)
    future::plan(future::cluster, workers = cl, gc = T)
    consensus = future_map(1:100, function(x) select_areas(x, clusters = clusters, model = model, condition = condition, model_id = model_id),
                           .progress = TRUE,
                           .options = furrr_options(seed = TRUE))
    consensus = do.call(rbind, consensus)
    parallel::stopCluster(cl)
    
    consensus_list[[i]] = consensus
    
  }
  
  consensus_df = dplyr::bind_rows(consensus_list)
  
  options(future.globals.maxSize = 500*1024^2)
  
  # save result
  saveRDS(consensus_df, "data/output/consensus.rds")
  
  # read consensus df
  consensus_df = readRDS("data/output/consensus.rds")
  
  ################################################################################
  # data consensus: repeat same analysis on actual data points
  ################################################################################
  # intersections data...
  # compute consensus
  consensus_data = clusters %>% filter(condition == "data")
  intersections_data = list()
  c=1
  for (j in unique(consensus_data$unique_id)){
    
    # subset polygon
    target = consensus_data %>% filter(unique_id == j)
    target = st_as_sf(target)
    
    # subset reference set
    reference = consensus_data %>% filter(day == target$day & unique_id != target$unique_id)
    reference = st_as_sf(reference)
    
    # sample points from polygon surface
    area = st_area(target)
    area = units::set_units(area, m^2)
    pts = st_sample(target, size = round(max(1, drop_units(area / 100))), type = "regular", replace = F)
    
    # how many intersections with other polygons for each point
    test = st_contains(reference, pts)
    consensus = length(unlist(test)) / (area %>% drop_units() / 100)
    
    intersections_data[[c]] = data.frame(unique_id = target$unique_id,
                                         consensus = consensus,
                                         day = target$day,
                                         camera_id = target$camera_id,
                                         confidence = target$confidence,
                                         model_id = "data")
    
    # print
    print(paste(c, "/", length(unique(consensus_data$unique_id)), sep = ""))
    c=c+1
  }
  
  consensus_data = dplyr::bind_rows(intersections_data)
  consensus_data$iteration = 0
  
  # add to full dataframe
  consensus_df = bind_rows(consensus_df, consensus_data)
  
  # save
  saveRDS(consensus_df, "data/output/consensus.rds")
  
  ################################################################################
  # illustrate
  ################################################################################
  # visualize area consensus
  select_areas_illustration <- function(clusters, model, condition, model_id){
    
    # extract model parameters
    if (condition == "random"){
      
    } else {
      main_effects = model %>% spread_draws(weights[feature]) %>% group_by(feature) %>% summarize(weights = mean(weights))
      lake_effects = model %>% spread_draws(v_lake[lake, feature]) %>% group_by(lake,feature) %>% summarize(weights = mean(v_lake))
      id_effects = model %>% spread_draws(v_id[id, feature]) %>% group_by(id, feature) %>% summarize(weights = mean(v_id))
      trip_effects = model %>% spread_draws(v_trip[trip, feature]) %>% group_by(trip, feature) %>% summarize(weights = mean(v_trip))
    }
    
    # loop through all clusters
    polygons = list()
    c = 1
    for (i in unique(clusters$unique_id)){
      
      # subset simulated clusters
      sub_clusters = clusters %>% filter(unique_id == i) %>% filter(condition == "simulated")
      
      # identifiers
      trip_tmp = sub_clusters$trip[1]
      lake_tmp = sub_clusters$lake[1]
      id_tmp = sub_clusters$id[1]
      
      # feature matrix
      f = sub_clusters[, c("displacement_scaled", "boundary_scaled", "depth_scaled", "shore_scaled", "roughness_scaled")] %>% st_drop_geometry()
      
      # compute feature weights
      if (condition == "random"){
        w = t(t(c(0, 0, 0, 0, 0)))
      } else {
        w = main_effects$weights + lake_effects[lake_effects$lake == lake_tmp,"weights"] + trip_effects[trip_effects$trip == trip_tmp,"weights"]  + id_effects[id_effects$id == id_tmp,"weights"]
      }
      
      # select polygon
      p = softmax(as.matrix(f)[, 1:nrow(w)] %*% as.matrix(w)) 
      
      # sample
      choice = sample(1:length(p), 1, prob = p)
      
      # select polygon
      polygons[[c]] = sub_clusters[choice,]
      
      # count
      c = c + 1
      
    }
    
    # to df
    polygons = dplyr::bind_rows(polygons)
    
    return(polygons)
    
  }
  
  save_locations = c("data/processed_data/model_fit/nomination_models/a_random_fit.rds",
                     "data/processed_data/model_fit/nomination_models/b_displacement_fit.rds",
                     "data/processed_data/model_fit/nomination_models/c_boundary_fit.rds",
                     "data/processed_data/model_fit/nomination_models/d_depth_fit.rds",
                     "data/processed_data/model_fit/nomination_models/e_shore_fit.rds",
                     "data/processed_data/model_fit/nomination_models/f_roughness_fit.rds")
  
  polygon_list = list()
  for (i in 1:length(save_locations)){
    
    # load model
    model = readRDS(save_locations[i])
    model_id = c("random", "+displacement", "+boundary", "+depth", "+shore", "+roughness")[i]
    
    if (i == 1){
      condition = "random"
    } else {
      condition = "not random"
    }
    
    # compute consensus
    polygons = select_areas_illustration(clusters = clusters, model = model, condition = condition, model_id = model_id)
    polygon_list[[i]] = polygons
    
  }
  
  # make visualizations
  competition_areas = readRDS(file = "data/study_lakes.rds")
  depth_profiles = readRDS("data/depth_profiles.rds")
  c=1
  plots_ill = list()
  for (i in 1:5){
    
    # take example lake
    example_lake = depth_profiles[[i]]
    example_area = competition_areas[[i]]
    
    example_lake <- st_transform(example_lake, st_crs(example_area))
    
    # 1: Data
    plots_ill[[c]] = ggplot() + 
        geom_sf(data = st_intersection(example_area,example_lake), color = "lightgrey", size = 0.5) + 
        geom_sf(data = st_as_sf(clusters %>% filter(day == i & condition == "data")), fill = "lightblue", alpha = .6) +
        scale_fill_manual(values = colors, guide = F)+
        xlab("Longitude") + 
        ylab("Latitude") + 
        theme(panel.background = element_rect(fill = NA),
              axis.title.x = element_text(),
              axis.title.y = element_text(),
              axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              title = element_text(size = 5, margin=margin(0,0,0,0)),
              plot.tag = element_text(size = 7, face = "bold", margin=margin(0,0,0,0)),
              axis.text.y = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA),
              legend.position = "none",
              plot.tag.position = "topleft",
              legend.key.size = unit(2, "lines"),
              legend.text = element_text(size = 3),
              legend.title = element_text(size = 3, margin=margin(0,0,0,0)),
              legend.background = element_rect(colour = "black", fill = "white"),
              plot.margin=unit(c(.1, .2, .1, .2), "cm"))+
        labs(tag = ifelse(c %% 7 == 1, LETTERS[round(c/7) + 1], " ")) + 
        ggtitle("Data")
    c = c + 1
    
    # 2: Random
    plots_ill[[c]] = ggplot() + 
        geom_sf(data = st_intersection(example_area,example_lake), color = "lightgrey", size = 0.5) + 
        geom_sf(data = st_as_sf(polygon_list[[1]] %>% filter(day == i)), fill = "lightblue", alpha = .6) +
        scale_fill_manual(values = colors, guide = F)+
        xlab("Longitude") + 
        ylab("Latitude") + 
        theme(panel.background = element_rect(fill = NA),
              axis.title.x = element_text(),
              axis.title.y = element_text(),
              axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              title = element_text(size = 5, margin=margin(0,0,0,0)),
              plot.tag = element_text(size = 7, face = "bold", margin=margin(0,0,0,0)),
              axis.text.y = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA),
              legend.position = "none",
              plot.tag.position = "topleft",
              legend.key.size = unit(2, "lines"),
              legend.text = element_text(size = 3),
              legend.title = element_text(size = 3, margin=margin(0,0,0,0)),
              legend.background = element_rect(colour = "black", fill = "white"),
              plot.margin=unit(c(.1, .2, .1, .2), "cm"))+
      theme(plot.tag = element_text("A", color = "white"))+
      ggtitle("Random")
    c = c + 1
    
    # 3: Displacement
    plots_ill[[c]] = ggplot() + 
        geom_sf(data = st_intersection(example_area,example_lake), color = "lightgrey", size = 0.5) + 
        geom_sf(data = st_as_sf(polygon_list[[2]] %>% filter(day == i)), fill = "lightblue", alpha = .6) +
        scale_fill_manual(values = colors, guide = F)+
        xlab("Longitude") + 
        ylab("Latitude") + 
        theme(panel.background = element_rect(fill = NA),
              axis.title.x = element_text(),
              axis.title.y = element_text(),
              axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              title = element_text(size = 5, margin=margin(0,0,0,0)),
              plot.tag = element_text(size = 7, face = "bold", margin=margin(0,0,0,0)),
              axis.text.y = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA),
              legend.position = "none",
              plot.tag.position = "topleft",
              legend.key.size = unit(2, "lines"),
              legend.text = element_text(size = 3),
              legend.title = element_text(size = 3, margin=margin(0,0,0,0)),
              legend.background = element_rect(colour = "black", fill = "white"),
              plot.margin=unit(c(.1, .2, .1, .2), "cm"))+
      theme(plot.tag = element_text("A", color = "white"))+
      ggtitle("+Displacement")
    c = c + 1
    
    # 4: Boundary
      plots_ill[[c]] = ggplot() + 
        geom_sf(data = st_intersection(example_area,example_lake), color = "lightgrey", size = 0.5) + 
        geom_sf(data = st_as_sf(polygon_list[[3]] %>% filter(day == i)), fill = "lightblue", alpha = .6) +
        scale_fill_manual(values = colors, guide = F)+
        xlab("Longitude") + 
        ylab("Latitude") + 
        theme(panel.background = element_rect(fill = NA),
              axis.title.x = element_text(),
              axis.title.y = element_text(),
              axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              title = element_text(size = 5, margin=margin(0,0,0,0)),
              plot.tag = element_text(size = 7, face = "bold", margin=margin(0,0,0,0)),
              axis.text.y = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA),
              legend.position = "none",
              plot.tag.position = "topleft",
              legend.key.size = unit(2, "lines"),
              legend.text = element_text(size = 3),
              legend.title = element_text(size = 3, margin=margin(0,0,0,0)),
              legend.background = element_rect(colour = "black", fill = "white"),
              plot.margin=unit(c(.1, .2, .1, .2), "cm"))+
        theme(plot.tag = element_text("A", color = "white"))+
        ggtitle("+Boundary")
    c = c + 1
    
    # 5: Depth
      plots_ill[[c]] = ggplot() + 
        geom_sf(data = st_intersection(example_area,example_lake), color = "lightgrey", size = .5) + 
        geom_sf(data = st_as_sf(polygon_list[[4]] %>% filter(day == i)), fill = "lightblue", alpha = .6) +
        scale_fill_manual(values = colors, guide = F)+
        xlab("Longitude") + 
        ylab("Latitude") + 
        theme(panel.background = element_rect(fill = NA),
              axis.title.x = element_text(),
              axis.title.y = element_text(),
              axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              title = element_text(size = 5, margin=margin(0,0,0,0)),
              plot.tag = element_text(size = 7, face = "bold", margin=margin(0,0,0,0)),
              axis.text.y = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA),
              legend.position = "none",
              plot.tag.position = "topleft",
              legend.key.size = unit(2, "lines"),
              legend.text = element_text(size = 3),
              legend.title = element_text(size = 3, margin=margin(0,0,0,0)),
              legend.background = element_rect(colour = "black", fill = "white"),
              plot.margin=unit(c(.1, .2, .1, .2), "cm"))+
        theme(plot.tag = element_text("A", color = "white"))+
        ggtitle("+Depth")
    c = c + 1
    
    # 6: Shore
    plots_ill[[c]] = ggplot() + 
        geom_sf(data = st_intersection(example_area,example_lake), color = "lightgrey", size = .5) + 
        geom_sf(data = st_as_sf(polygon_list[[5]] %>% filter(day == i)), fill = "lightblue", alpha = .6) +
        scale_fill_manual(values = colors, guide = F)+
        xlab("Longitude") + 
        ylab("Latitude") + 
        theme(panel.background = element_rect(fill = NA),
              axis.title.x = element_text(),
              axis.title.y = element_text(),
              axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              title = element_text(size = 5, margin=margin(0,0,0,0)),
              plot.tag = element_text(size = 7, face = "bold", margin=margin(0,0,0,0)),
              axis.text.y = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA),
              legend.position = "none",
              plot.tag.position = "topleft",
              legend.key.size = unit(2, "lines"),
              legend.text = element_text(size = 3),
              legend.title = element_text(size = 3, margin=margin(0,0,0,0)),
              legend.background = element_rect(colour = "black", fill = "white"),
              plot.margin=unit(c(.1, .2, .1, .2), "cm"))+
      theme(plot.tag = element_text("A", color = "white"))+
      ggtitle("+Shore")
    c = c + 1
    
    # 7: Roughness
    plots_ill[[c]] = ggplot() + 
        geom_sf(data = st_intersection(example_area,example_lake), color = "lightgrey", size = .5) + 
        geom_sf(data = st_as_sf(polygon_list[[6]] %>% filter(day == i)), fill = "lightblue", alpha = .6) +
        scale_fill_manual(values = colors, guide = F)+
        xlab("Longitude") + 
        ylab("Latitude") + 
        theme(panel.background = element_rect(fill = NA),
              axis.title.x = element_text(),
              axis.title.y = element_text(),
              axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              title = element_text(size = 5, margin=margin(0,0,0,0)),
              plot.tag = element_text(size = 7, face = "bold", margin=margin(0,0,0,0)),
              axis.text.y = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA),
              legend.position = "none",
              plot.tag.position = "topleft",
              legend.key.size = unit(2, "lines"),
              legend.text = element_text(size = 3),
              legend.title = element_text(size = 3, margin=margin(0,0,0,0)),
              legend.background = element_rect(colour = "black", fill = "white"),
              plot.margin=unit(c(.1, .2, .1, .2), "cm"))+
      theme(plot.tag = element_text("A", color = "white"))+
      ggtitle("+Roughness")
    c = c + 1
    
  }
  
  # arrange panel
  panel_simulations = gridExtra::arrangeGrob(
    grobs = plots_ill,
    heights = c(1, 1, 1, 1, 1),
    layout_matrix = t(matrix(1:35, ncol = 5))
  )
  
  ggsave(panel_simulations, file = "5_plots/panel_simulations.png", 
         width = 18, height = 15, bg = "white", dpi = 300, units = "cm")
  ggsave(panel_simulations, file = "5_plots/panel_simulations.pdf", 
         width = 18, height = 15, bg = "white", dpi = 300, units = "cm")
  ggsave(panel_simulations, file = "5_plots/panel_simulations.svg", 
         width = 18, height = 15, bg = "white", dpi = 300, units = "cm", device = svglite::svglite)

}
