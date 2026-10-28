################################################################################
#
# TITLE: Extract ecological knowledge from scanned maps
# 
# AUTHORS: Schakowski, A. (schakowski@mpib-berlin.mpg.de)
#
# DESCRIPTION: This script takes the processed lake maps with intensified colors 
# and returns dataframes where map coordinates are matched to the respective 
# colors used on the map. The script also creates images to check with the original
# map for any errors in the script. 
#
# LAST UPDATED (DD/MM/YYYY): 06/10/2026
#
################################################################################
extract_ecological_knowledge_data <- function(){

# We will use alpha hull computation to detect the marked areas. We need some
# technical helper functions for class conversion to the right format.
# These functions are taken from the credited source below. 
################################### HELPER #####################################
# SOURCE: https://babichmorrowc.github.io/post/2019-03-18-alpha-hull/

# function to convert an arc into line segments
# given the center of the arc, the radius, the vector, and the angle (radians)
arc2line <- function(center, r, vector, theta, npoints = 100) {
  # Get the angles at the extremes of the arcs
  angles <- anglesArc(vector, theta)
  # Generate sequence of angles along the arc to determine the points
  seqang <- seq(angles[1], angles[2], length = npoints)
  # Generate x coordinates for points along the arc
  x <- center[1] + r * cos(seqang)
  # Generate y coordinates for points along the arc
  y <- center[2] + r * sin(seqang)
  coords.xy <- cbind(x,y)
  line <- Line(coords = coords.xy)
  return(line)
}

ahull2lines <- function(hull){
  arclist <- hull$arcs
  lines <- list()
  for (i in 1:nrow(arclist)) {
    # Extract the attributes of arc i
    center_i <- arclist[i, 1:2]
    radius_i <- arclist[i, 3]
    vector_i <- arclist[i, 4:5]
    theta_i <- arclist[i, 6]
    # Convert arc i into a Line object
    line_i <- arc2line(center = center_i, r = radius_i, vector = vector_i, theta = theta_i)
    list_length <- length(lines)
    if(list_length > 0){
      # If a line has already been added to the list of lines
      # Define last_line_coords as the coordinates of the last line added to the list before the ith line
      last_line_coords <- lines[[list_length]]@coords
    }
    if(i == 1){
      # Add the first line to the list of lines
      lines[[i]] <- line_i
    } else if(isTRUE(all.equal(line_i@coords[1,], last_line_coords[nrow(last_line_coords),]))){
      # If the first coordinate in the ith line is equal to the last coordinate in the previous line
      # then those lines should be connected
      # Row bind the coordinates for the ith line to the coordinates of the previous line in the list
      lines[[list_length]]@coords <- rbind(last_line_coords, line_i@coords[2:nrow(line_i@coords),])
    } else {
      # If the first coordinate in the ith line does not match the last coordinate in the previous line
      # then the ith line represents a new line
      # Add the ith line to the list as a new element
      lines[[length(lines) + 1]] <- line_i
    }
  }
  # Convert the list of lines to a Line object
  lines <- Lines(lines, ID = 'l')
  # Convert the Line object to a SpatialLines object
  sp_lines <- SpatialLines(list(lines))
  return(sp_lines)
}


spLines2poly <- function(sp_lines){
  # Extract the lines slot
  lines_slot <- sp_lines@lines[[1]]
  # # Create a list of booleans indicating whether a given Line represents a polygon
  # poly_bool <- sapply(lines_slot@Lines, function(x){
  #   coords <- lines_slot@Lines[[1]]@coords
  #   # coords <- x@coords
  #   # Check if the first coordinate in the line is the same as the last
  #   all.equal(coords[1,], coords[nrow(coords),])
  # })
  # # Pull out the lines that form polygons
  # poly_lines <- sp_lines[poly_bool]
  # poly_lines_slot <- poly_lines@lines
  # Create SpatialPolygons
  sp_polys <- SpatialPolygons(list(Polygons(lapply(lines_slot@Lines, function(x) {
    Polygon(slot(x, "coords"))
  }
  ),
  ID = "1")))
  return(sp_polys)
}

ahull2poly <- function(hull){
  # Convert the alpha hull to SpatialLines
  hull2SpatialLines <- ahull2lines(hull)
  # Convert SpatialLines to SpatialPolygon
  SpatialLines2SpatialPolygon <- spLines2poly(hull2SpatialLines)
  return(SpatialLines2SpatialPolygon)
}

################################################################################
# The following function takes scanned pages and convert them to polygon objects
################################################################################
convert_scan_to_df <- function(scan, #file path of processed image 
                               tlc, #coordinates of top left map corner
                               brc, #coordinates of bottom right map corner
                               file, #file name of processed image
                               control_path, #path were control image is stored
                               plot, #should a control plot be saved
                               option
){
  
  # load scanned page
  map = readPNG(paste(scan))  
  
  # keep only 3 dimensions, R, G, and B
  map = map[, , 1:3]
  
  ## PART 1: Rotate and crop image
  
  # Rotate image if necessary
  # find the red dot outside the map
  
  # reshape image into a data frame
  df = data.frame(
    red = matrix(map[,,1], ncol=1),
    green = matrix(map[,,2], ncol=1),
    blue = matrix(map[,,3], ncol=1)
  )
  
  # create pixel id and add row and column id of pixel
  df$id = 1:nrow(df)
  df$col = rep(1:dim(map)[2], each = dim(map)[1])
  df$row = rep(1:dim(map)[1], dim(map)[2])
  
  # define noise term
  # modulates the color intensity of the blue color to be detected as blue
  epsilon = .6
  
  # filter blue points
  df <- df %>% 
    mutate(is_blue = blue > epsilon & blue > (green + red))
  
  # compute dbscan clustering to detect the blue points in the corners
  clusters = dbscan::dbscan(df[df$is_blue, c("col", "row")],eps = 10) 
  
  # assign cluster label to blue points
  clustered_data = df[df$is_blue, ]
  clustered_data$label = clusters$cluster
  clustered_data = clustered_data %>% filter(label != 0)
  
  # compute map extent
  map_extent = clustered_data %>% 
    group_by(label) %>% 
    summarize(r = mean(row), 
              c = mean(col)) %>% 
    summarize(top_left_y = round(.5*(r[rank(r)==2] + r[rank(r)==1])),
              top_left_x = round(.5*(c[rank(c)==2] + c[rank(c)==1])),
              bottom_right_y = round(.5*(r[rank(r)==3] + r[rank(r)==4])),
              bottom_right_x = round(.5*(c[rank(c)==3] + c[rank(c)==4])))
  
  # crop image
  df = df %>% 
    mutate(is_map = row > pull(map_extent[1]) & 
             row < pull(map_extent[3]) & 
             col > pull(map_extent[2]) & 
             col < pull(map_extent[4]))
  
  # define noise term: modulates color intensity threshold to be detected as red
  epsilon = .6
  df <- df %>% 
    mutate(is_red = !is_map & red > epsilon & red > (blue + green))
  
  # find position of red dot
  red_dot = c(mean(df[df$is_red, "row"]), mean(df[df$is_red, "col"]))
  
  # what is position of red dot relative to map extent
  check1 = map_extent[1] > red_dot[1] #if TRUE, then red dot on top
  check2 = map_extent[2] > red_dot[2] #if TRUE, then red dot on the left side
  check3 = map_extent[3] < red_dot[1] #if TRUE, then red dot on the bottom
  check4 = map_extent[4] < red_dot[2] #if TRUE, then red dot on the right
  
  # rotate to right orientation
  # orientation correct
  if (check1){
    
  } else if (check2){ #rotate by 90 degrees clockwise
    map = rotateFixed(map, 90)
    
  } else if (check3){ #rotate by 180 degrees
    map = rotateFixed(map, 180)
    
  } else if (check4){ #rotate by 90 degrees counterclockwise
    map = rotateFixed(map, 270)
    
  }
  
  # reshape image into a data frame
  df = data.frame(
    red = matrix(map[,,1], ncol=1),
    green = matrix(map[,,2], ncol=1),
    blue = matrix(map[,,3], ncol=1)
  )
  
  # create pixel id and add row and column id of pixel
  df$id = 1:nrow(df)
  df$col = rep(1:dim(map)[2], each = dim(map)[1])
  df$row = rep(1:dim(map)[1], dim(map)[2])
  
  # define noise term
  # modulates the color intensity of the blue color to be detected as blue
  epsilon = .6
  
  # filter blue points
  df <- df %>% 
    mutate(is_blue = blue > epsilon & blue > (green + red))
  
  # compute dbscan clustering to detect the blue points in the corners
  clusters = dbscan::dbscan(df[df$is_blue, c("col", "row")],eps = 10) 
  
  # assign cluster label to blue points
  clustered_data = df[df$is_blue, ]
  clustered_data$label = clusters$cluster
  clustered_data = clustered_data %>% filter(label != 0)
  
  # compute map extent
  map_extent = clustered_data %>% 
    group_by(label) %>% 
    summarize(r = mean(row), 
              c = mean(col)) %>% 
    summarize(top_left_y = round(.5*(r[rank(r)==2] + r[rank(r)==1])),
              top_left_x = round(.5*(c[rank(c)==2] + c[rank(c)==1])),
              bottom_right_y = round(.5*(r[rank(r)==3] + r[rank(r)==4])),
              bottom_right_x = round(.5*(c[rank(c)==3] + c[rank(c)==4])))
  
  # crop image
  df = df %>% 
    mutate(is_map = row > pull(map_extent[1]) & 
             row < pull(map_extent[3]) & 
             col > pull(map_extent[2]) & 
             col < pull(map_extent[4]))
  
  ## PART 2: Detect marked areas
  # subset map
  df_map = df %>% filter(is_map)
  
  # define noise
  epsilon = .8
  
  # exclude grey scale
  df_map = df_map %>% 
    mutate(exclude = abs(red - green) < epsilon & blue > .5*(red + green)-epsilon | (blue > .5 & (red < .2 | green < .2)))
  
  # subset areas
  df_areas = df_map %>% 
    filter(exclude == F)
  
  # break loop if missing data (i.e., no marked areas)
  if (nrow(df_areas)!=0){
    
    # compute the dbscan clustering to separate the colors
    eps = .3
    clusters = dbscan::dbscan(df_areas[, c("red", "green", "blue")], eps = eps)
    
    while((!(length(unique(clusters$cluster)) <= 3 & !0 %in% unique(clusters$cluster)) & !(length(unique(clusters$cluster)) == 4 & 0 %in% unique(clusters$cluster))) & eps > .01){
      clusters = dbscan::dbscan(df_areas[, c("red", "green", "blue")], eps = eps)
      eps = eps - .01
    }
    
    # each identified color gets a label. Labels are not yet linked to colors.
    df_areas$label = clusters$cluster
    
    # filter not assigned clusters
    df_areas = df_areas %>% filter(label != 0)
    
    # replacing colors with the centroids of the identified clusters. #
    # each label gets a distinct color value
    colors = df_areas %>% 
      group_by(label) %>% 
      summarize(R = mean(red),
                G = mean(green),
                B = mean(blue))
    
    # apply dbscan to detect spatially connected areas
    clusters = dbscan::dbscan(df_areas[, c("row", "col")], eps = 3)
    
    # assign clusters to data
    df_areas$clusters = clusters$cluster
    df_areas = df_areas %>% filter(clusters != 0)
    
    # get spatial extent of each cluster (compute alpha hull)
    hull <- df_areas %>%
      group_by(clusters, label) %>% 
      st_as_sf(coords = c("col", "row")) %>%
      summarize(geometry = st_union(geometry))
    
    hull_new = list()
    for (k in 1:nrow(hull)){
      hull_new[[k]] <- st_as_sf(ahull2poly(alphahull::ahull(st_coordinates(hull[k,])[,1:2], alpha = 3)))
    }
    
    hull_new = do.call(rbind, hull_new)
    hull_new$label = hull$label
    hull_new$clusters = hull$clusters
    hull = hull_new
    
    # for each pixel, we can now determine in which cluster they lie
    points = df_map[, c("col", "row")]
    cluster_id = st_intersects(hull, points %>% st_as_sf(coords = c("col", "row"))) 
    
    # we need a key dataframe that tells us the color label of each cluster- 
    # later colors are assigned to key
    key = cbind(hull$clusters, hull$label)
    
    # add cluster label to original df
    n_clusters = length(cluster_id)
    df_map$label_noise_reduction = NA
    df_map$cluster_id = NA
    for (i in 1:n_clusters){
      df_map$label_noise_reduction[cluster_id[[i]]] = key[i, 2]
      df_map$cluster_id[cluster_id[[i]]] = i
    }
    
    # match number of pixels in y,x dimension to coordinates, units are m
    y_diff = tlc[1] - brc[1]
    x_diff = brc[2] - tlc[2]
    
    # compute pixel size
    pixel_size_y = y_diff/length(unique(df_map$row))
    pixel_size_x = x_diff/length(unique(df_map$col))
    
    # if we did everything right, then the pixel size should be the same in both dimensions. Here the error is very low (~1cm)
    #pixel_size_y - pixel_size_x
    
    # obtain quadratic pixels, by averaging x/y pixel size
    pixel_size = .5 * (pixel_size_x + pixel_size_y)
    
    # reduce noise by averaging pixels computing the boundaries from all four starting points
    # add coordinates to dataframe
    df_map$pixel_center_y1 = tlc[1] - (abs(df_map$row - min(df_map$row)) + 1)*pixel_size - .5 * pixel_size
    df_map$pixel_center_x1 = tlc[2] + (abs(df_map$col - min(df_map$col)) + 1)*pixel_size + .5 * pixel_size
    
    # recompute coordinates starting from other corners
    # starting top right
    df_map$pixel_center_y2 = tlc[1] - (abs(df_map$row - min(df_map$row)) + 1)*pixel_size - .5 * pixel_size
    df_map$pixel_center_x2 = brc[2] - (abs(df_map$col - max(df_map$col)) + 1)*pixel_size - .5 * pixel_size
    
    # starting bottom right
    df_map$pixel_center_y3 = brc[1] + (abs(df_map$row - max(df_map$row)) + 1)*pixel_size + .5 * pixel_size
    df_map$pixel_center_x3 = brc[2] - (abs(df_map$col - max(df_map$col)) + 1)*pixel_size - .5 * pixel_size
    
    # starting bottom left
    df_map$pixel_center_y4 = brc[1] + (abs(df_map$row - max(df_map$row)) + 1)*pixel_size - .5 * pixel_size
    df_map$pixel_center_x4 = tlc[2] + (abs(df_map$col - min(df_map$col)) + 1)*pixel_size + .5 * pixel_size
    
    # average coordinates
    df_map = df_map %>% 
      mutate(pixel_center_x = 1/4 * (pixel_center_x1 + pixel_center_x2 + 
                                       pixel_center_x3 + pixel_center_x4),
             pixel_center_y = 1/4 * (pixel_center_y1 + pixel_center_y2 + 
                                       pixel_center_y3 + pixel_center_y4))
    # assign color label
    epsilon = .9
    colors$color = ifelse(colors$R > colors$G & colors$R > colors$B & abs(colors$R - colors$G) >= epsilon, "red",
                          ifelse(colors$G > colors$R & colors$G > colors$B & abs(colors$G - colors$R) >= epsilon, "green",
                                 ifelse(colors$R > colors$B & 
                                          colors$G > colors$B & abs(colors$R - colors$G) < (1-epsilon),"yellow", 
                                        NA)))
    colors$label_noise_reduction = colors$label
    
    # add to dataframe
    df_map = left_join(df_map, colors)
    
    # add camera id and day
    df_map$day = strsplit(strsplit(file, "_")[[1]][1], "D")[[1]][2]
    df_map$camera_id = strsplit(strsplit(strsplit(file, "_")[[1]][2], "[.]")[[1]][1], "C")[[1]][2]
    
    # convert hull to coordinates using same formula
    polygons = list()
    for (h in 1:nrow(hull)){
      
      coordinates_hull <- as.data.frame(st_coordinates(hull[h,])[, 1:2])
      
      coordinates_hull$y1 = tlc[1] - (abs(coordinates_hull$Y - min(df_map$row)) + 1)*pixel_size - .5 * pixel_size
      coordinates_hull$x1 = tlc[2] + (abs(coordinates_hull$X - min(df_map$col)) + 1)*pixel_size + .5 * pixel_size
      
      # recompute coordinates starting from other corners
      # starting top right
      coordinates_hull$y2 = tlc[1] - (abs(coordinates_hull$Y - min(df_map$row)) + 1)*pixel_size - .5 * pixel_size
      coordinates_hull$x2 = brc[2] - (abs(coordinates_hull$X - max(df_map$col)) + 1)*pixel_size - .5 * pixel_size
      
      # starting bottom right
      coordinates_hull$y3 = brc[1] + (abs(coordinates_hull$Y - max(df_map$row)) + 1)*pixel_size + .5 * pixel_size
      coordinates_hull$x3 = brc[2] - (abs(coordinates_hull$X - max(df_map$col)) + 1)*pixel_size - .5 * pixel_size
      
      # starting bottom left
      coordinates_hull$y4 = brc[1] + (abs(coordinates_hull$Y - max(df_map$row)) + 1)*pixel_size - .5 * pixel_size
      coordinates_hull$x4 = tlc[2] + (abs(coordinates_hull$X - min(df_map$col)) + 1)*pixel_size + .5 * pixel_size
      
      # average coordinates
      coordinates_hull = coordinates_hull %>% 
        mutate(X = 1/4 * (x1 + x2 + x3 + x4),
               Y = 1/4 * (y1 + y2 + y3 + y4))
      
      # remove 
      coordinates_hull = coordinates_hull %>% dplyr::select(X,Y)
      
      # convert to df
      polygons[[h]] = st_as_sf(coordinates_hull %>% st_as_sf(coords = c("X", "Y")) %>% st_combine() %>% st_cast("LINESTRING") %>% st_cast("POLYGON"))
      
    }
    
    # merge to df
    polygons = do.call(rbind, polygons)
    polygons$label = hull$label
    polygons$clusters = hull$clusters
    polygons$camera_id = df_map$camera_id[1]
    polygons$day = df_map$day[1]
    
    # cluster id: counts number of areas per id
    polygons$nomination_id = 1:nrow(polygons)
    
    # add color label
    polygons = left_join(polygons, colors %>% dplyr::select(label, color) %>% mutate(confidence = ifelse(color == "red", 3, ifelse(color == "yellow", 2, ifelse(color == "green", 1, NA)))))
    
    # only save relevant variables
    polygons = polygons %>% dplyr::select(-label, -clusters)
    
    # control plot
    if (plot){
      control_plot = df_map %>% 
        filter(!is.na(color)) %>% 
        ggplot(aes(x = pixel_center_x, y = pixel_center_y, color = as.factor(color))) + 
        scale_color_identity() + 
        geom_point(alpha = .4) +
        xlab("X") + 
        ylab("Y") + 
        coord_cartesian(xlim = c(tlc[2], brc[2]), ylim = c(brc[1], tlc[1]))+
        theme_minimal() + 
        theme(legend.position = c(.2, .2),
              aspect.ratio = 1)
      
      # save control image
      ggsave(file = paste(control_path, file, sep = "/"), control_plot, bg = "white")
    }
    
  } else {
    df_map = NULL
    polygons = NULL
  }
  # return
  if (option == "polygons"){
    return(polygons)
  } else {
    return(df_map)
  }
  
}
################################################################################
# define coordinates of map excerpts for all 5 days.
tlc_all = list(c(6945000, 637000),
               c(6920000, 650000),
               c(6918000, 652000),
               c(6948000, 633000),
               c(6973000, 639000))
brc_all = list(c(6942000, 640000),
               c(6917000, 653000),
               c(6916000, 654000),
               c(6944000, 636000),
               c(6970000, 643000))

################################################################################
# create polygons for each day
POLYGONS = list()
for (day in 1:5){
  
  # file were processed images are saved
  path = paste(getwd(),"/data/map_data/processed/", "day", day, sep = "")
  
  # file to store control images to check that algorithm is working properly
  control_path = paste(getwd(),"/data/map_data/control/", "day", day, sep = "")
  
  # reduce image sizes and store in reduced folders
  path_reduced = paste(getwd(),"/data/map_data/reduced/", "day", day, sep = "")
  
  # load files
  files = mixedsort(sort(list.files(path, pattern = ".png")))
  path_scan = paste(path, files, sep = "/")
  
  # rescale images
  for (i in 1:length(path_scan)){
    
    # resize image
    map_tmp <- image_read(paste(path, files[i], sep = "/"))
    map <- image_scale(map_tmp, "800")
    
    # save new file in new folder
    image_write(map, path = paste(path_reduced, files[i], sep = "/"), format = "png")
  }
  
  # start processing images
  files = mixedsort(sort(list.files(path_reduced, pattern = ".png")))
  
  polygons = list()
  for (i in 1:length(files)){
    polygons[[i]] <- convert_scan_to_df(scan = paste(path_reduced, files[i], sep = "/"), 
                                      tlc = tlc_all[[day]], 
                                      brc = brc_all[[day]], 
                                      file = files[i], 
                                      control_path = control_path,
                                      plot = T, 
                                      option = "polygons")
    
  }
  
  # convert to df
  polygons = do.call(rbind, polygons)
  
  # add to list
  POLYGONS[[day]] <- polygons
}

# gc
gc()

# convert to df
POLYGONS = do.call(rbind, POLYGONS)

# gc
gc()

# remove confusing label and cluster variables.
# nomination counts number of areas per id
#POLYGONS$nomination = POLYGONS$clusters
#POLYGONS = POLYGONS %>% dplyr::select(-label, -clusters)

# set coordinate reference system
st_crs(POLYGONS) <-"epsg:3067"
saveRDS(POLYGONS, "data/ecological_knowledge_polygons.rds")

}
###################################### END #####################################