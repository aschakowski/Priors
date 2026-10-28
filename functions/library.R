## Library

## If a package is installed, it will be loaded. If any 
## are not, the missing package(s) will be installed 
## from CRAN and then loaded.

## First specify the packages of interest
packages = c("readr",     #read data
             "tidyverse", #data processing
             "lubridate", #date transformations
             "eeptools",  #date transformations
             "ggplot2",   #visualizations
             "reshape2",  #data transformations
             "geosphere",
             "Rcpp",
             "devtools",
             "loo",
             "abind",
             "terra",
             "dbscan",
             "momentuHMM",
             "data.table",
             "igraph",
             "gridExtra",
             "tidybayes",
             "concaveman",
             "future",
             "furrr",
             "purrr",
             "ggforce",
             #"ggalt",
             "forcats",
             "readxl",
             "rstan",
             "brms",
             "cmdstanr",
             "RColorBrewer",
             "raster",
             "sf",
             "sp",
             "fasterize",
             "R.utils",
             "gtools",
             "fields",
             "sp",
             "units",
             "rasterize",
             "png",
             "alphahull"
)

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

library(devtools)
devtools::install_github("nschiett/fishualize", force = TRUE)
library(fishualize)

############################# helper functions #################################

# compute mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# inverse logit
inv_logit <- function(x){
  exp(x)/(1+exp(x))
}

# softmax
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

# compute x and y coordinates#Function: https://stackoverflow.com/questions/18639967/converting-latitude-and-longitude-points-to-utm
LongLatToUTM<-function(x,y,zone){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example# +datum=WGS84
  res <- spTransform(xy, CRS(paste("+proj=utm +datum=WGS84 +zone=",zone," ellps=WGS84",sep='')))
  return(as.data.frame(res))
}


# heading angles
heading_angle <- function(x, y){
  v <- c(y[1]-x[1],y[2]-x[2])
  angle <- atan2(v[2],v[1])
  
  while(angle<=(-pi))
    angle <- angle + 2*pi
  while(angle>pi)
    angle <- angle - 2*pi
  
  return(angle)
}


rm(list = c("package.check", "packages"))

# https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/


# key linking lake id 2022 and 2023
key = cbind(y2023 = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
            y2022 = c(7, 10, 6, 9, 8, 1, 2, 3, 5, 4))


#install.packages("OpenImageR")
#install.packages("magick")
#library(magick)
#library(OpenImageR)
