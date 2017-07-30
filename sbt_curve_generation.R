####################################    Space Beats Time Research  #######################################
################################ Generating datasets for examples of curves #######################################
#This script generates datasets for the Space Beats Time Framework.
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 07/28/2017 
#DATE MODIFIED: 07/30/2017
#Version: 1
#PROJECT:  with Marco Millones            
#
#COMMENTS: - Generate additional datasets
#
#TO DO:
# - 
#COMMIT:  
#
#################################################################################################

###Loading R library and packages                                                      

library(sp)
library(rgdal)
library(spdep)
library(BMS) #contains hex2bin and bin2hex
library(bitops)
library(gtools)
library(maptools)
library(parallel)
library(rasterVis)
library(raster)
library(zoo)  #original times series function functionalies and objects/classes
library(xts)  #extension of time series objects and functionalities
library(forecast) #arima and other time series methods
library(lubridate) #date and time handling tools
library(colorRamps) #contains matlab.like color palette
library(rgeos) #spatial analysis, topological and geometric operations e.g. interesect, union, contain etc.
library(sphet) #spatial analyis, regression eg.contains spreg for gmm estimation

###### Functions used in this script

function_sbt_curve_generation <- "sbt_curve_generation_functions_07302017.R"
script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts"
source(file.path(script_path,function_sbt_curve_generation))

#####  Parameters and argument set up ###########


###### PART 1: generate dataset ##########

n_row <- 3
n_col <- 3

r <- raster(nrows=n_row,ncols=n_col)

vtn1 <- c(0.8,0.8,0.2,
         0.8,0.2,0.2,
         0.2,0.8,0.8)

vt0 <- c(0.8,0.8,0.2,
         0.8,0.2,0.2,
         0.2,0.8,0.8)

vt1 <- c(0.3,0.3,0.2,
         0.3,0.2,0.3,
         0.2,0.3,0.3)

vt2 <- c(0.4,0.4,0.2,
         0.4,0.2,0.4,
         0.2,0.4,0.4)

vt3 <- c(0.6,0.6,0.2,
         0.6,0.2,0.6,
         0.2,0.6,0.6)

vt4 <- c(0.6,0.6,0.2,
         0.6,0.2,0.6,
         0.2,0.6,0.6)

list_vt <- list(vtn1,vt0,vt1,vt2,vt3,vt4)

list_r <-(lapply(list_vt, function(x,rast){rast[]<-x;return(rast)},rast=r))
r_var <- stack(list_r)

n_layers <- nlayers(r_var)
r_var_names <- paste0("time_",1:n_layers)
names(r_var) <- r_var_names

## generating additional variables
r_zonal <- subset(r_var,1)
r_zonal[] <- rep(1,ncell(r_zonal))
r_x <-init(r_var,v="x")
r_y <-init(r_var,v="y")

col_palette <- matlab.like(100)
#col_palette <- matlab.like(100)

levelplot(r_var,col.regions=col_palette,main="Time series")

#Call function for sbt?
writeRaster(r_var, names(r_var), bylayer=TRUE, format='GTiff')

########## Generate figures and metrics ##########


################# END OF SCRIPT ##################

