####################################    General utility functions   #######################################
############################          Processing data tile by tile        #######################################
# This script contains general function to process spatial data tile by tile
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 10/31/2017 
#DATE MODIFIED: 11/02/2017
#Version: 
#PROJECT: General

#TO DO:
#
#COMMIT: 
#
#################################################################################################

#This script currently contains 7 functions:

###Loading R library and packages                                                      

library(sp)
library(rgdal)
library(spdep)
library(gtools)
library(maptools)
library(parallel)
library(rasterVis)
library(raster)
library(forecast) #ARIMA forecasting
library(xts)
library(zoo)
library(lubridate)
library(colorRamps) #contains matlab.like color palette
library(rgeos)
library(sphet) #contains spreg
library(BMS) #contains hex2bin and bin2hex
library(bitops)
library(sf)

###### Functions used in this script

cropping_raster_stack <- function(i,tile_sp,r_stack,file_format,out_dir,out_suffix_s){
  ## This functions crops a stack of raster using a set of polygons.
  ## Polygons are mostly likely representing processing tiles
  tile_poly_sp <- tile_sp[i,]
  
  out_dir_tile <- paste("tile_",i,sep="")
  
  if(!dir.exists(file.path(out_dir,out_dir_tile))){
    dir.create(file.path(out_dir,out_dir_tile))
  }
  #out_dir_tile <- dir.create(paste("tile_",i,sep=""))
  r_s_crop <- crop(r_stack,
                   tile_poly_sp,
                   filename=file.path(out_dir,out_dir_tile,paste0("crop",file_format)),
                   bylayer=T,
                   suffix=paste0(names(r_stack),"_",paste("tile_",i,sep=""),out_suffix_s),
                   overwrite=T)
  
  lf <- list.files(path=file.path(out_dir,out_dir_tile),
                   pattern=paste0(file_format,"$"),
                   full.names=T)
  out_list_filename <- file.path(out_dir,out_dir_tile,
                                 paste0("raster_files_list_",out_dir_tile,".txt"))
  write.table(lf,out_list_filename,row.names=F,sep=",")
  
  return(out_list_filename)
}


##
#To make overlapping grid, buffer teh region outline by % of area 
##and then generate grid.
#Crop/clip from the original using intersect?

#Or select each new grid feature,
#buffer
#dissolve
#reassemle the features
generate_outline_and_grid <- function(reg_layer,n_grid=NULL,out_suffix="",out_dir="."){
  #This function generates a grid from an input layer
  if(class(reg_layer)=="RasterLayer"){
    r <- reg_layer
    r[] <- 1 #change this later
    r <- mask(r,reg_layer)
    reg_sp <- rasterToPolygons(r,dissolve=T)
    reg_layer <- as(reg_sp, "sf") #this makes a class of "sfc_polygon" "sfc"
  }
  if(class(reg_layer)=="SpatialPolygons"){
    reg_layer <- as(reg_layer, "sf")
  }
  ## generate grid: default is 2x2
  if(is.null(n_grid)){
    n_grid <- c(2,2)
  }
  #reg_grid <- st_make_grid(reg_layer,n=n_grid,offset = n_distance)
  reg_grid <- st_make_grid(reg_layer,n=n_grid)
  
  reg_sf <- st_sf(tile_id=1:length(reg_grid),reg_grid)
  
  ## save grid figures
  plot(reg_layer)
  plot(reg_grid,border="red",add=T)
  
  ### buffer: make overlapping area later!!!
  #n_distance <- 5000
  #reg_buffer <- st_buffer(reg_layer,dist=n_distance)
  
  out_filename <- file.path(out_dir,
                            paste0("reg_sf_",out_suffix,".shp"))
  ## overwrite if existing
  st_write(reg_sf,out_filename,delete_dsn =T)
  
  return(reg_sf)
}

#############################################

in_dir <- "/home/parmentier/Data/Space_beats_time/Data/data_Rita_NDVI/rev_project_output"
out_dir <- in_dir

lf <- list.files(path=in_dir,
                 pattern=".rst$",
                 full.names=T)
r_stack <- stack(lf) #input stack
r_val <- subset(r_stack,1)

file_format <- ".tif"
 
#debug(generate_outline_and_grid)
tile_sf <- generate_outline_and_grid(r_val,n_grid=c(2,1))

#tile_sf <- r_test 
tile_sp <- as(tile_sf,"Spatial")

num_cores <- 2
out_suffix <- ""

debug(cropping_raster_stack)
out_file_test <- cropping_raster_stack(2,
                     tile_sp = tile_sp,
                     r_stack = r_stack,
                     file_format =file_format,
                     out_dir = out_dir,
                     out_suffix_s = out_suffix)


out_file <- mclapply(1:nrow(tile_sp),
         FUN=cropping_raster_stack,
         tile_sp = tile_sp,
         r_stack = r_stack,
         file_format = ".tif",
         out_dir = out_dir,
         out_suffix_s = out_suffix,
         mc.cores = num_cores,
         mc.preschedule = FALSE)


############################# END OF SCRIPT ###########################
    