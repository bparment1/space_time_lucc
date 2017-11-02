####################################    Space Time Analyses Project   #######################################
############################  Yucatan case study: SAR, SARMA etc -PART 2  #######################################
#This script functions to produce predictions for the dates following the Hurricane Dean event.       
#The script uses spatial regression with weight matrix to predict NDVI values in the MOORE EDGY region. 
#Temporal predictions use OLS with the image of the previous time step rather than ARIMA.
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/09/2014 
#DATE MODIFIED: 08/08/2017
#Version: 2
#PROJECT: GLP Conference Berlin,YUCATAN CASE STUDY with Marco Millones            
#PROJECT: Workshop for William and Mary: an intro to spatial regression with R 
#PROJECT: Geocomputation and AAG 2015
#PROJECT: Space beats time paper

#TO DO:
# modify the rasterize_df_fun function to allow ref image
# add the ARIMA method to run more efficiently
#
#COMMIT: checking arima code
#
#################################################################################################

#This script currently contains 7 functions:

#[1] "calc_ac_stat_fun" : compute MAE for predictions wiht or without regions/zones
#[2] "create_dir_fun"   : create an output directory                    
#[3] "create_sp_poly_spatial_reg" : create a list of weights and polygon file for spatial regression
#[4] "create_uniform_mask_fun" : harmonize NA values given input layers with different valid values            
#[5] "load_obj" : load R object                          
#[6] "predict_spat_reg_fun" : function to perform spatial regresssion prediction
#[7] "predict_temp_reg_fun : function to preform temporal  predictions       
#[8] "rasterize_df_fun" : create raster from textfile with location information
#
#Add arima functions
#[9] "convert_arima_pred_to_raster"         
#[10] "extract_arima_mod_info"              
#[11] "pixel_ts_arima_predict"              
#[12] "raster_NA_image"                     
#[13] "raster_ts_arima"                      
#[14] "raster_ts_arima_predict"                         

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

###### Functions used in this script
#debug(generate_outline_and_grid)
r_val <- subset(s_raster,1)
r_test <- generate_outline_and_grid(r_val)
debug(generate_outline_and_grid)
r_test <- generate_outline_and_grid(r_val,n_grid=c(2,1))
tile_sf <- subset(r_test,1)
tile_sf <- r_test[1,]

tile_sf <- r_test 
tile_sp <- as(tile_sf,"Spatial")

num_cores <- 2
out_file <- mclapply(1:nrow(tile_sp),
         FUN=cropping_raster_stack,
         tile_sp = tile_sp,
         r_stack = s_raster,
         file_format = ".tif",
         out_dir = out_dir,
         out_suffix_s = out_suffix,
         mc.cores = num_cores,
         mc.preschedule = FALSE)

cropping_raster_stack <- function(i,tile_sp,r_stack,file_format,out_dir,out_suffix_s){
  ## This functions crops a stack of raster using a set of polygons.
  ## Polygons are mostly likely representing processing tiles
  tile_poly_sp <- tile_sp[i,]
  
  out_dir_tile <- paste("tile_",i,sep="")
  
  if(!dir.exists(out_dir_tile)){
    dir.create(out_dir_tile)
  }
  #out_dir_tile <- dir.create(paste("tile_",i,sep=""))
  r_s_crop <- crop(r_stack,
                   tile_poly_sp,
                   filename=file.path(out_dir_tile,paste0("crop",file_format)),
                   bylayer=T,
                   suffix=paste(names(r_stack),out_suffix_s),
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


    