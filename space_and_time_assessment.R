#######################################    Space Time Analyses Project   #######################################
############################ Assessing Space and Time Model for Space Beats Time framework #######################################
#This script contains functions to assess predictions for Space Beats Time Framework.       
#Temporal predictions use OLS with the image of the previous time step or ARIMA.
#Spatial predictions use spatial regression (lag error model) with different estimation methods (e.g. eigen, chebyshev etc.).
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/07/2017 
#DATE MODIFIED: 11/27/2017
#Version: 1

#PROJECT: Space beats time Framework
#TO DO:
#
#COMMIT: more changes and testing of code with tiles space and time predictions assessment
#

#################################################################################################

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

#Should use the data that is mosaiced!!

args<-commandArgs(TRUE)

args_table <- args[1]

#args_table <- "/home/bparmentier/Google Drive/Space_beats_time/Data/input_arguments_sbt_script_NDVI_Rita_10292017.csv"
#args_table <- "/home/parmentier/Data/Space_beats_time/Data/input_arguments_sbt_script_NDVI_Rita_11072017.csv"
args_table <- "/home/parmentier/Data/Space_beats_time/Data/input_arguments_sbt_assessment_script_NDVI_Rita_11152017.csv"
df_args <- read.table(args_table,sep=",",stringsAsFactors = FALSE)

### use column 2,3,4 etc.
#index_val <- 2 #this is set up for parallelization, if we have multiple regions/tiles, tile1
index_val <- 2 #this is set up for parallelization, if we have multiple regions/tiles, tile 2 rita

in_dir <- df_args[1,index_val]
out_dir <- df_args[2,index_val]
proj_str <- df_args[3,index_val]
file_format <- df_args[4,index_val]
NA_flag_val <- df_args[5,index_val]
out_suffix <- df_args[6,index_val]
create_out_dir_param <- df_args[7,index_val] 
data_fname <- df_args[8,index_val] 
coord_names <- df_args[9,index_val]  
zonal_colnames <- df_args[10,index_val] 
var_names <- df_args[11,index_val] 
num_cores <- df_args[12,index_val] 

n_time_event <- df_args[13,index_val]
time_window_selected <- df_args[14,index_val] 
#previous_step <- df_args[15,index_val] 
date_range <- df_args[15,index_val]  
r_ref <- df_args[16,index_val]  
r_temp_pred <- df_args[17,index_val]
r_spat_pred <- df_args[18,index_val]
method_space <- df_args[19,index_val] 
method_time <- df_args[20,index_val] 
pixel_index <- df_args[21,index_val] 

#P1
in_dir <- "/home/parmentier/Data/Space_beats_time/outputs/output_tile_2_NDVI_Rita_11062017"
r_temp_pred <- list.files(path=in_dir,
                          pattern="r_temp_pred_arima_arima_.*._tile_2_NDVI_Rita_11062017.tif",
                          full.names=T)
out_file <- file.path(in_dir,"raster_temp_files_list_tile_2.txt")
write.table(r_temp_pred,out_file)
#r_spat_pred <- list.files(path=in_dir,
#                                         pattern="r_spat_.*._tile_2_NDVI_Rita_11062017.tif",
#                                         full.names=T)
r_spat_pred <- list.files(path=in_dir,
                                         pattern="r_spat_pred_mle_eigen_no_previous_step_.*._tile_2_NDVI_Rita_11062017.tif",
                                         full.names=T)
out_file <- "/home/parmentier/Data/Space_beats_time/outputs/output_tile_2_NDVI_Rita_11062017/raster_spat_files_list_tile_2.txt"
#/home/parmentier/Data/Space_beats_time/outputs/output_tile_2_NDVI_Rita_11062017
write.table(r_spat_pred,out_file)
#r_spat_pred_mle_eigen_no_previous_step__t_113_tile_2_NDVI_Rita_11062017.tif
s_raster <- "/home/parmentier/Data/Space_beats_time/Data/data_Rita_NDVI/rev_project_output/tile_2/raster_files_list_tile_2.txt"
date_range <- "2001.01.01;2010.12.31;16"
time_window_predicted <- "105;114"
r_zonal <- "crop_r_zonal_rev_tile_2"
r_ref <- NULL
method_space <- "mle;eigen" #method for space and time used to predict space and time respectively
method_time <- "arima;arima;TRUE"
out_suffix <- "assessment_tile_2_NDVI_Rita_11062017"
out_dir <- "output_tile_1_2_combined_NDVI_Rita_11062017"
create_out_dir_param <- FALSE  

###### Functions used in this script

debug(accuracy_space_time_calc)

test <- accuracy_space_time_calc(r_temp_pred=r_temp_pred,
                                 r_spat_pred=r_spat_pred,
                                 s_raster=data_fname,
                                 time_window_predicted=time_window_predicted,
                                 r_zonal=zonal_colnames,
                                 method_space=method_space,
                                 method_time=method_time,
                                 r_ref=r_ref,
                                 out_suffix=out_suffix,
                                 var_names=var_names,
                                 NA_flag_val=NA_flag_val,
                                 file_format=file_format, 
                                 date_range=date_range,
                                 out_dir=out_dir,
                                 create_out_dir_param=create_out_dir_param)
  

#################### END OF SCRIPT ################################################################