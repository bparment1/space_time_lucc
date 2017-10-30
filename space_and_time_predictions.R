####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: SAR, SARMA etc -PART 2  #######################################
#This script produces predictions for the Space Beats Time Framework.
#Data used in the script are related to hurricane or other event  of any nature.
#The script uses spatial neighbours to predict values in the reference region. 
#Temporal predictions use OLS with the image of the previous time or the ARIMA method.
# Event type: Rita from 09/18 to 09/26
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/09/2014 
#DATE MODIFIED: 10/30/2017
#Version: 3
#PROJECT: GLP Conference Berlin,YUCATAN CASE STUDY with Marco Millones            
#PROJECT: Workshop for William and Mary: an intro to geoprocessing with R 
#PROJECT: AAG 2015 in Chicago, with Marco Millones
#PROJECT: Geocomputation conference in Dallas with Marco Millones
#
#COMMENTS: - Testing alternative methods to eigen for spatial predictions: "Chebyshev" on new light data
#         - clean up and organize code to be more general for any dataset
#TO DO:
# - add confidence interval around reg coef: this is important!!
# - add variance around the MAE values in the accuracy assessment
# - modify parallelization so that it works both on windows and linux/macos
# - automation to call from the terminal/shell
#
#
#COMMIT: modifications to arima outputs on Katrina NDVI 
#

#Rscript space_and_time_predictions_07292017b.R "/home/bparmentier/Google Drive/Space_beats_time/Data/input_arguments_sbt_script_NDVI_Katrina_09292017.csv"

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
library(sf)

###### Functions used in this script

#function_space_and_time_predictions <- "space_and_time_predictions_functions_08112017.R"
function_space_and_time_predictions <- "space_and_time_predictions_functions_10302017.R"
function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_functions_10302017.R" #PARAM 1
function_paper_figures_analyses <- "space_beats_time_sbt_paper_figures_functions_01092016.R" #PARAM 1
function_data_figures_reporting <- "spatial_analysis_data_figures_reporting_functions_08042017.R" #PARAM 1
script_path <- "/home/parmentier/Data/Space_beats_time/sbt_scripts" #path to script #PARAM 2
#script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts"
source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.
source(file.path(script_path,function_paper_figures_analyses)) #source all functions used in this script 1.
source(file.path(script_path,function_data_figures_reporting)) #source all functions used in this script 1.
source(file.path(script_path,function_space_and_time_predictions))

#Aggregation code
function_multilabel_fuzzy_analyses <- "classification_multilabel_processing_functions_03142017.R" #PARAM 1
#script_path <- "/home/bparmentier/Google Drive/LISER_Lux/R_scripts" #path to script #PARAM 2
source(file.path(script_path,function_multilabel_fuzzy_analyses)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

args<-commandArgs(TRUE)

args_table <- args[1]

#args_table <- "/home/bparmentier/Google Drive/Space_beats_time/Data/input_arguments_sbt_script_NDVI_Rita_10292017.csv"
args_table <- "/home/parmentier/Data/Space_beats_time/Data/input_arguments_sbt_script_NDVI_Rita_10292017.csv"

df_args <- read.table(args_table,sep=",",stringsAsFactors = FALSE)

### use column 2,3,4 etc.
index_val <- 2 #this is set up for parallelization, if we have multiple regions

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
previous_step <- df_args[15,index_val] 
date_range <- df_args[16,index_val]  
agg_fact <- df_args[17,index_val]  
agg_fun <- df_args[18,index_val]
use_majority <- df_args[19,index_val]
method_space <- df_args[20,index_val] 
re_initialize_arima <- df_args[21,index_val] 
method_time <- df_args[22,index_val]
pixel_index <- df_args[23,index_val] 

#### input values for Katrina NDVI data
# in_dir <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_Rita_NDVI/" #PARAM 1
# out_dir <- "/home/bparmentier/Google Drive/Space_beats_time/outputs" #PARAM 2
# proj_str<- "+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"  #PARAM 3
# file_format <- ".tif" #PARAM5 #PARAM 4
# NA_flag_val <- "-9999" #PARAM7 #PARAM5
# out_suffix <- "NDVI_Rita_10222017" # PARAM6, output suffix for the files and output folder
# create_out_dir_param <- TRUE #PARAM7
# data_fname <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_Rita_NDVI/project_output/raster_data_list_NDVI_houston.txt" #PARAM 8
# coord_names <- "x,y" #PARAM 9
# zonal_colnames <- "r_mask_rita" #PARAM 12
# var_names <- "1,230" #PARAM 10 #Data is stored in the columns 3 to 22
# #num_cores <- "4" #PARAM 11
# n_time_event <- "110" #PARAM 12 #this is the timestep corresponding to the event ie Hurricane Katrina (Aug 23- Aub 31 2005): 235-243 DOY, storm surge Aug 29 in New Orelans
# time_window_selected <- "100,116" #PARAM 13: use alll dates for now
# previous_step <- TRUE #PARAM 14
# date_range <- "2001.01.01,2010.12.31,16" #PARAM 15, NDVI Katrina
# agg_fact <- NULL #PARAM 16 , if NULL no aggregation is performed
# agg_fun <- "mean" #PARAM 17
# use_majority <- TRUE #PARAM 18
# method_space <- "mle,eigen" #PARAM 19, estimator <- "mle",estimation_method <- "eigen"
# re_initialize_arima <- TRUE #PARAM 20
# method_time <- "arima,arima,TRUE" #PARAM 21, estimator <- "arima",estimation_method <-"arima"
# pixel_index <- "800" #PARAM 22

################# START SCRIPT ###############################

### PART I READ AND PREPARE DATA FOR REGRESSIONS #######

### Check input argument data types:
NA_flag_val <- as.integer(NA_flag_val)
coord_names <- unlist(strsplit(coord_names,","))
num_cores <- as.integer(num_cores)
var_names <- as.integer(unlist(strsplit(var_names,",")))
var_names <- seq(var_names[1],var_names[2])
n_time_event <- as.integer(n_time_event)
time_window_selected <-  as.integer(unlist(strsplit(time_window_selected,",")))
time_window_selected <- seq(time_window_selected[1],time_window_selected[2])

method_space <- (unlist(strsplit(method_space,",")))
method_time <- (unlist(strsplit(method_time,",")))

date_range <- (unlist(strsplit(date_range,",")))

### Handle command line data type
if(agg_fact=="NULL"){
  agg_fact <- NULL
}else{
  agg_fact <- as.integer(agg_fact)
}
if(out_dir=="NULL"){
  out_dir <- NULL
}
  
#set up the working directory
#Create output directory

if(is.null(out_dir)){
  out_dir <- dirname(in_dir) #output will be created in the input dir
}

out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix_s)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

## Add check to see if raster tif or list of files!!

data_tb <- try(read.table(data_fname,sep=",",header=T,stringsAsFactors = F))
#if(inherits(data_tb)=="try-error"){
  #<- stack(data_fname)
#}

#dates1 <- generate_dates_by_step(date_range1[1],date_range1[2],16)$dates
#dates2 <- unique(year(generate_dates_by_step(date_range2[1],date_range2[2],1)$dates)) #extract year
dates_val <- generate_dates_by_step(date_range[1],date_range[2],as.integer(date_range[3]))$dates #NDVI Katrina

#Transform table text file into a raster image

###################################
### Aggregate data if necessary
#debug(rasterize_df_fun)
#This function is very slow and inefficienct, needs improvement (add parallelization)

### Read in dataset
### If ncol is greater than 1 then, it is a textfile with data rather 
#than a list of raster
if(ncol(data_tb)>1){
  l_rast <- rasterize_df_fun(data_tb,
                             coord_names,
                             proj_str,
                             out_suffix,
                             out_dir=out_dir,
                             file_format=file_format,
                             NA_flag_val,
                             tolerance_val=0.000120005)
  zonal_colnames <- paste0("r_",zonal_colnames,"_",out_suffix)
  
}else{
  l_rast <- data_tb[,1]
}

##### Aggregate data if not NULL

if(!is.null(agg_fact)){
  #debug(aggregate_raster_fun)
  obj <- aggregate_raster_fun(l_rast,
                              zonal_colnames,
                              use_majority,
                              agg_fact,
                              agg_fun,
                              file_format,
                              rast_ref,
                              num_cores,
                              out_suffix,
                              out_dir)

  l_rast <- obj$l_rast  
  zonal_colnames <- obj$zonal_colnames
  l_rast_original <- obj$l_rast_original
}

###########################
#### PART II: data reporting: description

### Generate data description and figures: make this a markdown output later on!
#debug(explore_and_summarize_data)
#test <- explore_and_summarize_data(l_rast,zonal_colnames, var_names,n_time_event)
#debug(explore_and_summarize_data)
explore_obj <- explore_and_summarize_data(l_rast,
                                   zonal_colnames, 
                                   var_names,
                                   n_time_event,
                                   proj_str =proj_str,
                                   pixel_index = pixel_index,
                                   out_dir = out_dir,
                                   out_suffix=out_suffix)

#l_rast,zonal_colnames,var_names,n_time_event,pixel_index=800,out_dir=NULL,out_suffix=NULL
s_raster <- stack(l_rast)
#names(s_raster) <- names(data_tb)

###########################
#### PART III: run space and time model

#num_cores <- 4
#previous_step <- T
#method_space <- c("mle","eigen") #estimator <- "mle",estimation_method <- "eigen"
#re_initialize_arima <- T
#method_time <- c("arima","arima",re_initialize_arima) # estimator <- "arima",estimation_method <-"arima"
#method_time <- c(estimator,estimation_method,F)
#method_time <- c("lm","ols",FALSE)

#debug(run_space_and_time_models)
function_space_and_time_predictions <- "space_and_time_predictions_functions_10302017.R"
script_path <- "/home/parmentier/Data/Space_beats_time/sbt_scripts"
source(file.path(script_path,function_space_and_time_predictions))

run_space_and_time_models(s_raster,
                          n_time_event,
                          time_window_selected,
                          method_space=method_space,
                          method_time=method_time, 
                          NA_flag_val=NA_flag_val,
                          file_format=file_format,
                          rast_ref=NULL,
                          zonal_colnames,
                          num_cores=num_cores,
                          out_dir=out_dir, 
                          out_suffix=out_suffix)


###########################  END OF SCRIPT #########################################