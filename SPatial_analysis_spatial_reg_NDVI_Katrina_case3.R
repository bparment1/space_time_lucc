####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: SAR, SARMA etc -PART 2  #######################################
#This script produces predictions for the Space Beats Time Framework.
#Data used in the script are related to hurricane Katrina.
#The script uses spatial neighbours to predict Light data values in the Hurricane Katrina New Orleans region. 
#Temporal predictions use OLS with the image of the previous time or the ARIMA method.
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/09/2014 
#DATE MODIFIED: 07/29/2017
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
#COMMIT:  moving aggreate raster function for SBT spatial resolution test
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

function_space_and_time_predictions <- "space_and_time_predictions_functions_07282017.R"
function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_functions_07282017.R" #PARAM 1
function_paper_figures_analyses <- "space_beats_time_sbt_paper_figures_functions_01092016.R" #PARAM 1
function_data_figures_reporting <- "spatial_analysis_data_figures_reporting_functions_06172017.R" #PARAM 1
#script_path <- "/home/parmentier/Data/Space_beats_time/sbt_scripts" #path to script #PARAM 2
script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts"
source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.
source(file.path(script_path,function_paper_figures_analyses)) #source all functions used in this script 1.
source(file.path(script_path,function_data_figures_reporting)) #source all functions used in this script 1.
source(file.path(script_path,function_space_and_time_predictions))

#Aggregation code
function_multilabel_fuzzy_analyses <- "classification_multilabel_processing_functions_03142017.R" #PARAM 1
script_path <- "/home/bparmentier/Google Drive/LISER_Lux/R_scripts" #path to script #PARAM 2
source(file.path(script_path,function_multilabel_fuzzy_analyses)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

#in_dir <- "~/Data/Space_beats_time/case3data/" #lights/table" #PARAM3
#in_dir <- "/home/parmentier/Data/Space_beats_time/Case2_data_NDVI/"
#in_dir <- "~/Data/Space_beats_time/case3data/lights/table"
#in_dir <- "~/Data/Space_beats_time/Case1a_data"
in_dir <- "/home/bparmentier/Google Drive/Space_beats_time/Case2_data_NDVI/" #PARAM 1

#out_dir <- "/home/parmentier/Data/Space_beats_time/outputs"
out_dir <- "/home/bparmentier/Google Drive/Space_beats_time/outputs" #PARAM 2
proj_str<- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0"  #PARAM 3

## Constant
file_format <- ".tif" #PARAM5 #PARAM 4
NA_flag_val <- -9999 #PARAM7 #PARAM5

out_suffix <-"NDVI_Katrina_07282017" # PARAM6, output suffix for the files and output folder 
create_out_dir_param=TRUE #PARAM7

#data_fname <- file.path("/home/parmentier/Data/Space_beats_time/R_Workshop_April2014","Katrina_Output_CSV - Katrina_pop.csv")
#data_fname <- file.path(in_dir,"lights/table","Kat_lights.txt") #PARAM 10
#data_fname <- file.path(in_dir,"output_Katrina_04082015","dat_reg_var_list_NDVI_Katrina_04082015.txt")
data_fname <- file.path(in_dir,"dat_reg2_var_list_NDVI_NDVI_Katrina_04102015.txt") #PARAM 8

coord_names <- c("x","y") #PARAM 9
#coord_names <- c("Long","Lat") 
#coord_names <- c("XCoord","YCoord")
#coord_names <- c("POINT_X1","POINT_Y1")

zonal_colnames <- "r_srtm_Katrina_rec2" #PARAM 12

var_names <- 1:230 #PARAM 10 #Data is stored in the columns 3 to 22
#num_cores <- 11 #PARAM 
num_cores <- 4 #PARAM 11

n_time_event <- 108 #PARAM 12 #this is the timestep corresponding to the event ie Hurricane Katrina (Aug 23- Aub 31 2005): 235-243 DOY, storm surge Aug 29 in New Orelans
time_window_selected <- 100:116 #PARAM 13: use alll dates for now

previous_step <- T #PARAM 14

#date_range1 <- c("2001.01.01","2012.12.31") #EDGY DEAN
#date_range2 <- c("1992.01.01","2013.12.31") #Light Katrina: annual
date_range3 <- c("2001.01.01","2010.12.31") #PARAM 15, NDVI Katrina

agg_fact <- 5 #PARAM 16
agg_fun <- "mean" #PARAM 17
use_majority <- T #PARAM 18

## Constant?
method_space <- c("mle","eigen") #PARAM 19, estimator <- "mle",estimation_method <- "eigen"
re_initialize_arima <- T #PARAM 20
method_time <- c("arima","arima",re_initialize_arima) #PARAM 21, estimator <- "arima",estimation_method <-"arima"

#pixel index
pixel_index <- 800 #PARAM 22

################# START SCRIPT ###############################

### PART I READ AND PREPARE DATA FOR REGRESSIONS #######

#set up the working directory
#Create output directory
create_sp_poly_spatial_reg
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

data_tb <-read.table(data_fname,sep=",",header=T)
#dates1 <- generate_dates_by_step(date_range1[1],date_range1[2],16)$dates
#dates2 <- unique(year(generate_dates_by_step(date_range2[1],date_range2[2],1)$dates)) #extract year
dates3 <- generate_dates_by_step(date_range3[1],date_range3[2],16)$dates #NDVI Katrina

#Transform table text file into a raster image

###################################
### Aggregate data if necessary
#debug(rasterize_df_fun)
#This function is very slow and inefficienct, needs improvement (add parallelization)
l_rast <- rasterize_df_fun(data_tb,coord_names,proj_str,out_suffix,out_dir=".",file_format=file_format,NA_flag_val,tolerance_val=0.000120005)

if(!is.null(agg_fact)){
  #debug(aggregate_raster_fun)
  obj <- aggregate_raster_fun(l_rast,zonal_colnames,use_majority,agg_fact,agg_fun,file_format,rast_ref,num_cores,out_suffix, out_dir)

  l_rast <- obj$l_rast  
  zonal_colnames <- obj$zonal_colnames
  l_rast_original <- obj$l_rast_original
  
}

###########################
#### PART II: data reporting: description

### Generate data description and figures: make this a markdown output later on!
#debug(explore_and_summarize_data)
#test <- explore_and_summarize_data(l_rast,zonal_colnames, var_names,n_time_event)

test <- explore_and_summarize_data(l_rast,
                                   zonal_colnames, 
                                   var_names,
                                   n_time_event,
                                   proj_str =proj_str,
                                   pixel_index = 800,
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

#debug(run_space_and_time_models)

run_space_and_time_models(s_raster,
                          n_time_event,
                          time_window_selected,
                          method_space=c("mle","eigen"),
                          method_time=c("arima","arima",T), 
                          NA_flag_val=NA_flag_val,
                          file_format=".tif",
                          rast_ref=NULL,
                          zonal_colnames,
                          num_cores=num_cores,
                          out_dir=out_dir, 
                          out_suffix=out_suffix)


###########################  END OF SCRIPT #########################################