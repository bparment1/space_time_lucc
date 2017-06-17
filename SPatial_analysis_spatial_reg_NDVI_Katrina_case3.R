####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: SAR, SARMA etc -PART 2  #######################################
#This script produces predictions for the Space Beats Time Framework.
#Data used in the script are related to hurricane Katrina.
#The script uses spatial neighbours to predict Light data values in the Hurricane Katrina New Orleans region. 
#Temporal predictions use OLS with the image of the previous time or the ARIMA method.
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/09/2014 
#DATE MODIFIED: 06/17/2017
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
#COMMIT: testing and debugging function run_space_and_time_models
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

function_space_and_time_predictions <- "space_and_time_predictions_functions_06172017c.R"
function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_functions_06072017.R" #PARAM 1
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
in_dir <- "/home/bparmentier/Google Drive/Space_beats_time/Case2_data_NDVI/"

#out_dir <- "/home/parmentier/Data/Space_beats_time/outputs"
out_dir <- "/home/bparmentier/Google Drive/Space_beats_time/outputs"
proj_str<- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0"  

file_format <- ".tif" #PARAM5
NA_flag_val <- -9999 #PARAM7
out_suffix <-"NDVI_Katrina_06172017" #output suffix for the files and output folder #PARAM 8
create_out_dir_param=TRUE #PARAM9

#data_fname <- file.path("/home/parmentier/Data/Space_beats_time/R_Workshop_April2014","Katrina_Output_CSV - Katrina_pop.csv")
#data_fname <- file.path(in_dir,"lights/table","Kat_lights.txt") #PARAM 10
#data_fname <- file.path(in_dir,"output_Katrina_04082015","dat_reg_var_list_NDVI_Katrina_04082015.txt")
data_fname <- file.path(in_dir,"dat_reg2_var_list_NDVI_NDVI_Katrina_04102015.txt")

coord_names <- c("x","y") #PARAM 11
#coord_names <- c("Long","Lat") #PARAM 11
#coord_names <- c("XCoord","YCoord")
#coord_names <- c("POINT_X1","POINT_Y1")

zonal_colnames <- "r_srtm_Katrina_rec2" #PARAM 12

var_names <- 1:230 #PARAM 13 #Data is stored in the columns 3 to 22
#num_cores <- 11 #PARAM 14
num_cores <- 4 #PARAM 14

n_time_event <- 108 #PARAM 15 #this is the timestep corresponding to the event ie Hurricane Katrina (Aug 23- Aub 31 2005): 235-243 DOY, storm surge Aug 29 in New Orelans
time_window_selected <- var_names #PARAM 16: use alll dates for now
time_window_selected <- 100:116 #PARAM 16: use alll dates for now

re_initialize_arima <- T #PARAM 17, use re-initialization ie apply arima model with one step forward at each time step
previous_step <- T #PARAM 18

#date_range1 <- c("2001.01.01","2012.12.31") #EDGY DEAN
#date_range2 <- c("1992.01.01","2013.12.31") #Light Katrina: annual
date_range3 <- c("2001.01.01","2010.12.31") #NDVI Katrina

agg_fact <- 5
agg_fun <- "mean" 

################# START SCRIPT ###############################

### PART I READ AND PREPARE DATA FOR REGRESSIONS #######

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

  lf_agg <- mclapply(l_rast,
                     FUN=aggregate_raster,
                     #r_in=raster(lf_layerized_bool[1]),
                     agg_fact=agg_fact,
                     reg_ref_rast=NULL,
                     #agg_fun="mean",
                     agg_fun=agg_fun,
                     out_suffix=NULL,
                     file_format=file_format,
                     out_dir=out_dir,
                     out_rast_name = NULL,
                     mc.preschedule=FALSE,
                     mc.cores = num_cores) 
  
  l_rast_original <- l_rast
  l_rast <- unlist(lf_agg) 
  
  ###Break out and get mean per class and do majority rule!
  
  #l_rast_original
  #r_r_srtm_Katrina_rec2_NDVI_Katrina_03162017.rst"
  #r <- raster(paste0("r_",zonal_colnames,"_",out_suffix,file_format,sep=""))
  raster_name <- (paste0("r_",zonal_colnames,"_",out_suffix,file_format,sep=""))
  out_suffix_str <- paste0("agg5_zonal","_",out_suffix)
  #debug(generate_soft_cat_aggregated_raster_fun)
  lf_agg_soft <- generate_soft_cat_aggregated_raster_fun(raster_name,
                                                         reg_ref_rast=NULL,
                                                         agg_fact,
                                                         agg_fun,
                                                         num_cores,
                                                         NA_flag_val=NA_flag_val,
                                                         file_format,
                                                         out_dir,
                                                         out_suffix_str)

  reclass_val <- unique(raster(raster_name)) #unique zonal values to reassign
  #reclass_val <- c(0,1,2) # value for the elevation reclassified
  
  #function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_functions_04072017b.R" #PARAM 1
  #script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts"
  #source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.
  #debug(reclass_in_majority)

  r_reclass_obj <- reclass_in_majority(r_stack=stack(lf_agg_soft),
                                       threshold_val=NULL,
                                       max_aggregation = TRUE,
                                       reclass_val = reclass_val)
  
  plot(r_reclass_obj$r_rec)
  rast_zonal <- r_reclass_obj$r_rec
  #zonal_colnames
  raster_name <- paste0("agg_",agg_fact,"_","r_",zonal_colnames,"_",out_suffix,file_format)
  
  writeRaster(rast_zonal,
               filename=file.path(out_dir,raster_name),
               overwrite=TRUE)  
  
  #r_srtm_Katrina_rec2
  #-rw-rw-r-- 1 bparmentier bparmentier 1894 Apr  7 12:33 r_r_srtm_Katrina_rec2_NDVI_Katrina_04062017.tif
  #-rw-rw-r-- 1 bparmentier bparmentier 1016 Apr  7 12:34 agg_5_r_r_srtm_Katrina_rec2_NDVI_Katrina_04062017.tif
  
  ###
  zonal_colnames <- gsub(extension(raster_name),"",raster_name)
  ##
  
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

num_cores <- 4
previous_step <- T
method_space <- c("mle","eigen") #estimator <- "mle",estimation_method <- "eigen"
re_initialize_arima <- T
method_time <- c("arima","arima",re_initialize_arima) # estimator <- "arima",estimation_method <-"arima"

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