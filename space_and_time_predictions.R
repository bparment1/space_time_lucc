####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: SAR, SARMA etc -PART 2  #######################################
#This script produces predictions for the Space Beats Time Framework.
#Data used in the script are related to hurricane or other event  of any nature.
#The script uses spatial neighbours to predict values in the reference region. 
#Temporal predictions use OLS with the image of the previous time or the ARIMA method.
# Event type: Hurricanes, Urbanization etc.
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/09/2014 
#DATE MODIFIED: 06/07/2018
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
#COMMIT: setting up inputs for Lagos
#

#### RUNNNG script
#Script can be called from the shell using "Rscript" command:
#Rscript space_and_time_predictions_07292017b.R "/home/bparmentier/Google Drive/Space_beats_time/Data/input_arguments_sbt_script_NDVI_Katrina_09292017.csv"

# The input argument is a csv file with a minimum of two columns:
# - the first column describes the name of input parameters
# - the second column contains values of the input parameters
# - additional column may be added for additionals region to run in parallel

## The content of the csv input file parameter is described below:
#### There are 23 input parameters:
#in_dir,"/media/dan/Space_beats_time/Space_beats_time/Case2_data_NDVI/"
#out_dir,"/media/dan/Space_beats_time/Space_beats_time/outputs"
#proj_str,"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0"
#file_format,".tif"
#NA_flag_val,-9999
#out_suffix,"NDVI_Katrina_05252018"
#create_out_dir_param,"TRUE"
#data_fname,"/media/dan/Space_beats_time/Space_beats_time/Case2_data_NDVI/dat_reg2_var_list_NDVI_NDVI_Katrina_04102015.txt"
#coord_names,"x;y"
#zonal_colnames,"r_srtm_Katrina_rec2"
#var_names,"1;230"
#num_cores,"4"
#n_time_event,"108"
#time_window_selected,"100;116"
#previous_step,"TRUE" 
#date_range,"2001.01.01;2010.12.31;16"
#agg_fact,"5"
#agg_fun,"mean"
#use_majority,"TRUE"
#method_space,"mle;eigen"
#re_initialize_arima,"TRUE" 
#method_time,"arima;arima;TRUE"
#pixel_index,"800"

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
library(snow)

###### Functions used in this script

## space beats time predictions run on specific dataset
function_space_and_time_predictions <- "space_and_time_predictions_functions_06072018.R"
function_space_and_time_assessment <- "space_and_time_assessment_functions_06082018.R"
function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_functions_06072018.R" #PARAM 1
function_paper_figures_analyses <- "space_beats_time_sbt_paper_figures_functions_06042018.R" #PARAM 1
function_data_figures_reporting <- "spatial_analysis_data_figures_reporting_functions_06042018.R" #PARAM 1
#Aggregation code
function_aggregation_raster <- "aggregation_raster_functions_06032018.R" #PARAM 1

script_path <- "/media/dan/Space_beats_time/Space_beats_time/sbt_scripts" #path to script #PARAM 2
#script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts"
source(file.path(script_path,function_space_and_time_predictions))
source(file.path(script_path,function_space_and_time_assessment))
source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.
source(file.path(script_path,function_paper_figures_analyses)) #source all functions used in this script 1.
source(file.path(script_path,function_data_figures_reporting)) #source all functions used in this script 1.
source(file.path(script_path,function_aggregation_raster)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

args<-commandArgs(TRUE)

args_table <- args[1]

###Comment this out if run from shell script
#args_table <- "/home/bparmentier/Google Drive/Space_beats_time/Data/input_arguments_sbt_script_NDVI_Rita_10292017.csv"
#args_table <- "/media/dan/Space_beats_time/Space_beats_time/Data/input_arguments_sbt_script_REACT_Lagos_NDVI_mod13_05242018.csv"
args_table <- "/media/dan/Space_beats_time/Space_beats_time/Data/input_arguments_sbt_script_NDVI_Katrina_06092018.csv"
df_args <- read.table(args_table,sep=",",stringsAsFactors = FALSE)

### use column 2,3,4 etc.
#index_val <- 2 #this is set up for parallelization, if we have multiple regions/tiles, tile1
index_val <- 2 #this is set up for parallelization, if we have multiple regions/tiles, tile 2 rita

#### There are 23 input parameters:
in_dir <- df_args[1,index_val] #input directory, path to data
out_dir <- df_args[2,index_val] #output directory for analyses
proj_str <- df_args[3,index_val] #projection for the region 
file_format <- df_args[4,index_val] #image/raster format e.g. tif
NA_flag_val <- df_args[5,index_val] #flag value for NA, e.g. -9999
out_suffix <- df_args[6,index_val] # 
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
rast_ref <- df_args[24,index_val] 

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
coord_names <- unlist(strsplit(coord_names,";"))
num_cores <- as.integer(num_cores)
var_names <- as.integer(unlist(strsplit(var_names,";")))
var_names <- seq(var_names[1],var_names[2])
n_time_event <- as.integer(n_time_event)
time_window_selected <-  as.integer(unlist(strsplit(time_window_selected,";")))
time_window_selected <- seq(time_window_selected[1],time_window_selected[2])

method_space <- (unlist(strsplit(method_space,";")))
method_time <- (unlist(strsplit(method_time,";")))

date_range <- (unlist(strsplit(date_range,";")))

### Handle command line data type
if(agg_fact=="NULL"){
  agg_fact <- NULL
}else{
  agg_fact <- as.integer(agg_fact)
}
if(out_dir=="NULL"){
  out_dir <- NULL
}
 
if(zonal_colnames=="NULL"){
  zonal_colnames <- NULL
}
### reference image
if(rast_ref=="NULL"){
  rast_ref <- NULL
}

### reference image
if(use_majority=="TRUE"){
  agg_method_cat <- "majority"
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

dates_val <- generate_dates_by_step(date_range[1],date_range[2],date_range[3])$dates #NDVI Katrina

## Add check to see if raster tif or list of files!!
## Assumes that it is text file of csv type
input_data_fname <- file.path(in_dir,
                              data_fname)
#sep=",",header=T,stringsAsFactors = F)
data_tb <- try(read.table(input_data_fname,sep=",",header=T,stringsAsFactors = F))
## if this is not a text file: check for stack or 
if(inherits(data_tb,"try-error")){
  ###  Load data if needed:
  data_tb <-brick(file.path(in_dir,data_fname))
  #dim(data_tb)
  #generate list of files and write out separate files to run predictions
  raster_name_tmp <- basename(filename(data_tb))
  bylayer_val=TRUE
  suffix_str <- paste0(1:length(dates_val),"_",gsub("-","_",dates_val))
  writeRaster(data_tb,
              filename=file.path(out_dir,raster_name_tmp),
              bylayer=bylayer_val,
              #suffix=paste(names(r),"_",out_suffix,sep=""),
              #format=format_raster,
              suffix=suffix_str,
              overwrite=TRUE,
              NAflag=NA_flag_val,
              #datatype=data_type_str,
              options=c("COMPRESS=LZW"))
  pattern_str <- gsub(extension(raster_name_tmp),"",raster_name_tmp)
  pattern_raster_str <- paste0(pattern_str,".*",file_format) 
  lf <- mixedsort(list.files(pattern=pattern_raster_str,path=out_dir))
  write.table(lf,file=paste0("list_raster_files_",pattern_str,".txt"))
  data_tb <- read.table(paste0("list_raster_files_",pattern_str,".txt"),
                        stringsAsFactors = F)
}

#Transform table text file into a raster image

###################################
### Aggregate data if necessary
#debug(rasterize_df_fun)
#This function is very slow and inefficienct, needs improvement (add parallelization)

### Read in dataset
### If ncol is greater than 1 then, it is a textfile with data rather 
#than a list of raster files (e.g. tif) or brick
if(ncol(data_tb)>1 && class(data_tb)!="raster"){
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

#### set zonal stat if NULL
if(is.null(zonal_colnames)){
  #this does not load everything in memory:
  r_stack <- stack(l_rast[var_names])

  beginCluster(num_cores)
  r_mean <- clusterR(r_stack, calc, args=list(mean, na.rm=T))
  endCluster()
  
  raster_name <- file.path(out_dir,"r_mean.tif")
  writeRaster(r_mean,file=raster_name,overwrite=T)

  r_zonal <- !is.na(r_mean)
  plot(r_zonal)
  
  #set1f <- function(x){rep(1, x)}
  zonal_colnames <- "r_zonal"
  raster_name <- file.path(out_dir,"r_zonal.tif")
  #r_zonal <- init(r_FID, fun=set1f, filename=raster_name, overwrite=TRUE)
  #s_raster <- addLayer(s_raster,r_zonal)
  writeRaster(r_zonal,file=raster_name,overwrite=T)
  
  l_rast <- c(l_rast,raster_name)
}

if(is.null(rast_ref)){
  ## if null then use mean as reference image
  raster_name_ref <- file.path(out_dir,"r_mean.tif")
  if(!file.exists(raster_name_ref)){
    #this does not load everything in memory:
    r_stack <- stack(l_rast[var_names])
    
    beginCluster(num_cores)
    r_mean <- clusterR(r_stack, calc, args=list(mean, na.rm=T))
    endCluster()
    
    #raster_name <- file.path(out_dir,"r_mean.tif")
    writeRaster(r_mean,file=raster_name_ref,overwrite=T)
    rast_ref <- raster_name_ref
  }
  l_rast <- c(l_rast,raster_name_ref)
}

##### Aggregate data if not NULL

if(!is.null(agg_fact)){
  #undebug(aggregate_raster_fun)
  #list_l_rast <- list(l_rast)
  names(l_rast) <- sub(extension(l_rast),"",basename(l_rast))
  #undebug(aggregate_raster_fun)

  obj <- aggregate_raster_fun(l_rast=l_rast,
                              cat_names=zonal_colnames,
                              agg_method_cat=agg_method_cat,
                              agg_fact=agg_fact,
                              agg_fun=agg_fun,
                              file_format=file_format,
                              rast_ref=rast_ref,
                              num_cores=num_cores,
                              out_suffix=out_suffix,
                              out_dir=out_dir)

  l_rast <- obj$l_rast  
  zonal_colnames <- obj$l_rast_cat
  zonal_colnames<- sub(extension(zonal_colnames),"",
      basename(zonal_colnames))
  l_rast_original <- obj$l_rast_original

  #rast_ref_tmp <- sub(extension(rast_ref[[1]]),"",
  #                    basename(rast_ref[[1]]))
  rast_ref_tmp <- sub(extension(rast_ref),"",
                      basename(rast_ref))
  
  no_cat <- which(rast_ref_tmp==sub(extension(l_rast_original),"",
                                      basename(l_rast_original)))
  
  rast_ref <- l_rast[no_cat]
  
}

###########################
#### PART II: data reporting: description

### Generate data description and figures: make this a markdown output later on!
#debug(explore_and_summarize_data)
#test <- explore_and_summarize_data(l_rast,zonal_colnames, var_names,n_time_event)
#debug(explore_and_summarize_data)
#if(is.null(pixel_index)){
#  gCentroid()
#}


explore_obj <- explore_and_summarize_data(as.character(l_rast), #list o
                                   zonal_colnames, 
                                   var_names,
                                   n_time_event,
                                   proj_str =proj_str,
                                   pixel_index = pixel_index,
                                   animation = TRUE, #will geneate animation for the TS
                                   out_dir = out_dir,
                                   out_suffix=out_suffix)

#l_rast,zonal_colnames,var_names,n_time_event,pixel_index=800,out_dir=NULL,out_suffix=NULL
s_raster <- stack(l_rast)
#names(s_raster) <- names(data_tb)

###########################
#### PART III: run space and time model

#function_space_and_time_predictions <- "space_and_time_predictions_functions_10302017.R"
names(s_raster)
#undebug(run_space_and_time_models)
space_and_time_prediction_obj <- run_space_and_time_models(s_raster,
                                 n_time_event,
                                 time_window_selected,
                                 method_space=method_space,
                                 method_time=method_time, 
                                 NA_flag_val=NA_flag_val,
                                 file_format=file_format,
                                 rast_ref=as.character(rast_ref), # was a lit before
                                 zonal_colnames,
                                 var_names, ## this is the relevant columns to use in s_raster for space and time modeling
                                 num_cores=num_cores,
                                 out_dir=out_dir, 
                                 out_suffix=out_suffix)

  
###################################
#### PART IV: Assessment space and time model


#### 
out_suffix_tmp <- paste("assessment_no_previous_",out_suffix,sep="")
r_spat_pred_no_previous <- stack(space_and_time_prediction_obj$r_spat_pred_no_previous)
r_temp_pred <- stack(space_and_time_prediction_obj$r_temp_pred)

#undebug(plot_map_predictions)

r_plot_no_previous <- plot_map_predictions(n_time_event=n_time_event,
                               r_var=subset(s_raster,var_names),
                               r_spat_pred=r_spat_pred_no_previous,
                               r_temp_pred=r_temp_pred,
                               subset(s_raster,zonal_colnames),
                               z_lim_range=NULL,
                               variable_name=NULL,
                               out_suffix=out_suffix_tmp,
                               out_dir=out_dir)

var_names_tmp <- df_args[11,index_val]
time_window_selected_tmp <- df_args[14,index_val]

#var_names <- df_args[11,index_val] 
#num_cores <- df_args[12,index_val] 
#n_time_event <- df_args[13,index_val]
#time_window_selected <- df_args[14,index_val] 

#debug(accuracy_space_time_calc)
accuracy_space_and_time_obj_no_previous <- accuracy_space_time_calc(
  r_temp_pred=r_temp_pred,
  #r_spat_pred=r_spat_pred,
  r_spat_pred=r_spat_pred_no_previous,
  #s_raster = data_fname,
  s_raster= s_raster,#observed stack
  proj_str = proj_str,
  time_window_selected =time_window_selected_tmp,
  n_time_event = n_time_event,
  r_zonal = zonal_colnames,
  method_space = method_space,
  method_time = method_time,
  #r_ref = r_ref,
  r_ref = rast_ref,
  #out_suffix = out_suffix_test,
  out_suffix = out_suffix_tmp,
  var_names = var_names_tmp,
  NA_flag_val = NA_flag_val,
  file_format =file_format,
  date_range = date_range,
  out_dir = out_dir,
  create_out_dir_param =create_out_dir_param)


r_temp_pred <- stack(space_and_time_prediction_obj$r_temp_pred)
r_spat_pred_with_previous <- stack(space_and_time_prediction_obj$r_spat_pred_with_previous)
#s_raster

var_names_tmp <- paste(var_names[1],var_names[length(var_names)],sep=";")
time_window_selected_tmp <- paste(time_window_selected[1],time_window_selected[length(time_window_selected)],sep=";")
out_suffix_tmp <- paste("assessment_with_previous_",out_suffix,sep="")
#debug(accuracy_space_time_calc)
accuracy_space_and_time_obj_with_previous <- accuracy_space_time_calc(
  r_temp_pred=r_temp_pred,
  #r_spat_pred=r_spat_pred,
  r_spat_pred=r_spat_pred_with_previous,
  #s_raster = data_fname,
  s_raster= s_raster,#observed stack
  proj_str = proj_str,
  time_window_selected =time_window_selected_tmp,
  n_time_event = n_time_event,
  r_zonal = zonal_colnames,
  method_space = method_space,
  method_time = method_time,
  #r_ref = r_ref,
  r_ref = rast_ref,
  #out_suffix = out_suffix_test,
  out_suffix = out_suffix_tmp,
  var_names = var_names_tmp,
  NA_flag_val = NA_flag_val,
  file_format =file_format,
  date_range = date_range,
  out_dir = out_dir,
  create_out_dir_param =create_out_dir_param)

r_plot_with_previous <- plot_map_predictions(n_time_event=n_time_event,
                                           r_var=subset(s_raster,var_names),
                                           r_spat_pred=r_spat_pred_with_previous,
                                           r_temp_pred=r_temp_pred,
                                           subset(s_raster,zonal_colnames),
                                           z_lim_range=NULL,
                                           variable_name=NULL,
                                           out_suffix=out_suffix_tmp,
                                           out_dir=out_dir)

accuracy_space_and_time_obj_no_previous$mae_zones_tb
accuracy_space_and_time_obj_with_previous$mae_zones_tb


###########################  END OF SCRIPT #########################################
