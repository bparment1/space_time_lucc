####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: SAR, SARMA etc -PART 2  #######################################
#This script produces a prediction for the dates following the Hurricane event.       
#The script uses spatial neighbours to predict Light data values in the Hurricane Katrina New Orleans region. 
#Temporal predictions use OLS with the image of the previous time or the ARIMA method.
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/09/2014 
#DATE MODIFIED: 06/08/2017
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
#COMMIT: modifying summary figures and checking all outputs
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

function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_functions_06072017.R" #PARAM 1
function_paper_figures_analyses <- "space_beats_time_sbt_paper_figures_functions_01092016.R" #PARAM 1
function_data_figures_reporting <- "spatial_analysis_data_figures_reporting_functions_06082017.R" #PARAM 1
#script_path <- "/home/parmentier/Data/Space_beats_time/sbt_scripts" #path to script #PARAM 2
script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts"
source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.
source(file.path(script_path,function_paper_figures_analyses)) #source all functions used in this script 1.
source(file.path(script_path,function_data_figures_reporting)) #source all functions used in this script 1.

#Aggregation code
function_multilabel_fuzzy_analyses <- "classification_multilabel_processing_functions_03142017.R" #PARAM 1
script_path <- "/home/bparmentier/Google Drive/LISER_Lux/R_scripts" #path to script #PARAM 2
source(file.path(script_path,function_multilabel_fuzzy_analyses)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

#in_dir <- "~/Data/Space_beats_time/case3data/" #lights/table" #PARAM3
#in_dir <- "/home/parmentier/Data/Space_beats_time/Case2_data_NDVI/"
in_dir <- "/home/bparmentier/Google Drive/Space_beats_time/Case2_data_NDVI/"

#in_dir <- "~/Data/Space_beats_time/case3data/lights/table"
#in_dir <- "~/Data/Space_beats_time/Case1a_data"
#out_dir <- "/home/parmentier/Data/Space_beats_time/outputs"
out_dir <- "/home/bparmentier/Google Drive/Space_beats_time/outputs"

#moore_window <- file.path(in_dir,"window_test4.rst")
#winds_zones_fname <- file.path(in_dir,"00_windzones_moore_sin.rst")

proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" #CONST 1
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 # CONST 2
proj_str<- CRS_WGS84 
CRS_reg <- CRS_WGS84 # PARAM 4

file_format <- ".tif" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"NDVI_Katrina_06082017" #output suffix for the files and output folder #PARAM 8
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

#dates1 <- generate_dates_by_step(date_range1[1],date_range1[2],16)$dates
#dates2 <- unique(year(generate_dates_by_step(date_range2[1],date_range2[2],1)$dates)) #extract year
dates3 <- generate_dates_by_step(date_range3[1],date_range3[2],16)$dates #NDVI Katrina
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

#Transform table text file into a raster image

### 
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
                                   pixel_index = 800,
                                   out_dir = out_dir,
                                   out_suffix=out_suffix)

#l_rast,zonal_colnames,var_names,n_time_event,pixel_index=800,out_dir=NULL,out_suffix=NULL
s_raster <- stack(l_rast)
#names(s_raster) <- names(data_tb)

########################## RUN SPACE AND TIME MODEL

num_cores <- 4
previous_step <- T
method_time <- c("arima","arima") # estimator <- "arima",estimation_method <-"arima"
method_space <- c("mle","eigen") #estimator <- "mle",estimation_method <- "eigen"


run_space_and_time_models <- function(s_raster,n_time_event,time_window_selected,
                                      method_space=c("mle","eigen"),
                                      method_time=c("arima","arima"), 
                                      file_format=".tif",
                                      rast_ref=NULL,
                                      zonal_colnames,
                                      num_cores=1,out_dir=".", out_suffix=NULL){
  #Function to run space and time model given several inputs
  #
  
  if(is.null(out_suffix)){
    out_suffix=""
  }
  ###########################################################################
  ############## PART III PREDICT MODELS  SPAT REGRESSION OVER MULTIPLE time steps ####
  
  ##This times will we use an automated function to generate predictions over multiple dates
  
  ##########################
  #### RUN FOR SELECTED DATES and the three methods..
  
  ###########SPATIAL METHODS
  ## Now predict for four dates using "mle": this does not work for such large area such as EDGY!!
  
  ### Predict using spatial regression: this should be a master function...
  #r_spat_var <- subset(r_stack,139:161) #predict before (step 152) and three dates after (step 153)
  ### Important:
  #As input!!
  #r_spat_var should contain one image before the on being predicted, i.e. to predict step 100, you need at least 
  #a stack of raster 99 and raster 100.
  #num_cores <- 4
  #num_cores_tmp <- 11
  time_step_start <- n_time_event - 8 #this is the time step for which to start the arima model with, start at 99
  time_step_end <- n_time_event + 8
  time_step_subset <- time_step_start -1 #use 99
  time_window_selected <- time_step_subset:time_step_end
  
  time_window_predicted <- time_step_start:time_step_end #100 to 116
  r_spat_var <- subset(s_raster,time_window_selected) #predict before and after event
  
  if(is.null(rast_ref)){
    rast_ref <- subset(s_raster,1)
  }
  
  ########
  ### TEST SPATIAL Prediction on one date....
  
  ### CLEAN OUT AND SCREEN NA and list of neighbours
  #Let's use our function we created to clean out neighbours
  #Prepare dataset 1 for function: date t-2 (mt2) and t-1 (mt1) before hurricane
  #this is for prediction t -1 
  #r_var <- subset(r_stack,n_time_event) #this is the date we want to use to run the spatial regression
  r_clip <- rast_ref #this is the image defining the study area
  
  #############################################
  ### MLE EIGEN Using previous step first
  
  list_models <- NULL
  proj_str <- NULL #if null the raster images are not reprojected
  #the ouput suffix was wrong, needs to be 153!!!
  #Use 100 to 116
  #out_suffix_s <- paste("t_",100:length(time_window_selected),"_",out_suffix,sep="")#this should really be automated!!!
  
  #estimator <- "mle"
  #estimator <- "mle"
  #estimation_method <- "eigen"
  #estimation_method <- "LU"
  #estimation_method <- "Chebyshev"
  #estimator <- gmm
  
  #function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_03022017_functions.R" #PARAM 1
  #source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.
  
  previous_step <- T #Use t-1 neighbour and t-1 value to predict t0
  
  out_suffix_s <- paste("t_",time_window_predicted,"_",out_suffix,sep="")#this should really be automated!!!
  
  #list_param_spat_reg <- list(out_dir,r_spat_var,r_clip,proj_str,list_models,out_suffix_s,file_format,estimator)
  list_param_spat_reg <- list(out_dir,r_spat_var,r_clip,proj_str,list_models,out_suffix_s,
                              file_format,estimator,estimation_method,previous_step)
  names(list_param_spat_reg) <- c("out_dir","r_var_spat","r_clip","proj_str","list_models",
                                  "out_suffix","file_format","estimator","estimation_method","previous_step")
  n_pred <- nlayers(r_spat_var) - 1 # minus one because the first date is not predicted
  #debug(predict_spat_reg_fun)
  #predict_spat_reg_fun(1,list_param=list_param_spat_reg)
  
  pred_spat_mle_eigen_with_previous  <- mclapply(1:n_pred,FUN=predict_spat_reg_fun,
                                                 list_param=list_param_spat_reg,
                                                 mc.preschedule=FALSE,
                                                 mc.cores = num_cores)
  
  save(pred_spat_mle_eigen_with_previous,file=file.path(out_dir,paste("pred_spat_",estimator,"_",estimation_method,
                                                                      "_","with_previous_",out_suffix,".RData",sep="")))
  
  #pred_spat_mle_eigen: extract raster images from object
  spat_pred_rast_mle_eigen_with_previous <- stack(lapply(pred_spat_mle_eigen_with_previous,FUN=function(x){x$raster_pred})) #get stack of predicted images
  spat_res_rast_mle_eigen_with_previous <- stack(lapply(pred_spat_mle_eigen_with_previous,FUN=function(x){x$raster_res})) #get stack of predicted images
  levelplot(spat_pred_rast_mle_eigen_with_previous,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
  levelplot(spat_res_rast_mle_eigen_with_previous,col.regions=matlab.like(25)) #view the four predictions using mle spatial reg.
  
  
  #############################################
  ### STEP 2:
  ### MLE EIGEN not using previous step first
  
  previous_step <- F
  
  #list_param_spat_reg <- list(out_dir,r_spat_var,r_clip,proj_str,list_models,out_suffix_s,file_format,estimator)
  list_param_spat_reg<- list(out_dir,r_spat_var,r_clip,proj_str,list_models,out_suffix_s,file_format,estimator,estimation_method,previous_step)
  names(list_param_spat_reg) <- c("out_dir","r_var_spat","r_clip","proj_str","list_models","out_suffix","file_format","estimator","estimation_method","previous_step")
  n_pred <- nlayers(r_spat_var) -1
  
  #debug(predict_spat_reg_fun)
  #predict_spat_reg_fun(1,list_param=list_param_spat_reg)
  
  #debug(predict_spat_reg_fun)
  #predict_spat_reg_fun(1,list_param=list_param_spat_reg)
  
  pred_spat_mle_eigen_no_previous <- mclapply(1:n_pred,FUN=predict_spat_reg_fun,
                                              list_param=list_param_spat_reg,
                                              mc.preschedule=FALSE,
                                              mc.cores = num_cores)
  
  save(pred_spat_mle_eigen_no_previous,file=file.path(out_dir,paste("pred_spat_",estimator,"_",estimation_method,"_",
                                                                    "_","no_previous_",out_suffix,".RData",sep="")))
  
  #pred_spat_mle_eigen: extract raster images from object
  spat_pred_rast_mle_eigen_no_previous <- stack(lapply(pred_spat_mle_eigen_no_previous,FUN=function(x){x$raster_pred})) #get stack of predicted images
  spat_res_rast_mle_eigen_no_previous <- stack(lapply(pred_spat_mle_eigen_no_previous,FUN=function(x){x$raster_res})) #get stack of predicted images
  levelplot(spat_pred_rast_mle_eigen_no_previous,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
  levelplot(spat_res_rast_mle_eigen_no_previous,col.regions=matlab.like(25)) #view the four predictions using mle spatial reg.
  
  
  ##############################################################################################
  ############## PART III TEMPORAL METHODS: PREDICT MODELS  FOR USING TEMP REGRESSION OVER MULTIPLE time steps ####
  
  ## Predict using temporal info: time steps 
  
  ### Use ARIMA!!! with time before
  
  #out_dir
  #r_clip
  #proj_str
  #file_format
  
  r_temp_var <- subset(s_raster,time_window_selected) # relies on the previous date in contrast to spat reg, this is r_var param
  list_models <-NULL
  
  #the ouput suffix was wrong, needs to be 153!!!
  #Use 100 to 116
  time_step_start <- n_time_event - 8 #this is the time step for which to start the arima model with, start at 99
  #100 for NDVI
  out_suffix_s <- paste("t_",time_step_start:length(time_window_selected),"_",out_suffix,sep="")#this should really be automated!!!
  estimator <- "lm"
  estimation_method <-"ols"
  
  #ARIMA specific
  
  num_cores_tmp <- num_cores
  time_step <- n_time_event - 8 #this is the time step for which to start the arima model with, start at 99
  n_pred_ahead <- 16
  rast_ref <- subset(s_raster,1) #first image ID
  r_stack_arima <- mask(s_raster,rast_ref)
  
  #r_stack <- r_stack_arima
  arima_order <- NULL
  
  #list_param_temp_reg <- list(out_dir,r_temp_var,r_clip,proj_str,list_models,out_suffix_s,file_format,estimator,estimation_method,
  #                            num_cores_tmp,time_step,n_pred_ahead,r_stack_arima,arima_order,NA_flag_val)
  #names(list_param_temp_reg) <- c("out_dir","r_var","r_clip","proj_str","list_models","out_suffix_s","file_format","estimator","estimation_method",
  #                            "num_cores","time_step","n_pred_ahead","r_stack","arima_order","NA_flag_val")
  #n_pred <- nlayers(r_temp_var) -1
  #undebug(predict_temp_reg_fun)
  #test_temp <- predict_temp_reg_fun(14,list_param_temp_reg)
  #plot(raster(test_temp$raster_pred),main=basename(test_temp$raster_pred))
  
  #source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.
  
  #pred_temp_lm <- lapply(1:n_pred,FUN=predict_temp_reg_fun,list_param=list_param_temp_reg) 
  
  ##ARIMA
  
  #not run here...
  
  r_temp_var <- subset(s_raster,time_window_selected) # relies on the previous date in contrast to spat reg, this is r_var param
  list_models <-NULL
  
  #the ouput suffix was wrong, needs to be 153!!!
  #Use 100 to 116
  time_step_start <- n_time_event - 8 #this is the time step for which to start the arima model with, start at 99
  
  out_suffix_s <- paste("t_",time_step_start:length(time_window_selected),"_",out_suffix,sep="")#this should really be automated!!!
  estimator <- "lm"
  estimation_method <-"ols"
  
  #ARIMA specific
  
  num_cores_tmp <- num_cores
  time_step <- n_time_event - 8 #this is the time step for which to start the arima model with, start at 99
  n_pred_ahead <- 16
  rast_ref <- subset(s_raster,1) #first image ID
  r_stack_arima <- mask(s_raster,rast_ref)
  
  #r_stack <- r_stack_arima
  arima_order <- NULL
  
  estimator <- "arima"
  estimation_method <-"arima"
  r_clip_tmp <- NULL
  r_clip_tmp <- rast_ref
  
  num_cores_tmp <- num_cores
  
  list_param_temp_reg <- list(out_dir,r_temp_var,r_clip_tmp,proj_str,list_models,out_suffix_s,file_format,estimator,estimation_method,
                              num_cores_tmp,time_step,n_pred_ahead,r_stack_arima,arima_order,NA_flag_val)
  names(list_param_temp_reg) <- c("out_dir","r_var","r_clip","proj_str","list_models","out_suffix_s","file_format","estimator","estimation_method",
                                  "num_cores","time_step","n_pred_ahead","r_stack","arima_order","NA_flag_val")
  #n_pred <- nlayers(r_temp_var) -1
  n_pred <- 16
  #debug(predict_temp_reg_fun)
  
  if(re_initialize_arima==T){
    l_pred_temp_arima <- vector("list",length=length(time_window_selected))
    
    for(i in 1:length(time_window_selected)){
      n_pred <- 1
      n_pred_ahead <- n_pred
      #time_step <- n_time_event - 8 #this is the time step for which to start the arima model with, start at 99
      
      #n_time_pred_start <- 99 + i
      n_time_pred_start <- time_window_selected[1] - 1 + i
      time_step <-n_time_pred_start
      n_time_pred_end <- n_time_pred_start + length(time_window_selected)
      
      out_suffix_s <- paste("t_",n_time_pred_start:n_time_pred_end,"_",out_suffix,sep="")#this should really be automated!!!
      list_param_temp_reg <- list(out_dir,r_temp_var,r_clip_tmp,proj_str,list_models,out_suffix_s,file_format,estimator,estimation_method,
                                  num_cores_tmp,time_step,n_pred_ahead,r_stack_arima,arima_order,NA_flag_val)
      names(list_param_temp_reg) <- c("out_dir","r_var","r_clip","proj_str","list_models","out_suffix_s","file_format","estimator","estimation_method",
                                      "num_cores","time_step","n_pred_ahead","r_stack","arima_order","NA_flag_val")
      #undebug(predict_temp_reg_fun)
      pred_temp_arima <- predict_temp_reg_fun(i,list_param_temp_reg) #only one date predicted...four step ahead
      l_pred_temp_arima[[i]] <- pred_temp_arima  #only one date predicted...one step ahead
    }
  }else{
    time_step <- n_time_event - 8 #this is the time step for which to start the arima model with, start at 99
    #n_pred <- 16
    n_pred <- length(time_window_selected) -1
    n_pred_ahead <- n_pred
    
    #n_time_pred_start <- 100 
    n_time_pred_start <- time_window_selected[1] 
    
    n_time_pred_end <- n_time_pred_start + length(time_window_selected)
    
    out_suffix_s <- paste("t_",n_time_pred_start:n_time_pred_end,"_",out_suffix,sep="")#this should really be automated!!!
    list_param_temp_reg <- list(out_dir,r_temp_var,r_clip_tmp,proj_str,list_models,out_suffix_s,file_format,estimator,estimation_method,
                                num_cores_tmp,time_step,n_pred_ahead,r_stack_arima,arima_order,NA_flag_val)
    names(list_param_temp_reg) <- c("out_dir","r_var","r_clip","proj_str","list_models","out_suffix_s","file_format","estimator","estimation_method",
                                    "num_cores","time_step","n_pred_ahead","r_stack","arima_order","NA_flag_val")
    
    #debug(predict_temp_reg_fun)
    pred_temp_arima <- predict_temp_reg_fun(1,list_param_temp_reg) #only one date predicted...four step ahead
    
  }
  

  ############ PART V COMPARE MODELS IN PREDICTION ACCURACY #################
  
  #source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.
  
  #using 11 cores
  #pred_temp_arima <- lapply(1:n_pred,FUN=predict_temp_reg_fun,list_param=list_param_temp_reg)
  
  #pred_temp_arima <- load_obj("temp_reg_obj_arima_arima_t_153_EDGY_predictions_03182015.RData")
  #pred_temp_arima <- load_obj("temp_reg_obj_arima_arima_t_100_NDVI_Katrina_04102015.RData")
  
  #extract_files <- function(i,x){obj<-x[[i]];obj$raster_pred}
  #debug(extract_files)
  #extract_files(1,l_pred_temp_arima)
  r_temp_pred_rast_arima <- stack(unlist(lapply(1:length(l_pred_temp_arima),FUN=function(i,x){obj<-x[[i]];obj$raster_pred},x=l_pred_temp_arima)))
  r_temp_res_rast_arima <- stack(unlist(lapply(1:length(l_pred_temp_arima),FUN=function(i,x){obj<-x[[i]];obj$raster_res},x=l_pred_temp_arima)))
  
  levelplot(r_temp_pred_rast_arima,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
  levelplot(r_temp_res_rast_arima,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
  
  #r_temp_pred_rast_lm <- stack(lapply(pred_temp_lm,FUN=function(x){x$raster_pred}))
  #r_temp_res_rast_lm <- stack(lapply(pred_temp_lm,FUN=function(x){x$raster_res}))
  #levelplot(r_temp_pred_rast_lm,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
  #levelplot(r_temp_res_rast_lm,col.regions=rev(terrain.colors(255)),main="Var residuals after hurricane")
  
  levelplot(spat_pred_rast_mle_eigen_no_previous,col.regions=rev(matlab.like(255))) #view the four predictions using mle spatial reg.
  levelplot(spat_res_rast_mle_eigen_no_previous,col.regions=rev(matlab.like(255))) #view the four predictions using mle spatial reg.
  
  levelplot(spat_pred_rast_mle_eigen_with_previous,col.regions=matlab.like(255)) #view the four predictions using mle spatial reg.
  levelplot(spat_res_rast_mle_eigen_with_previous,col.regions=matlab.like(255)) #view the four predictions using mle spatial reg.
  
  #spat_pred_rast_mle_eigen_no_previous <- stack(lapply(pred_spat_mle_eigen_no_previous,FUN=function(x){x$raster_pred})) #get stack of predicted images
  #spat_res_rast_mle_eigen_no_previous <- stack(lapply(pred_spat_mle_eigen_no_previous,FUN=function(x){x$raster_res})) #get stack of predicted images
  #levelplot(spat_pred_rast_mle_eigen_no_previous,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
  #levelplot(spat_res_rast_mle_eigen_no_previous,col.regions=matlab.like(25)) #view the four predictions using mle spatial reg.
  
  projection(r_temp_pred_rast_arima) <- proj_str
  projection(r_temp_res_rast_arima) <- proj_str
  
  #projection(r_temp_res_rast_lm) <- CRS_reg
  #projection(r_temp_pred_rast_lm) <- CRS_reg
  
  # Problem with rast_zonal: when aggregated categorical value become continuous
  #
  #
  
  projection(spat_pred_rast_mle_eigen_no_previous) <- proj_str
  projection(spat_res_rast_mle_eigen_no_previous) <- proj_str
  projection(spat_pred_rast_mle_eigen_with_previous) <- proj_str
  projection(spat_res_rast_mle_eigen_with_previous) <- proj_str
  
  projection(r_temp_pred_rast_arima) <- proj_str
  #projection(r_temp_pred_rast_lm) <- CRS_reg
  
  rast_zonal <- subset(s_raster,zonal_colnames)
  
  projection(rast_ref) <- proj_str
  projection(rast_zonal) <- proj_str
  projection(s_raster) <- proj_str
  
  #r_huric_obs <- subset(s_raster,time_window_selected[-1]) #remove 99 because not considred in the prediction!!
  r_huric_obs <- subset(s_raster,time_window_predicted) #remove 99 because not considred in the prediction!!
  
  #plot(r_huric_obs)
  #r_huric_w <- crop(r_huric_w,rast_ref)
  levelplot(r_huric_obs,col.regions=matlab.like(25))
  levelplot(r_temp_pred_rast_arima,col.regions=matlab.like(25))
  
  #r_winds_m <- crop(winds_wgs84,res_temp_s) #small test window
  #res_temp_s_lm <- temp_pred_rast_lm - r_huric_obs
  temp_pred_rast_arima <- subset(r_temp_pred_rast_arima,1:length(time_window_predicted))
  
  res_temp_s_arima <- temp_pred_rast_arima - r_huric_obs
  res_spat_s <- spat_pred_rast_mle_eigen_no_previous - r_huric_obs
  #res_spat_s <- spat_pred_rast_mle_eigen_no_previous - r_huric_obs
  
  #names(res_temp_s_lm) <- sub("pred","res",names(res_temp_s_lm))
  names(res_temp_s_arima) <- sub("pred","res",names(res_temp_s_arima))
  
  #names(res_spat_mle_s) <- sub("pred","res",names(res_spat_mle_s))
  names(res_spat_s) <- sub("pred","res",names(res_spat_s))
  #names(res_temp_s) <- paste("r_res_s_",1:nlayers(res_temp_s),"_",out_suffix,sep="")
  
  #r_results <- stack(s_raster,rast_zonal,temp_pred_rast_arima,
  #                   spat_pred_rast_mle_eigen,spat_pred_rast_mle_Chebyshev,res_temp_s_arima,res_temp_s_lm,res_spat_s)
  
  r_results <- stack(s_raster,rast_zonal,temp_pred_rast_arima,
                     spat_pred_rast_mle_eigen_no_previous,spat_pred_rast_mle_eigen_with_previous,res_temp_s_arima,res_spat_s)
  
  dat_out <- as.data.frame(r_results)
  dat_out <- na.omit(dat_out)
  write.table(dat_out,file=paste("dat_out_",out_suffix,".txt",sep=""),
              row.names=F,sep=",",col.names=T)
  
  ### Now accuracy assessment using MAE
  
  out_suffix_s <- paste("temp_arima_",out_suffix,sep="_")
  
  #undebug(calc_ac_stat_fun)
  ac_temp_arima_obj <- calc_ac_stat_fun(r_pred_s=temp_pred_rast_arima,
                                        r_var_s=r_huric_obs,
                                        r_zones=rast_zonal,
                                        file_format=file_format,
                                        out_suffix=out_suffix_s)
  
  #out_suffix_s <- paste("temp_lm_",out_suffix,sep="_")
  
  #ac_temp_lm_obj <- calc_ac_stat_fun(r_pred_s=temp_pred_rast_lm,r_var_s=r_huric_obs,r_zones=rast_zonal,
  #                                file_format=file_format,out_suffix=out_suffix_s)
  
  out_suffix_s <- paste("spat_mle_eigen_with_previous",out_suffix,sep="_")  
  ac_spat_mle_eigen_with_previous_obj <- calc_ac_stat_fun(r_pred_s=spat_pred_rast_mle_eigen_with_previous,
                                                          r_var_s=r_huric_obs,
                                                          r_zones=rast_zonal,
                                                          file_format=file_format,
                                                          out_suffix=out_suffix_s)
  
  out_suffix_s <- paste("spat_mle_eigen_no_previous",out_suffix,sep="_")  
  ac_spat_mle_eigen_no_previous_obj <- calc_ac_stat_fun(r_pred_s=spat_pred_rast_mle_eigen_no_previous,
                                                        r_var_s=r_huric_obs,
                                                        r_zones=rast_zonal,
                                                        file_format=file_format,
                                                        out_suffix=out_suffix_s)
  
  #mae_tot_tb <- t(rbind(ac_spat_obj$mae_tb,ac_temp_obj$mae_tb))
  #mae_tot_tb <- (cbind(ac_spat_obj$mae_tb,ac_temp_obj$mae_tb))
  #mae_tot_tb <- (cbind(ac_spat_mle_eigen_no_previous_obj$mae_tb,
  #ac_spat_mle_eigen_with_previous_obj$mae_tb,
  #ac_temp_arima_obj$mae_tb
  #,ac_temp_lm_obj$mae_tb))
  
  mae_tot_tb <- (cbind(ac_spat_mle_eigen_no_previous_obj$mae_tb,
                       ac_spat_mle_eigen_with_previous_obj$mae_tb,
                       ac_temp_arima_obj$mae_tb
  ))
  
  mae_tot_tb <- as.data.frame(mae_tot_tb)
  row.names(mae_tot_tb) <- NULL
  names(mae_tot_tb)<- c("spat_reg_no_previous","spat_reg_with_previous","temp_arima")#,"temp_lm")
  #mae_tot_tb$time <- 2:nrow(mae_tot_tb)
  #mae_tot_tb$time <- 2:n_pred
  mae_tot_tb$time <- 1:nlayers(r_huric_obs)
  
  y_range<- range(cbind(mae_tot_tb$spat_reg_no_previous,mae_tot_tb$spat_reg_with_previous,mae_tot_tb$temp_arima))
  
  res_pix<-960
  col_mfrow<-1
  row_mfrow<-1
  png(filename=paste("Figure_temporal_profiles_time_1_",out_suffix,".png",sep=""),
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  plot(spat_reg_no_previous ~ time, type="b",col="cyan",data=mae_tot_tb,ylim=y_range)
  lines(temp_arima ~ time, type="b",col="magenta",data=mae_tot_tb)
  lines(spat_reg_with_previous ~ time, type="b",col="blue",data=mae_tot_tb)
  
  #lines(temp_lm ~ time, type="b",col="red",data=mae_tot_tb)
  legend("topleft",legend=c("spat_no","temp_arima","spat_with"),col=c("cyan","magenta","blue"),lty=1)
  title("Overall MAE for spatial and temporal models") #Note that the results are different than for ARIMA!!!
  
  dev.off()
  
  write.table(mae_tot_tb,file=paste("mae_tot_tb","_",out_suffix,".txt",sep=""))
  
  #### BY ZONES ASSESSMENT: Ok now it is general so it should be part of the function...
  
  #mae_zones_tb <- rbind(ac_spat_mle_obj$mae_zones_tb[1:3,],
  #                      ac_temp_obj$mae_zones_tb[1:3,])
  mae_zones_tb <- rbind(ac_spat_mle_eigen_no_previous_obj$mae_zones_tb,
                        ac_spat_mle_eigen_with_previous_obj$mae_zones_tb,
                        ac_temp_arima_obj$mae_zones_tb)
  
  mae_zones_tb <- as.data.frame(mae_zones_tb)
  
  n_zones <- length(unique(mae_zones_tb$zone))
  
  mae_zones_tb$method <- c(rep("spat_reg_no",n_zones),rep("temp_with_reg",n_zones),rep("temp_arima_reg",n_zones))
  
  n_time <- ncol(mae_zones_tb) -1
  pred_names <- c("zone",paste("t",2:n_time,sep="_"),"method")
  #names(mae_zones_tb) <- c("zones","pred1","pred2","pred3","pred4","method")
  names(mae_zones_tb) <- pred_names
  
  write.table(mae_zones_tb,file=paste("mae_zones_tb","_",out_suffix,".txt",sep=""))
  
  mydata<- mae_zones_tb
  dd <- do.call(make.groups, mydata[,-ncol(mydata)]) 
  #dd$lag <- mydata$lag 
  #drop first few rows that contain no data but zones...
  n_start <-n_zones*3 +1 #3 because we have 3 methods...
  #n_start <-n_zones*2 +1
  dd <- dd[n_start:nrow(dd),]
  
  dd$zones <- mydata$zones
  dd$method <- mydata$method
  
  dd$zones <- mydata$zone #use recycle rule
  
  #Note that to get the title correct, one needs to add
  res_pix<-960
  col_mfrow<-1
  row_mfrow<-1
  png(filename=paste("Figure1_ref_layer_time_1_",out_suffix,".png",sep=""),
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  xyplot(data~which | as.factor(zones) ,group=method,data=dd,type="b",xlab="time",ylab="VAR",
         strip = strip.custom(factor.levels=unique(as.character(dd$zones))), #fix this!!!
         auto.key = list("topright", corner = c(0,1),# col=c("black","red"),
                         border = FALSE, lines = TRUE,cex=1.2)
  )
  dev.off()
  #histogram(r_zonal)
  
  ############## Prepare return object #############
  
  ###
  
  ##
  return()
}


################### END OF SCRIPT ##################