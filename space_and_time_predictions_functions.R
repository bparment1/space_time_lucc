#######################################    Space Time Analyses Project   #######################################
############################  Running Space and Time Model for Space Beats Time framework #######################################
#This script contains functions to produce predictions for Space Beats Time Framework.       
#The script uses functions described additional sript for the project.
#Temporal predictions use OLS with the image of the previous time step or ARIMA.
#Spatial predictions use spatial regression (lag error model) with different estimation methods (e.g. eigen, chebyshev etc.).
#A model with space and time is implemented using neighbours from the previous time step.
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 06/23/2017 
#DATE MODIFIED: 08/07/2017
#Version: 1
#PROJECT: GLP Conference Berlin,YUCATAN CASE STUDY with Marco Millones            
#PROJECT: Workshop for William and Mary: an intro to spatial regression with R 
#PROJECT: Geocomputation and AAG 2015
#PROJECT: Space beats time paper

#TO DO:
# Simplify and clean up code
#
#COMMIT: adding aggregation raster function to test SBT at different spatial resolutions
#
#################################################################################################

#This script currently contains one functions.

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

aggregate_raster_fun <- function(l_rast,zonal_colnames,use_majority,agg_fact,agg_fun,file_format,rast_ref,num_cores,out_suffix, out_dir){
  #
  #Function to aggregate input raster stack
  #if use majority then the zonal layer is aggregated and then reclassfied based by the majority rule
  
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
  
  if(use_majority==TRUE){
    
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
    
  }
  
  if(use_majority==FALSE){
    #sure SRTM and reclass based on threshold values?
    
  }
  
  #r_srtm_Katrina_rec2
  #-rw-rw-r-- 1 bparmentier bparmentier 1894 Apr  7 12:33 r_r_srtm_Katrina_rec2_NDVI_Katrina_04062017.tif
  #-rw-rw-r-- 1 bparmentier bparmentier 1016 Apr  7 12:34 agg_5_r_r_srtm_Katrina_rec2_NDVI_Katrina_04062017.tif
  
  ###
  zonal_colnames <- gsub(extension(raster_name),"",raster_name)
  ##
  
  ##########################
  #### prepare return object
  
  obj <- list(zonal_colnames,l_rast,l_rast_original)
  names(obj) <- c("zonal_colnames","l_rast","l_rast_original")
  
  return(obj)
}


run_space_and_time_models <- function(s_raster,n_time_event,time_window_selected,
                                      method_space=c("mle","eigen"),
                                      method_time=c("arima","arima",T),
                                      NA_flag_val=-9999,
                                      file_format=".tif",
                                      rast_ref=NULL,
                                      zonal_colnames,
                                      num_cores=1,out_dir=".", 
                                      out_suffix=NULL){
  #Function to run space and time model given several inputs
  #
  
  if(is.null(out_suffix)){
    out_suffix=""
  }
  
  projection(s_raster) <- proj_str
  
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
  time_step_start <- time_window_selected[1]
  time_step_end <- time_window_selected[length(time_window_selected)]
  
  #time_step_start <- n_time_event - 8 #this is the time step for which to start the arima model with, start at 99
  #time_step_end <- n_time_event + 8
  time_step_subset <- time_step_start +1 # this is because we miss the first date of pred!!
  #time_window_selected <- time_step_subset:time_step_end
  
  time_window_predicted <- time_step_subset:time_step_end #100 to 116
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
  #proj_str <- NULL #if null the raster images are not reprojected
  
  #Use 100 to 116
  #out_suffix_s <- paste("t_",100:length(time_window_selected),"_",out_suffix,sep="")#this should really be automated!!!
  
  estimator <- method_space[1]
  estimation_method <- method_space[2]
  
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
  #browser()
  
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
  
  estimator <- method_time[1]
  estimation_method <- method_time[2]
  re_initialize_arima <- method_time[3] # set to TRUE or FALSE for ARIMA method
  
  #rast_ref <- subset(s_raster,1) #first image ID
  r_stack <- mask(s_raster,rast_ref)
  
  r_temp_var <- subset(r_stack,time_window_selected) # relies on the previous date in contrast to spat reg, this is r_var param
  #r_temp_var <- subset(s_raster,time_window_selected) # relies on the previous date in contrast to spat reg, this is r_var param
  list_models <-NULL
  num_cores_tmp <- num_cores
  
  browser()
  ###############################
  ### Predict using OLS 
  if(estimator=="lm"){
    
    #Use 100 to 116
    #time_step_start <- n_time_event - 8 #this is the time step for which to start the arima model with, start at 99
    #100 for NDVI
    #out_suffix_s <- paste("t_",time_step_start:length(time_window_selected),"_",out_suffix,sep="")#this should really be automated!!!
    #Note the first date is dropped!!!
    out_suffix_s <- paste("t_",time_window_predicted,"_",out_suffix,sep="")#this should really be automated!!!
    
    #estimator <- "lm"
    #estimation_method <-"ols"
    
    #ARIMA specific
    #time_step <- time_step_start - 1
    time_step <- time_step_start + 1 # We are
    
    #time_step <- n_time_event - 8 #this is the time step for which to start the arima model with, start at 99
    #n_pred_ahead <- 16
    arima_order <- NULL
    r_clip_tmp <- rast_ref
    
    n_pred <- length(time_window_predicted)
    
    #r_temp_var <- subset(r_stack,c(time_step,time_window_selected)) # relies on the previous date in contrast to spat reg, this is r_var param
    r_temp_var <- subset(r_stack,time_window_selected) # relies on the previous date in contrast to spat reg, this is r_var param
    
    list_param_temp_reg <- list(out_dir,r_temp_var,r_clip_tmp,proj_str,list_models,out_suffix_s,file_format,estimator,estimation_method,
                                num_cores_tmp,time_step,n_pred,r_stack,arima_order,NA_flag_val)
    names(list_param_temp_reg) <- c("out_dir","r_var","r_clip","proj_str","list_models","out_suffix_s","file_format","estimator","estimation_method",
                                    "num_cores","time_step","n_pred_ahead","r_stack","arima_order","NA_flag_val")
    
    #list_param_temp_reg <- list(out_dir,r_temp_var,r_clip,proj_str,list_models,out_suffix_s,
    #                            file_format,estimator,estimation_method,
    #                          num_cores_tmp,time_step,n_pred_ahead,r_stack,arima_order,NA_flag_val)
    #names(list_param_temp_reg) <- c("out_dir","r_var","r_clip","proj_str","list_models","out_suffix_s","file_format","estimator","estimation_method",
    #                            "num_cores","time_step","n_pred_ahead","r_stack","arima_order","NA_flag_val")
    #n_pred <- nlayers(r_temp_var) -1
    #undebug(predict_temp_reg_fun)
    #test_temp <- predict_temp_reg_fun(14,list_param_temp_reg)
    #plot(raster(test_temp$raster_pred),main=basename(test_temp$raster_pred))
    
    #source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.
    
    pred_temp_lm <- lapply(1:n_pred,FUN=predict_temp_reg_fun,list_param=list_param_temp_reg) 
    
  }
  
  ###############################
  ### Predict using ARIMA
  browser()
  
  if(estimator=="arima"){

    #time_step_start <- time_window_selected[1]
    #time_step_end <- time_window_selected[length(time_window_selected)]
    
    #time_step_start <- n_time_event - 8 #this is the time step for which to start the arima model with, start at 99
    #time_step_end <- n_time_event + 8
    #time_step_subset <- time_step_start +1 # this is because we miss the first date of pred!!
    
    #time_window_selected <- time_step_subset:time_step_end
    #Use 100 to 116
    #time_step_start <- n_time_event - 8 #this is the time step for which to start the arima model with, start at 99
    #out_suffix_s <- paste("t_",time_step_start:length(time_window_selected),"_",out_suffix,sep="")#this should really be automated!!!
    out_suffix_s <- paste("t_",time_window_predicted,"_",out_suffix,sep="")#this should really be automated!!!
    
    #ARIMA specific
    #num_cores_tmp <- num_cores
    #time_step <- n_time_event - 8 #this is the time step for which to start the arima model with, start at 99
    time_step <- time_step_start + 1 # We are
    
    n_pred_ahead <- nlayers(r_temp_var) -1 
    
    arima_order <- NULL
    r_clip_tmp <- rast_ref
    
    list_param_temp_reg <- list(out_dir,r_temp_var,r_clip_tmp,proj_str,list_models,out_suffix_s,file_format,estimator,estimation_method,
                                num_cores_tmp,time_step,n_pred_ahead,r_stack,arima_order,NA_flag_val)
    names(list_param_temp_reg) <- c("out_dir","r_var","r_clip","proj_str","list_models","out_suffix_s","file_format","estimator","estimation_method",
                                    "num_cores","time_step","n_pred_ahead","r_stack","arima_order","NA_flag_val")
    #n_pred <- nlayers(r_temp_var) -1
    n_pred <- n_pred_ahead
    #debug(predict_temp_reg_fun)
    ## problem to resolve the arima model from re_initialization are overwritten
    
    if(re_initialize_arima==T){
      #l_pred_temp_arima <- vector("list",length=length(time_window_selected))
      l_pred_temp_arima <- vector("list",length=length(time_window_selected))
      
      #To predict the last date, we need the previous date
      for(i in 1:length(time_window_selected)){
      #for(i in 1:length(time_window_predicted)){
          
        n_pred <- 1
        n_pred_ahead <- n_pred
        #time_step <- n_time_event - 8 #this is the time step for which to start the arima model with, start at 99
        
        #n_time_pred_start <- 99 + i
        n_time_pred_start <- time_window_selected[1] - 1 + i
        
        time_step <-n_time_pred_start
        n_time_pred_end <- n_time_pred_start + length(time_window_selected)
        
        #out_suffix_s <- paste("t_",n_time_pred_start:n_time_pred_end,"_",out_suffix,sep="")#this should really be automated!!!
        
        list_param_temp_reg <- list(out_dir,r_temp_var,r_clip_tmp,proj_str,list_models,out_suffix_s,file_format,estimator,estimation_method,
                                    num_cores_tmp,time_step,n_pred_ahead,r_stack,arima_order,NA_flag_val)
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
                                  num_cores_tmp,time_step,n_pred_ahead,r_stack,arima_order,NA_flag_val)
      names(list_param_temp_reg) <- c("out_dir","r_var","r_clip","proj_str","list_models","out_suffix_s","file_format","estimator","estimation_method",
                                      "num_cores","time_step","n_pred_ahead","r_stack","arima_order","NA_flag_val")
      
      #debug(predict_temp_reg_fun)
      pred_temp_arima <- predict_temp_reg_fun(1,list_param_temp_reg) #only one date predicted...four step ahead
      
    }
    
    
  }
  
  ############ PART V COMPARE MODELS IN PREDICTION ACCURACY #################
  browser()
  #source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.
  
  #using 11 cores
  #pred_temp_arima <- lapply(1:n_pred,FUN=predict_temp_reg_fun,list_param=list_param_temp_reg)
  
  #pred_temp_arima <- load_obj("temp_reg_obj_arima_arima_t_153_EDGY_predictions_03182015.RData")
  #pred_temp_arima <- load_obj("temp_reg_obj_arima_arima_t_100_NDVI_Katrina_04102015.RData")
  
  #extract_files <- function(i,x){obj<-x[[i]];obj$raster_pred}
  #debug(extract_files)
  #extract_files(1,l_pred_temp_arima)
  if(method_time[1]=="arima"){
    r_temp_pred_rast_arima <- stack(unlist(lapply(1:length(l_pred_temp_arima),FUN=function(i,x){obj<-x[[i]];obj$raster_pred},x=l_pred_temp_arima)))
    r_temp_res_rast_arima <- stack(unlist(lapply(1:length(l_pred_temp_arima),FUN=function(i,x){obj<-x[[i]];obj$raster_res},x=l_pred_temp_arima)))
    levelplot(r_temp_pred_rast_arima,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
    levelplot(r_temp_res_rast_arima,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
    projection(r_temp_pred_rast_arima) <- proj_str
    projection(r_temp_res_rast_arima) <- proj_str
    temp_pred_rast <- subset(r_temp_pred_rast_arima,1:length(time_window_predicted))
    
  }
  
  if(method_time[1]=="lm"){
    r_temp_pred_rast_lm <- stack(lapply(pred_temp_lm,FUN=function(x){x$raster_pred}))
    r_temp_res_rast_lm <- stack(lapply(pred_temp_lm,FUN=function(x){x$raster_res}))
    levelplot(r_temp_pred_rast_lm,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
    levelplot(r_temp_res_rast_lm,col.regions=rev(terrain.colors(255)),main="Var residuals after hurricane")
    projection(r_temp_res_rast_lm) <- proj_str
    projection(r_temp_pred_rast_lm) <- proj_str
    temp_pred_rast <- subset(r_temp_pred_rast_lm,1:length(time_window_predicted))
    
  }
  
  levelplot(spat_pred_rast_mle_eigen_no_previous,col.regions=rev(matlab.like(255))) #view the four predictions using mle spatial reg.
  levelplot(spat_res_rast_mle_eigen_no_previous,col.regions=rev(matlab.like(255))) #view the four predictions using mle spatial reg.
  
  levelplot(spat_pred_rast_mle_eigen_with_previous,col.regions=matlab.like(255)) #view the four predictions using mle spatial reg.
  levelplot(spat_res_rast_mle_eigen_with_previous,col.regions=matlab.like(255)) #view the four predictions using mle spatial reg.
  
  #spat_pred_rast_mle_eigen_no_previous <- stack(lapply(pred_spat_mle_eigen_no_previous,FUN=function(x){x$raster_pred})) #get stack of predicted images
  #spat_res_rast_mle_eigen_no_previous <- stack(lapply(pred_spat_mle_eigen_no_previous,FUN=function(x){x$raster_res})) #get stack of predicted images
  #levelplot(spat_pred_rast_mle_eigen_no_previous,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
  #levelplot(spat_res_rast_mle_eigen_no_previous,col.regions=matlab.like(25)) #view the four predictions using mle spatial reg.
  
  # Problem with rast_zonal: when aggregated categorical value become continuous
  #
  #
  
  projection(spat_pred_rast_mle_eigen_no_previous) <- proj_str
  projection(spat_res_rast_mle_eigen_no_previous) <- proj_str
  projection(spat_pred_rast_mle_eigen_with_previous) <- proj_str
  projection(spat_res_rast_mle_eigen_with_previous) <- proj_str
  
  #r_huric_obs <- subset(s_raster,time_window_selected[-1]) #remove 99 because not considred in the prediction!!
  r_huric_obs <- subset(s_raster,time_window_predicted) #remove 99 because not considred in the prediction!!
  
  #plot(r_huric_obs)
  #r_huric_w <- crop(r_huric_w,rast_ref)
  levelplot(r_huric_obs,col.regions=matlab.like(25))
  #levelplot(r_temp_pred_rast_arima,col.regions=matlab.like(25))
  
  #r_winds_m <- crop(winds_wgs84,res_temp_s) #small test window
  #res_temp_s_lm <- temp_pred_rast_lm - r_huric_obs
  #res_temp_s_arima <- temp_pred_rast_arima - r_huric_obs
  
  res_temp_s <- temp_pred_rast - r_huric_obs
  res_spat_s <- spat_pred_rast_mle_eigen_no_previous - r_huric_obs
  #res_spat_s <- spat_pred_rast_mle_eigen_no_previous - r_huric_obs
  
  #names(res_temp_s_lm) <- sub("pred","res",names(res_temp_s_lm))
  #names(res_temp_s_arima) <- sub("pred","res",names(res_temp_s_arima))
  names(res_temp_s) <- sub("pred","res",names(res_temp_s))
  
  #names(res_spat_mle_s) <- sub("pred","res",names(res_spat_mle_s))
  names(res_spat_s) <- sub("pred","res",names(res_spat_s))
  #names(res_temp_s) <- paste("r_res_s_",1:nlayers(res_temp_s),"_",out_suffix,sep="")
  
  #r_results <- stack(s_raster,rast_zonal,temp_pred_rast_arima,
  #                   spat_pred_rast_mle_eigen,spat_pred_rast_mle_Chebyshev,res_temp_s_arima,res_temp_s_lm,res_spat_s)
  
  #browser()
  rast_zonal <- subset(s_raster,match(zonal_colnames,names(s_raster)))
  
  #r_results <- stack(s_raster,rast_zonal,temp_pred_rast_arima,
  #                   spat_pred_rast_mle_eigen_no_previous,
  #                   spat_pred_rast_mle_eigen_with_previous,
  #                   res_temp_s_arima,res_spat_s)
  r_results <- stack(s_raster,rast_zonal,temp_pred_rast,
                     spat_pred_rast_mle_eigen_no_previous,
                     spat_pred_rast_mle_eigen_with_previous,
                     res_temp_s,res_spat_s)
  
  dat_out <- as.data.frame(r_results)
  dat_out <- na.omit(dat_out)
  
  filename_dat_out <- file.path(out_dir,paste("dat_out_",out_suffix,".txt",sep=""))
  write.table(dat_out,file=filename_dat_out,
              row.names=F,sep=",",col.names=T)
  
  ### Now accuracy assessment using MAE
  
  method_time
  out_suffix_s <- paste("temp_",method_time[1],"_",out_suffix,sep="")
  

  #undebug(calc_ac_stat_fun)
  #ac_temp_arima_obj <- calc_ac_stat_fun(r_pred_s=temp_pred_rast_arima,
  #                                      r_var_s=r_huric_obs,
  #                                      r_zones=rast_zonal,
  #                                      file_format=file_format,
  #                                      out_suffix=out_suffix_s)
  browser()
  ac_temp_obj <- calc_ac_stat_fun(r_pred_s=temp_pred_rast,
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
  
  mae_tot_tb <- cbind(ac_spat_mle_eigen_no_previous_obj$mae_tb,
                       ac_spat_mle_eigen_with_previous_obj$mae_tb,
                       ac_temp_obj$mae_tb)
  
  name_method_time <- paste0("temp_",method_time[1])
  
  mae_tot_tb <- as.data.frame(mae_tot_tb)
  row.names(mae_tot_tb) <- NULL
  names(mae_tot_tb)<- c("spat_reg_no_previous","spat_reg_with_previous",name_method_time)#,"temp_lm")
  #mae_tot_tb$time <- 2:nrow(mae_tot_tb)
  #mae_tot_tb$time <- 2:n_pred
  mae_tot_tb$time <- 1:nlayers(r_huric_obs)
  y_range<- range(cbind(mae_tot_tb$spat_reg_no_previous,mae_tot_tb$spat_reg_with_previous,mae_tot_tb[[name_method_time]]))
  
  res_pix<-960
  col_mfrow<-1
  row_mfrow<-1
  png(filename=paste("Figure_temporal_profiles_time_1_",out_suffix,".png",sep=""),
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  temp_formula_str <- paste0(paste0("temp_",method_time[1])," ~ ","time")
  plot(spat_reg_no_previous ~ time, type="b",col="cyan",data=mae_tot_tb,ylim=y_range)
  lines(as.formula(temp_formula_str), type="b",col="magenta",data=mae_tot_tb)
  lines(spat_reg_with_previous ~ time, type="b",col="blue",data=mae_tot_tb)
  legend("topleft",
         legend=c("spat_no",name_method_time,"spat_with"),
         col=c("cyan","magenta","blue"),
         lty=1,
         cex=0.8)
  title("Overall MAE for spatial and temporal models") #Note that the results are different than for ARIMA!!!
  
  dev.off()
  
  write.table(mae_tot_tb,file=paste("mae_tot_tb","_",out_suffix,".txt",sep=""))
  
  #### BY ZONES ASSESSMENT: Ok now it is general so it should be part of the function...
  
  #mae_zones_tb <- rbind(ac_spat_mle_obj$mae_zones_tb[1:3,],
  #                      ac_temp_obj$mae_zones_tb[1:3,])
  mae_zones_tb <- rbind(ac_spat_mle_eigen_no_previous_obj$mae_zones_tb,
                        ac_spat_mle_eigen_with_previous_obj$mae_zones_tb,
                        ac_temp_obj$mae_zones_tb)
  
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
  
  p <- xyplot(data~which | as.factor(zones) ,group=method,data=dd,type="b",xlab="time",ylab="VAR",
              strip = strip.custom(factor.levels=unique(as.character(dd$zones))), #fix this!!!
              auto.key = list("topright", corner = c(0,1),# col=c("black","red"),
                              border = FALSE, lines = TRUE,cex=1.2))
  try(print(p))
  
  dev.off()
  #histogram(r_zonal)
  
  ############## Prepare return object #############
  
  space_and_time_prediction_obj <-  list(filename_dat_out,s_raster,mae_tot_tb,mae_zones_tb)
  names(space_and_time_prediction_obj) <- c("filename_dat_out","s_raster","mae_tot_tb","mae_zones_tb")
  return(space_and_time_prediction_obj)
  
}
