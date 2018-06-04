#######################################    Space Time Analyses Project   #######################################
############################  Running Space and Time Model for Space Beats Time framework #######################################
#This script contains functions to produce predictions for Space Beats Time Framework.       
#The script uses functions described additional sript for the project.
#Temporal predictions use OLS with the image of the previous time step or ARIMA.
#Spatial predictions use spatial regression (lag error model) with different estimation methods (e.g. eigen, chebyshev etc.).
#A model with space and time is implemented using neighbours from the previous time step.
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 06/23/2017 
#DATE MODIFIED: 06/04/2018
#Version: 1
#PROJECT: GLP Conference Berlin,YUCATAN CASE STUDY with Marco Millones            
#PROJECT: Workshop for William and Mary: an intro to spatial regression with R 
#PROJECT: Geocomputation and AAG 2015
#PROJECT: Space beats time paper

#TO DO:
# Simplify and clean up code
#
#COMMIT: fixing aggregation for categorical variable
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
library(sf)

###### Functions used in this script

run_space_and_time_models <- function(s_raster,n_time_event,time_window_selected,
                                      method_space=c("mle","eigen"),
                                      method_time=c("arima","arima",T),
                                      NA_flag_val=-9999,
                                      file_format=".tif",
                                      rast_ref=NULL,
                                      zonal_colnames,
                                      num_cores=1,
                                      out_dir=".", 
                                      out_suffix=NULL){
  #
  #This function  runs space and time models given several inputs.
  #Default models are ARIMA and spatial lag models.
  #Spatial models included: lagsarlm (spdep package) with "mle" and 
  # estimation_method: "eigen", "LU", "Chebyshev"
  # Other spatial model method include: gmm and lm.
  # Note that the lm method generates a lag variable and uses OLS. This is not
  # advised as there will be a biased and different standard errors in the estimates.
  #Temporal models included: ARIMA and LM/OLS.
  #
  #AUTHORS: Benoit Parmentier
  #CREATED: 06/23/2017
  #MODIFIED: 06/04/2018
  #
  ##INPUTS
  #1) n_time_event: time step number of the event
  #2) time_window_selected: subset window used to run the models
  #3) method_space: method and algorithm used for spatial models
  #                  defaults are c("mle","eigen"): maximum likelihood estimator and eigen
  #4) method_time: method and algorithm used for temporal models
  #                 defauts are c("arima","arima",T):
  #                 - arima model
  #                 - algorithm arima
  #                 - reinitialization is TRUE
  #5) NA_flag_val: No data value, default is -9999
  #6) file_format: output raster file format, default is ".tif",
  #7) rast_ref: reference raster used as mask, default is NULL
  #             then the first image of the time series is used
  #8) zonal_colnames: variable name for the zonal/strata variable
  #9) num_cores: number of cores used, with default value 1
  #10) out_dir: output directory used, default is to use the current working dir "." 
  #11) out_suffix: output suffix added to the names of all outputs
  #                default value is "NULL"
  #
  ##OUTPUTS
  # Output is a list with the following items":
  #1) filename_dat_out:
  #2) list_output_filename:
  #3) s_raster:
  #4) r_temp_pred
  #5) r_spat_pred_no_previous
  #6) r_spat_pred_with_previous
  #7) res_spat_s_no_previous
  #8) res_spat_s_with_previous
  
  ####################################################
  
  ###### Start script ######
  
  if(is.null(out_suffix)){
    out_suffix=""
  }
  
  projection(s_raster) <- proj_str
  
  ###########################################################################
  ############## PART I PREDICT MODELS  SPAT REGRESSION OVER MULTIPLE time steps ####
  
  ### Important:
  #As input!!
  #r_spat_var should contain one image before the on being predicted, i.e. to predict step 100, you need at least 
  #a stack of raster 99 and raster 100.
  #num_cores <- 4
  time_step_start <- time_window_selected[1]
  time_step_end <- time_window_selected[length(time_window_selected)]
  
  # this is because we miss the first date of pred!!
  time_step_subset <- time_step_start + 1 
  time_window_predicted <- time_step_subset:time_step_end #100 to 116
  #### Subset raster stack for spatial model:
  r_spat_var <- subset(s_raster,
                       time_window_selected) #predict before and after event
  
  ##### Use a raster reference image to define the study area:
  if(is.null(rast_ref)){
    rast_ref <- subset(s_raster,1) #use default first image
  }else{
    rast_ref <- raster(rast_ref)
  }
  
  ########
  ### TEST SPATIAL Prediction on one date....
  
  ### CLEAN OUT AND SCREEN NA and list of neighbours
  #Let's use our function we created to clean out neighbours
  r_clip <- rast_ref #this is the image defining the study area
  
  #############################################
  ### MLE EIGEN Using previous step first
  
  list_models <- NULL

  estimator <- method_space[1] #spatial method used for space model
  estimation_method <- method_space[2] # algorithm used for the space model
  
  ## 
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
  #browser()
  #test <- predict_spat_reg_fun(1,list_param=list_param_spat_reg)
  
  
  pred_spat_mle_eigen_with_previous  <- mclapply(1:n_pred,
                                                 FUN=predict_spat_reg_fun,
                                                 list_param=list_param_spat_reg,
                                                 mc.preschedule=FALSE,
                                                 mc.cores = num_cores)
  
  save(pred_spat_mle_eigen_with_previous,file=file.path(out_dir,paste("pred_spat_",estimator,"_",estimation_method,
                                                                      "_","with_previous_",out_suffix,".RData",sep="")))
  
  #pred_spat_mle_eigen: extract raster images from object
  #r_spat_pred_with_previous <- stack(lapply(pred_spat_mle_eigen_with_previous,FUN=function(x){x$raster_pred})) #get stack of predicted images
  #r_spat_res_with_previous <- stack(lapply(pred_spat_mle_eigen_with_previous,FUN=function(x){x$raster_res})) #get stack of predicted images
  #levelplot(spat_pred_rast_mle_eigen_with_previous,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
  #levelplot(spat_res_rast_mle_eigen_with_previous,col.regions=matlab.like(25)) #view the four predictions using mle spatial reg.
  
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
  ## Generate image raster stack in PART IV
  #r_spat_pred_no_previous <- stack(lapply(pred_spat_mle_eigen_no_previous,FUN=function(x){x$raster_pred})) #get stack of predicted images
  #r_spat_res_no_previous <- stack(lapply(pred_spat_mle_eigen_no_previous,FUN=function(x){x$raster_res})) #get stack of predicted images
  
  
  ##############################################################################################
  ############## PART II TEMPORAL METHODS: PREDICT MODELS  FOR USING TEMP REGRESSION OVER MULTIPLE time steps ####
  
  estimator <- method_time[1] #temporal method used
  estimation_method <- method_time[2] #algorithm used for the temporal method
  re_initialize_arima <- method_time[3] # set to TRUE or FALSE for ARIMA method
  
  #rast_ref <- subset(s_raster,1) #first image ID
  r_stack <- mask(s_raster,rast_ref)
  
  r_temp_var <- subset(r_stack,time_window_selected) # relies on the previous date in contrast to spat reg, this is r_var param
  #r_temp_var <- subset(s_raster,time_window_selected) # relies on the previous date in contrast to spat reg, this is r_var param
  list_models <-NULL
  num_cores_tmp <- num_cores
  
  #browser()
  ###############################
  ### Predict using OLS 
  if(estimator=="lm"){
    #estimator <- "lm"
    #estimation_method <-"ols" #add other options later
    
    #Note the first date is dropped!!!
    out_suffix_s <- paste("t_",time_window_predicted,"_",out_suffix,sep="")#this should really be automated!!!
    
    #ARIMA specific
    #time_step <- time_step_start - 1
    time_step <- time_step_start + 1 # We are
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
  #browser()
  
  if(estimator=="arima"){

    out_suffix_s <- paste("t_",time_window_predicted,"_",out_suffix,sep="")#this should really be automated!!!
    
    #ARIMA specific
    #num_cores_tmp <- num_cores
    #time_step <- n_time_event - 8 #this is the time step for which to start the arima model with, start at 99
    time_step <- time_step_start + 1 # We are
    
    n_pred_ahead <- nlayers(r_temp_var) -1 
    
    arima_order <- NULL
    r_clip_tmp <- rast_ref
    
    #list_param_temp_reg <- list(out_dir,r_temp_var,r_clip_tmp,proj_str,list_models,out_suffix_s,file_format,estimator,estimation_method,
    #                            num_cores_tmp,time_step,n_pred_ahead,r_stack,arima_order,NA_flag_val)
    #names(list_param_temp_reg) <- c("out_dir","r_var","r_clip","proj_str","list_models","out_suffix_s","file_format","estimator","estimation_method",
    #                                "num_cores","time_step","n_pred_ahead","r_stack","arima_order","NA_flag_val")
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
        
        #fix the output dir
        list_param_temp_reg <- list(out_dir,r_temp_var,r_clip_tmp,proj_str,list_models,
                                    out_suffix_s[i],file_format,estimator,estimation_method,
                                    num_cores_tmp,time_step,n_pred_ahead,
                                    r_stack,arima_order,NA_flag_val)
        names(list_param_temp_reg) <- c("out_dir","r_var","r_clip","proj_str","list_models",
                                        "out_suffix_s","file_format","estimator","estimation_method",
                                        "num_cores","time_step","n_pred_ahead","r_stack","arima_order","NA_flag_val")
        #undebug(predict_temp_reg_fun)
        pred_temp_arima <- predict_temp_reg_fun(i,list_param_temp_reg) #only one date predicted...four step ahead
        l_pred_temp_arima[[i]] <- pred_temp_arima  #only one date predicted...one step ahead
        
        #r_temp_pred <- stack(unlist(lapply(1:length(l_pred_temp_arima),FUN=function(i,x){obj<-x[[i]];obj$raster_pred},x=l_pred_temp_arima)))
        #r_temp_res <- stack(unlist(lapply(1:length(l_pred_temp_arima),FUN=function(i,x){obj<-x[[i]];obj$raster_res},x=l_pred_temp_arima)))
        
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
      
      
      #r_temp_pred <- stack(unlist(lapply(1:length(l_pred_temp_arima),FUN=function(i,x){obj<-x[[i]];obj$raster_pred},x=l_pred_temp_arima)))
      #r_temp_res <- stack(unlist(lapply(1:length(l_pred_temp_arima),FUN=function(i,x){obj<-x[[i]];obj$raster_res},x=l_pred_temp_arima)))
      
    }
    
    
  }
  
  ############ PART V COMPARE MODELS IN PREDICTION ACCURACY #################
  #browser()
  #using 11 cores
  #pred_temp_arima <- lapply(1:n_pred,FUN=predict_temp_reg_fun,list_param=list_param_temp_reg)
  
  if(method_time[1]=="arima"){
    r_temp_pred_rast_arima <- stack(unlist(lapply(1:length(l_pred_temp_arima),FUN=function(i,x){obj<-x[[i]];obj$raster_pred},x=l_pred_temp_arima)))
    r_temp_res_rast_arima <- stack(unlist(lapply(1:length(l_pred_temp_arima),FUN=function(i,x){obj<-x[[i]];obj$raster_res},x=l_pred_temp_arima)))
    #note we need to subset!!
    r_temp_pred <- subset(r_temp_pred_rast_arima,1:length(time_window_predicted))
  }
  
  if(method_time[1]=="lm"){
    r_temp_pred_rast_lm <- stack(lapply(pred_temp_lm,FUN=function(x){x$raster_pred}))
    r_temp_res_rast_lm <- stack(lapply(pred_temp_lm,FUN=function(x){x$raster_res}))
    r_temp_pred <- subset(r_temp_pred_rast_lm,1:length(time_window_predicted))
    
  }
  
  if(method_space[1]=="mle"){
   # note that estimator may vary for maximum likelihood method
    
    #pred_spat_mle_eigen: extract raster images from object
    #spat_pred_rast_mle_eigen_no_previous <- stack(lapply(pred_spat_mle_eigen_no_previous,FUN=function(x){x$raster_pred})) #get stack of predicted images
    #spat_res_rast_mle_eigen_no_previous <- stack(lapply(pred_spat_mle_eigen_no_previous,FUN=function(x){x$raster_res})) #get stack of predicted images
    
    r_spat_pred_no_previous <- stack(lapply(pred_spat_mle_eigen_no_previous,FUN=function(x){x$raster_pred})) #get stack of predicted images
    r_spat_res_no_previous <- stack(lapply(pred_spat_mle_eigen_no_previous,FUN=function(x){x$raster_res})) #get stack of predicted images
   
    r_spat_pred_with_previous <- stack(lapply(pred_spat_mle_eigen_with_previous,FUN=function(x){x$raster_pred})) #get stack of predicted images
    r_spat_res_with_previous <- stack(lapply(pred_spat_mle_eigen_with_previous,FUN=function(x){x$raster_res})) #get stack of predicted images
    
  }
  
  #### Now subset observed data:
  r_obs <- subset(s_raster,time_window_predicted) #remove 99 because not considred in the prediction!!
  
  # Problem with rast_zonal: when aggregated categorical value become continuous

  projection(r_spat_pred_with_previous) <- proj_str
  projection(r_spat_pred_no_previous) <- proj_str
  projection(r_temp_pred) <- proj_str
  projection(r_obs) <- proj_str
  
  #### Make quick pannel plots?
  #levelplot(r_obs,col.regions=matlab.like(25))

  #### Step 2: compute residuals
  #browser()
  #this is better: use this to not store temporary files
  #res_temp_s <- overlay(r_temp_pred,r_obs,fun=function(x,y){x-y})
  res_temp_s <- r_temp_pred - r_obs
  res_spat_s_no_previous <- r_spat_pred_no_previous - r_obs
  res_spat_s_with_previous <- r_spat_pred_with_previous - r_obs
  
  ### Assign names for new layers
  names(res_temp_s) <- sub("pred","res",names(res_temp_s))
  names(res_spat_s_no_previous) <- sub("pred","res",names(res_spat_s_no_previous))
  names(res_spat_s_with_previous) <- sub("pred","res",names(res_spat_s_with_previous))
  
  #writeRaster(res_spat_s_with_previous,file=names(res_spat_s_with_previous),file_format)
  
  writeRaster(res_temp_s,
              filename=paste0(file.path(out_dir,names(res_temp_s)),file_format),
              bylayer=T,
              options=c("COMPRESS=LZW"),
              overwrite=TRUE)

  writeRaster(res_spat_s_no_previous,
              filename=paste0(file.path(out_dir,names(res_spat_s_no_previous)),file_format),
              bylayer=T,
              options=c("COMPRESS=LZW"),
              overwrite=TRUE)
  
  writeRaster(res_spat_s_with_previous,
              filename=paste0(file.path(out_dir,names(res_spat_s_with_previous)),file_format),
              bylayer=T,
              options=c("COMPRESS=LZW"),
              overwrite=TRUE)
  
  ### Compute residuals by zones
  rast_zonal <- subset(s_raster,match(zonal_colnames,names(s_raster)))

  r_results <- stack(s_raster,
                     r_temp_pred,
                     r_spat_pred_no_previous,
                     r_spat_pred_with_previous,
                     res_temp_s,res_spat_s_no_previous,
                     res_spat_s_with_previous)
  
  dat_out <- as.data.frame(r_results)
  dat_out <- na.omit(dat_out)
  
  filename_dat_out <- file.path(out_dir,paste("dat_out_",out_suffix,".txt",sep=""))
  write.table(dat_out,file=filename_dat_out,
              row.names=F,sep=",",col.names=T)
  
  #browser()

  #### getfilenames for all outputs
  ## Also write out list of files with path for each!
  s_raster_filenames <- file.path(out_dir,paste0(names(s_raster),file_format))
  r_temp_pred_filenames <- file.path(out_dir,paste0(names(r_temp_pred),file_format))
  r_spat_pred_no_previous_filenames <- file.path(out_dir,paste0(names(r_spat_pred_no_previous),file_format))
  r_spat_pred_with_previous_filenames <- file.path(out_dir,paste0(names(r_spat_pred_with_previous),file_format))
  res_spat_s_no_previous_filenames <- file.path(out_dir,paste0(names(res_spat_s_no_previous),file_format))
  res_spat_s_with_previous_filenames <- file.path(out_dir,paste0(names(res_spat_s_with_previous),file_format))
  
  ## Write out separe and common files of outputs
  write.table(s_raster_filenames,
              file= paste0("list_s_raster_files_",out_suffix,".txt"))
  write.table(r_temp_pred_filenames,
              file= paste0("list_r_temp_pred_files_",out_suffix,".txt"))
  write.table(r_spat_pred_no_previous_filenames,
              file= paste0("list_r_spat_pred_no_previous_files_",out_suffix,".txt"))
  write.table(r_spat_pred_with_previous_filenames,
              file= paste0("list_r_spat_pred_with_previous_files_",out_suffix,".txt"))
  write.table(res_spat_s_no_previous_filenames,
              file= paste0("list_res_spat_s_no_previous_files_",out_suffix,".txt"))
  write.table(res_spat_s_with_previous_filenames,
              file= paste0("list_res_spat_s_with_previous_files_",out_suffix,".txt"))
  
  #### now all lists combined

  list_files_df <- data.frame(
  s_raster=paste0("list_s_raster_files_",out_suffix,".txt"), 
  r_temp_pred=paste0("list_r_temp_pred_files_",out_suffix,".txt"),
  r_spat_pred_no_previous=paste0("list_r_spat_pred_no_previous_files_",out_suffix,".txt"),
  r_spat_pred_with_previous=paste0("list_r_spat_pred_with_previous_files_",out_suffix,".txt"),
  res_spat_s_no_previous=paste0("list_res_spat_s_no_previous_files_",out_suffix,".txt"),
  res_spat_s_with_previous=paste0("list_res_spat_s_with_previous_files_",out_suffix,".txt"))

  list_output_files_df <- data.frame(
    name=c("s_raster","r_temp_pred","r_spat_pred_no_previous",
      "r_spat_pred_with_previous","res_spat_s_no_previous",
      "res_spat_s_with_previous"),
    file=c(file.path(out_dir,paste0("list_s_raster_files_",out_suffix,".txt")), 
    file.path(out_dir,paste0("list_r_temp_pred_files_",out_suffix,".txt")),
    file.path(out_dir,paste0("list_r_spat_pred_no_previous_files_",out_suffix,".txt")),
    file.path(out_dir,paste0("list_r_spat_pred_with_previous_files_",out_suffix,".txt")),
    file.path(out_dir,paste0("list_res_spat_s_no_previous_files_",out_suffix,".txt")),
    file.path(out_dir,paste0("list_res_spat_s_with_previous_files_",out_suffix,".txt")))
  )
  
  #write.table()
  head(list_output_files_df)
  #View(list_files_df)
  #paste(method_space,collpase="_")
  ### Need to add column with type and method used

  list_output_filename <- file.path(out_dir,paste0("list_output_files_df_",out_suffix,".txt"))
  write.table(list_output_files_df,
              file= list_output_filename)
  
  ############## Prepare return object #############

  #space_and_time_prediction_obj <-  list(filename_dat_out,s_raster,mae_tot_tb,mae_zones_tb)
  #names(space_and_time_prediction_obj) <- c("filename_dat_out","s_raster","mae_tot_tb","mae_zones_tb")
  
  space_and_time_prediction_obj <-  list(filename_dat_out,
                                         list_output_filename,
                                         s_raster_filenames,
                                         r_temp_pred_filenames,
                                         r_spat_pred_no_previous_filenames,
                                         r_spat_pred_with_previous_filenames,
                                         res_spat_s_no_previous_filenames,
                                         res_spat_s_with_previous_filenames)
                                         
  names(space_and_time_prediction_obj) <- c("filename_dat_out",
                                             "list_output_filename",
                                             "s_raster",
                                             "r_temp_pred",
                                             "r_spat_pred_no_previous",
                                             "r_spat_pred_with_previous",
                                             "res_spat_s_no_previous",
                                             "res_spat_s_with_previous")
  
  #names(space_and_time_prediction_obj) <- c("filename_dat_out","s_raster","mae_tot_tb","mae_zones_tb")
  
  save(space_and_time_prediction_obj,
       file=paste("space_and_time_prediction_obj_",out_suffix,".RData",sep=""))             
  
  return(space_and_time_prediction_obj)
  
}
