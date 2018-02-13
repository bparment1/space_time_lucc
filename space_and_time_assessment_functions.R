#######################################    Space Time Analyses Project   #######################################
############################ Assessing Space and Time Model for Space Beats Time framework #######################################
#This script contains functions to assess predictions for Space Beats Time Framework.       
#Temporal predictions use OLS with the image of the previous time step or ARIMA.
#Spatial predictions use spatial regression (lag error model) with different estimation methods (e.g. eigen, chebyshev etc.).
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/07/2017 
#DATE MODIFIED: 02/13/2017
#Version: 1

#PROJECT: Space beats time Framework
#TO DO:
#
#COMMIT: splitting function and main script for the assessment
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

## Function to mosaic modis or other raster images

mosaic_m_raster_list<-function(j,list_param){
  #This functions returns a subset of tiles from the modis grid.
  #Arguments: modies grid tile,list of tiles
  #Output: spatial grid data frame of the subset of tiles
  #Note that rasters are assumed to be in the same projection system!!
  
  #rast_list<-vector("list",length(mosaic_list))
  #for (i in 1:length(mosaic_list)){  
  # read the individual rasters into a list of RasterLayer objects
  # this may be changed so that it is not read in the memory!!!
  
  #parse output...
  
  #j<-list_param$j
  mosaic_list<-list_param$mosaic_list
  out_path<-list_param$out_path
  out_names<-list_param$out_rastnames
  file_format <- list_param$file_format
  NA_flag_val <- list_param$NA_flag_val
  ## Start
  
  input.rasters <- lapply(as.character(mosaic_list[[j]]), raster)
  mosaiced_rast<-input.rasters[[1]]
  
  for (k in 2:length(input.rasters)){
    mosaiced_rast<-mosaic(mosaiced_rast,input.rasters[[k]], fun=mean)
    #mosaiced_rast<-mosaic(mosaiced_rast,raster(input.rasters[[k]]), fun=mean)
  }
  
  data_name<-paste("mosaiced_",sep="") #can add more later...
  #raster_name<-paste(data_name,out_names[j],".tif", sep="")
  raster_name<-paste(data_name,out_names[j],file_format, sep="")
  
  writeRaster(mosaiced_rast, NAflag=NA_flag_val,filename=file.path(out_path,raster_name),overwrite=TRUE)  
  #Writing the data in a raster file format...  
  rast_list<-file.path(out_path,raster_name)
  
  ## The Raster and rgdal packages write temporary files on the disk when memory is an issue. This can potential build up
  ## in long  loops and can fill up hard drives resulting in errors. The following  sections removes these files 
  ## as they are created in the loop. This code section  can be transformed into a "clean-up function later on
  ## Start remove
  tempfiles<-list.files(tempdir(),full.names=T) #GDAL transient files are not removed
  files_to_remove<-grep(out_suffix,tempfiles,value=T) #list files to remove
  if(length(files_to_remove)>0){
    file.remove(files_to_remove)
  }
  #now remove temp files from raster package located in rasterTmpDir
  removeTmpFiles(h=0) #did not work if h is not set to 0
  ## end of remove section
  
  return(rast_list)
}

calc_ac_stat_fun <- function(r_pred_s,r_var_s,r_zones,file_format=".tif",out_suffix){
  #Purpose: Calculate accuracy statistics for given regions/zones of the study area
  #Statistics are MAE (Mean Absolute Error) and RMSE(Root Mean Square Error)
  #Parameters:
  #Input:
  #r_pred_s: raster stack of layers predictions (for example predicted NDVI)
  #r_var_s: raster stack of layers actual values (for example observed NDVI)
  #r_zones: raster defining zones of relevance to accuracy (e.g.hurricane winds zones)
  #Output:
  #
  #
  #
  
  ##Functions used
  rmse_fun <-function(x){sqrt(mean(x^2,na.rm=TRUE))}
  mae_fun <-function(x){mean(abs(x),na.rm=TRUE)}
  #rmse_fun <-function(x){sqrt(mean(x^2,na.rm=TRUE))}  
  sd_rmse_fun <-function(x){(sd(x^2,na.rm=TRUE))}  
  #mae_fun <-function(x){sd(abs(x),na.rm=TRUE)}
  sd_mae_fun <- function(x){sd(abs(x))} #sd Absolute Error give a residuals vector
  
  ###
  
  ##Start script
  
  #Accuracy/errors by zones
  r_res_s <- r_pred_s - r_var_s #residuals, stack of raster layers
  mse_zones_tb <- zonal(r_res_s^2,r_zones,fun="mean") #mean square error
  mae_zones_tb <- zonal(abs(r_res_s),r_zones,fun="mean") #absolute error
  sd_mae_zones_tb <- zonal(abs(r_res_s),r_zones,fun="sd") #absolute error
  sd_rmse_tb <- zonal(r_res_s^2,r_zones,fun="sd") #
  #sd_mse_zones_tb <- zonal(r_res_s^2,r_zones,fun="sd") #mean square error
  
  rmse_zones_tb <- cbind(mse_zones_tb[,1],
                         sqrt(mse_zones_tb[,2:dim(mse_zones_tb)[2]])) #root mean square error
  if(!is.null(colnames(rmse_zones_tb)[1])){
    colnames(rmse_zones_tb)[1] <- c("zone")
  }
  
  #Overall Accuracy/errors 
  
  mae_tb <- cellStats(abs(r_res_s),mean) #calculate MAE for layer stack
  rmse_tb <- sqrt(cellStats((r_res_s)^2,mean)) #calculate rmse overall
  sd_mae_tb <- cellStats(abs(r_res_s),sd) #sd for absolute error
  
  #write out residuals rasters
  r_res_s <- writeRaster(r_res_s,filename=file.path(out_dir,"r_res_s.tif"),bylayer=TRUE,
                         suffix=paste(1:nlayers(r_res_s),out_suffix,sep="_"),overwrite=TRUE)
  ac_obj <- list(mae_tb,rmse_tb,sd_mae_tb,mae_zones_tb,rmse_zones_tb,sd_mae_zones_tb,sd_rmse_tb)
  names(ac_obj) <- c("mae_tb","rmse_tb","sd_mae_tb","mae_zones_tb","rmse_zones_tb","sd_mae_zones_tb","sd_rmse_tb")
  
  return(ac_obj)
}


accuracy_space_time_calc <- function(r_temp_pred,r_spat_pred,s_raster,proj_str,n_time_event,time_window_selected,
                                     r_zonal,method_space,method_time,r_ref,out_suffix,
                                     var_names,NA_flag_val,file_format,date_range,
                                     out_dir,create_out_dir_param){
  
  ##INPUTS
  #1)r_temp_pred
  #2) r_spat_pred
  #3) s_raster
  #4) proj_str
  #5) n_time
  #6) time_window_selected
  #7) r_zonal
  #8) methods_space
  #9) method_time
  #10) r_ref
  #11) out_suffix
  #12) var_nams
  #13) NA_flag_val
  #14) file_format
  #15) date_range
  #16) out_dir
  #17) out_dir_param
  #
  ##OUTPUTS
  #1)
  #2)
  #
  
  ######## Start script #####
  
  ###  Load data if needed:
  if(class(r_temp_pred)=="character"){
    #/home/parmentier/Data/Space_beats_time/Data/data_Rita_NDVI/rev_project_output/tile_2
    r_temp_pred_df <- read.table(r_temp_pred,stringsAsFactors = F)
    lf <-r_temp_pred_df[,1] 
    r_temp_pred <- stack(lf)
  }
  if(class(r_spat_pred)=="character"){
    r_spat_pred_df <- read.table(r_spat_pred,stringsAsFactors = F)
    lf <-r_spat_pred_df[,1] 
    r_spat_pred <- stack(lf)
  }
  if(class(s_raster)=="character"){
    df_lf <- read.table(s_raster,sep=",",stringsAsFactors = F,header=T) 
    s_raster <- stack(df_lf[,1])
  }
  
  date_range <- (unlist(strsplit(date_range,";")))
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

  #if same length, don't remove
  if(length(time_window_selected)==nlayers(r_temp_pred)){
    r_obs <- subset(s_raster,time_window_selected) #remove 99 because not considred in the prediction!!
  }else{
    r_obs <- subset(s_raster,time_window_selected[-1]) #remove 99 because not considred in the prediction!!
  }
  
  projection(r_obs) <- proj_str
  projection(r_temp_pred) <- proj_str
  projection(r_spat_pred) <- proj_str
  
  #### Step 2: compute residuals
  #browser()
  #this is better: use this to not store temporary files
  #res_temp_s <- overlay(r_temp_pred,r_obs,fun=function(x,y){x-y})
  res_temp_s <- r_temp_pred - r_obs
  res_spat_s <- r_spat_pred - r_obs

  ### Assign names for new layers
  names(res_temp_s) <- sub("pred","res",names(res_temp_s))
  names(res_spat_s) <- sub("pred","res",names(res_spat_s))

  ### Compute residuals by zones
  rast_zonal <- subset(s_raster,match(zonal_colnames,names(s_raster)))
  r_x <- init(rast_zonal,"x")
  r_y <- init(rast_zonal,"y")
  
  ### Generate output in text format
  r_results <- stack(r_x,r_y,s_raster,r_temp_pred,r_spat_pred,
                                          res_temp_s,res_spat_s)
  
  dat_out <- as.data.frame(r_results)
  dim(dat_out)
  #dat_out <- na.omit(dat_out) #don't do this
  
  freq_tb_zonal <- as.data.frame(freq(rast_zonal))
  
  filename_freq_tb_zonal <- file.path(out_dir,paste("freq_tb_zonal_",out_suffix,".txt",sep=""))
  write.table(freq_tb_zonal,
              file=filename_freq_tb_zonal,
              row.names=F,sep=",",col.names=T)
  
  
  filename_dat_out <- file.path(out_dir,paste("dat_out_",out_suffix,".txt",sep=""))
  write.table(dat_out,file=filename_dat_out,
              row.names=F,sep=",",col.names=T)
  
  ### Now accuracy assessment using MAE
  
  out_suffix_s <- paste("temp_",method_time[1],"_",method_time[2],"_",out_suffix,sep="")

  #undebug(calc_ac_stat_fun)
  #Find where that function is!!!
  
  #function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_functions_11072017.R" #PARAM 1

  ac_temp_obj <- calc_ac_stat_fun(r_pred_s=r_temp_pred,
                                  r_var_s=r_obs,
                                  r_zones=rast_zonal,
                                  file_format=file_format,
                                  out_suffix=out_suffix_s)  
  
  out_suffix_s <- paste("spat_",method_space[1],"_",method_space[2],"_",out_suffix,sep="")
  ac_spat_obj <- calc_ac_stat_fun(r_pred_s=r_spat_pred,
                              r_var_s=r_obs,
                              r_zones=rast_zonal,
                              file_format=file_format,
                              out_suffix=out_suffix_s)
  
  mae_tot_tb <- cbind(ac_spat_obj$mae_tb,
                      ac_temp_obj$mae_tb)
  
  name_method_time <- paste0("temp_",method_time[1],"_",method_time[2])
  name_method_space <- paste0("spat_",method_space[1],"_",method_space[2])
  mae_tot_tb <- as.data.frame(mae_tot_tb)
  row.names(mae_tot_tb) <- NULL
  names(mae_tot_tb)<- c(name_method_space,name_method_time)
  
  mae_tot_tb$time <- 1:nlayers(r_obs)
  y_range<- range(cbind(mae_tot_tb[[name_method_space]],mae_tot_tb[[name_method_time]]))
  
  res_pix<-960
  col_mfrow<-1
  row_mfrow<-1
  png(filename=paste("Figure_temporal_profiles_time_1_",out_suffix,".png",sep=""),
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  temp_formula_str <- as.formula(paste0(paste0(name_method_time," ~ ","time")))
  spat_formula_str <- as.formula(paste0(paste0(name_method_space," ~ ","time")))
  
  plot(temp_formula_str,type="b",
       col="cyan",
       data=mae_tot_tb,
       ylim=y_range,
       ylab="MAE")
  lines(spat_formula_str, type="b",col="magenta",data=mae_tot_tb)
  
  legend("topleft",
         legend=c(name_method_time,name_method_space),
         col=c("cyan","magenta"),
         lty=1,
         cex=0.8)
  title("Overall MAE for spatial and temporal models") #Note that the results are different than for ARIMA!!!
  
  dev.off()
  
  write.table(mae_tot_tb,file=paste("mae_tot_tb","_",out_suffix,".txt",sep=""))
  
  #### BY ZONES ASSESSMENT: Ok now it is general so it should be part of the function...
  
  #mae_zones_tb <- rbind(ac_spat_mle_obj$mae_zones_tb[1:3,],
  #                      ac_temp_obj$mae_zones_tb[1:3,])
  #mae_zones_tb <- rbind(ac_spat_mle_eigen_no_previous_obj$mae_zones_tb,
  #                      ac_spat_mle_eigen_with_previous_obj$mae_zones_tb,
  #                      ac_temp_obj$mae_zones_tb)
  mae_zones_tb <- rbind(ac_spat_obj$mae_zones_tb,
                        ac_temp_obj$mae_zones_tb)
  
  mae_zones_tb <- as.data.frame(mae_zones_tb)
  
  n_zones <- length(unique(mae_zones_tb$zone))
  
  #mae_zones_tb$method <- c(rep("spat_reg_no",n_zones),rep("temp_with_reg",n_zones),rep("temp_arima_reg",n_zones))
  mae_zones_tb$method <- c(rep(name_method_space,n_zones),rep(name_method_time,n_zones))
  
  n_time <- ncol(mae_zones_tb) -1
  pred_names <- c("zone",paste("t",2:n_time,sep="_"),"method")
  #names(mae_zones_tb) <- c("zones","pred1","pred2","pred3","pred4","method")
  names(mae_zones_tb) <- pred_names
  
  write.table(mae_zones_tb,file=paste("mae_zones_tb","_",out_suffix,".txt",sep=""))
  
  mydata<- mae_zones_tb
  dd <- do.call(make.groups, mydata[,-ncol(mydata)]) 
  #dd$lag <- mydata$lag 
  #drop first few rows that contain no data but zones...
  n_methods <- 2 #2 because we have 2 methods... (one for space and one for time)
  n_start <-n_zones*2 +1 
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
