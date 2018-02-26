#######################################    Space Time Analyses Project   #######################################
############################ Assessing Space and Time Model for Space Beats Time framework #######################################
#This script contains functions to assess predictions for Space Beats Time Framework.       
#Temporal predictions use OLS with the image of the previous time step or ARIMA.
#Spatial predictions use spatial regression (lag error model) with different estimation methods (e.g. eigen, chebyshev etc.).
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/07/2017 
#DATE MODIFIED: 02/25/2018
#Version: 1

#PROJECT: Space beats time Framework
#TO DO: Add movie generation: one for raster images and one for sbt sequence + time profile
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

format_mae_zones_tb <- function(df_zone){
  #Reformat data from mae_zone computation in easy format
  #
  val_zone_mae <- t(df_zone)[-1,]
  head(val_zone_mae)
  val_mae <- as.vector(val_zone_mae)
  df_zone_mae <- data.frame(mae=val_mae)
  n_times_step <- nrow(df_zone_mae)/2
  n_zones <- ncol(val_zone_mae)
  df_zone_mae$zone <- unlist(lapply(1:n_zones,function(x)(rep(x,n_times_step))))
  #mae_zones_tb <- mae_zones_tb[-1,]
  #names(df_zone_mae_tb) <- paste0("zone_",1:ncol(mae_zones_tb))
  df_zone_mae$time <- 1:n_times_step 
  #mae_zones_tb$time <- 1:n_times_step 
  return(df_zone_mae)
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
  #View(mae_tot_tb)
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
  png(filename=paste("Figure_temporal_profiles_MAE_overall_",out_suffix,".png",sep=""),
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
  

  df_zone <- ac_spat_obj$mae_zones_tb
  mae_zones_tb_spat <- format_mae_zones_tb(df_zone)
  df_zone <- ac_temp_obj$mae_zones_tb
  mae_zones_tb_temp <- format_mae_zones_tb(df_zone)
  
  #View(mae_zones_tb_spat)
  mae_zones_tb_spat$method <- name_method_space
  mae_zones_tb_temp$method <- name_method_time
  
  mae_zones_tb_spat$type <- "spatial"
  mae_zones_tb_temp$type <- "temporal"
  
  mae_zones_tb <- rbind(mae_zones_tb_spat,mae_zones_tb_temp)
  #method_time <- unlist(strsplit(method_time,";"))
  #method_space <- unlist(strsplit(method_space,";"))
  #name_method_time <- paste0("temp_",method_time[1],"_",method_time[2])
  #name_method_space <- paste0("spat_",method_space[1],"_",method_space[2])
  
  #View(mae_zones_tb)

  write.table(mae_zones_tb,
              file=file.path(out_dir,paste("mae_zones_tb","_",out_suffix,".txt",sep="")))
  
  ##### Let's set up the figure production now
  #View(mae_zones_tb)


  ### Now generate plot
  #y_range <- range(cbind(mae_zones_tb[[name_method_space]],mae_zones_tb[[name_method_time]]))
  
  y_range <- range(mae_zones_tb$mae)
  legend_val <- c("tempor model","spatial model")
  
  i <- 1
  list_figures_mae_profiles_png <- vector("list",length=n_zones)
    
  for(i in 1:n_zones){
    
    zone_val <- i
    data_mae <- subset(mae_zones_tb,mae_zones_tb$zone==zone_val)
    data_mae_space <- subset(data_mae,method==name_method_space)
    data_mae_time <- subset(data_mae,method==name_method_time)
    
    #View(data_mae)
    
    res_pix <- 960
    col_mfrow<- 1
    row_mfrow<- 0.7
    
    png_filename <- paste("Figure_temporal_profiles_MAE_zones_",zone_val,"_",out_suffix,".png",sep="")
    png(filename=png_filename,
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    
    formula_str <- as.formula(paste0(paste0("mae"," ~ ","time")))
    #spat_formula_str <- as.formula(paste0(paste0(name_method_space," ~ ","time")))
    
    #par(mgp=c(0,1,0))
    ## margin for side 2 is 7 lines in size
    #op <- par(mar=c(5,7,4,2) +0.1) ## default is c(5,4,4,2) + 0.1
    op <- par(mar=c(5,7,6,2) +0.1) ## default is c(5,4,4,2) + 0.1, 4 is 4+32 for top
    
    plot(formula_str,type="b",
         col="cyan",
         pch=16, #filled in circles...
         lwd=3,
         data= data_mae_time,
         ylim=y_range,
         ylab="Mean Absolute Error (MAE)",
         xaxt="n",
         xlab="Time",
         cex.lab=2,
         cex.axis=1.5,
         cex=1.5)
    lines(formula_str, 
          type="b",
          lwd=3,
          pch=16, #filled in circles
          col="magenta",
          data=data_mae_space)
    
    x_labels <- c("T-5","T-4","T-3","T-2","T-1","T+1","T+2","T+3","T+4","T+5")
    axis(1, at = 1:10, labels = x_labels , cex.axis = 1.5)
    #axis(2, cex.axis = 2)
    
    legend("topleft",
           legend=legend_val,
           col=c("cyan","magenta"),
           lty=1,
           lwd=3,
           cex=1.5,
           bty="n")
    title_str <- paste0("Zone ",zone_val," MAE for spatial and temporal models ")
    title(title_str,cex.main=2.3,font=2)
    #title("Overall MAE for spatial and temporal models",cex.main=2.3,font=2) #Note that the results are different than for ARIMA!!!
    par(op)
    dev.off()
    
    list_figures_mae_profiles_png[[i]] <- png_filename 
  }
  
  #### Add thing to convert using 
  
  #Use ImageMagick
  #convert Figure_4a_agri_hinterland_cat_sum_flow_10182016.png Figure_4b_meat_hinterland_cat_sum_flow_10182016.png Figure_4c_livestock_hinterland_cat_sum_flow_10182016.png -append test.png
  #use convert fig1.png fig2.png fig3.png -append to join vertically
  #use convert fig1.png fig2.png fig3.png +append to join left to right
  
  #paste(unlist(list_figures_mae_profiles_png),collapse=" ")
  
  png_filename_mae_combined <- file.path(out_dir,paste("Figure_3_mae_profiles_by_zones_","_",out_suffix,".png",sep=""))
  png_filename_mae_combined <- paste("Figure_3_mae_profiles_by_zones_","_",out_suffix,".png",sep="")
  
  cmd_str <- paste("convert",
                   #png_filename4a,
                   #png_filename4b,
                   #png_filename4c,
                   paste(unlist(list_figures_mae_profiles_png),collapse=" "),
                   "+append",
                   png_filename_mae_combined,sep=" ")
  system(cmd_str)
  
  #To look up later to change resolution with ImageMagick
  #http://www.imagemagick.org/discourse-server/viewtopic.php?t=18241
  #convert -units PixelsPerInch rose1.png -resample 300 rose2d.jpg
  
  #http://superuser.com/questions/479197/i-want-to-change-dpi-with-imagemagick-without-changing-the-actual-byte-size-of-t
  #convert -units PixelsPerInch input.png -density 300 output.png
  
  ####### Generate animation:
  
  browser()
  
  x_labels <- c("T-5","T-4","T-3","T-2","T-1","T+1","T+2","T+3","T+4","T+5")
  
  i<-1
  l_dates <- x_labels
  r_mosaiced_scaled <- r_obs
  NA_flag_val
  region_name <- "RITA"
  variable_name <- "NDVI"
  zlim_val <- NULL
  stat_opt <- T
  
  list_param <- list(i,l_dates,r_mosaiced_scaled,NA_flag_val, 
  out_dir,out_suffix,region_name,variable_name,zlim_val,stat_opt)
  
  names(list_param) <- c("i","l_dates","r_mosaiced_scaled","NA_flag_val", 
                         "out_dir","out_suffix","region_name","variable_name","zlim_val","stat_opt")
  
  debug(plot_raster_mosaic)
  test <- plot_raster_mosaic(i,list_param)

  #list_figures <- lapply(1:length(l_dates),
  #       FUN=plot_raster_mosaic,
  #       list_param=list_param)
  
  list_plot_fig_obj <- mclapply(1:length(l_dates),
                                  FUN = plot_raster_mosaic,
                                  list_param = list_param,
                                  mc.preschedule = FALSE,
                                  mc.cores = num_cores)
  #list_plot_fig_obj[[1]]$min_max_df
  #list_plot_fig_obj[[1]]$stat_df
  #list_plot_fig_obj[[1]]$png_filename
  
  lf_plot_fig <- lapply(list_plot_fig_obj,function(x){x$png_filename})
  lf_raster <- filename(r_obs)
  
  out_suffix_str <- out_suffix #need to change this
  if(stat_opt==TRUE){
    l_stat_df <- lapply(list_plot_fig_obj,function(x){x$stat_df})
    stat_df <- do.call(rbind,l_stat_df)
    stat_df$date <- l_dates
    stat_df$files <- lf_raster
    
    ### Write out information
    stat_df_fname <- file.path(out_dir, paste0("stat_df_", out_suffix_str, ".txt"))
    write.table(stat_df,stat_df_fname,sep=",",row.names = F)
    
  }else{
    stat_df <- NULL
  }
  l_min_max_df <- lapply(list_plot_fig_obj,function(x){x$stat_df})
  min_max_df <- do.call(rbind,l_min_max_df)
  min_max_df$date <- l_dates
  min_max_df$files <- lf_raster
  ### Write out information
  min_max_df_fname <- file.path(out_dir, paste0("min_max_df_", out_suffix_str, ".txt"))
  write.table(min_max_df,file = min_max_df_fname,sep = ",",row.names = F)
  
  if(is.null(zlim_val)){
    out_suffix_movie <- paste("min_max_", out_suffix_str, sep = "")
  } else{
    zlim_val_str <- paste(zlim_val, sep = "_", collapse = "_")
    out_suffix_movie <- paste(zlim_val_str, "_", out_suffix, sep = "")
  }
  filenames_figures_mosaic <- paste0("list_figures_animation_", out_suffix_movie, ".txt")
  
  write.table(unlist(lf_plot_fig),
              filenames_figures_mosaic,row.names = F,col.names = F,quote = F)
  
  #frame_speed <- 25
  frame_speed <- 50
  animation_format <- ".mp4"
  
  debug(generate_animation_from_figures_fun)
  out_dir <- "." #problem with file path that includes Google Drive (space), resolves this issue later
  out_filename_figure_animation <- generate_animation_from_figures_fun(filenames_figures = filenames_figures_mosaic,
                                                                       frame_speed = frame_speed,
                                                                       format_file = animation_format,
                                                                       out_suffix = out_suffix_movie,
                                                                       out_dir = out_dir,
                                                                       out_filename_figure_animation = NULL)
  
  
  ####### NOW DO AVERAGE PROFILES ########
  
  #debug(compute_avg_by_zones)
  #test <- compute_avg_by_zones(r_stack,r_zonal,out_suffix_str="",out_dir=".")
    
  
  
  ############## Prepare return object #############
  
  space_and_time_prediction_obj <-  list(filename_dat_out,s_raster,mae_tot_tb,mae_zones_tb)
  names(space_and_time_prediction_obj) <- c("filename_dat_out","s_raster","mae_tot_tb","mae_zones_tb")
  
  return(space_and_time_prediction_obj)
}


#create animation from figures:
generate_animation_from_figures_fun <- function(filenames_figures,frame_speed=50,format_file=".gif",in_dir=NULL,out_suffix="",out_dir=".",out_filename_figure_animation=NULL){
  #This function generates an animation given a list of files or textfile. The default format is .gif.
  #The function requires ImageMagick to produce animation.
  #INPUTS:
  #1) filenames_figures: list of files as "list" or "character, or file name of a text file containing the list of figures.
  #2) frame_speed: delay option in constructing the animation, default is 50,
  #the unit is 1/100th of second so 50 is 2 frame per second
  #3) format_file=".gif", ".mp4" or ".avi
  #4) in_dir : location where figures for animations are stored
  #4) out_suffix: ouput string added as suffix, default is ""
  #5) out_dir: output directory, default is current directory
  #6) out_filename_figure_animation: output filename if NULL (default)
  #                                   then generated from input suffix
  
  #OUTPUTS:
  #1) out_filename_figure_animation: file containing the movie
  
  ######## Beging script #############
  
  if(is.null(out_filename_figure_animation)){
    #out_filename_figure_movies <- file.path(out_dir,paste("mosaic_movie_",out_suffix_movie,".gif",sep=""))
    out_filename_figure_animation <- file.path(out_dir,paste("animation_frame_",frame_speed,"_",out_suffix,format_file,sep=""))
  }
  
  if(class(filenames_figures)=="list" | (length(filenames_figures)>1)){ #if length of object greater than 1 then assume a list of files
    #filename_figures_mosaic <- file.path(out_dir,"mosaic_plot_fig.txt")
    out_filenames_figures <- file.path(out_dir,paste("list_figures_animation_",out_suffix,".txt",sep=""))
    write.table(unlist(filenames_figures),out_filenames_figures,row.names = F,col.names = F,quote = F)
    filenames_figures <- out_filenames_figures
  }
  
  #now generate movie with imageMagick
  
  if(format_file==".gif"){ #file format for the animation
    if(file.exists(out_filename_figure_animation)){
      file.remove(out_filename_figure_animation)
    }
    
    #-delay 20
    #delay_option <- 60
    delay_option <- frame_speed #this might change if format is changed
    
    cmd_str <- paste("convert",
                     paste("-delay",delay_option),
                     paste0("@",filenames_figures),
                     out_filename_figure_animation)
    #convert @myimages.txt mymovie.gif
    #save cmd_str in text file!!!
    
  }
  
  #if format is .mp4
  if(format_file==".mp4" | format_file==".avi"){
    
    if(file.exists(out_filename_figure_animation)){
      file.remove(out_filename_figure_animation)
    }
    
    #ffmpeg -f image2 -pattern_type glob -i '*.png' out.mp4
    #ffmpeg -f image2 -r 1 -pattern_type glob -i '*.png' out.mp4
    
    #-r1 one second by frame
    #rate of one per second:
    #ffmpeg -f image2 -r 1 -pattern_type glob -i '*.png' out.mp4
    #crf: used for compression count rate factor is between 18 to 24, the lowest is the highest quality
    #ffmpeg -f image2 -r 1 -vcodec libx264 -crf 24 -pattern_type glob -i '*.png' out.mp4
    
    #out_fig <- read.table(filenames_figures)
    #out_fig
    if(is.null(in_dir)){
      in_dir <- "./fig_animation" #Assume that the figures are here...
    }
    
    frame_rate <- frame_speed/100 # convert frame rate in second
    cmd_str <- paste("ffmpeg",
                     paste("-f","image2",sep=" "),
                     paste("-r",frame_rate,sep=" "),
                     #paste("-pattern_type glob -i","'*.png'",sep=" "),
                     paste("-pattern_type glob -i","'./fig_animation/*.png'",sep=" "),
                     #paste("-pattern_type glob -i ",filenames_figures),
                     out_filename_figure_animation,sep=" ")
    
  }
  system(cmd_str)
  
  #####
  
  return(out_filename_figure_animation)
  
}


plot_raster_mosaic <- function(i,list_param){
  #Function to plot raster image
  #
  #INPUTS
  #1) l_dates
  #2) r_stack
  #3) NA_flag_val
  #4) out_dir,
  #5) out_suffix_str
  #6) region_name
  #7) variable_name
  #8) zlim_val
  #
  
  ############# Start script #########
  
  l_dates <- list_param$l_dates
  r_mosaiced_scaled <- list_param$r_mosaiced_scaled
  NA_flag_val <- list_param$NA_flag_val
  out_dir <- list_param$out_dir
  out_suffix <- list_param$out_suffix
  region_name <- list_param$region_name
  variable_name <- list_param$variable_name
  zlim_val <- list_param$zlim_val
  stat_opt <- list_param$stat_opt #if TRUE computer stats for the image
  
  #for (i in 1:length(nlayers(r_mosaic_scaled))){
  
  date_proc <- l_dates[i]
  r_pred <- subset(r_mosaiced_scaled,i)
  NAvalue(r_pred)<- NA_flag_val 
  
  raster_name <- filename(r_pred)
  extension_str <- extension(raster_name)
  raster_name_tmp <- gsub(extension_str,"",basename(raster_name))
  
  date_proc <- l_dates[i]
  if(class(date_proc)!="Date"){
    date_val <- as.Date(strptime(date_proc,"%Y%m%d"))
    #month_name <- month.name(date_val)
  }else{
    date_val <- date_proc
  }
  
  if(is.na(date_val)){
    date_str <- date_proc #it means it is a label e.g. T-4
  }else{
    month_str <- format(date_val, "%b") ## Month, char, abbreviated
    year_str <- format(date_val, "%Y") ## Year with century
    day_str <- as.numeric(format(date_val, "%d")) ## numeric month
    date_str <- paste(month_str," ",day_str,", ",year_str,sep="")
  }
  
  res_pix <- 1200
  #res_pix <- 480
  col_mfrow <- 1
  row_mfrow <- 1
  
  
  if(is.null(zlim_val)){
    
    if(is.na(minValue(r_pred))){ #if min is NA
      r_pred <- setMinMax(r_pred) #then set min and max values
    }
    
    min_max_val <- c(minValue(r_pred),maxValue(r_pred))
    min_max_df <- data.frame(min=min_max_val[1],max=min_max_val[2])
    
    zlim_val_str <- paste(min_max_val,sep="_",collapse="_")
    #png_filename <-  file.path(out_dir,paste("Figure4_clim_mosaics_day_","_",date_proc,"_",region_name,"_zlim_",zlim_val_str,"_",out_suffix,".png",sep =""))
    #raster_name_tmp
    png_filename <-  file.path(out_dir,paste("Figure_",raster_name_tmp,"_zlim_",zlim_val_str,"_",out_suffix,".png",sep =""))
    
    title_str <-  paste("Predicted ",variable_name, " on ",date_str , " ", sep = "")
    #browser()
    
    png(filename=png_filename,width = col_mfrow * res_pix,height = row_mfrow * res_pix)
    
    plot(r_pred,main =title_str,cex.main =1.5,col=matlab.like(255),
         legend.shrink=0.8,legend.width=0.8)
    #axis.args = list(cex.axis = 1.6), #control size of legend z
    #legend.args=list(text='dNBR', side=4, line=2.5, cex=2.2))
    #legend.args=list(text='dNBR', side=4, line=2.49, cex=1.6))
    dev.off()
  }else{
    min_max_val <- NULL
    zlim_val_str <- paste(zlim_val,sep="_",collapse="_")
    #png_filename <-  file.path(out_dir,paste("Figure_mosaics_day_","_",date_proc,"_",region_name,"_",zlim_val_str,"_",out_suffix,".png",sep =""))
    
    png_filename <-  file.path(out_dir,paste("Figure_",raster_name_tmp,"_zlim_",zlim_val_str,"_",out_suffix,".png",sep =""))
    
    title_str <-  paste("Predicted ",variable_name, " on ",date_str , " ", sep = "")
    
    png(filename=png_filename,width = col_mfrow * res_pix,height = row_mfrow * res_pix)
    
    plot(r_pred,main =title_str,cex.main =1.5,col=matlab.like(255),
         zlim=zlim_val, #this time we set the range limit for visualization
         legend.shrink=0.8,legend.width=0.8)
    #axis.args = list(cex.axis = 1.6), #control size of legend z
    #legend.args=list(text='dNBR', side=4, line=2.5, cex=2.2))
    #legend.args=list(text='dNBR', side=4, line=2.49, cex=1.6))
    dev.off()
  }
  
  ### Option to compute stats for the raster image
  if(stat_opt==TRUE){
    if(is.null(min_max_val)){
      min_val <- cellStats(r_pred,stat="min",na.rm=T)
      max_val <- cellStats(r_pred,stat="max",na.rm=T)
    }else{
      min_val <- min_max_val[1]
      max_val <- min_max_val[2]
    }
    
    NA_no <- freq(r_pred,value=NA)
    n_cell <- ncell(r_pred)
    perc_NA <- 100*(NA_no/n_cell)
    if(dataType(r_pred)=="INT2S" | dataType(r_pred)=="INT4S"){
      #integer overflow for sum
      #filename(r_pred)
      r_pred_tmp <- r_pred
      raster_name_tmp <- paste(raster_name_tmp,"_tmp",file_format,sep="")
      writeRaster(r_pred_tmp,raster_name_tmp,datatype="FLT8S",overwrite=T) #this very slow
      r_pred_tmp <-raster(raster_name_tmp) #read in the converted image
      #dataType(r_pred_tmp)
      mean_val <- cellStats(r_pred_tmp,stat="mean",na.rm=T)
      sd_val <- cellStats(r_pred_tmp,stat="sd",na.rm=T)
      file.remove(raster_name_tmp)
    }else{
      mean_val <- cellStats(r_pred,stat="mean",na.rm=T)
      sd_val <- cellStats(r_pred,stat="sd",na.rm=T)
    }
    stat_df <- data.frame(min=min_val,
                          max=max_val,
                          sd=sd_val,
                          mean=mean_val,
                          NA_no=NA_no,
                          n_cell=n_cell,
                          perc_NA)
  }else{
    stat_df <- NULL
  }
  
  fig_obj <- list(png_filename,min_max_df,stat_df)
  names(fig_obj) <- c("png_filename","min_max_df","stat_df")
  
  #return(png_filename)
  return(fig_obj)
}
  
###################### END OF SCRIPT ##########################  
  
