####################################    Space Time Analyses Project   #######################################
############################  Yucatan case study: SAR, SARMA etc -PART 2  #######################################
#This script functions produce predictions for the dates following the Hurricane Dean event.       
#The script uses spatial regression with weight matrix to predict NDVI values in the MOORE EDGY region. 
#Temporal predictions use OLS with the image of the previous time step rather than ARIMA.
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/15/2014 
#DATE MODIFIED: 08/24/2017
#Version: 2
#PROJECT: GLP Conference Berlin,YUCATAN CASE STUDY with Marco Millones            
#PROJECT: Workshop for William and Mary: an intro to spatial regression with R 
#PROJECT: Geocomputation and AAG 2015
#PROJECT: Space beats time paper

#TO DO:
# modify the rasterize_df_fun function to allow ref image
# add the ARIMA method to run more efficiently
#
#COMMIT: debugging explore_data function to account for raster stack rather than data table
#
#################################################################################################

#This script currently contains 1 function:

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

#log_file_create_list_param_fun <- function(list_param,function_name,out_file=TRUE){
#  num_parameters <- length(list_param)
#  df_list_param <- as.data.frame(name)
#  for(i in 1:length(list_param)){
#    df_list_param$argument <- names(list_param[[i]])
#   df_list_param$class <- class(list_param[[i]])
#    if(class(list_param[[i]])==char){
#     df_list_param$value <- list_param[[i]]
#    }else{
#      df_list_param$value <- NA
#    }
#      
#    out_dir  <- list_param$out_dir
#    r_ref_s    <- list_param$r_var #if NULL, no image is created, this is the reference image
#    #list_param$ <- rast_ref
#    r_clip     <- list_param$r_clip
#    proj_str <- list_param$proj_str
#    list_models <- list_param$list_models
#    file_format <- list_param$file_format
#    estimator <- list_param$estimator
#    estimation_method <- list_param$estimation_method #currently used only for mle from errorsarlm
#    NA_flag_val <- list_param$NA_flag_val
#    #ARIMA specific
#    num_cores <- list_param$num_cores #paraallelization in space.this should be done by row or til enot by pixel!!!!
#    time_step <- list_param$time_step #this is the time step for which to start the arima model with
#    n_pred_ahead <- list_param$n_pred_ahead
#    r_stack <- list_param$r_stack
#    arima_order <- list_param$arima_order
#
#}
#    return()
#  }
#}

### Make a function to uniformize NA accros given dates!!!

explore_and_summarize_data <- function(l_rast,zonal_colnames,var_names,n_time_event,proj_str,pixel_index=800,out_dir=NULL,out_suffix=NULL){
  ## Function that generates a set of figures and summary of basics information on data provided for SBT.
  # This include time profile, raster time steps and before after event figures.
  #
  #INPUTS
  #1) l_rast: list of raster time series and variables (zones, etc.)
  #2) zonal_colnames: name of the column matching the zone
  #3) var_names: index vector matching the time series
  #4) n_time_event: time step matching the event studied
  #5) proj_str: projection string in PROJ4 format
  #6) pixel_index: pixel index corresponding to time series to be plotted
  #7) out_dir: if null use the current directory
  #8) out_suffix: if null use "" (maybe process id later on)
  #OUTPUTS
  ## ojbect with three items:
  #1) zones_avg: data.frame with average by zone
  #2) zones_vag_long_tb: data.frame with average by zone in long format
  #3) out_files: list of outputs files generated (mostly figures)

  ############# BEGIN SCRIPT ###########
  
  if(is.null(out_dir)){
    out_dir = "."
  }
  
  if(is.null(out_suffix)){
    out_suffix <- ""
  }
  
  s_raster <- stack(l_rast) #stack with all the variables
  projection(s_raster) <- proj_str
  #names(s_raster) <- names(data_tb)                
  r_FID <- subset(s_raster,1) #Assumes ID or reference image is the first image of the stack
  
  list_fig_filename <- vector("list",length=8)
  
  #############
  ##Figure 1: reference layer
  
  out_fig_filename <- file.path(out_dir,paste("Figure1_ref_layer_time_1_",out_suffix,".png",sep=""))
  
  res_pix<-960
  col_mfrow<-1
  row_mfrow<-1
  png(filename= out_fig_filename,
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  plot(r_FID,main=l_rast[1])
  dev.off()
  
  list_fig_filename[[1]] <- out_fig_filename
  
  freq_tb <- (freq(r_FID))
  #writeRaster()
  
  ###############
  ##Figure 2: zonal layer
  
  out_fig_filename <- file.path(out_dir,paste("Figure2_layer_zonal_areas_",out_suffix,".png",sep=""))
  
  res_pix<-960
  col_mfrow<-1
  row_mfrow<-1
  png(filename=out_fig_filename,
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  plot(subset(s_raster,zonal_colnames),main=zonal_colnames)
  
  dev.off()
  
  list_fig_filename[[2]] <-  out_fig_filename
  
  reg_var_list <- l_rast[var_names] #only select population raster
  r_stack <- stack(reg_var_list)
  projection(r_stack) <- proj_str
  #names(r_stack) <- names(data_tb)[var_names]
  
  ##### 
  #Report this basic information in  a text file 
  filename_stat <- paste0("raster_time_series_basics_info_",out_suffix,".txt")
  
  sink(file=filename_stat)
  "Raster stack dimension for variable:"
  dim(r_stack) #34x49x230 
  "Number of cells in a raster:"
  ncell(s_raster) #1666
  "Number of NA values:"
  freq(r_FID,value=NA) #122
  "Number of valid (not NA):"
  ncell(s_raster) - freq(r_FID,value=NA) #1544
  "Spatial resolution in x and y:"
  res(r_stack) #about 1km
  "Projection information:"
  projection(r_stack)
  sink()
  
  #improve this later!
  
  ##################
  ## Figure 3: visualization of time series
  
  out_fig_filename <- file.path(out_dir,paste("Figure3_visualization_time_series_time_step_1_to_4_",out_suffix,".png",sep=""))
  
  #projection(r_stack) <- CRS_WGS84 #making sure the same  CRS format is used
  res_pix<-960
  col_mfrow<-1
  row_mfrow<-1
  png(filename=out_fig_filename,
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  p <- levelplot(r_stack,layers=1:4,col.regions=matlab.like(125)) #show first four images (half a year more or less)
  #plot(r_stack,y=1:4)
  print(p)
  dev.off()
  
  list_fig_filename[[3]] <- out_fig_filename
  
  ## Figure 4: visualization of event in time series
  
  #n_time_event #108
  #dates3[108]
  #[1] "2005-08-29"
  
  ## plottting after and before event
  n_before <- n_time_event - 2
  n_after <- n_time_event + 2
  
  out_fig_filename <- file.path(out_dir,paste("Figure4_visualization_time_series_time_event_before_after_",out_suffix,".png",sep=""))
  
  res_pix<-960
  col_mfrow<-1
  row_mfrow<-1
  png(filename=out_fig_filename,
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  p<- levelplot(r_stack,
                layers=n_before:n_after,
                col.regions=matlab.like(125))

  print(p)
  
  dev.off()
  
  list_fig_filename[[4]] <- out_fig_filename
  
  ## Figure 5: histogram before and after event

  out_fig_filename <- file.path(out_dir,paste("Figure5_histogram_time_series_time_event_before_after_",out_suffix,".png",sep=""))
  
  res_pix<-960
  col_mfrow<-1
  row_mfrow<-1
  png(filename=out_fig_filename,
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  p <- histogram(subset(r_stack,n_before:n_after))
  print(p)
  dev.off()
  
  list_fig_filename[[5]] <- out_fig_filename
  
  ##############
  ## Figure 6: Profile for the time series

  out_fig_filename <- file.path(out_dir,paste("Figure6_time_series_time_profiles_pixel_",pixel_index,"_",out_suffix,".png",sep=""))
  
  #mean_vals <- colMeans(data_tb[,var_names],na.rm=T)
  mean_vals <- cellStats(r_stack,"mean")
  
  pix_val <- r_stack[index_val]
  #pixval <- data_tb[pixel_index,var_names] #should choose earlier wich pixel
  
  #pix300 <- data_tb[300,var_names]
  
  res_pix<-960
  col_mfrow<-1
  row_mfrow<-1
  png(filename=out_fig_filename ,
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  plot(1:length(var_names),mean_vals,type="b",ylab="var",xlab="time step",
       main="Average variable and pixel 800 profile")
  lines(1:length(var_names),pix_val,type="b",ylab="var",xlab="time step",col=c("red"),
        main="Average variable and pixel 800 profile")
  legend("bottomleft",legend=c("Overall average VAR ","PIX 800 "),
         col=c("black","red"),lty=c(1,2))
  abline(v=n_time_event,col="blue")
  
  dev.off()
  
  list_fig_filename[[6]] <- out_fig_filename
  
  ##################
  ## Figure 7: Average pop per elevation zones (observed data)
  ## By zone/strata
  
  r_zonal <- subset(s_raster,zonal_colnames)
  zones_tb_avg<- zonal(r_stack,r_zonal,fun='mean')
  
  zones_avg_df <- as.data.frame(zones_tb_avg)
  n_zones <- length(unique(zones_avg_df$zone))
  
  n_time <- ncol(zones_avg_df) -1
  #pred_names <- c("zone",paste("t",2:n_time,sep="_"),"method")
  
  mydata<- zones_avg_df
  #dd <- do.call(make.groups, mydata[,-ncol(mydata)]) 
  dd <- do.call(make.groups, mydata) 
  #dd$lag <- mydata$lag 
  dd$zones <- mydata$zones
  #dd$method <- mydata$method
  #drop first few rows that contain no data but zones...
  n_start <-n_zones +1
  dd <- dd[n_start:nrow(dd),]
  tmp_time <-unlist(lapply(1:n_time,FUN=function(i){rep(i,n_zones)}))
  dd$time  <- tmp_time
  dd$zones <- mydata$zone #use recycle rule
  
  out_fig_filename <- file.path(out_dir,paste("Figure7a_average_by_zonal_areas_time_series_time_",out_suffix,".png",sep=""))
  
  res_pix<-960
  col_mfrow<-1
  row_mfrow<-1
  png(filename=out_fig_filename,
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  p <- xyplot(data~time |zones,#group=method,
         data=dd,type="b",xlab="time",ylab="VAR",
         #strip = strip.custom(factor.levels=c("z3","z4","z5")), #fix this!!!
         auto.key = list("topright", corner = c(0,1),# col=c("black","red"),
                         border = FALSE, lines = TRUE,cex=1.2))
         
  print(p)
  
  dev.off()
  
  list_fig_filename[[7]] <- out_fig_filename 
  
  ###############
  
  out_fig_filename <- file.path(out_dir,paste("Figure7b_average_by_zonal_areas_time_series_time_xyplot_",out_suffix,".png",sep=""))
  
  res_pix<-960
  col_mfrow<-1
  row_mfrow<-1
  png(filename= out_fig_filename,
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  p<- xyplot(data~time,group=zones,
         data=dd,type="b",xlab="time",ylab="VAR",
         #strip = strip.custom(factor.levels=c("z3","z4","z5")), #fix this!!!
         auto.key = list("topright", corner = c(0,1),# col=c("black","red"),
                         border = FALSE, lines = TRUE,cex=1.2),
         main="Average by zones for VAR"
  )
  print(p)
  dev.off()
  
  list_fig_filename[[8]] <- out_fig_filename 
  
  ##### Prepare object to  return
  
  list_out_filename <- c(list_fig_filename,filename_stat)
  
  write.table(zones_avg_df,file=paste("zones_avg_df","_",out_suffix,".txt",sep=""))
  write.table(dd,file=paste("zones_avg_df_long_table","_",out_suffix,".txt",sep=""))
  
  explore_obj <- list(zones_avg_df,dd,list_out_filename)
  names(explore_obj) <- c("zones_avg","zones_avg_long_tb","out_files")
  
  return(explore_obj)
}

################### END OF SCRIPT ##################