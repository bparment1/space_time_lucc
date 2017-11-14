#######################################    Space Time Analyses Project   #######################################
############################ Assessing Space and Time Model for Space Beats Time framework #######################################
#This script contains functions to assess predictions for Space Beats Time Framework.       
#Temporal predictions use OLS with the image of the previous time step or ARIMA.
#Spatial predictions use spatial regression (lag error model) with different estimation methods (e.g. eigen, chebyshev etc.).
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/07/2017 
#DATE MODIFIED: 11/14/2017
#Version: 1

#PROJECT: Space beats time Framework
#TO DO:
#
#COMMIT: more changes and testing of code with tiles space and time predictions assessment
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

args<-commandArgs(TRUE)

args_table <- args[1]

#args_table <- "/home/bparmentier/Google Drive/Space_beats_time/Data/input_arguments_sbt_script_NDVI_Rita_10292017.csv"
args_table <- "/home/parmentier/Data/Space_beats_time/Data/input_arguments_sbt_script_NDVI_Rita_11072017.csv"

df_args <- read.table(args_table,sep=",",stringsAsFactors = FALSE)

### use column 2,3,4 etc.
#index_val <- 2 #this is set up for parallelization, if we have multiple regions/tiles, tile1
index_val <- 3 #this is set up for parallelization, if we have multiple regions/tiles, tile 2 rita

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
#previous_step <- df_args[15,index_val] 
date_range <- df_args[15,index_val]  
r_ref <- df_args[16,index_val]  
temp_fname <- df_args[17,index_val]
spat_fname <- df_args[18,index_val]
method_space <- df_args[19,index_val] 
method_time <- df_args[20,index_val] 
pixel_index <- df_args[21,index_val] 

#P1
in_dir <- "/home/parmentier/Data/Space_beats_time/outputs/output_tile_2_NDVI_Rita_11062017"
r_temp_pred <- list.files(path=in_dir,
                          pattern="r_temp_pred_arima_arima_.*._tile_2_NDVI_Rita_11062017.tif",
                          full.names=T)
out_file <- "/home/parmentier/Data/Space_beats_time/Data/data_Rita_NDVI/rev_project_output/tile_1/raster_spat_files_list_tile_2.txt"
write.table(r_temp_pred,out_file)
#r_spat_pred <- list.files(path=in_dir,
#                                         pattern="r_spat_.*._tile_2_NDVI_Rita_11062017.tif",
#                                         full.names=T)
r_spat_pred <- list.files(path=in_dir,
                                         pattern="r_spat_pred_mle_eigen_no_previous_step_.*._tile_2_NDVI_Rita_11062017.tif",
                                         full.names=T)
out_file <- "/home/parmentier/Data/Space_beats_time/Data/data_Rita_NDVI/rev_project_output/tile_2/raster_spat_files_list_tile_2.txt"
#r_spat_pred_mle_eigen_no_previous_step__t_113_tile_2_NDVI_Rita_11062017.tif
s_raster <- "/home/parmentier/Data/Space_beats_time/Data/data_Rita_NDVI/rev_project_output/tile_2/raster_files_list_tile_2.txt"
date_range <- "2001.01.01;2010.12.31;16"
time_window_predicted <- "105;114"
r_zonal <- "crop_r_zonal_rev_tile_2"
r_ref <- NULL
method_space <- "mle;eigen" #method for space and time used to predict space and time respectively
method_time <- "arima;arima;TRUE"
out_suffix <- "assessment_tile_2_NDVI_Rita_11062017"
out_dir <- "output_tile_1_2_combined_NDVI_Rita_11062017"
create_out_dir <- FALSE  

###### Functions used in this script

debug(accuracy_space_time_calc)
test <- accuracy_space_time_calc(r_temp_pred,r_spat_pred,s_raster,time_window_predicted,
                                     r_zonal,methods_space,method_time,r_ref,out_suffix,
                                     var_names,NA_flag_val,file_format, date_range,
                                   out_dir,create_out_dir)
  
accuracy_space_time_calc <- function(r_temp_pred,r_spat_pred,s_raster,time_window_predicted,
                                     r_zonal,methods_space,method_time,r_ref,out_suffix,
                                     var_names,NA_flag_val,file_format,
                                     out_dir,create_out_dir){
  
  ##INPUTS
  #1)r_temp_pred
  #2) r_spat_pred
  #3) s_raster
  #4) time_window_predicted
  #5) r_zonal
  #6) r_ref
  #7) methods_name
  #7) out_suffix
  #8) out_dir
  ##OUTPUTS
  #1)
  #2)
  #
  
  ######## Start script #####
  
  ###  Load data if needed:
  if(class(r_temp_pred)=="character"){
    lf <-r_temp_pred 
    r_temp_pred <- stack(lf)
  }
  if(class(r_spat_pred)=="character"){
    lf <-r_spat_pred 
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
  
  date_range <- (unlist(strsplit(date_range,";")))
  
  r_obs <- subset(s_raster,time_window_predicted[1]:time_window_predicted[2]) #remove 99 because not considred in the prediction!!
  
  #plot(r_huric_obs)
  #r_huric_w <- crop(r_huric_w,rast_ref)
  levelplot(r_obs,col.regions=matlab.like(25))
  levelplot(r_temp_pred,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
  projection(r_temp_pred) <- proj_str

  projection(spat_pred_rast_mle_eigen_with_previous) <- proj_str

  #### Step 2: compute residuals
  
  res_temp_s <- r_temp_pred - r_obs
  res_spat_s <- r_spat_pred - r_obs

  names(res_temp_s) <- sub("pred","res",names(res_temp_s))
  names(res_spat_s) <- sub("pred","res",names(res_spat_s))

  ### Compute residuals by zones
  rast_zonal <- subset(s_raster,match(zonal_colnames,names(s_raster)))
  
  r_results <- stack(s_raster,rast_zonal,temp_pred_rast,
                                          res_temp_s,res_spat_s)
  
  dat_out <- as.data.frame(r_results)
  dat_out <- na.omit(dat_out)
  
  filename_dat_out <- file.path(out_dir,paste("dat_out_",out_suffix,".txt",sep=""))
  write.table(dat_out,file=filename_dat_out,
              row.names=F,sep=",",col.names=T)
  
  ### Now accuracy assessment using MAE
  
  out_suffix_s <- paste("temp_",method_time[1],"_",out_suffix,sep="")

  #undebug(calc_ac_stat_fun)

  ac_temp_obj <- calc_ac_stat_fun(r_pred_s=r_temp_pred,
                                  r_var_s=r_obs,
                                  r_zones=rast_zonal,
                                  file_format=file_format,
                                  out_suffix=out_suffix_s)  
  
  #out_suffix_s <- paste("spat_mle_eigen_with_previous",out_suffix,sep="_")  
  out_suffix_s <- paste("spat_",method_spatial,out_suffix,sep="_")
  ac_spat <- calc_ac_stat_fun(r_pred_s=r_spat_pred,
                              r_var_s=r_huric_obs,
                              r_zones=rast_zonal,
                              file_format=file_format,
                              out_suffix=out_suffix_s)
  
  mae_tot_tb <- cbind(ac_spat_obj$mae_tb,
                      ac_temp_obj$mae_tb)
  
  name_method_time <- paste0("temp_",method_time[1],"_",method_time[2])
  name_method_space <- paste0("spat_",method_time[1],"_",method_space[2])
  mae_tot_tb <- as.data.frame(mae_tot_tb)
  row.names(mae_tot_tb) <- NULL
  #names(mae_tot_tb)<- c("spat_reg_no_previous","spat_reg_with_previous",name_method_time)#,"temp_lm")
  #mae_tot_tb$time <- 2:nrow(mae_tot_tb)
  #mae_tot_tb$time <- 2:n_pred
  mae_tot_tb$time <- 1:nlayers(r_obs)
  y_range<- range(cbind(mae_tot_tb$spat_reg_no_previous,mae_tot_tb$spat_reg_with_previous,mae_tot_tb[[name_method_time]]))
  
  res_pix<-960
  col_mfrow<-1
  row_mfrow<-1
  png(filename=paste("Figure_temporal_profiles_time_1_",out_suffix,".png",sep=""),
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  temp_formula_str <- paste0(paste0("temp_",method_time[1])," ~ ","time")
  plot(spat_reg_no_previous ~ time, type="b",col="cyan",data=mae_tot_tb,ylim=y_range)
  #lines(as.formula(temp_formula_str), type="b",col="magenta",data=mae_tot_tb)
  lines(spat_reg_with_previous ~ time, type="b",col="blue",data=mae_tot_tb)
  legend("topleft",
         legend=c(namee_method_space,name_method_time),
         col=c("cyan","blue"),
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
