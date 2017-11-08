#######################################    Space Time Analyses Project   #######################################
############################ Assessing Space and Time Model for Space Beats Time framework #######################################
#This script contains functions to assess predictions for Space Beats Time Framework.       
#Temporal predictions use OLS with the image of the previous time step or ARIMA.
#Spatial predictions use spatial regression (lag error model) with different estimation methods (e.g. eigen, chebyshev etc.).
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/07/2017 
#DATE MODIFIED: 11/08/2017
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



r_temp_pred <-
r_spat_pred <-
s_raster <-
time_window_predicted <-
r_zonal <-
r_ref <-
methods_name <- c("mle_eigen","arima") #method for space and time used to predict space and time respectively
out_suffix <-
out_dir <-
  
###### Functions used in this script

accuracy_space_time_calc <- function(r_temp_pred,r_spat_pred,s_raster,time_window_predicted,
                                     r_zonal,methods_name,r_ref,out_suffix,out_dir){
  
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
  
  r_obs <- subset(s_raster,time_window_predicted) #remove 99 because not considred in the prediction!!
  
  #plot(r_huric_obs)
  #r_huric_w <- crop(r_huric_w,rast_ref)
  levelplot(r_huric_obs,col.regions=matlab.like(25))
  #levelplot(r_temp_pred_rast_arima,col.regions=matlab.like(25))
  
  method_time[1]=="arima"
  r_temp_pred_rast <- stack(unlist(lapply(1:length(l_pred_temp_arima),
                                                FUN=function(i,x){obj<-x[[i]];
                                                obj$raster_pred},
                                                x=l_pred_temp_arima)))
  
  #r_temp_res_rast <- stack(unlist(lapply(1:length(l_pred_temp_arima),
  #                                       FUN=function(i,x){obj<-x[[i]];
  #                                       obj$raster_res},
  #                                       x=l_pred_temp_arima)))
  levelplot(r_temp_pred_rast_arima,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
  levelplot(r_temp_res_rast_arima,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
  projection(r_temp_pred_rast_arima) <- proj_str
  projection(r_temp_res_rast_arima) <- proj_str
  temp_pred_rast <- subset(r_temp_pred_rast_arima,1:length(time_window_predicted))
    
  
  


  levelplot(spat_pred_rast_mle_eigen_no_previous,col.regions=rev(matlab.like(255))) #view the four predictions using mle spatial reg.
  levelplot(spat_res_rast_mle_eigen_no_previous,col.regions=rev(matlab.like(255))) #view the four predictions using mle spatial reg.
  
  levelplot(spat_pred_rast_mle_eigen_with_previous,col.regions=matlab.like(255)) #view the four predictions using mle spatial reg.
  levelplot(spat_res_rast_mle_eigen_with_previous,col.regions=matlab.like(255)) #view the four predictions using mle spatial reg.
  
  projection(spat_pred_rast_mle_eigen_no_previous) <- proj_str
  projection(spat_res_rast_mle_eigen_no_previous) <- proj_str
  projection(spat_pred_rast_mle_eigen_with_previous) <- proj_str
  projection(spat_res_rast_mle_eigen_with_previous) <- proj_str
  
  
  #r_winds_m <- crop(winds_wgs84,res_temp_s) #small test window
  #res_temp_s_lm <- temp_pred_rast_lm - r_huric_obs
  #res_temp_s_arima <- temp_pred_rast_arima - r_huric_obs
  
  #### Step 2: compute residuals
  
  res_temp_s <- temp_pred_rast - r_obs
  res_spat_s <- spat_pred_rast_mle_eigen_no_previous - r_obs

  names(res_temp_s) <- sub("pred","res",names(res_temp_s))
  names(res_spat_s) <- sub("pred","res",names(res_spat_s))

  #browser()
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
  
  method_time
  out_suffix_s <- paste("temp_",method_time[1],"_",out_suffix,sep="")

  #undebug(calc_ac_stat_fun)

  ac_temp_obj <- calc_ac_stat_fun(r_pred_s=temp_pred_rast,
                                  r_var_s=r_obs,
                                  r_zones=rast_zonal,
                                  file_format=file_format,
                                  out_suffix=out_suffix_s)  
  
  method_spatial
  #out_suffix_s <- paste("spat_mle_eigen_with_previous",out_suffix,sep="_")  
  out_suffix_s <- paste("spat_",method_spatial,out_suffix,sep="_")
  ac_spat <- calc_ac_stat_fun(r_pred_s=spat_pred_rast,
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
  
  mae_tot_tb <- cbind(ac_spat_obj$mae_tb,
                      ac_spat_obj$mae_tb,
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
