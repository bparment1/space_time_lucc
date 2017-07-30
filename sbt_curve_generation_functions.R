####################################    Space Beats Time Research  #######################################
################################ Generating datasets for examples of curves #######################################
#This script generates datasets for the Space Beats Time Framework.
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 07/28/2017 
#DATE MODIFIED: 07/30/2017
#Version: 1
#PROJECT:  with Marco Millones            
#
#COMMENTS: - Generate additional datasets
#
#TO DO:
# - 
#COMMIT:  
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

###################### Functions used ############


t_corr_fun <- function(i,list_rast){
  r_subset <- stack(list_rast[[i-1]],list_rast[[i]]) 
  cor_val<-layerStats(r_subset,'pearson')
  cor_val <- cor_val[[1]][1,2]
  return(cor_val)
}

compute_change_shape_metrics <- function(var_mean,method="differencing"){
  ## TO find A, B, C: use differencing
  var_diff <- diff(var_mean)
  
  b_range <- c(0,as.numeric(var_diff!=0))
  b_range <- grep(1,b_range)
  
  b_metric <- length(b_range)
  which.min(var_diff)
  
  #class(var_mean)
  test <- as.numeric((names(var_mean)==names(which.min(var_diff))))
  a_metric_index <- grep(1,test) 
  
  a_metric <- var_mean[b_range[1]-1] - var_mean[a_metric_index]  #drop
  c_metric <- var_mean[b_range[length(b_range)]]- var_mean[a_metric_index]
  
  df_shape_metrics <- data.frame(a=a_metric,b=b_metric,c=c_metric)
  
  obj_metrics <- list(df_shape_metrics,b_range)
  names(obj_metrics) <- c("df_shape_metrics","b_range")
  
  return(obj_metrics)
}

generate_plots_table <- function(r_var,mae_tot_tb,moran_type="queen",out_suffix="",out_dir="."){

  ###
  #?Moran
  
  if(moran_type=="queen"){
  
    f <- matrix(c(1,1,1,
                  1,0,1,
                  1,1,1), nrow=3)
    
  }
  
  #undebug(t_corr_fun)
  t_corr_fun(2,list_rast=list_r)
  
  var_mean <- cellStats(r_var,mean)
  s_corr <- unlist(lapply(list_r, function(rast){Moran(rast,w=f)}))
  t_corr <- unlist(lapply(2:length(list_r),FUN=t_corr_fun,list_rast=list_r))
  t_corr <- c(NA,t_corr)
  
  plot(var_mean,type="l")
  
  obj_metrics <- compute_change_shape_metrics(var_mean,method = "differencing")
  
  plot(diff(var_mean),type="l")
  
  plot(t_corr,type="b",col="pink",ylim=c(-1,1))
  lines(s_corr,type="b",col="blue")
  
  par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
  plot(t_corr,type="l",col="pink")
  par(new = TRUE)
  plot(s_corr,type="l",col="blue",axes = FALSE, bty="n",xlab = "", ylab = "")
  axis(side=4, at = pretty(range(s_corr)))
  mtext("spat corr", side=4, line=3)
  
  ### Generate tables:
  
  df_ts <- data.frame(s_corr=s_corr,t_corr=t_corr,var=var_mean)
  df_ts$time_step <- 1:nrow(df_ts)
  #View(df_ts)compute_change_shape_metrics
  plot(t_corr,s_corr,type="b")
  text(t_corr,s_corr,labels=df_ts$time_step,cex=2)
  #text(c(2,2),c(37,35),labels=c("Non-case","Case"))
  
  return(df_ts)
}


################# END OF SCRIPT ##################

