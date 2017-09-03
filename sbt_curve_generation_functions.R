####################################    Space Beats Time Research  #######################################
################################ Function to generate and analyze curves #######################################
#This script contains function to analyze datasets for the Space Beats Time Framework.
#This includes:  time series correlation, spatial correlation for steps of time series stack.
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 08/17/2017 
#DATE MODIFIED: 09/03/2017
#Version: 1
#PROJECT:  Space Beats Time             
#
#COMMENTS: - Generate additional datasets
#
#TO DO:
# - 
#COMMIT: adding Moran's I calculation and documentation
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
library(sf)

############## start of script #################

t_corr_fun <- function(i,list_rast){
  #
  #
  
  if(class(list_rast)=="character"){
    r_subset <- stack(list_rast[[i-1]],list_rast[[i]]) 
  }else{
    r_subset <- stack(subset(list_rast,i-1),subset(list_rast,i))
  }
  cor_val<-layerStats(r_subset,'pearson',na.rm=T)
  cor_val <- cor_val[[1]][1,2]
  return(cor_val)
}

Moran_run<- function (i,r_stack,f=NULL){
  #function to compute Moran's I on a raster stack
  if(is.null(f)){
    f <- matrix(c(1,1,1,
                  1,0,1,
                  1,1,1), nrow=3)
  }
  moran_val <- Moran(subset(r_stack,i),w=f)
  return(moran_val)
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
  a_metric_index <- grep(1,test) +1 #because diff is shifter left
  
  a_metric <- var_mean[b_range[1]-1] - var_mean[a_metric_index]  #drop
  c_metric <- var_mean[b_range[length(b_range)]]- var_mean[a_metric_index]
  
  df_shape_metrics <- data.frame(a=a_metric,b=b_metric,c=c_metric)
  
  obj_metrics <- list(df_shape_metrics,b_range)
  names(obj_metrics) <- c("df_shape_metrics","b_range")
  
  return(obj_metrics)
}

generate_dates_by_step <-function(start_date,end_date,step_date){
  #library(xts) declare out of this function
  #library(zoo)
  #library(lubridate)
  
  st <- as.Date(start_date,format="%Y.%m.%d")
  en <- as.Date(end_date,format="%Y.%m.%d")
  #year_list <-seq(format(st,"%Y"),format(en,"%Y")) #extract year
  year_list <- seq(as.numeric(strftime(st,"%Y")),as.numeric(strftime(en,"%Y"))) #extract year
  
  ll_list <- vector("list",length=length(year_list))
  for (i in 1:length(year_list)){
    if(i==1){
      first_date <-st
    }else{
      first_date<-paste(year_list[[i]],"-01","-01",sep="")
    }
    if(i==length(year_list)){
      last_date <-en
    }else{
      last_date<-paste(year_list[[i]],"-12","-31",sep="")
    }
    #ll <- seq.Date(st, en, by=step)
    ll <- seq.Date(as.Date(first_date), as.Date(last_date), by=step_date)
    ll_list[[i]]<-as.character(ll)
    #paste(yday(ll,)
  }
  
  #
  dates_modis <-as.Date(unlist((ll_list))) 
  #wiht 001
  dates_DOY_modis <- paste(year(dates_modis),sprintf("%03d", yday(dates_modis)),sep="")
  dates_obj <- list(dates_modis,dates_DOY_modis)
  names(dates_obj) <- c("dates","doy")  
  return(dates_obj)
}

generate_plots_table <- function(r_var,mae_tot_tb,moran_type="queen",out_suffix="",out_dir="."){
  ## This function generates plots and tables for mean profile, correlaion and shape metrics.
  #
  #
  #
  
  ##########
  ### PART 1: compute means, metrics and correlations
  
  #?Moran
  
  if(moran_type=="queen"){
  
    f <- matrix(c(1,1,1,
                  1,0,1,
                  1,1,1), nrow=3)
    
  }
  
  var_mean <- cellStats(r_var,mean,na.rm=T)

  t_corr <- unlist(lapply(2:n_layers,FUN=t_corr_fun,list_rast=r_var))

  #undebug(Moran_run)
  #Moran_run(1,r_var,f=f)
  
  s_corr <- unlist(lapply(2:nlayers(r_var),FUN=Moran_run,r_stack=r_var,f=f))

  obj_metrics <- compute_change_shape_metrics(var_mean,method = "differencing")
  
  #### PART 2: Plots of mean profile and shape metrics 
  #debug(compute_change_shape_metrics)
  
  res_pix<- 500
  col_mfrow<-1
  row_mfrow<-1
  
  png_filename_var_mean <- paste("Figure_mean_plot_var_with_shape_metrics_",
                        out_suffix,".png",sep="")
  png(png_filename_var_mean,
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  ### Now mean plot?
  plot(var_mean,type="l",
       ylim=c(2500,6500),
       ylab="Mean NDVI",
       xlab="Time Steps")
  
  ###  A metric 
  
  a_metric_coords <- list(c(9.5,5500),c(9.5,3500))
  #segments(a_metric_coords[[1]][1],
  #        a_metric_coords[[1]][2],
  #        a_metric_coords[[2]][1],
  #        a_metric_coords[[2]][2])
  arrows(a_metric_coords[[1]][1],
         a_metric_coords[[1]][2],
         a_metric_coords[[2]][1],
         a_metric_coords[[2]][2],
         code=3,
         cex=0.7,
         col="green")
  #text(a_metric_coords[[1]][1]/2,a_metric_coords[[1]][2]/2,"A")
  n_time_event <- 9 # for 108
  
  text(11.5,4500,pos=2,
       labels="A",
       cex=1)
  
  #### Add B metric
  #do B, 200 below the lowest point
  b_metrics <- 13
  b_metric_coords <- list(c(9,2800),c(22,2800))
  
  arrows(b_metric_coords[[1]][1],
         b_metric_coords[[1]][2],
         b_metric_coords[[2]][1],
         b_metric_coords[[2]][2],
         code=3,
         cex=0.4,
         col="green")
  
  text(15,2850,pos=1,
       labels="B",
       cex=1)
  
  #### Add C metric
  #do C, 100 to the right
  c_metrics <- 13
  c_metric_coords <- list(c(23,2900),c(23,3500))
  
  arrows(c_metric_coords[[1]][1],
         c_metric_coords[[1]][2],
         c_metric_coords[[2]][1],
         c_metric_coords[[2]][2],
         code=3,
         cex=0.4,
         lwd=0.7,
         col="green")
  
  text(24,3300,pos=1,
       labels="C",
       cex=1)

  legend("topright",
         legend=c("A: Loss","B: Length of event","C: Recovery"),
         #col=c("magenta","blue"),
         #lty=1,
         cex=1,
         bty="n")
  title("Mean value and event shape metrics") 
  
  dev.off()
  
  ### Part 4: Space and Time predictions
  ### Now plot space beats time

  res_pix<- 500
  col_mfrow<-1
  row_mfrow<-1
  
  png_filename_sbt <- paste("Figure_space_and_time_predictions_",
                        out_suffix,".png",sep="")
  png(png_filename_sbt,
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  plot(2:24,mae_tot_tb$temp_arima,
       #2:23,
       type="l",
       col="magenta",
       xlab="Time Steps",
       ylab="Mean Absolute Error (MAE)")
  lines(mae_tot_tb$spat_reg_no_previous,
        type="l",
        col="blue")
  legend("topright",
         legend=c("temporal","spatial"),
         col=c("magenta","blue"),
         lty=1,
         bty="n",
         cex=1)
  title("Space and Time prediction errors")
  
  #### Add D metric
  #do D, 100 to the right
  
  d_metrics <- 13
  #d_metric_coords <- list(c(8,620),c(8,2080))
  d_metric_coords <- list(c(12,630),c(12,2080))
  
  arrows(d_metric_coords[[1]][1],
         d_metric_coords[[1]][2],
         d_metric_coords[[2]][1],
         d_metric_coords[[2]][2],
         code=3,
         cex=0.4,
         lwd=0.7,
         col="green")
  
  text(12.5,1500,pos=1,
       labels="D",
       cex=1)

  
  #### Add E metric
  #do e, 100 to the right
  
  e_metrics <- 13
  #d_metric_coords <- list(c(8,620),c(8,2080))
  e_metric_coords <- list(c(8,500),c(14,500))
  
  arrows(e_metric_coords[[1]][1],
         e_metric_coords[[1]][2],
         e_metric_coords[[2]][1],
         e_metric_coords[[2]][2],
         code=3,
         cex=0.4,
         lwd=0.7,
         col="green")
  
  text(11,490,pos=1,
       labels="E",
       cex=1)
  
  legend("topleft",
         legend=c("D: Strength","E: Length of event"),
         #col=c("magenta","blue"),
         #lty=1,
         cex=1,
         bty="n")

  dev.off()  
  
  ###### Part 5:
  ## Plot temporal and spatial correlation:
  res_pix<- 500
  col_mfrow<-1
  row_mfrow<-1
  
  png_filename_temp_spat_correlation <- paste("Figure_temporal_spatial_correlation_",
                                 out_suffix,".png",sep="")
  png(png_filename_temp_spat_correlation,
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  plot(t_corr,type="l",col="magenta",
       ylim=c(-1,1),
       ylab="correlation",
       xlab="Time steps")
  lines(s_corr,type="l",col="blue")
  legend("bottomleft",
         legend=c("temporal","spatial"),
         col=c("magenta","blue"),
         lty=1,
         cex=1,
         bty="n")
  title("Spatial and temporal correlation") 
  
  dev.off()
  
  
  ### Generate tables:
  
  df_info <- data.frame(s_corr=s_corr,t_corr=t_corr,var=var_mean[2:length(var_mean)])
  #df_ts$time_step <- 1:nrow(df_ts)
  #View(df_ts)
  
  ######
  
  obj_info <- list(df_info,
                   list(png_filename_var_mean,
                   png_filename_var_mean,png_filename_sbt,
                   png_filename_temp_spat_correlation))
  names(obj_info) <- c("df_info","list_png")
  
  #####
  
  return(obj_info)
}


################# END OF SCRIPT ##################

