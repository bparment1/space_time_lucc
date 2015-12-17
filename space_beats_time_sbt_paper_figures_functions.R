####################################    Space Time Analyses PAPER   #######################################
###########################################  Figures production  #######################################
#This script produces figures based on outputs from earlier analyses for the space beats time project.
#The script uses spatial regression and temporal model (ARIMA and temporal OLS) values predicted for various. 
#Current case studies include:

#-raster NDVI MODIS data in the Yucatan after hurricane Dean
#-raster DSMP ligth data in New Orleans after hurricane Katrina
#-raster NDVI MODIS data in New Orleans after hurricane Katrina

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 04/20/2015 
#DATE MODIFIED: 12/17/2015
#Version: 1
#PROJECT: GLP Conference Berlin,YUCATAN CASE STUDY with Marco Millones            
#PROJECT: Workshop for William and Mary: an intro to geoprocessing with R 
#PROJECT: Chicago AAG 2015 with Marco Millones
#PROJECT: Dallas Geocomputation 2015 with Marco Millones, Springer book
#
#COMMENTS: - 
#         - 
#TO DO:
#
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

###### Functions used in this script sourced from other files

#function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_05172015_functions.R" #PARAM 1
#script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts" #path to script #PARAM 2
#source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.

##### Functions used in this script 

generate_dates_modis <-function(start_date,end_date,step_date){
  #library(xts) declare out of this function
  #library(zoo)
  #library(lubridate)
  
  st <- as.Date(start_date,format="%Y.%m.%d")
  en <- as.Date(end_date,format="%Y.%m.%d")
  #year_list <-seq(format(st,"%Y"),format(en,"%Y")) #extract year
  year_list <- seq(strftime(st,"%Y"),strftime(en,"%Y")) #extract year
  
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
  }
  
  dates_modis <-as.Date(unlist((ll_list))) 
  return(dates_modis)
}

plot_by_zones_and_timestep_fun <- function(plot_filename,var_name,event_timestep,pix_res,input_data_df,layout_m){
  #plot_filename: name of the png file storing the figure
  #event_timestep: time step for the event on the plot e.g. 1.5 for Dean (if four timesteps are displayed)
  #var_name : name of the variable plotted
  #pix_res : resolution of the plot
  #
  #
  
  ## STARTE SCRIPT ##
  
  png(plot_filename,
    height=pix_res*layout_m[1],width=pix_res*layout_m[2])

  #par(mfrow=layout_m)
  par(mfrow=c(layout_m[1],layout_m[2]))

  #Get range, flatten first the data.frame
  zones_range <- range(input_data_df$zone)
  list_zones_range <- unique(input_data_df$zone)
  df_tmp <- subset(input_data_df,select=c(names(input_data_df)!=c("zones")))
  vect_tmp <- as.vector(unlist(subset(df_tmp,select=c(names(df_tmp)!="method"))))
  y_range_all <- range(vect_tmp)

  
  for ( i in 1:length(list_zones_range)){
    val_zone <- list_zones_range[i]
    col_names <- as.character(unique(input_data_df$method))
    input_data <- subset(input_data_df,zones==val_zone)
    input_data <- subset(input_data,select=c(names(input_data)!="zones"))
    input_data <- subset(input_data,select=c(names(input_data)!="method"))
    input_data <- as.data.frame(t(input_data))
    names(input_data) <- col_names #make spat, time column
    input_data$time <- 1:nrow(input_data) #add a time colum
    
    y_range_offset <- (y_range_all[2] - y_range_all[1])*0.1 #this can be an additional input argument as well
    y_range<- c(y_range_all[1],y_range_all[2]+y_range_offset)
    #xlab_tick <- mae_tot_tb$time
    xlab_tick <- c("T-1","T+1","T+2","T+3")
    x_tick_position <- input_data$time
  
    plot(spat_reg ~ time, type="l",col="cyan",data=input_data,ylim=y_range,
         ylab= var_name,xlab="Time step",xaxt="n",lwd=3,pch=16,cex.lab=1.6)
    points(spat_reg ~ time, col="cyan",data=input_data,ylim=y_range,
         ylab= var_name,xlab="Time step",xaxt="n",lwd=3,pch=16,cex.lab=1.6,cex.pch=1.6)
    lines(temp ~ time, col="magenta",lwd=3,pch=16,data=input_data)
    points(temp ~ time, col="magenta",lwd=3,pch=16,data=input_data,cex.pch=1.6)
    axis(1,at= x_tick_position,labels=xlab_tick,cex=1.2)
  
    abline(v=event_timestep,lty="dashed")
    #abline(v=1.5,lty="dashed")
    #legend("topleft",legend=c("hurricane event"),lty="dashed",bty="n")
    legend("top",legend=c("hurricane event"),lty="dashed",bty="n",cex=1.3)
    #legend("topleft",legend=c("hurricane event"),lty="dashed",bty="n",cex=1.4)
  
    #legend("topright",legend=c("spatial","temporal"),col=c("cyan","magenta"),lty=1,lwd=3,bty="n",cex=1.6)
    #title("Overall MAE for Spatial and Temporal models", cex.main=2) #Note that the results are different than for ARIMA!!!
    legend("topright",legend=c("spatial","temporal"),col=c("cyan","magenta"),lty=1,lwd=3,bty="n",cex=1.3)
    title(paste("Zone ",val_zone,sep=""), cex.main=2) #Note that the results are different than for ARIMA!!!
  }

  dev.off()
  
  ## Prepare function object
  plot_obj <- list(input_data,plot_filename)
  names(plot_obj) <- c("input_data","plot_filename")
  return(plot_obj)
}

plot_by_tot_and_timestep_fun <- function(plot_filename,var_name,event_timestep,pix_res,input_data_df,layout_m){
  #plot_filename: name of the png file storing the figure
  #event_timestep: time step for the event on the plot e.g. 1.5 for Dean (if four timesteps are displayed)
  #var_name : name of the variable plotted
  #pix_res : resolution of the plot
  #
  #
  
  ## STARTE SCRIPT ##
  
  png(plot_filename,
    height=pix_res*layout_m[1],width=pix_res*layout_m[2])

  #par(mfrow=layout_m)

  #Get range, flatten first the data.frame
  #zones_range <- range(input_data_df$zone)
  #list_zones_range <- unique(input_data_df$zone)
  df_tmp <- subset(input_data_df,select=c(names(input_data_df)!=c("time"))) #remove time column before getting the range
  vect_tmp <- as.vector(unlist(subset(df_tmp,select=c(names(df_tmp)!="method"))))
  y_range_all <- range(vect_tmp)

  input_data <- input_data_df
  y_range_offset <- (y_range_all[2] - y_range_all[1])*0.1 #this can be an additional input argument as well
  y_range<- c(y_range_all[1],y_range_all[2]+y_range_offset)
  #xlab_tick <- mae_tot_tb$time
  xlab_tick <- c("T-1","T+1","T+2","T+3")
  x_tick_position <- input_data$time
  
  plot(spat_reg ~ time, type="l",col="cyan",data=input_data,ylim=y_range,
         ylab= var_name,xlab="Time step",xaxt="n",lwd=3,pch=16,cex.lab=1.6)
  points(spat_reg ~ time, col="cyan",data=input_data,ylim=y_range,
         ylab= var_name,xlab="Time step",xaxt="n",lwd=3,pch=16,cex.lab=1.6,cex.pch=1.6)
  lines(temp ~ time, col="magenta",lwd=3,pch=16,data=input_data)
  points(temp ~ time, col="magenta",lwd=3,pch=16,data=input_data,cex.pch=1.6)
  axis(1,at= x_tick_position,labels=xlab_tick,cex=1.2)
  
  abline(v=event_timestep,lty="dashed")
  #abline(v=1.5,lty="dashed")
  #legend("topleft",legend=c("hurricane event"),lty="dashed",bty="n")
  legend("top",legend=c("hurricane event"),lty="dashed",bty="n",cex=1.3)
  #legend("topleft",legend=c("hurricane event"),lty="dashed",bty="n",cex=1.4)
  
  legend("topright",legend=c("spatial","temporal"),col=c("cyan","magenta"),lty=1,lwd=3,bty="n",cex=1.3)
  title("Overall MAE for Spatial and Temporal models", cex.main=1.6) #Note that the results are different than for ARIMA!!!  legend("topright",legend=c("spatial","temporal"),col=c("cyan","magenta"),lty=1,lwd=3,bty="n",cex=1.3)

  #print(p)
  dev.off()

  ## Prepare function object
  plot_obj <- list(input_data,plot_filename)
  names(plot_obj) <- c("input_data","plot_filename")
  return(plot_obj)
}



## Function for spatial patterns to come here...

################### END OF SCRIPT ##################