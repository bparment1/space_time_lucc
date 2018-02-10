####################################    Space Time Analyses PAPER   #######################################
###########################################  Figures production  #######################################
#This script produces figures based on outputs from earlier analyses for the space beats time framework.

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 04/20/2015 
#DATE MODIFIED: 02/10/2018
#Version: 1
#PROJECT: SBT RITA book chapter
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

compute_avg_by_zones <- function(r_stack,r_zonal,out_suffix_str="",out_dir="."){
  #Quick function to compute average by zone using a raster zone and a raster stack
  #
  #
  library(raster)
  zones_tb_avg<- raster::zonal(r_stack,r_zonal,fun='mean')
  
  zones_avg_df <- as.data.frame(zones_tb_avg)
  n_zones <- length(unique(zones_avg_df$zone))
  #unique(r_zonal)
  
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
  
  #xyplot(data~time |zones,#group=method,
  #       data=dd,type="b",xlab="time",ylab="VAR",
  #       #strip = strip.custom(factor.levels=c("z3","z4","z5")), #fix this!!!
  #       auto.key = list("topright", corner = c(0,1),# col=c("black","red"),
  #                       border = FALSE, lines = TRUE,cex=1.2)
  #)
  
  #xyplot(data~time,group=zones,
  #       data=dd,type="b",xlab="time",ylab="VAR",
  #       #strip = strip.custom(factor.levels=c("z3","z4","z5")), #fix this!!!
  #       auto.key = list("topright", corner = c(0,1),# col=c("black","red"),
  #                      border = FALSE, lines = TRUE,cex=1.2),
  #       main="Average by zones for VAR"
  #)
  
  outfile1 <- file.path(out_dir,paste("zones_avg_df","_",out_suffix_str,".txt",sep=""))
  outfile2 <- file.path(out_dir,paste("zones_avg_df_long_table","_",out_suffix_str,".txt",sep=""))
  write.table(zones_avg_df,file=outfile1)
  write.table(dd,file=outfile2)
  zones_obj <- list(zones_avg_df,dd)
  names(zones_obj) <- c("zones_avg_df","zones_avg_df_long_table")
  return(zones_obj)
}

plot_temporal_time_series_profile_by_zones <- function(start_date,end_date,dates,date_event,n_time_event,data_tb,r_var,r_zonal,var_name,y_range,x_label,title_str,out_dir,out_suffix_str){
  ##This is assuming a maximum of three regions... change this later...
  #
  #
  #
  #
  date_event <- as.Date(date_event)
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)
  index_dates_selected <- dates >= start_date & dates <= end_date
  dates_selected <- dates[index_dates_selected]
  #dates3[n_time_event3]
  n_time_event_selected <- which(dates_selected==dates[n_time_event]) #time of the event...
  
  #debug(compute_avg_by_zones)
  zonal_obj <- compute_avg_by_zones(r_stack=r_var,r_zonal=r_zonal,out_suffix_str=out_suffix_str,
                                    out_dir=out_dir)
  zones_avg_df <- zonal_obj$zones_avg_df
  n_zones <- nrow(zones_avg_df)
  #find out which date is 107!!!
  df <- as.data.frame(t(zones_avg_df))
  n_zones_labels <- zones_avg_df$zone
  #mean_vals <- colMeans(data_tb[,index_dates_selected],na.rm=T)
  mean_vals <- cellStats(subset(r_var,which(index_dates_selected)),mean)
  #df$mean_vals<-mean_vals
  df_ts <- zoo(mean_vals,dates_selected)
  #pixval <- data_tb[800,var_names]
  #pix300 <- data_tb[300,var_names]
  n_step_selected <- length(mean_vals)
  layout_m <- c(1.5,1.1)
  
  png(paste("Figure","_1b_","average_temporal_profiles_by_zones_subset",var_name,"_data_",out_suffix_str,".png", sep=""),
      height=480*layout_m[2],width=480*layout_m[1])
  # Set margins to make room for x axis labels
  #par(mar = c(7, 4, 4, 2) + 0.3)
  plot(df_ts,type="b",col="red",ylim=y_range,ylab=var_name,
       xlab=x_label)
  col_pal <- c("red","black","green","blue")
  pch_type <- c(1,2,3,4)
  for(i in 1:n_zones){
    par(new=TRUE)
    #lines(1:n_step_selected,zones_avg_df[1,index_dates_selected3],type="b",pch=2,col="black") #zone 4
    if(ncol(zones_avg_df)>length(mean_vals)){
      zones_avg_df_tmp <- zones_avg_df[,-1] #drop first column with zones
    }
    plot(1:n_step_selected,zones_avg_df_tmp[i,index_dates_selected],type="b",pch=pch_type[1+i],
         col=col_pal[1+i],
         ylim=y_range,ylab="",xlab="",axes=F) #zone 4
  }
  ##Add doted vertical line at the time of the event
  if(is.null(date_event)){
    #abline(v= n_time_event_selected+0.5,lty="dashed")
    abline(v= n_time_event_selected,lty="dashed")
  }else{
    n_time_event_selected_tmp <- dates_selected==date_event
    
    if(sum(n_time_event_selected_tmp)==1){ #if one, it means that there is a date matching!
      n_time_event_selected <- which(dates_selected==date_event) #time of the event...
      abline(v= n_time_event_selected,lty="dashed")
    }
    if(sum(n_time_event_selected_tmp)==0){ #if one, it means that there is no date matching!
      #n_time_event_selected <- which(dates_selected==date_event) #time of the event...
      diff_date_event <- dates_selected - date_event #find the closest matching date
      diff_min_index <- which(diff_date_event==min(as.numeric(abs(diff_date_event)))) #min distance in date term
     if(dates_selected[diff_min_index] > date_event){
       n_time_event_selected <- diff_min_index - 0.5
     }
      abline(v= n_time_event_selected,lty="dashed")
    }
  }
  
  legend("topleft",legend=c("Overall",paste("zone ",n_zones_labels,sep="")),
         cex=1, col=col_pal,bty="n",
         lty=1,pch=1:4)
  legend("topright",legend=c("hurricane event"),cex=1,lty="dashed",bty="n")
  title(title_str,cex=1.6, font=2)
  dev.off()
}

## Function for spatial patterns to come here...

################### END OF SCRIPT ##################