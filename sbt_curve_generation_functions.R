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

  ###
  #?Moran
  
  if(moran_type=="queen"){
  
    f <- matrix(c(1,1,1,
                  1,0,1,
                  1,1,1), nrow=3)
    
  }
  
  #undebug(t_corr_fun)
  #t_corr_fun(2,list_rast=list_r)
  
  #var_mean_stack <- cellStats(r_stack,mean,na.rm=T)
  var_mean <- cellStats(r_var,mean,na.rm=T)
  #Make this a time series object?
  
  t_corr <- unlist(lapply(2:n_layers,FUN=t_corr_fun,list_rast=r_var))
  #t_corr_val <- t_corr_fun(2,list_rast=r_var)
  
  plot(t_corr,type="b")
  plot(var_mean,type="b")
  
  #undebug(Moran_run)
  #Moran_run(1,r_var,f=f)
  
  s_corr <- unlist(lapply(2:nlayers(r_var),FUN=Moran_run,r_stack=r_var,f=f))
  
  ## Plot temporal and spatial correlation:
  
  plot(t_corr,type="b",col="magenta",
       ylim=c(-1,1),
       ylab="correlation",
       xlab="Time steps")
  lines(s_corr,type="b",col="blue")
  legend("bottomleft",
         legend=c("temporal","spatial"),
         col=c("magenta","blue"),
         lty=1,
         cex=0.8)
  title("Spatial and temporal correlation") 
  
  ### Plot on different scale
  par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
  plot(t_corr,type="l",
       col="magenta",
       ylab="temporal correlation",
       xlab="Time steps")
  par(new = TRUE)
  plot(s_corr,
       type="l",
       col="blue",
       axes = FALSE, 
       bty="n",xlab = "", ylab = "")
  axis(side=4, at = pretty(range(s_corr)))
  mtext("spatial correlation", side=4, line=3)
  legend("bottomleft",
         legend=c("temporal","spatial"),
         col=c("magenta","blue"),
         lty=1,
         cex=0.8)
  title("Spatial and temporal correlation") 
  
  #curve(var_mean)
  
  #### Now examine shape of inputs:
  #debug(compute_change_shape_metrics)
  obj_metrics <- compute_change_shape_metrics(var_mean,method = "differencing")
  
  #plots()
  ### Now mean plot?
  plot(var_mean,type="l",
       ylim=c(2500,6500))
  
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
  
  text(n_time_event,4500,pos=1,
       labels="A",
       cex=1)

  #annotate(c("A"),X=a_metric_coords[[1]][1],
  #         x=a_metric_coords[[1]][1],cex=c(1.4,3),
  #         arrow.lwd=c(2,4),
  #         arrow.length=c(0.1,0.25),
  #         border.lwd=c(1,3), 
  #        fill=c("pink","white"), col=c("blue","green"),
  #         adj=c(0,1),border.col=c("orange","purple"))
  
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

  text(15,2900,pos=1,
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
  
  text(24,3200,pos=1,
       labels="C",
       cex=1)
  
  ### Now plot space beats time
  plot(mae_tot_tb$temp_arima, 2:23,
       type="l",
       col="magenta")
  lines(mae_tot_tb$spat_reg_no_previous,
        type="l",
        col="blue")
  legend("topright",
         legend=c("temporal","spatial"),
         col=c("magenta","blue"),
         lty=1,
         cex=0.8)
  
  ### Generate tables:
  
  df_ts <- data.frame(s_corr=s_corr,t_corr=t_corr,var=var_mean)
  df_ts$time_step <- 1:nrow(df_ts)
  View(df_ts)
  
  #compute_change_shape_metrics
  plot(t_corr,s_corr,type="b")
  text(t_corr,s_corr,pos=1,labels=1:22,cex=1)
  #text(c(2,2),c(37,35),labels=c("Non-case","Case"))
  
  return()
}


################# END OF SCRIPT ##################

