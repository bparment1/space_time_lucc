####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: Comparison of results  #######################################
#This script compares predictions results obtain from spatial regression and ARIAM.                        
#AUTHORS: Marco Millones and Benoit Parmentier                                             
#DATE CREATED: 03/10/2014 
#DATE MODIFIED: 03/11/2014
#Version: 1
#PROJECT: GLP Conference Berlin,YUCATAN CASE STUDY with Marco Millones             
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
library(forecast)
library(xts)
library(zoo)
library(lubridate)
library(colorRamps) #contains matlab.like color palette
library(rgeos)

###### Functions used in this script

### Make a function to uniformize NA accros given dates!!!

create_uniform_mask_fun <- function(r_var,proj_str=NULL,r_clip=NULL){
  #Function to create uniform mask accross raster layer stack
  
  if(!is.null(r_clip)){
      r_var_w <- crop(r_var,r_clip) #clip (or "window") raster image to Moore subset area
      r_var_w <- mask(r_var_w,r_clip) #mask areas
  }else{
    r_var_w <- r_var
  }
  
  r_NA_mask <-calc(r_var_w,sum,na.rm=FALSE)
  if(!is.null(proj_str)){
    r_NA_mask <- projectRaster(r_NA_mask,crs=proj_str) #project to latlong
  }
  return(r_NA_mask)
}

calc_ac_stat_fun <- function(r_pred_s,r_var_s,r_zones){
  
  ##Functions used
  rmse_fun <-function(x){sqrt(mean(x^2,na.rm=TRUE))}
  mae_fun <-function(x){mean(abs(x),na.rm=TRUE)}

  ##Start script
  
  #Accuracy/errors by zones
  r_res_s <- r_pred_s - r_var_s
  mse_zones_tb <- zonal(r_res_s^2,r_zones,fun="mean")
  mae_zones_tb <- zonal(abs(r_res_s),r_zones,fun="mean")
  rmse_zones_tb <- sqrt(mse_zones_tb)
  
  #Overall Accuracy/errors 
  
  mae_tb <- cellStats(abs(r_res_s),mean) #calculate MAE for layer stack
  rmse_tb <- sqrt(cellStats((r_res_s)^2,mean))

  ac_obj <- list(mae_tb,rmse_tb,mae_zones_tb,rmse_zones_tb)
  names(ac_obj) <- c("mae_tb","rmse_tb","mae_zones_tb","rmse_zones_tb")
  
  return(ac_obj)
}


#####  Parameters

#in_dir<-"C:/Users/mmmillones/Dropbox/Space_Time"
in_dir <-"/home/parmentier/Data/Space_Time"
#in_dir <- "/Users/benoitparmentier/Dropbox/Data/Space_Time/"

#set up the working directory
setwd(in_dir)

moore_window <- file.path(in_dir,"moore_window.rst")
test_shp_path <- file.path(in_dir,"GEODA_Analysis") #path to get to the files for spatial pred
#  "/Users/benoitparmentier/Dropbox/Data/Space_Time/GEODA_Analysis/TEST_SPATIAL_ONE.shp"
proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
proj_str<- CRS_WGS84
CRS_interp <- proj_modis_str

out_suffix <-"03102014" #output suffix for the files that are masked for quality and for 

########### START SCRIPT #################

## Read in data
moore_w <- raster(moore_window) 
projection(moore_w) <- CRS_interp

#reads input NDVI time series images
reg_var_list <- list.files(path=in_dir,pattern="reg_mosaiced_MOD13A2_A.*.__005_1_km_16_days_NDVI_09242013_09242013.rst$")
#reg_mosaiced_MOD13A2_A2012353__005_1_km_16_days_NDVI_09242013_09242013.rst
#reg_var_list <- list.files(pattern="reg_mosaiced_MOD13A2_A.*.__005_1_km_16_days_NDVI_09242013_09242013.rst$")

r_stack <- stack(reg_var_list) #276 images from 2001 to 2012 included
#levelplot(r_stack,layers=1:2) #show first 2 images (half a year more or less)

projection(r_stack) <- CRS_interp #this is modis sinusoidal
r_var <- subset(r_stack,153:154) #date before hurricane is 153
r_clip <- moore_w
NAvalue(r_clip) <- 0


#list_spat_reg_obj <- lapply(1:length(list_models),FUN=predict_var_spatial,list_param=list_param_spat_reg)

r_pred_spat <- stack(list.files(path=out_dir,pattern=paste(out_suffix,file_format,sep="")))
plot(r_pred_s,colNA="black")

r_pred_spat <- subset(r_pred_spat,1:4)
#r_pred_s <- subset(r_pred_s,1:4)

#### NOW TAKE A LOOK AT THE ACCURACY OF PREDICTION!!!

mae_val <- c(mae1,mae2,mae3,mae4)
ref_sejidos_name <- "~/Data/Space_Time/00_sejidos_group_sel5_ids.rst"
ref_winds_name <- "~/Data/Space_Time/00_windzones_moore_sin.rst"

#r_sejidos <- raster(ref_sejidos_name)
r_winds <- raster(ref_winds_name)
r_winds<- crop(r_winds,r_clip)
projection(r_winds) <- CRS_interp

r_winds_m <- projectRaster(r_winds,crs=CRS_WGS84)

#Now calculate the mae

mse_winds <- zonal(r_pred_s^2,r_winds_m,fun="mean")
mae_winds <- zonal(abs(r_pred_s),r_winds_m,fun="mean")
sqrt(mse_winds)
mae_winds

## Quick look at temporal predictions...

pattern_str <-"test_mooore_fire_auto_.*._03072014.*"
r_pred_temp_all <- stack(list.files(path=out_dir,pattern=paste(pattern_str,file_format,sep="")))
r_pred_temp <- subset(r_pred_temp_all,c(1,3,5,7))
#plot(r_pred_temp,colNA="black")

#r_pred_temp <- projectRaster(r_pred_temp,crs=CRS_WGS84)
r_pred_temp <- projectRaster(r_pred_temp,r_winds_m)

r_huric_w <- subset(r_stack,c(153,154,155,156))
r_huric_w <- projectRaster(r_huric_w,r_winds_m)
r_huric_w <- mask(r_huric_w,r_pred_temp)

#debug(calc_ac_stat_fun)
ac_temp_obj <- calc_ac_stat_fun(r_pred_s=r_pred_temp,r_var_s=r_huric_w,r_zones=r_winds_m)
  
ac_spat_obj <- calc_ac_stat_fun(r_pred_s=r_pred_spat,r_var_s=r_huric_w,r_zones=r_winds_m)

##OVERALL ASSESSMENT

#mae_tot_tb <- t(rbind(ac_spat_obj$mae_tb,ac_temp_obj$mae_tb))
mae_tot_tb <- (cbind(ac_spat_obj$mae_tb,ac_temp_obj$mae_tb))
mae_tot_tb <- as.data.frame(mae_tot_tb)
row.names(mae_tot_tb) <- NULL
names(mae_tot_tb)<- c("spat_reg","arima")
mae_tot_tb$time <- as.integer(1:4)

layout_m<-c(1,1.5) #one row two columns

png(paste("Figure_accuracy_","mae","_","total_by_model","_",out_suffix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])

plot(spat_reg ~ time, type="b",col="magenta",pch=1,data=mae_tot_tb,ylim=c(800,1900),
     ylab="MAE for NDVI")
lines(arima ~ time, type="b",col="cyan",pch=2,data=mae_tot_tb)
#axis(side=1,las=1,tick=TRUE,
#       at=breaks_lab,labels=tick_lab, cex.axis=1.2,font=2) #reduce number of labels to Jan and June
#  #text(tick_lab, par(\u201cusr\u201d)[3], labels = tick_lab, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
#  axis(2,cex.axis=1.4,font=2)
#  box()
legend("topleft",legend=names(mae_tot_tb)[1:2], 
        cex=1.4, col=c("magenta","cyan"),bty="n",
        lty=1,pch=1:2)
title("Average Mean Absolute Error (MAE) by method",cex=1.6, font=2)

dev.off()

write.table(mae_tot_tb,file=paste("mae_tot_tb","_",out_suffix,".txt",sep=""))

#### BY ZONES ASSESSMENT
mae_zones_tb <- rbind(ac_spat_obj$mae_zones_tb[4:6,],
                      ac_temp_obj$mae_zones_tb[4:6,])
mae_zones_tb <- as.data.frame(mae_zones_tb)
mae_zones_tb$method <- c("spat_reg","spat_reg","spat_reg","arima","arima","arima")
names(mae_zones_tb) <- c("zones","pred1","pred2","pred3","pred4","method")

write.table(mae_zones_tb,file="mae_zones_tb","_",out_suffix,".txt",sep=""))

mae_val <- (as.vector(as.matrix(mae_zones_tb[,2:5])))
avg_ac_tb <- as.data.frame(mae_val)

avg_ac_tb$metric <- rep("mae",length(mae_val))
avg_ac_tb$zones <- rep(c(3,4,5),4)
avg_ac_tb$time <- c(rep(1,6),rep(2,6),rep(3,6),rep(4,6))
avg_ac_tb$method <- rep(c("spat_reg","spat_reg","spat_reg","arima","arima","arima"),4)
names(avg_ac_tb)[1]<- "ac_val"
names_panel_plot <- c("time -1","time +1","time +2","time +3")
layout_m<-c(1,1.5) #one row two columns

png(paste("Figure_accuracy_","mae","_","by_winds_zones","_",out_suffix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])

p <- xyplot(ac_val~zones|time, # |set up pannels using method_interp
            group=method,data=avg_ac_tb, #group by model (covariates)
            main="Average MAE by winds zones and time step ",
            type="b",as.table=TRUE,
            #strip = strip.custom(factor.levels=c("time -1","time +1","time +2","time +3")),
            strip=strip.custom(factor.levels=names_panel_plot),
            index.cond=list(c(1,2,3,4)),    #this provides the order of the panels)
            pch=1:length(avg_ac_tb$method),
            par.settings=list(superpose.symbol = list(
            pch=1:length(avg_ac_tb$method))),
            auto.key=list(columns=1,space="right",title="Model",cex=1),
            #auto.key=list(columns=5),
            xlab="Winds zones",
            ylab="MAE for NDVI")
print(p)
dev.off()
## Other plots

#Compare predictions for 153 (t-1), before hurricane

plot(mae_zones_tb$zones,mae_zones_tb$spat_pred1,type="b",col="red",ylim=c(800,2000))
lines(mae_zones_tb$zones,mae_zones_tb$temp_pred1,type="b",col="black",ylim=c(800,2000))

#Compare predictions for 154 (t+1), one time step after hurricane

plot(mae_zones_tb$zones,mae_zones_tb$spat_pred2,type="b",col="red",ylim=c(800,2000))
lines(mae_zones_tb$zones,mae_zones_tb$temp_pred2,type="b",col="black",ylim=c(800,2000))

#Compare predictions for 155 (t+2), one time step after hurricane

plot(mae_zones_tb$zones,mae_zones_tb$spat_pred3,type="b",col="red",ylim=c(800,2000))
lines(mae_zones_tb$zones,mae_zones_tb$temp_pred3,type="b",col="black",ylim=c(800,2000))

#Compare predictions for 156 (t+2), one time step after hurricane

plot(mae_zones_tb$zones,mae_zones_tb$spat_pred4,type="b",col="red",ylim=c(700,2000))
lines(mae_zones_tb$zones,mae_zones_tb$temp_pred4,type="b",col="black",ylim=c(700,2000))
#titles()



