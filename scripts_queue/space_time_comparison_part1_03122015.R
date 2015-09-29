####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: Comparison of results  #######################################
#This script compares predictions results obtained from spatial regression and ARIMA.                        
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/10/2014 
#DATE MODIFIED: 03/12/2015
#Version: 2
#PROJECT: AAG and geocomputation with Marco Millones             
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

calc_ac_stat_fun <- function(r_pred_s,r_var_s,r_zones,file_format=".tif",out_suffix){
  #Parameters:
  #Input:
  #r_pred_s: raster stack of layers predictions (for example predicted NDVI)
  #r_var_s: raster stack of layers actual values (for example observed NDVI)
  #r_zones: raster defining zones of relevance to accuracy (e.g.hurricane winds zones)
  #
  
  ##Functions used
  rmse_fun <-function(x){sqrt(mean(x^2,na.rm=TRUE))}
  mae_fun <-function(x){mean(abs(x),na.rm=TRUE)}

  ##Start script
  
  #Accuracy/errors by zones
  r_res_s <- r_pred_s - r_var_s #residuals, stack of raster layers
  mse_zones_tb <- zonal(r_res_s^2,r_zones,fun="mean") #mean square error
  mae_zones_tb <- zonal(abs(r_res_s),r_zones,fun="mean") #absolute error
  rmse_zones_tb <- sqrt(mse_zones_tb) #root mean square error
  
  #Overall Accuracy/errors 
  
  mae_tb <- cellStats(abs(r_res_s),mean) #calculate MAE for layer stack
  rmse_tb <- sqrt(cellStats((r_res_s)^2,mean)) #calculate rmse overall
  #write out residuals rasters
  r_res_s <- writeRaster(r_res_s,filename=file.path(out_dir,"r_res_s.tif"),bylayer=TRUE,
                         suffix=paste(1:nlayers(r_res_s),out_suffix,sep="_"),overwrite=TRUE)
  ac_obj <- list(mae_tb,rmse_tb,mae_zones_tb,rmse_zones_tb)
  names(ac_obj) <- c("mae_tb","rmse_tb","mae_zones_tb","rmse_zones_tb")
  
  return(ac_obj)
}


#####  Parameters

#in_dir<-"C:/Users/mmmillones/Dropbox/Space_Time"
in_dir <-"/home/parmentier/Data/Space_Time"
#in_dir <- "/Users/benoitparmentier/Dropbox/Data/Space_Time/"
out_dir <- in_dir
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

out_suffix <-"03172014" #output suffix for the files that are masked for quality and for 
file_format <- ".tif"
ref_sejidos_name <- "~/Data/Space_Time/00_sejidos_group_sel5_ids.rst"
ref_winds_name <- "~/Data/Space_Time/00_windzones_moore_sin.rst"

########### START SCRIPT #################

####   Read in data  ####

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

#### Read in spatial predictions  ####

#out_suffix_s <- "03102014"
out_suffix_s <- paste(c("_spat_lag_.*"),out_suffix,sep="_")
r_pred_spat <- stack(list.files(path=out_dir,pattern=paste(out_suffix_s,
                                                           file_format,sep="")))
#plot(r_pred_spat,colNA="black")
r_pred_spat <- subset(r_pred_spat,1:4) #select predictions layers

## Read in temporal predictions...

pattern_str <-"test_mooore_auto_.*._03042014.*"
r_pred_temp_all <- stack(list.files(path=out_dir,pattern=paste(pattern_str,file_format,"$",sep="")))
r_pred_temp <- subset(r_pred_temp_all,c(1,3,5,7)) #get predictions images
plot(r_pred_temp,colNA="black") #this is a in modis sinusoidal projection
projection(r_pred_temp) <- CRS_interp

#### NOW TAKE A LOOK AT THE ACCURACY OF PREDICTION!!!

#ref_sejidos_name <- "~/Data/Space_Time/00_sejidos_group_sel5_ids.rst"
#ref_winds_name <- "~/Data/Space_Time/00_windzones_moore_sin.rst"

## Information used to break down accuracy by regions
r_winds <- raster(ref_winds_name)
r_winds<- crop(r_winds,r_clip)
projection(r_winds) <- CRS_interp

r_winds_m <- projectRaster(r_winds,crs=CRS_WGS84) #Check that it is using ngb

### Process predctions data to match i.e. reproject and clip
#r_pred_temp <- projectRaster(r_pred_temp,crs=CRS_WGS84)
r_pred_temp <- projectRaster(r_pred_temp,r_winds_m)

r_huric_w <- subset(r_stack,c(153,154,155,156)) #select relevant actual observations
r_huric_w <- projectRaster(r_huric_w,r_winds_m)## check use nearest neighbour option!!! we have negative NDVI

#r_huric_w2 <- subset(r_stack,c(152,153,154,155,156)) #select relevant actual observations
#r_huric_w2 <- projectRaster(r_huric_w2,r_winds_m)## check use nearest neighbour option!!! we have negative NDVI
#moore_m <- projectRaster(moore_w,r_winds_m)
#r_huric_w2 <- mask(r_huric_w2,moore_m)

#names(r_huric_wy) <- 
#names(r_huric_w2)<-c("t_step_152","t_step_153","t_step_154","t_step_155","t_step_156")

names(r_huric_w)<-c("t_step_153","t_step_154","t_step_155","t_step_156")
#writeRaster(r_huric_w,filename="r_huric_w.tif",bylayer=TRUE,suffix=names(r_huric_w))
r_huric_w <- writeRaster(r_huric_w,filename=file.path(out_dir,"r_huric_w.tif"),
                           bylayer=TRUE,suffix=names(r_huric_w),overwrite=TRUE)

## year 2007 profile of NDVI values

r_huric_wy <- subset(r_stack,139:161) #select relevant actual observations
r_huric_wy <- projectRaster(r_huric_wy,r_winds_m)## check use nearest neighbour option!!! we have negative NDVI
moore_m <- projectRaster(moore_w,r_winds_m)
r_huric_wy <- mask(r_huric_wy,moore_m)
#mean_year_2007 <- cellStats(r_huric_wy,mean)


mean_year_2007 <- cellStats(r_huric_wy,mean,na.rm=TRUE)
mean_year_2007$zone <- 0 #zone 0 will be the overall area!!!
mean_year_2007 <- do.call(cbind,mean_year_2007)
mean_year_2007 <- as.data.frame(mean_year_2007)
mean_year_2007 <- mean_year_2007[,c(24,1:23)]

zonal_mean_2007 <- zonal(r_huric_wy,r_winds_m,mean)
zonal_mean_2007 <- zonal_mean_2007[4:6,]
zonal_mean_2007 <- as.data.frame(zonal_mean_2007)
zonal_mean_2007 <- rbind(zonal_mean_2007,mean_year_2007) #NOTE THAT ZONE=0 is used to store overall mean NDVI
#x<- merge(zonal_mean_2007,mean_year_2007,by="zone")

write.table(zonal_mean_2007,file=paste("NDVI_2007_avg_profile_by_zones_and_overall","_",out_suffix,".txt",sep=""))

layout_m <- c(1,1)
png(paste("Figure_ndvi_","profile","_","mean","_",out_suffix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])

plot(1:23,zonal_mean_2007[4,2:24],type="b",pch=1,col="red",ylim=c(2000,10000)) #overall
lines(1:23,zonal_mean_2007[1,2:24],type="b",pch=2,col="black") #zone 4
lines(1:23,zonal_mean_2007[2,2:24],type="b",pch=3,col="green") #zone 5
lines(1:23,zonal_mean_2007[3,2:24],type="b",pch=4,col="blue") #zone 6
abline(v=15.5,lty="dashed")
legend("topleft",legend=c("Overall","zone 4","zone 5", "zone 6"),
        cex=0.8, col=c("red","black","green","blue"),bty="n",
        lty=1,pch=1:4)
legend("topright",legend=c("hurricane event"),cex=0.8,lty="dashed",bty="n")
title("Average NDVI in the study area and by zone",cex=1.6, font=2)

dev.off()

#axis(side=1,las=1,tick=TRUE,
#       at=breaks_lab,labels=tick_lab, cex.axis=1.2,font=2) #reduce number of labels to Jan and June
#  #text(tick_lab, par(\u201cusr\u201d)[3], labels = tick_lab, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
#  axis(2,cex.axis=1.4,font=2)
#  box()

#Now calculate accuracy metrics and residuals images
#undebug(calc_ac_stat_fun)
out_suffix_s <- paste("temp_",out_suffix,sep="_")
ac_temp_obj <- calc_ac_stat_fun(r_pred_s=r_pred_temp,r_var_s=r_huric_w,r_zones=r_winds_m,
                                file_format=file_format,out_suffix=out_suffix_s)
out_suffix_s <- paste("spat_",out_suffix,sep="_")  
ac_spat_obj <- calc_ac_stat_fun(r_pred_s=r_pred_spat,r_var_s=r_huric_w,r_zones=r_winds_m,
                                file_format=file_format,out_suffix=out_suffix_s)

## quick look at residuals
#r_res_s <- stack(list.files(pattern="r_res_s_.*__03112014.tif"))

#Create correlation matrix among 
cor_mat <- layerStats(r_huric_w,stat="pearson",na.rm=TRUE)
cor_mat2 <- layerStats(r_huric_w2,stat="pearson",na.rm=TRUE)

cor_mat <- cor_mat[[1]]
write.table(cor_mat,file=paste("correlation_matrix","_",out_suffix,".txt",sep=""))

## COMPARE METHODS OF PREDICTIONS : OVERALL ASSESSMENT

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

#### COMPARE METHODS OF PREDICTIONS BY ZONES ASSESSMENT

mae_zones_tb <- rbind(ac_spat_obj$mae_zones_tb[4:6,],
                      ac_temp_obj$mae_zones_tb[4:6,])
mae_zones_tb <- as.data.frame(mae_zones_tb)
mae_zones_tb$method <- c("spat_reg","spat_reg","spat_reg","arima","arima","arima")
names(mae_zones_tb) <- c("zones","pred1","pred2","pred3","pred4","method")

write.table(mae_zones_tb,file=paste("mae_zones_tb","_",out_suffix,".txt",sep=""))

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


########### Combining all accuracy: spat_reg, spat_lag, reg_ols and arima for four dates...

#spat_reg :  v2= W1_v1 + v1
#spat_lag :  v2= W1_v1
#reg_ols  :  v2= v1
#arima :  

#accuracy combined by the total area of study: overall MAE

mae_tot_tb1 <- read.table("mae_tot_tb_03112014.txt") #spat reg 
mae_tot_tb2 <- read.table("mae_tot_tb_03132014.txt") #reg ols
mae_tot_tb3 <- read.table("mae_tot_tb_03172014.txt") # spat lag

#mae_tot_tb2$method <- c("ols_reg","ols_reg","ols_reg","arima","arima","arima")
mae_tot_tb_combined <- cbind(mae_tot_tb1[,1],mae_tot_tb2[,1],mae_tot_tb3)
names(mae_tot_tb_combined) <- c("spat_reg","ols_reg","spat_lag","arima","time")
#mae_tot_tb_combined$time <- s:4 
write.table(mae_tot_tb_combined,file=paste("mae_tot_tb_combined","_",out_suffix,".txt",sep=""))

layout_m<-c(1,1.5) #one row two columns

png(paste("Figure_accuracy_","mae","_","total_combined_by_model","_",out_suffix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])

plot(spat_reg ~ time, type="b",col="magenta",pch=1,data=mae_tot_tb_combined,ylim=c(800,1900),
     ylab="MAE for NDVI")
lines(ols_reg  ~ time, type="b",col="green",pch=2,data=mae_tot_tb_combined)
lines(spat_lag  ~ time, type="b",col="red",pch=3,data=mae_tot_tb_combined)
lines(arima ~ time, type="b",col="cyan",pch=4,data=mae_tot_tb_combined)

#axis(side=1,las=1,tick=TRUE,
#       at=breaks_lab,labels=tick_lab, cex.axis=1.2,font=2) #reduce number of labels to Jan and June
#  #text(tick_lab, par(\u201cusr\u201d)[3], labels = tick_lab, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
#  axis(2,cex.axis=1.4,font=2)
#  box()
legend("topleft",legend=names(mae_tot_tb_combined)[1:4], 
        cex=1.4, col=c("magenta","green","red","cyan"),bty="n",
        lty=1,pch=1:4)
title("Average Mean Absolute Error (MAE) by method",cex=1.6, font=2)

dev.off()

#accuracy combined by zones

mae_zones_tb1 <- read.table("mae_zones_tb_03112014.txt")
mae_zones_tb2 <- read.table("mae_zones_tb_03132014.txt")
mae_zones_tb3 <- read.table("mae_zones_tb_03172014.txt")


mae_zones_tb2$method <- c("ols_reg","ols_reg","ols_reg","arima","arima","arima")
mae_zones_tb3$method <- c("spat_lag","spat_lag","spat_lag","arima","arima","arima")

mae_zones_tb_combined <- rbind(mae_zones_tb1[1:3,],mae_zones_tb2[1:3,],mae_zones_tb3)

write.table(mae_zones_tb_combined,file=paste("mae_zones_tb_combined","_",out_suffix,".txt",sep=""))

mae_val <- (as.vector(as.matrix(mae_zones_tb_combined[,2:5])))
avg_ac_tb <- as.data.frame(mae_val)

avg_ac_tb$metric <- rep("mae",length(mae_val))
avg_ac_tb$zones <- rep(c(3,4,5),4)
avg_ac_tb$time <- c(rep(1,12),rep(2,12),rep(3,12),rep(4,12))
avg_ac_tb$method <- rep(c("spat_reg","spat_reg","spat_reg",
                          "ols_reg","ols_reg","ols_reg",
                          "spat_lag","spat_lag","spat_lag",
                          "arima","arima","arima"),4)
names(avg_ac_tb)[1]<- "ac_val"
names_panel_plot <- c("time -1","time +1","time +2","time +3")
layout_m<-c(1,1.5) #one row two columns

png(paste("Figure_accuracy_","mae","_combined","by_winds_zones","_",out_suffix,".png", sep=""),
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

############# END OF SCRIPT ###################