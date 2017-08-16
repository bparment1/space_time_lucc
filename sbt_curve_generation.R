####################################    Space Beats Time Research  #######################################
################################ Generating datasets for examples of curves #######################################
#This script generates datasets for the Space Beats Time Framework.
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 07/28/2017 
#DATE MODIFIED: 08/15/2017
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
library(sf)

###### Functions used in this script

function_sbt_curve_generation <- "sbt_curve_generation_functions_08022017.R"
script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts"
source(file.path(script_path,function_sbt_curve_generation))

#####  Parameters and argument set up ###########

in_dir <- "/home/bparmentier/Google Drive/Space_beats_time/Data" #PARAM 1

out_dir <- "/home/bparmentier/Google Drive/Space_beats_time/outputs" #PARAM 2
proj_str<- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0"  #PARAM 3

## Constant
file_format <- ".tif" #PARAM5 #PARAM 4
NA_flag_val <- -9999 #PARAM7 #PARAM5

out_suffix <-"Katrina_08152017" # PARAM6, output suffix for the files and output folder 
create_out_dir_param=TRUE #PARAM7

#coord_names <- c("x","y") #PARAM 9
#zonal_colnames <- "r_srtm_Katrina_rec2" #PARAM 12
num_cores <- 4 #PARAM 11

#n_time_event <- 3
n_time_event <- "108" #PARAM 12 #this is the timestep corresponding to the event ie Hurricane Katrina (Aug 23- Aub 31 2005): 235-243 DOY, storm surge Aug 29 in New Orelans

#time_window_selected <- 2:5
time_window_selected <- 100:116 #PARAM 13: use alll dates for now

#n_time_event <- 108 #PARAM 12 #this is the timestep corresponding to the event ie Hurricane Katrina (Aug 23- Aub 31 2005): 235-243 DOY, storm surge Aug 29 in New Orelans
#time_window_selected <- 100:116 #PARAM 13: use alll dates for now
date_range <- c("2001.01.01","2010.12.31",16) #PARAM 15, NDVI Katrina

################# START SCRIPT ###############################

### PART I READ AND PREPARE DATA FOR REGRESSIONS #######

#set up the working directory
#Create output directory

if(is.null(out_dir)){
  out_dir <- dirname(in_dir) #output will be created in the input dir
}

out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix_s)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

###### PART 1: generate dataset ##########

n_row <- 3
n_col <- 3

r <- raster(nrows=n_row,ncols=n_col)

vtn1 <- c(0.8,0.8,0.2,
         0.8,0.2,0.2,
         0.2,0.8,0.8)

vt0 <- c(0.8,0.8,0.2,
         0.8,0.2,0.2,
         0.2,0.8,0.8)

vt1 <- c(0.3,0.3,0.2,
         0.3,0.2,0.3,
         0.2,0.3,0.3)

vt2 <- c(0.4,0.4,0.2,
         0.4,0.2,0.4,
         0.2,0.4,0.4)

vt3 <- c(0.6,0.6,0.2,
         0.6,0.2,0.6,
         0.2,0.6,0.6)

vt4 <- c(0.6,0.6,0.2,
         0.6,0.2,0.6,
         0.2,0.6,0.6)

list_vt <- list(vtn1,vt0,vt1,vt2,vt3,vt4)

list_r <-(lapply(list_vt, function(x,rast){rast[]<-x;return(rast)},rast=r))
r_var <- stack(list_r)

n_layers <- nlayers(r_var)
r_var_names <- paste0("time_",1:n_layers)
names(r_var) <- r_var_names

## generating additional variables
r_zonal <- subset(r_var,1)
r_zonal[] <- rep(1,ncell(r_zonal))
r_x <-init(r_var,v="x")
r_y <-init(r_var,v="y")
r_covar <- stack(r_zonal,r_x,r_y)
names(r_covar) <- c("rec_zones","x","y")

r_stack <- stack(r_var,r_covar)

#Call function for sbt?
writeRaster(r_stack, names(r_stack), bylayer=TRUE, format='GTiff',overwrite=T)


in_dir_rast <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_NDVI_Katrina"
lf <-list.files(pattern="NDVI_Katrina_08092017.tif",in_dir_rast,full.names=T)
lf <- mixedsort(lf)
r_stack <- stack(lf)

r_var <- subset(r_stack,time_window_selected)

n_layers <- nlayers(r_var)
r_var_names <- paste0("time_",1:n_layers)
names(r_var) <- r_var_names

n_time_event <- 9 # for 108
plot(r_var)

col_palette <- matlab.like(100)
#col_palette <- matlab.like(100)

levelplot(r_var,col.regions=col_palette,main="Time series")


### Write text file with names of input layers

lf <- file.path(out_dir,paste0(names(r_stack),".tif"))
lf
  
list_raster_filename <- paste("time_raster_",out_suffix,".txt",sep="")
write.table(lf,file=list_raster_filename,sep=",")
#df_rast <- read.table(list_raster_filename,sep=",")

########## Generate figures and metrics ##########

var_mean_stack <- cellStats(r_stack,mean)
var_mean <- cellStats(r_var,mean)

plot(var_mean,type="b")
plot(var_mean_stack,type="b")

t_corr <- unlist(lapply(2:n_layers,FUN=t_corr_fun,list_rast=r_var))
t_corr_val <- t_corr_fun(2,list_rast=r_var)

plot(t_corr_val,type="b")


f <- matrix(c(1,1,1,
              1,0,1,
              1,1,1), nrow=3)

undebug(Moran_run)
Moran_run(1,r_var,f=f)
Moran_run<- function (i,r_stack,f=NULL){
  if(is.null(f)){
    f <- matrix(c(1,1,1,
                  1,0,1,
                  1,1,1), nrow=3)
  }
  moran_val <- Moran(subset(r_stack,i),w=f)
  return(moran_val)
}

s_corr <- unlist(lapply(1:nlayers(r_var),FUN=Moran_run,r_stack=r_var,f=f))

plot(s_corr,type="b")

plot(t_corr,type="b",col="pink",ylim=c(-1,1))
lines(s_corr,type="b",col="blue")

par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
plot(t_corr,type="l",col="pink")
par(new = TRUE)
plot(s_corr,type="l",col="blue",axes = FALSE, bty="n",xlab = "", ylab = "")
axis(side=4, at = pretty(range(s_corr)))
mtext("spat corr", side=4, line=3)
################# END OF SCRIPT ##################

