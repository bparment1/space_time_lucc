####################################    Space Beats Time Research  #######################################
################################ Generating datasets for examples of curves #######################################
#This script generates datasets for the Space Beats Time Framework.
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 07/28/2017 
#DATE MODIFIED: 08/18/2017
#Version: 1
#PROJECT:  with Marco Millones            
#
#COMMENTS: - Generate additional datasets
#
#TO DO:
# - 
#COMMIT: initial commit curve analysis
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

create_dir_fun <- function(outDir,out_suffix=NULL){
  #if out_suffix is not null then append out_suffix string
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    outDir <- file.path(outDir,out_name)
  }
  #create if does not exists
  if(!file.exists(outDir)){
    dir.create(outDir)
  }
  return(outDir)
}

#Used to load RData object saved within the functions produced.
load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

function_sbt_curve_generation <- "sbt_curve_generation_functions_08182017.R"
script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts"
source(file.path(script_path,function_sbt_curve_generation))

#####  Parameters and argument set up ###########

in_dir <- "/home/bparmentier/Google Drive/Space_beats_time/Data" #PARAM 1

out_dir <- "/home/bparmentier/Google Drive/Space_beats_time/outputs" #PARAM 2
proj_str<- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0"  #PARAM 3

## Constant
file_format <- ".tif" #PARAM5 #PARAM 4
NA_flag_val <- -9999 #PARAM7 #PARAM5

out_suffix <-"curve_analysis_Katrina_08172017" # PARAM6, output suffix for the files and output folder 
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

###### Part 1: load in data ###################

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

#lf <- file.path(out_dir,paste0(names(r_stack),".tif"))
#lf
  
#list_raster_filename <- paste("time_raster_",out_suffix,".txt",sep="")
#write.table(lf,file=list_raster_filename,sep=",")
#df_rast <- read.table(list_raster_filename,sep=",")

########## Generate figures and metrics ##########

var_mean_stack <- cellStats(r_stack,mean,na.rm=T)
var_mean <- cellStats(r_var,mean,na.rm=T)

plot(var_mean,type="b")
plot(var_mean_stack,type="b")

#Make this a time series object?

t_corr <- unlist(lapply(2:n_layers,FUN=t_corr_fun,list_rast=r_var))
t_corr_val <- t_corr_fun(2,list_rast=r_var)

plot(t_corr_val,type="b")


f <- matrix(c(1,1,1,
              1,0,1,
              1,1,1), nrow=3)

#undebug(Moran_run)
Moran_run(1,r_var,f=f)

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

