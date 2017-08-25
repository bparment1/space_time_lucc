####################################    Space Beats Time Research  #######################################
################################ Generating datasets for examples of curves #######################################
#This script performs analyses of curves correlation and shape.
#It use datasets for the Space Beats Time Framework.
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 07/28/2017 
#DATE MODIFIED: 08/25/2017
#Version: 1
#PROJECT:  with Marco Millones            
#
#COMMENTS: - Generate additional datasets
#
#TO DO:
# - 
#COMMIT: changes to generate plot and tables function
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

function_sbt_curve_generation <- "sbt_curve_generation_functions_08252017.R"
script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts"
source(file.path(script_path,function_sbt_curve_generation))

#####  Parameters and argument set up ###########

in_dir <- "/home/bparmentier/Google Drive/Space_beats_time/Data" #PARAM 1

out_dir <- "/home/bparmentier/Google Drive/Space_beats_time/outputs" #PARAM 2
proj_str<- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0"  #PARAM 3

## Constant
file_format <- ".tif" #PARAM5 #PARAM 4
NA_flag_val <- -9999 #PARAM7 #PARAM5

out_suffix <-"curve_analysis_Katrina_08252017" # PARAM6, output suffix for the files and output folder 
create_out_dir_param=TRUE #PARAM7

#coord_names <- c("x","y") #PARAM 9
#zonal_colnames <- "r_srtm_Katrina_rec2" #PARAM 12
num_cores <- 4 #PARAM 11

#n_time_event <- 3
n_time_event <- "108" #PARAM 12 #this is the timestep corresponding to the event ie Hurricane Katrina (Aug 23- Aub 31 2005): 235-243 DOY, storm surge Aug 29 in New Orelans

#time_window_selected <- 2:5
#time_window_selected <- 100:116 #PARAM 13: use alll dates for now
time_window_selected <- 100:123

#n_time_event <- 108 #PARAM 12 #this is the timestep corresponding to the event ie Hurricane Katrina (Aug 23- Aub 31 2005): 235-243 DOY, storm surge Aug 29 in New Orelans
#time_window_selected <- 100:116 #PARAM 13: use alll dates for now
date_range <- c("2001.01.01","2010.12.31",16) #PARAM 15, NDVI Katrina

sbt_results_filename <- "/home/bparmentier/Google Drive/Space_beats_time/outputs/output_NDVI_Katrina_08242017/mae_tot_tb_NDVI_Katrina_08242017.txt"

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



f <- matrix(c(1,1,1,
              1,0,1,
              1,1,1), nrow=3)


##### remove seasonality

date_range <- c("2001.01.01","2010.12.31") #PARAM 15, NDVI Katrina
range_dates <- generate_dates_by_step(date_range[1],date_range[2],16)$dates #NDVI Katrina
class(range_dates)

#s <- setZ(s, as.Date('2000-1-1') + 0:2)
#r_ts <- setZ(subset(r_stack,1:230), range_dates)

#r_tmp <- rollmean(r_ts,k=12)

# x <- zApply(s, by=as.yearqtr, fun=mean, name='quarters')
#r_tmp <- zApply(r_ts, by=12, fun=rollmean)

data_df <- as.data.frame(r_stack)
data_df <- na.omit(data_df)

df_ts <- (t(data_df))
dim(df_ts)

df_ts <- zoo(df_ts,range_dates)

(df_ts[1:10,])

dim(df_ts)

df_ts_smoothed <- rollmean(df_ts,k=23)
plot(df_ts[,1])
lines(df_ts_smoothed[,1],col="red")

#time_window_selected <- 100:116 #PARAM 13: use alll dates for now

range_dates[100:123]
df_ts_w <- df_ts[100:123,]
plot(df_ts_w[,1])

## Run shape code on the window and also get the space beats time predictions for 100:123
#/home/bparmentier/Google Drive/Space_beats_time/outputs

sbt_results_filename <- "/home/bparmentier/Google Drive/Space_beats_time/outputs/output_NDVI_Katrina_08242017/mae_tot_tb_NDVI_Katrina_08242017.txt"
mae_tot_sbt <- read.table(sbt_results_filename,sep=" ",header=T,stringsAsFactors = F)

View(mae_tot_sbt)

debug(generate_plots_table)
lf_file <- generate_plots_table(r_var,
                        mae_tot_tb=mae_tot_sbt,
                        moran_type="queen",
                        out_suffix=out_suffix,
                        out_dir=out_dir)
  
  
################# END OF SCRIPT ##################

