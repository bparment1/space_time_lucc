########################################  MODIS AND RASTER PROCESSING #######################################
########################################### Read, project, crop and process rasters #####################################
#This script contains general functions to processs raster images, raster time series as well MODIS specific functions.
#This script will form the basis of a library of functions for raster processing of for GIS and Remote Sensing applications.
#AUTHOR: Benoit Parmentier                                                                       
#CREATED ON: 07/23/2018
#MODIFIED ON: 07/23/2018
#PROJECT: None, general utility functions for raster (GIS) processing. 
#COMMIT: multiband option changes for mosaic of MOD09A1
#
#TODO:
#1)Modify generation of CRS for additional projected system (only LCC, Lambert Conformal at this stage)

### List of 22 functions currently available:
# Add documentation later

#[1] "assign_projection_crs"            
#[2] "change_names_file_list"           
#[23] "screening_val_r_stack_fun" 

   
###Loading R library and packages                                                      
#library(gtools)    # loading some useful tools 

library(sp)
library(raster)
library(rgdal)
require(rgeos)
library(BMS) #contains hex2bin and bin2hex
library(bitops)
require(RCurl)
require(stringr)
require(XML)
library(lubridate)

###Loading R library and packages                                                      
#library(gtools)    # loading some useful tools 

library(sp)
library(raster)
library(BMS) #contains hex2bin and bin2hex
library(bitops)
library(gtools)
library(parallel)
library(rasterVis)
library(spdep)                               # Spatial pacakge with methods and spatial stat. by Bivand et al.
library(rgdal)                               # GDAL wrapper for R, spatial utilities
library(gstat)                               # Kriging and co-kriging by Pebesma et al.
library(fields)                              # NCAR Spatial Interpolation methods such as kriging, splines
library(raster)                              # Hijmans et al. package for raster processing
library(rasterVis)
library(spgwr)
library(reshape)
library(sf)
library(dplyr)
library(lubridate)

#################################
### Functions used in the script:

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

function_raster_processing <- "MODIS_and_raster_processing_functions_07182018.R"
function_processing_modis_data <- "processing_MODIS_data_functions_07182018.R"
function_qc_modis_processing <-"QC_layers_modis_processing_functions_07182018.R"

#script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts"  #path to script functions
script_path <- "/media/dan/Space_beats_time/Space_beats_time/sbt_scripts"  #path to script functions

source(file.path(script_path,function_raster_processing)) #source all functions used in this script.
source(file.path(script_path,function_processing_modis_data)) #source all functions used in this script.
source(file.path(script_path,function_qc_modis_processing)) #source all functions used in this script.

################################
###### Parameters and arguments

#ARG1
#in_dir <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_Harvey_NDVI" #param1
#in_dir <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_RITA_reflectance" #param1
in_dir <- "/media/dan/Space_beats_time/Space_beats_time/Data/data_modis_NDVI_maria"
#ARG2
#out_dir <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_Harvey_NDVI" #param2
#out_dir <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_RITA_reflectance" #param1
out_dir <- "/media/dan/Space_beats_time/Space_beats_time/Data/data_modis_NDVI_maria"
function_sbt_curve_generation <- "sbt_curve_generation_functions_09042017.R"
script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts"
source(file.path(script_path,function_sbt_curve_generation))

#####  Parameters and argument set up ###########

in_dir <- "/home/bparmentier/Google Drive/Space_beats_time/Data" #PARAM 1

out_dir <- "/home/bparmentier/Google Drive/Space_beats_time/outputs" #PARAM 2
proj_str<- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0"  #PARAM 3

## Constant
file_format <- ".tif" #PARAM5 #PARAM 4
NA_flag_val <- -9999 #PARAM7 #PARAM5

out_suffix <-"temporal_analyses_Katrina_09032017" # PARAM6, output suffix for the files and output folder 
create_out_dir_param=TRUE #PARAM7

#coord_names <- c("x","y") #PARAM 9
#zonal_colnames <- "r_srtm_Katrina_rec2" #PARAM 12
num_cores <- 4 #PARAM 11

#n_time_event <- "108" #PARAM 12 #this is the timestep corresponding to the event ie Hurricane Katrina (Aug 23- Aub 31 2005): 235-243 DOY, storm surge Aug 29 in New Orelans
time_window_selected <- 100:123
n_time_event <- 9 # for 108

date_range <- c("2001.01.01","2010.12.31",16) #PARAM 15, NDVI Katrina

#sbt_results_filename <- "/home/bparmentier/Google Drive/Space_beats_time/outputs/output_NDVI_Katrina_08242017/mae_tot_tb_NDVI_Katrina_08242017.txt"

#debug(processing_modis_data)

#function_raster_processing <-"MODIS_and_raster_processing_functions_02262018.R"
#function_processing_modis_data <-"processing_MODIS_data_functions_02262018.R"

#script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts"  #path to script functions

#source(file.path(script_path,function_raster_processing)) #source all functions used in this script.
#source(file.path(script_path,function_processing_modis_data)) #source all functions used in this script.

#https://gis.stackexchange.com/questions/237272/mean-by-month-on-r-stacked-raster

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

#in_dir_rast <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_NDVI_Katrina"
in_dir_rast <- "/media/dan/Space_beats_time/Space_beats_time/Data/data_NDVI_Katrina"

lf <-list.files(pattern="NDVI_Katrina_08092017.tif",in_dir_rast,full.names=T)
lf <- mixedsort(lf)
r_stack <- stack(lf)

r_var <- subset(r_stack,time_window_selected)

n_layers <- nlayers(r_var)
r_var_names <- paste0("time_",1:n_layers)
names(r_var) <- r_var_names

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

https://rpubs.com/mohammadshadan/288218

########## Generate figures and metrics ##########

# Create the object data using 5 random numbers
data <- rnorm(5)

# Create dates as a Date class object starting from 2016-01-01
dates <- seq(as.Date("2016-01-01"), length = 5, by = "days")

# Use xts() to create smith
smith <- xts(x = data, order.by = dates)

# Create bday (1899-05-08) using a POSIXct date class object
bday <- as.POSIXct("1899-05-08")

# Create hayek and add a new attribute called born
hayek <- xts(x = data, order.by = dates, born = bday)

range_dates <- generate_dates_by_step(date_range[1],date_range[2],16)$dates #NDVI Katrina

r_ts <- subset(r_stack,1:230)
setZ(r_ts,range_dates)
length(range_dates)

r_monthly_mean <- zApply(r_ts, by=month,fun=mean)

################################ END OF SCRIPT ##############################

