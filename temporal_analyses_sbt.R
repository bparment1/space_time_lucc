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
library(zoo)
library(xts)

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
setZ(r_ts,range_dates,"16-day")
length(range_dates)
nlayers(r_ts)
nlayers(r_ts)
names(r_ts)
r_monthly_mean <- zApply(r_ts, by=as.yearmon,fun=mean)

class(r_ts)

library(MODIS)
interpolation/global_product_assessment_part0_functions.R

time_period <- "month"

lf_var <- lf[1:230]

generate_aggregate_raster <- function(i,lf_periods_vals,fun_val= "mean",file_format=".itf",out_rastname=NULL, out_dir="NULL",out_suffix=NULL){
  
  lf_in <- lf_periods_vals[[i]]
  r_s <- stack(lf_in)
  
  #out_rastname <- "test.tif"
  if(!is.null(out_rastname)){
    rast_name <- out_rastname[i]
  }else{
    rast_name <- paste0("agg_",i,"_",out_suffix,file_format)
  }
  
  if(is.null(out_dir)){
    out_dir <- "."
  }
  
  rast_name <- file.path(out_dir,rast_name)
  if(is.character(fun_Val)){
    fun_val <- get(fun_val)
  }
  
  r_temp_agg <- raster::calc(r_s,
                fun=fun_val,
                filename=rast_name,
                overwrite=T)
  
  return(rast_name)
}


get_time_interval <- function(date_vals,time_period,lf_var){
  ##takes time series objects and returns locations of observation in the last observations (function from xts)
  #endpoints(range_dates, on="months")
  #endpoints(range_dates, on="weeks",k=2)
  
  index_periods <- endpoints(date_vals,on=time_period,lf_var)
  
  ## generate tables of date with columns: dates year, quarter,month,week, doy

  list_date_vals <- as.Date(strptime(date_vals, "%Y%m%d"))
  month_str <- format(list_date_vals, "%b") ## Month, char, abbreviated
  year_str <- format(list_date_vals, "%Y") ## Year with century
  day_str <- as.numeric(format(list_date_vals, "%d")) ## numeric month
  
  df_files <- data.frame(lf = lf_var,
                         date = list_date_vals,
                         month_str = month_str,
                         year = year_str,
                         day = day_str,
                         dir = dirname(lf_var))
  
  lf_period_vals <- vector("list",length=length(index_periods)/2)
  
  for(i in 2:length(lf_var)){
    index_selected_t1 <- index_periods[i-1] + 1
    index_selected_t2 <- index_periods[i]
    
    lf_period_vals[[i-1]] <- c(lf_var[index_selected_t1],lf_var[index_selected_t2])
  }
  
  #class(range_dates)
  df_dates <- zoo(1:length(range_dates),range_dates)
  
  if(time_period=="month"){
    #dates_agg <- apply.monthly(df_dates,mean)
    #generate dates here:
    as.yearmon(date_vals,"%b-%Y")
  }
  
  if(time_period=="week"){
    dates_agg <- apply.weekly(range_dates,mean)
  }
  
  return(lf_period_vals,dates_agg)
}

#}else{ #generate empty data.frame
#  #declare data.frame with specific type, with zero row initialization
#  df_files <- data.frame(lf = character(),
#                         date = as.Date(character()), 
#                         month_str = character(), 
#                         year = character(),
#                         day = character(),
#                         dir = character()) 
#  list_dates_produced_date_val <- as.Date(character()) #declare date of raster predictions, empty
#}

temporal_aggregation <- function(date_vals,time_period,r_stack){
  # 1) range_dates
  # 2) time_period: month, week, year,doy
  # 3) r_stack
  
  #### Start script ########
  
  if(inherits(r_stack,"raster")){
    lf_var <- names(r_stack)
  }else{
    lf_var <- r_stack
  }

  test <- get_time_interval(date_vals,time_period,lf_var)
  
  ### now aggregate
  
  i <- 1
  browser()
  debug(generate_aggregate_raster)
  generate_aggregate_raster(i,
                            lf_period_vals,
                            fun_val= "mean",
                            out_rastname=NULL, 
                            out_dir="NULL",
                            out_suffix=NULL)
  
  mclapply(1:lenght(list_period_vals),
           Fun)
  
  selected_periods <- 
  time_window_selected
  
}



#: separate in month and day range_dates 
################################ END OF SCRIPT ##############################

