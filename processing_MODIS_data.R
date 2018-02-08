##################################################  PROCESSING MODIS DATA  #######################################
########################################### TILE PROCESSING OF MODIS PRODUCTS #####################################
#The current version implements a workflow for processing MODIS products.
#This script downloads and processes MODIS tiles using Quality Flag. 
#Tiles are mosaiced and reprojected for a specific study region.
#MODIS currently stores information in HDF4 format. Layers must be extracted and must be listed first
#using for example gdalinfo to identify the relevant dataset and QC flag options. 
#Note that QC flags are store bitpacks of 8bits (byte) in big endian!!!
#A data frame matching flag values is created to facilate the processing.            
#Inspiration and some code for the MODIS flag function originates from Steve Mosher:
#http://stevemosher.wordpress.com/2012/12/05/modis-qc-bits/
## MODIS WORKFLOW
# Processing of MODIS HDF files is done in 5 steps:
# Step 1: download modis tiles for specified product and version (e.g. version 5)
# Step 2: import modis tiles for specified file format (e.g. ".tif",".rst)
# Step 3: deal with modis flags (multiple levels)
# Step 4: mosaic tiles for every time step
# Step 5: reproject and crop extent to study region
#
#AUTHOR: Benoit Parmentier                                                                       
#CREATED ON : 09/16/2013  
#MODIFIED ON : 02/08/2018
#PROJECT: General MODIS processing of all projects
#COMMIT: debugging mosaicing step
#
#TODO: 
#1)Test additional Quality Flag levels for ALBEDO and other products (MOD09)
#2) Add function to report statistics: missing files
#3) Currently 20 input arguments (param), reduce to 15 or less
#4) Make this script a function callable from shell!!
#5) This script can be transformed to process other datasets using the https://lpdaac.usgs.gov/data_access/data_pool
#   e.g."https://e4ftl01.cr.usgs.gov/WELD/" for WELD datasets.
#6) Update mosaicing to use gdal for large areas


###################################################################################################

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

function_raster_processing <-"MODIS_and_raster_processing_functions_02072018.R"
function_processing_modis_data <-"processing_MODIS_data_functions_02082018e.R"

script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts"  #path to script functions

source(file.path(script_path,function_raster_processing)) #source all functions used in this script.
source(file.path(script_path,function_processing_modis_data)) #source all functions used in this script.

################################
###### Parameters and arguments

#ARG1
in_dir <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_Harvey_NDVI" #param1
#ARG2
out_dir <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_Harvey_NDVI" #param2
#ARG3
#http://spatialreference.org/ref/epsg/nad83-texas-state-mapping-system/proj4/
CRS_reg <- "+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" 
#ARG4
file_format <- ".tif" #raster format used #param4
#ARG5
NA_flag_val <- -9999
#ARG6
out_suffix <- "harvey_02082018"
#ARG7
create_out_dir_param=FALSE #param7
#ARG8
infile_reg_outline <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_Harvey_NDVI/revised_area_Rita/new_strata_rita_10282017.shp"
#ARG9
#local raster name defining resolution, extent
ref_rast_name <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_Harvey_NDVI/rev_project_output/mosaiced_MOD13A2_A2010353__006_1_km_16_days_NDVI_houston_10222017_crop_proj_reg_rev.rst"
#ref_rast_name <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_Harvey_NDVI/rita_outline_reg/Study_Area_Rita_New.shp"
#ARG10
MODIS_product <- "MOD13A2.006" #NDVI/EVI 1km product (monthly) #param12
#ARG11
date_param <- "2017.01.01;2017.12.31;16" #start date, end date, time_step
#ARG12
list_tiles_modis <- NULL #if NULL determine the tiles to download
#ARGS13
#scaling_factors <- c(1,-273.15) #set up as slope (a) and intercept (b), if NULL, no scaling done, setting for LST 
scaling_factors <- c(0.0001,0) #set up as slope (a) and intercept (b), if NULL, no scaling done, setting for NDVI 
#ARGS14
product_type = c("NDVI") #can be LST, ALBEDO etc.#this can be set from the modis product!! #param 19
#ARGS15
var_name <- "LST_Night_1km" #can be LST_Day_1km
#ARGS16
qc_name <- "QC_Night"
#ARGS17
num_cores <- 4 #param 20
#ARGS18
selected_flags <- list(QA_word1 ="VI Good Quality",QA_word1 ="VI Produced,check QA") #if NULL use default
#Select level 2:
#qc_product_l2_valid <- list(x=qc_lst_valid,QA_word2 %in% unique(QC_data_ndvi$QA_word2)[1:8]) #"Highest quality, 1","Lower quality, 2","Decreasing quality, 3",...,"Decreasing quality, 8" 
#ARGS19
agg_param <- c(FALSE,NULL,"mean") #False means there is no aggregation!!! #param 11
#ARG20
steps_to_run <- list(download=TRUE,        #1rst step
                     import=TRUE,          #2nd  step
                     apply_QC_flag=TRUE,   #3rd  step
                     mosaic=TRUE,          #4th  step
                     reproject=TRUE)       #5th  step 
### Constants

#Constant1:
proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
#Constant2:
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 
#Constatn3:
file_format_download <- "hdf"  
#Constant4:
infile_modis_grid <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_Harvey_NDVI/modis_reference_grid/modis_sinusoidal_grid_world.shp" #param11
#Constant5:
save_textfile <- TRUE

## Maybe add in_dir for all the outputs
#Constant6:
###Maybe have a text file??
download_dir <- NULL #step 1, multiple folders
import_dir <- NULL # step 2, multiple folders
mask_qc_dir <- NULL # step 3, this should be a list?, better as a text input
mosaic_dir <- NULL # step 4
project_dir <- NULL # step 5

out_dir_processing_steps <- list(download_dir,import_dir,mask_qc_dir,mosaic_dir,project_dir)
names(out_dir_processing_steps) <- c("download_dir","import_dir","mask_qc_dir","mosaic_dir","project_dir")

#debug(processing_modis_data)

processing_modis_data(in_dir,
                      out_dir,
                      CRS_reg,
                      file_format, 
                      NA_flag_val,
                      out_suffix,
                      create_out_dir_param,
                      infile_reg_outline,
                      ref_rast_name, 
                      MODIS_product,
                      date_param,
                      list_tiles_modis,
                      scaling_factors,
                      product_type,
                      var_name,
                      qc_name,
                      num_cores,
                      selected_flags, 
                      agg_param,
                      steps_to_run,
                      proj_modis_str,
                      CRS_WGS84,
                      file_format_download,
                      infile_modis_grid,
                      save_textfile,
                      out_dir_processing_steps)


################################ END OF SCRIPT ##############################

