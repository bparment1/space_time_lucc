########################################  TILES PROCESSING #######################################
########################################### Read QC flags and mosaic tiles in R #####################################
#The current version generate data forthe Arizona region.
#This script download and processes MODIS tiles using Quality Flag. 
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
#MODIFIED ON : 10/22/2017
#PROJECT: General MODIS processing of all projects
#COMMIT: modification of dowloading function
#
#TODO: 
#1)Test additional Quality Flag levels for ALBEDO and other product
#2) Add function to report statistics: missing files
#3) Currently 20 input arguments (param), reduce to 15 or less
#4) Make this script a function callable from shell!!
#5) This script can be transform to process other datasets using the https://lpdaac.usgs.gov/data_access/data_pool
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

function_analyses_paper <-"MODIS_and_raster_processing_functions_10222017.R"
script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts"  #path to script functions
source(file.path(script_path,function_analyses_paper)) #source all functions used in this script.

################################
###### Parameters and arguments

#in_dir <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_AZ_jacob" #param1
in_dir <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_Rita_NDVI"
#data_fname <- file.path("~/Data/Space_beats_time/stu/Katrina/run2/csv","Katrina2.csv")
out_dir <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_Rita_NDVI" #param2

proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" 
#CRS_reg <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 

proj_str<- CRS_WGS84 #check if in use somewhere?
#http://spatialreference.org/ref/epsg/nad83-texas-state-mapping-system/proj4/
CRS_reg <- "+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" 
file_format <- ".rst" #raster format used #param4
NA_value <- -9999 #param5
NA_flag_val <- NA_value
out_suffix <-"houston_10222017" #output suffix for the files that are masked for quality and for ...param6
#out_suffix <- "arizona_10182017"
create_out_dir_param=FALSE #param7

#in_dir <- "/data/project/layers/commons/modis/MOD11A1_tiles" #ATLAS SERVER 
#out_dir <- "/data/project/layers/commons/Oregon_interpolation/MODIS_processing_07072014/" #ATLAS SERVER 

#infile_reg_outline=""  #input region outline defined by polygon: none Katrina
infile_reg_outline=NULL #param9
infile_reg_outline<- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_Rita_NDVI/rita_outline_reg/Study_Area_Rita_New.shp"
#This is the shape file of outline of the study area                                                      #It is an input/output of the covariate script
#infile_reg_outline <- "/data/project/layers/commons/Oregon_interpolation/MODIS_processing_07072014/region_outlines_ref_files/OR83M_state_outline.shp" #input region outline defined by polygon: Oregon

#local raster name defining resolution, exent: oregon

#ref_rast_name <- "~/Data/Space_beats_time/Case2_data_NDVI/ref_rast_New_Orleans.rst"
ref_rast_name <- NULL #if null use the first image to define projection area
#ref_rast_name <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_Rita_Houston/rita_outline_reg/Study_Area_Rita_New.shp"
#ref_rast_name <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_AZ_jacob/reference_region_study_area_AZ.rst"
#ref_rast_name <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_AZ_jacob/Arizona_Outline_State_Plane/Arizona_Outlline_State_Plane.shp"
#ref_rast_name<-"/home/parmentier/Data/IPLANT_project/MODIS_processing_0970720134/region_outlines_ref_files/mean_day244_rescaled.rst" #local raster name defining resolution, exent: oregon
infile_modis_grid <- "/home/bparmentier/Google Drive/Space_beats_time/Data/modis_reference_grid/modis_sinusoidal_grid_world.shp" #param11

## Other specific parameters

MODIS_product <- "MOD13A2.006" #NDVI/EVI 1km product (monthly) #param12
#MODIS_product <- "MOD11A1.006"
#MODIS_product <- "MOD11A2.006" #should be product name
start_date <- "2001.01.01"  #param13
end_date <- "2010.12.31"  #param14
#end_date <- "2001.01.10"

#/home/bparmentier/Google Drive/Space_beats_time/Data/modis_reference_grid
#list_tiles_modis<- NULL #if NULL, determine tiles using the raster ref or reg_outline file  #param14
#list_tiles_modis<- c("h10v05,h10v06") # for Rita
#list_tiles_modis <- c("h08v05")
list_tiles_modis <- NULL #if NULL determine the tiles to download

file_format_download <- "hdf"  #param15
product_version <- 6 #param16
#temporal_granularity <- "Daily" #deal with options( 16 day, 8 day and monthly) #param17
#temporal_granularity <- "16 Day" #deal with options( 16 day, 8 day and monthly), unused at this stage...

#scaling_factors <- c(1,-273.15) #set up as slope (a) and intercept (b), if NULL, no scaling done, setting for LST 
scaling_factors <- c(0.0001,0) #set up as slope (a) and intercept (b), if NULL, no scaling done, setting for NDVI 
#scaling_factors <- NULL #set up as slope (a) and intercept (b), if NULL, no scaling done #param18 

product_type = c("NDVI") #can be LST, ALBEDO etc.#this can be set from the modis product!! #param 19
#product_type = c("LST") #can be LST, ALBEDO etc.

num_cores <- 4 #param 20

#Parameters/arguments not in use yet...
#num_cores <- 10 #number of cores used in parallel processing...
selected_flags <- list(QA_word1 ="VI Good Quality",QA_word1 ="VI Produced,check QA")
#Select level 2:
#qc_product_l2_valid <- list(x=qc_lst_valid,QA_word2 %in% unique(QC_data_ndvi$QA_word2)[1:8]) #"Highest quality, 1","Lower quality, 2","Decreasing quality, 3",...,"Decreasing quality, 8" 

agg_param <- c(FALSE,NULL,"mean") #False means there is no aggregation!!! #param 11

#run all processing steps
#param21
save_textfile <- TRUE

steps_to_run <- list(download=TRUE,       #1rst step
                     import=TRUE,          #2nd step
                     apply_QC_flag=TRUE,   #3rd step
                     mosaic=TRUE,          #4th step
                     reproject=TRUE)       #5th step 
### Constants

## Maybe add in_dir for all the outputs

###Maybe have a text file??
download_dir <- NULL #step 1, multiple folders
import_dir <- NULL # step 2, multiple folders
mask_qc_dir <- NULL # step 3, this should be a list?, better as a text input
mosaic_dir <- NULL # step 4
#mosaic_dir <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_AZ_jacob/mask_qc_h08v05"

project_dir <- NULL # step 5

######################################################
########################  BEGIN SCRIPT  #############

#Create output directory

if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

#####################################
#### STEP 1:  DOWNLOAD MODIS PRODUCT  ####

#Note that files are downloaded in the ouput directory in subdirectory with tile_name (e.g. h08v04)
#if(is.null(infile_reg_outline))
#  reg_outline<- create_polygon_from_extent(ref_rast_name,out_suffix)
#  #debug(get_modis_tiles_list)
#}
#list_tiles_modis <- get_modis_tiles_list(modis_grid,reg_outline,CRS_reg)

if(is.null(list_tiles_modis)){
  #
  #debug(get_modis_tiles_list)
  list_tiles_modis <- get_modis_tiles_list(modis_grid,
                                           reg_outline=infile_reg_outline,
                                           CRS_reg)
  
}
list_tiles_modis <- unlist(strsplit(list_tiles_modis,","))  # transform string into separate element in char vector

#debug(modis_product_download)
if(steps_to_run$download==TRUE){
  #debug(modis_product_download)
  download_modis_obj <- modis_product_download(MODIS_product,product_version,start_date,end_date,list_tiles_modis,file_format_download,out_dir,temporal_granularity)
  out_dir_tiles <- (file.path(in_dir,list_tiles_modis))
  list_files_by_tiles <-download_modis_obj$list_files_by_tiles #Use mapply to pass multiple arguments
  colnames(list_files_by_tiles) <- list_tiles_modis #note that the output of mapply is a matrix
}else{
  
  ### Add download_dir here
  out_dir_tiles <- (file.path(in_dir,list_tiles_modis))
  #list_files_by_tiles <- mapply(1:length(out_dir_tiles),FUN=list.files,MoreArgs=list(pattern="*.hdf$",path=out_dir_tiles,full.names=T),SIMPLIFY=T) #Use mapply to pass multiple arguments 
  list_files_by_tiles <-mapply(1:length(out_dir_tiles),FUN=function(i,x){list.files(path=x[[i]],pattern="*.hdf$",full.names=T)},MoreArgs=(list(x=out_dir_tiles)),SIMPLIFY=T) #Use mapply to pass multiple arguments
  colnames(list_files_by_tiles) <- list_tiles_modis #note that the output of mapply is a matrix
}

####################################
##### STEP 2: IMPORT MODIS LAYERS ###

##Modify this section into a function to extract names of product and quality flag automatically!!
infile_var <- list_files_by_tiles[[1]]
GDALinfo_hdf<-GDALinfo(infile_var[1],returnScaleOffset=FALSE)
str(GDALinfo_hdf)
modis_subdataset <- attributes(GDALinfo_hdf)$subdsmdata
print(modis_subdataset)

#Select automatically QC flag!!
if(product_type=="NDVI"){
  modis_layer_str1 <- unlist(strsplit(modis_subdataset[1],"\""))[3] #Get day NDVI layer
  modis_layer_str2 <- unlist(strsplit(modis_subdataset[5],"\""))[3] #Get day VI QC layer
}
if(product_type=="LST"){
  modis_layer_str1 <- unlist(strsplit(modis_subdataset[1],"\""))[3] #Get day LST layer
  modis_layer_str2 <- unlist(strsplit(modis_subdataset[3],"\""))[3] #Get day QC layer
}

### import list of files before mosaicing

file_format_import <- file_format
#NA_flag_val=-9999
var_modis_name <- unlist(strsplit(modis_layer_str1,":"))[3]
qc_modis_name <- unlist(strsplit(modis_layer_str2,":"))[3]

var_modis_name <- gsub(" ","_",var_modis_name) #suffix name for product, may contain white space so replace with "_"
qc_modis_name <- gsub(" ","_",qc_modis_name)

##loop over tiles:
#Took 10 minutes for 506 files and one tile
if(steps_to_run$import==TRUE){
  list_imported_files <- vector("list",length=length(list_tiles_modis))
  for(j in 1:length(list_tiles_modis)){
    #infile_var <- download_modis_obj$list_files_by_tiles[,j] 
    infile_var <-list_files_by_tiles[,j] #note can be any variable even thought LST presented  here
    #infile_var <- list_files_by_tiles[[j]]
    out_dir_tmp <- paste0("import_",list_tiles_modis[j])
    #out_dir_s <- file.path(out_dir,list_tiles_modis[j])
    out_dir_s <- file.path(out_dir,out_dir_tmp)
    
    out_suffix_s <- var_modis_name
    list_param_import_modis <- list(i=1,hdf_file=infile_var,subdataset=modis_layer_str1,NA_flag_val=NA_flag_val,out_dir=out_dir_s,
                                    out_suffix=out_suffix_s,file_format=file_format_import,scaling_factors=scaling_factors)
    #undebug(import_list_modis_layers_fun)
    #r_var_s_filename <- import_list_modis_layers_fun(1,list_param_import_modis)    
    #r_var_s <- raster(r_var_s_filename)
    #r_var_s <- mclapply(1:12,
    #                    FUN=import_list_modis_layers_fun,
    #                    list_param=list_param_import_modis,
    #                    mc.preschedule=FALSE,
    #                    mc.cores = num_cores) #This is the end bracket from mclapply(...) statement

    list_r_var_s <- mclapply(1:length(infile_var),
                        FUN=import_list_modis_layers_fun,
                        list_param=list_param_import_modis,
                        mc.preschedule=FALSE,
                        mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
    
    ##### Now do the qc flags
    out_suffix_s <- qc_modis_name
    list_param_import_modis <- list(i=1,
                                    hdf_file=infile_var,subdataset=modis_layer_str2,NA_flag_val=NA_flag_val,out_dir=out_dir_s,
                                    out_suffix=out_suffix_s,file_format=file_format_import,scaling_factors=NULL)
    #r1<-import_list_modis_layers_fun(1,list_param_import_modis)
    #r_qc_s <-mclapply(1:12,
    #                  FUN=import_list_modis_layers_fun,
    #                  list_param=list_param_import_modis,
    #                 mc.preschedule=FALSE,
    #                  mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
    
    list_r_qc_s <-mclapply(1:length(infile_var),
                      FUN=import_list_modis_layers_fun,
                      list_param=list_param_import_modis,
                      mc.preschedule=FALSE,
                      mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
    l_files <- list(var=list_r_var_s,qc=list_r_qc_s)
    list_imported_files[[j]] <- l_files
  }
}


if(steps_to_run$import==FALSE){
  ## Need to deal with multiple tiles
  list_imported_files <- vector("list",length=length(list_tiles_modis))
  for(j in 1:length(list_tiles_modis)){
    
    #j <-1
    #for(j in 1:length(list_tiles_modis)){
    if(is.null(import_dir)){
      out_dir_tmp <- paste0("import_",list_tiles_modis[j])
      #out_dir_s <- file.path(out_dir,list_tiles_modis[j])
      out_dir_s <- file.path(out_dir,out_dir_tmp)
    }else{
      out_dir_s <- import_dir
    }
  
    #out_dir_s <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_AZ_jacob/import_h08v05"
    #drwxrwxr-x 2 bparmentier bparmentier 401408 Oct  9 17:27 mask_qc_h08v05
  
    #LST_Day_1km_arizona_10092017.rst
    #MOD11A2_A2012361_h08v05_006_LST_Day_1km.rst
  
    #strsplit(x=MODIS_product, split="[.]")
    #list_r_var_s <- list.files(path=out_dir_s,
    #                     pattern=".*.LST_Day_1km.rst$")
  
  
    #file_pattern <- paste0(sub("[.]","_",MODIS_product),".*.",product_type,".*",file_format,"$")
    file_pattern <- paste0(".*.",product_type,".*",file_format,"$")
  
    list_r_var_s <- list.files(path=out_dir_s,
                             pattern=file_pattern,
                             full.names=T)
    file_pattern <- paste0(".*.","QC",".*",file_format,"$")
    #MOD11A2_A2002001_h08v05_006_QC_Day.rst
  
    list_r_qc_s <- list.files(path=out_dir_s,
                            pattern=file_pattern,
                            full.names=T)
    
    l_files <- list(var=list_r_var_s,qc=list_r_qc_s)
    list_imported_files[[j]] <- l_files
  }
}

### Should report on errors:
#import_error_list <- unlist(lapply(1:length(r_qc_s),FUN=function(x){class(x)=="try-error"}))
#import_error_list[as.numeric(import_error_list)==1]
#sum(as.numeric(import_error_list))

list_r_var_s <- unlist(list_r_var_s) #list of files as character vector
list_r_qc_s <- unlist(list_r_qc_s) #list of files as character vector
list_r_var_s1 <- list_imported_files[[1]]$var #first
list_r_var_s2 <- list_imported_files[[2]]$var #first

plot(raster(list_r_qc_s[1]))
plot(raster(list_r_var_s1[[2]]))
plot(raster(list_r_var_s2[[2]]))

#print(r_var_s)
#print(r_qc_s)

#################################
##### STEP 3: APPLY/DEAL WITH QC FLAG AND SCREEN VALUE FOR VALID RANGE ###

## Get QC information for lST/NDVI and mask values: imporove and automate this later
if(product_type=="NDVI"){
  QC_obj <- create_MODIS_QC_table(LST=FALSE, NDVI=TRUE) #Get table corresponding to QC for LST
  names(QC_obj)
  QC_data_ndvi <- QC_obj$NDVI
  #For NDVI: use this section to process. This is the default processing quality, this should go top in the parameters!!!
  #Select level 1:
  qc_lst_valid <- subset(x=QC_data_ndvi,QA_word1 == "VI Good Quality" | QA_word1 =="VI Produced,check QA")
  #Select level 2:
  qc_lst_valid <- subset(x=qc_lst_valid,QA_word2 %in% unique(QC_data_ndvi$QA_word2)[1:8]) #"Highest quality, 1","Lower quality, 2","Decreasing quality, 3",...,"Decreasing quality, 8" 

  #Select level 3:
  #...Not implemented at this stage

}

## Get QC information for lST/NDVI and mask values: imporove and automate this later
if(product_type=="LST"){
  #undebug(create_MODIS_QC_table)
  QC_obj <- create_MODIS_QC_table(LST=TRUE, NDVI=FALSE) #Get table corresponding to QC for LST
  names(QC_obj)
  #For LST: Use default value:
  QC_data_lst <- QC_obj$LST
  #Select level 1:
  qc_lst_valid <- subset(x=QC_data_lst,QA_word1 == "LST Good Quality" | QA_word1 =="LST Produced,Check QA")
  #Select level 2:
  qc_lst_valid <- subset(x=qc_lst_valid,QA_word2 == "Good Data" | QA_word2 =="Other Quality")
  #Select level 3:
  #...

}
#QC_obj <- create_MODIS_QC_table(LST=TRUE, NDVI=TRUE) #Get table corresponding to QC for LST

names(qc_lst_valid)
qc_valid<-qc_lst_valid$Integer_Value #value integer values
#NA_flag_val <- -9999

#r_lst_LST <- download_modis_obj$list_files_by_tiles[,j]
list_r_lst <- vector("list",length(list_tiles_modis)) #to contain image
list_r_qc <- vector("list",length(list_tiles_modis)) #to contain qc mask image
list_r_stack <- vector("list",length(list_tiles_modis)) #to contain results

#26 minutes for 230 files to apply NDVI mask
if(steps_to_run$apply_QC_flag==TRUE){
  for(j in 1:length(list_tiles_modis)){
    #list all files first
    #out_dir_s <- file.path(dirname(out_dir),list_tiles_modis)[j]
    #
    #out_dir_s <- file.path(out_dir,list_tiles_modis)[j]
    
    #out_dir_tmp <- paste0("mask_qc_",list_tiles_modis[j])
    out_dir_tmp <- paste0("import_",list_tiles_modis[j])
    #out_dir_s <- file.path(out_dir,list_tiles_modis[j])
    out_dir_s <- file.path(out_dir,out_dir_tmp) #input dir is import out dir
    
    out_suffix_s <- paste(var_modis_name,sep="") #for MODIS product (var)
    file_format_s <-file_format
    #undebug(create_raster_list_from_file_pat)
    list_r_lst[[j]] <- create_raster_list_from_file_pat(out_suffix_s,file_pat="",
                                                        in_dir=out_dir_s,
                                                        out_prefix="",
                                                        file_format=file_format_s)
    
    #out_dir_s <- file.path(dirname(out_dir),list_tiles_modis)[j]
    #out_dir_s <- file.path(out_dir,list_tiles_modis)[j] #same as above
    #debug(create_raster_list_from_file_pat)
    out_suffix_s <- paste(qc_modis_name,sep="") #for qc flags
    list_r_qc[[j]] <- create_raster_list_from_file_pat(out_suffix_s,
                                                       file_pat="",
                                                       in_dir=out_dir_s,
                                                       out_prefix="",
                                                       file_format=file_format_s)
    
    ##### Now prepare the input for screening using QC flag values
    
    ##Set output dir
    out_dir_tmp <- paste0("mask_qc_",list_tiles_modis[j])
    #out_dir_s <- file.path(out_dir,list_tiles_modis[j])
    out_dir_s <- file.path(out_dir,out_dir_tmp) #input dir is import out dir
    
    list_param_screen_qc <- list(qc_valid,
                                 list_r_qc[[j]], 
                                 list_r_lst[[j]],
                                 rast_mask=TRUE,
                                 NA_flag_val,out_dir_s,out_suffix) 
    names(list_param_screen_qc) <- c("qc_valid",
                                     "rast_qc", 
                                     "rast_var","rast_mask",
                                     "NA_flag_val","out_dir","out_suffix") 
    #undebug(screen_for_qc_valid_fun)
    #r_stack[[1]] <- screen_for_qc_valid_fun(1,list_param=list_param_screen_qc)
    #r_stack[[j]] <- lapply(1:length(list_r_qc[[j]]),FUN=screen_for_qc_valid_fun,list_param=list_param_screen_qc)
    #r_test <-mclapply(1:11,FUN=screen_for_qc_valid_fun,list_param=list_param_screen_qc,mc.preschedule=FALSE,mc.cores = 11) #This is the end bracket from mclapply(...) statement
    
    list_r_stack[[j]] <-mclapply(1:length(list_r_qc[[j]]),
                            FUN=screen_for_qc_valid_fun,
                            list_param=list_param_screen_qc,
                            mc.preschedule=FALSE,
                            mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
    
  }
}
#33 minutes for 4 tiles of NDVI 230 images
#26 minutes for 230 for NDVI
#19minutes for 505 files
#r_lst_by_tiles <-mapply(1:length(list_tiles_modis),FUN=list.files,pattern=paste".*.day_LST.*.rst$",path=out_dir_s,full.names=T) #Use mapply to pass multiple arguments
#r_lst <- mapply(1:length(out_suffix_s),FUN=create_raster_list_from_file_pat,
 #               file_pat="",in_dir=out_dir_s,out_prefix="",file_format=file_format_s)

### Now loop over a series of files...
#extract_list_from_list_obj(unlist(r_stack),"var")
  
#################################
##### STEP 4: MOSAIC TILES  ###

#debug(create_raster_list_from_file_pat)
if(steps_to_run$mosaic==TRUE){
  
  if(is.null(mosaic_dir)){

    #out_dir_tmp <- paste0("mosaic_",out_suffix)
    out_dir_tmp <- paste0("mosaic_output")
    
    #out_dir_s <- file.path(out_dir,list_tiles_modis[j])
    out_dir_s <- file.path(out_dir,out_dir_tmp)
  }else{
    out_dir_s <- mosaic_dir
  }
  
  list_m_var <- vector("list",length(list_tiles_modis))  
  l_df_raster_name <- vector("list",length(list_tiles_modis))  

  names(list_m_var)<- list_tiles_modis
  list_m_qc <- vector("list",length(list_tiles_modis))  
  names(list_m_qc)<- list_tiles_modis

  for (j in 1:length(list_tiles_modis)){
    #out_suffix_s <- paste(list_tiles_modis[j],"_",sprintf( "%03d", product_version),"_",var_modis_name,"_",out_suffix,sep="")
    #file_format_s <-file_format
    #out_suffix_s <- paste(list_tiles_modis[j],"_",sprintf( "%03d", product_version),"_",var_modis_name,"_",out_suffix,file_format_s,sep="")
    #out_dir_s <- file.path(dirname(out_dir),list_tiles_modis)[j]
    #out_dir_s <- file.path(out_dir,list_tiles_modis)[j]
    #list_m_var[[j]]<-list.files(pattern=paste(out_suffix_s,"$",sep=""),
    #                            path=out_dir_s,
    #                            full.names=TRUE) #inputs for moasics
    file_pattern <- paste0(".*.",product_type,".*.",
                           out_suffix,file_format,"$")
    
    in_dir_tmp <- paste0("mask_qc_",list_tiles_modis[j])
    #out_dir_s <- file.path(out_dir,list_tiles_modis[j])
    in_dir_s <- file.path(out_dir,in_dir_tmp) #input dir is import out dir
    
    #list_r_var_s <- list.files(path=out_dir_s,
    #                           pattern=file_pattern,
    #                           full.names=T)
    
    #MOD11A2_A2012353_h08v05_006_LST_Day_1km_arizona_10092017.rst
    #list_var_mosaiced <-list.files(pattern=".*.LST_Day_1km_arizona_10092017.rst$",
    #                               path=out_dir_s,
    #                               full.names=TRUE) #inputs for moasics
    
    list_r_var_s <-list.files(pattern=file_pattern,
                                   #patter=".*.NDVI.*.arizona_10182017.rst$",
                                   path=in_dir_s,
                                   full.names=TRUE) #inputs for moasics
    list_m_var[[j]] <- list_r_var_s
    #list_m_var[[j]] <-create_raster_list_from_file_pat(out_suffix_s,file_pat="",in_dir=out_dir_s,out_prefix="",file_format=file_format_s)
    df_m_var <- lapply(1:length(list_m_var[[j]]),
                       FUN=extract_dates_from_raster_name,
                       list_files=list_m_var[[j]])
    #df_m_var <- lapply(1:length(list_m_var),
    #                   FUN=extract_dates_from_raster_name,
    #                   list_files=list_m_var)
    df_m_var <- do.call(rbind,df_m_var)
    names(df_m_var) <- c(paste("raster_name",j,sep="_"),"date")
    
    df_m_var[,1] <- as.character(df_m_var[,1])
    l_df_raster_name[[j]] <- df_m_var
    
  }
  
  #test <- merge_all(l_df_raster_name,by="date") #,all.x=T,all.y=T) #does not work properly since df have to be order from most to least complete!!!
  #test <- merge(l_df_raster_name[[1]],l_df_raster_name[[2]],by="date",all.x=T,all.y=T)
  df_m_mosaics <- merge_multiple_df(l_df_raster_name,"date")
  x <- subset(df_m_mosaics,select= -c(date))# drop column with name "date"
    
  #report on missing dates:
  #st <- as.Date(start_date,format="%Y.%m.%d")
  #en <- as.Date(end_date,format="%Y.%m.%d")
  #ll <- seq.Date(st, en, by="1 day")
  #dates_queried <- format(ll,"%Y.%m.%d")
  #mosaic_list_var <-mapply(FUN="c",list_m_var,SIMPLIFY=T)
  #x <-mapply(FUN=as.matrix,list_m_var,SIMPLIFY=T)
    
  #x <-mapply(FUN="c",x,SIMPLIFY=T)
  #MODIS_product <- "MOD13A2.005" #NDVI/EVI 1km product (monthly) #param12
  #strsplit(MODIS_product,"[.]")
  #MODIS_product <- "MOD11A1.005"
  MODIS_product_name <- gsub("[.]","_",MODIS_product)
  date_str <- df_m_mosaics$date
  df_m_mosaics$out_rastnames_var <- paste(MODIS_product_name,date_str,"mosaic",product_type,out_suffix,sep="_") #no file format added!
  mosaic_list_var <- lapply(seq_len(nrow(x)), function(i){x[i,]}) #list of tiles by batch to mosaic
  #Prepare list of output names without extension
  out_rastnames_var <- (basename(gsub(list_tiles_modis[1],"",list_m_var[[1]])))
  out_rastnames_var <- gsub(extension(out_rastnames_var),"",out_rastnames_var)
  #out_rastnames_var <- df_m_mosaics$out_rastnames_var
    
  j <- 1
  
  #out_dir_mosaic <-     
  #out_dir_tmp <- paste0("mosaic_output_",out_suffix)
  out_dir_tmp <- paste0("mosaic_output")
  #out_dir_s <- file.path(out_dir,list_tiles_modis[j])
  out_dir_s <- file.path(out_dir,out_dir_tmp)
  if(!file.exists(out_dir_s)){
    dir.create(out_dir_s)
  }
  
  out_dir_mosaic <- out_dir_s
  #list_param_mosaic<-list(j,mosaic_list_var,out_rastnames_var,out_dir,file_format,NA_flag_val)
  list_param_mosaic<-list(j,
                          mosaic_list_var,
                          out_rastnames_var,
                          out_dir_mosaic,
                          file_format,
                          NA_flag_val)
    
  names(list_param_mosaic)<-c("j",
                              "mosaic_list",
                              "out_rastnames",
                              "out_path",
                              "file_format",
                              "NA_flag_val")
  #debug(mosaic_m_raster_list)
  #list_var_mosaiced <- mosaic_m_raster_list(1,list_param_mosaic)
  #list_var_mosaiced <-mclapply(1:11, 
  #                             list_param=list_param_mosaic, 
  #                             mosaic_m_raster_list,
  #                             mc.preschedule=FALSE,
  #                             mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
  
  #Parallelization
  #started at 22:07 -22:11
  list_var_mosaiced <-mclapply(1:length(mosaic_list_var), 
                               list_param=list_param_mosaic, 
                               mosaic_m_raster_list,
                               mc.preschedule=FALSE,
                               mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
    
  #list_var_mosaiced <- lapply(1:length(mosaic_list_var), list_param=list_param_mosaic, mosaic_m_raster_list) #This is the end bracket from mclapply(...) statement
    
  #r_test <- stack(list_var_mosaiced)
  #plot(r_test,y=1:2)
    

}

if(steps_to_run$mosaic==FALSE){
  
  if(is.null(mosaic_dir)){
    out_dir_tmp <- paste0("mosaic_output_",out_suffix)
    #out_dir_s <- file.path(out_dir,list_tiles_modis[j])
    out_dir_s <- file.path(out_dir,out_dir_tmp)
  }else{
    out_dir_s <- mosaic_dir
    #out_dir_s <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_AZ_jacob/import_h08v05"
    
  }
  
  #out_dir_s <- out_dir_mosaic
  #list_var_mosaiced <-    
  #list_m_var[[j]]<-list.files(pattern=paste(out_suffix_s,"$",sep=""),path=out_dir_s,full.names=TRUE) #inputs for moasics
  file_pattern <- paste0(".*.",product_type,".*.",
                         out_suffix,file_format,"$")
  
  #list_r_var_s <- list.files(path=out_dir_s,
  #                           pattern=file_pattern,
  #                           full.names=T)
  
  #MOD11A2_A2012353_h08v05_006_LST_Day_1km_arizona_10092017.rst
  #list_var_mosaiced <-list.files(pattern=".*.LST_Day_1km_arizona_10092017.rst$",
  #                               path=out_dir_s,
  #                               full.names=TRUE) #inputs for moasics
  
  list_var_mosaiced <-list.files(pattern=file_pattern,
                                 path=out_dir_s,
                                 full.names=TRUE) #inputs for moasics
  
}

#################################
##### STEP 5: REPROJECT AND CROP TO STUDY REGION  ###

if(steps_to_run$reproject==TRUE){
  
  # FIRST SET UP STUDY AREA ####
  
  # NOW PROJECT AND CROP WIHT REF REGION ####
  
  if(is.null(project_dir)){
    #out_dir_tmp <- paste0("project_output_",out_suffix)
    out_dir_tmp <- paste0("project_output")
    
    #out_dir_s <- file.path(out_dir,list_tiles_modis[j])
    out_dir_s <- file.path(out_dir,out_dir_tmp)
    if(!file.exists(out_dir_s)){
      dir.create(out_dir_s)
    }
  }else{
    out_dir_s <- project_dir
    #out_dir_s <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_AZ_jacob/import_h08v05"
    if(!file.exists(out_dir_s)){
      dir.create(out_dir_s)
    }
  }
  
  if(is.null(ref_rast_name)){
    
    #infile_reg_outline<- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_Rita_Houston/rita_outline_reg/Study_Area_Rita_New.shp"
    #Use one mosaiced modis tile as reference image...We will need to add a function 
    ref_rast_tmp <-raster(list_var_mosaiced[[1]]) 
    ref_rast_prj <-projectRaster(from=ref_rast_tmp,
                             res=res(ref_rast_tmp), #set resolution to the same as input
                             crs=CRS_reg,
                             method="ngb")
    #to define a local reference system and reproject later!!
    #Assign new projection system here in the argument CRS_reg (!it is used later)
    if(!is.null(infile_reg_outline)){
      reg_sf <- st_read(infile_reg_outline)
      reg_sf <- st_transform(reg_sf,crs=CRS_reg)
      reg_sp <-as(reg_sf, "Spatial") 
      ref_rast <- crop(ref_rast_prj,reg_sp)    
    }

  #Use the reference raster
  if(!is.null(ref_rast_name)){
    ref_rast<-raster(ref_rast_name) #This is the reference image used to define the study/processing area
    projection(ref_rast) <- CRS_reg #Assign given reference system from master script...
  }
  
  ##Create output names for region
  list_var_mosaiced_tmp <- remove_from_list_fun(list_var_mosaiced,condition_class ="try-error")$list
  
  out_suffix_var <-paste(out_suffix,file_format,sep="")          
  var_list_outnames <- change_names_file_list(list_var_mosaiced_tmp,out_suffix_var,"reg_",file_format,out_path=out_dir)     
  
  #list_param_create_region<-list(j,raster_name=list_var_mosaiced,reg_ref_rast=ref_rast,out_rast_name=var_list_outnames)
  #j<-1
  #list_param_create_region<-list(j,list_var_mosaiced,ref_rast,var_list_outnames,NA_flag_val)
  #names(list_param_create_region) <-c("j","raster_name","reg_ref_rast","out_rast_name","NA_flag_val")
  
  #undebug(create__m_raster_region)
  
  #### pasted from other file:
  out_rast_name <- NULL
  lf_r <- list_var_mosaiced
  list_param_create_region <- list(as.list(lf_r),
                                   ref_rast, out_rast_name,agg_param,
                                   file_format,NA_flag_val,
                                   input_proj_str=NULL,out_suffix="",out_dir_s)
  names(list_param_create_region) <- c("raster_name",
                                       "reg_ref_rast", "out_rast_name","agg_param",
                                       "file_format","NA_flag_val",
                                       "input_proj_str","out_suffix","out_dir")
  #debug(create__m_raster_region)
  #r_filename <- create__m_raster_region(1,list_param=list_param_create_region)
  #r <- raster(r_filename)
  reg_var_list <- mclapply(1:length(lf_r),
                       FUN=create__m_raster_region,
                       list_param=list_param_create_region,
                       mc.preschedule=FALSE,
                       mc.cores = num_cores)

  #reg_var_list <- mclapply(1:12,
  #                         FUN=create__m_raster_region,
  #                         list_param=list_param_create_region,
  #                         mc.preschedule=FALSE,
  #                         mc.cores = num_cores)
  
  ###Note: add x,y raster and mask defining the study area for the stack below!!
  
  lf <- list.files(path=out_dir_s,
                   pattern=paste0("*.",file_format,"$"),
                   full.names = T)
  write.table(lf,
              "raster_list_data.txt",
              sep=",",
              rownames=F)
  
  if(save_textfile==TRUE){
    
    r_reg_var<- stack(reg_var_list)
    
    dat_reg_var_spdf <- as(r_reg_var,"SpatialPointsDataFrame")
    ##66,806
    #drop NA
    dat_reg_var <- as.data.frame(dat_reg_var_spdf) 
    #Now
    
    #dat_out <- as.data.frame(r_reg_var)
    #dat_out <- na.omit(dat_out)
    out_filename <- paste0("dat_reg_var_list_",product_type,"_",out_suffix,".txt")
    write.table(dat_reg_var,
                file.path(out_dir_s,out_filename),
                row.names=F,sep=",",col.names=T)
    
    #write.table(dat_reg_var)
  }
  
  }
  
}

########################## END OF SCRIPT ##############################

### need to change this:
#r_test1 <- raster(create__m_raster_region(1,list_param_create_region))
#r_tmp <- raster(unlist(list_var_mosaiced)[1])
# r_test2 <- projectRaster(r_tmp,to=ref_rast,method="ngb")
# ### need to set the resolution
# r_test2 <- projectRaster(r_tmp,crs=CRS_reg, res=res(r_tmp),method="ngb")
# test_sp <- as(test_sf, "Spatial") 
# r_crop <- crop(r_test2,test_sp)
# ## Now turn the polygon file to mask using the reprojected and crop area.
# r_az <- rasterize(test_sp,r_crop)
# compareCRS(r_az,r_crop)
# projection(r_az) <- projection(r_crop)
# compareCRS(r_az,CRS_reg) #TRUE
# projection(r_az) <- CRS_reg
# #https://www.nceas.ucsb.edu/scicomp/recipes/projections
# NAvalue(r_az) <- -9999
# writeRaster(r_az,file.path("reference_region_study_area_AZ.rst"),
#             overwrite=T)
# writeRaster(r_az,file.path("reference_region_study_area_AZ.tif"),
#             overwrite=T)

#Currently we use only this:


###
#"gdalwarp -overwrite -t_srs '+proj=tmerc +lat_0=31 +lon_0=-111.9166666666667 +k=0.9999 +x_0=213360 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048 +no_defs' -of RST /home/bparmentier/Google Drive/Space_beats_time/Data/data_AZ_jacob/mask_qc_h08v05/MOD11A2_A2002001_h08v05_006_LST_Day_1km_arizona_10092017.rst /home/bparmentier/Google Drive/Space_beats_time/Data/data_AZ_jacob/test_projection_az.rst"
