########################################  MODIS TILES PROCESSING #######################################
########################################### Read QC flags and mosaic tiles in R #####################################
#The current version generat data forthe Ecuador region.
#This script processes MODIS tiles using Quality Flag. Tiles are mosaiced and reprojected for a specific study region.
#MODIS currently stores information in HDF4 format. Layers must be extracted and must be listed first
#using for example gdalinfo to identify the relevant dataset and QC flag options. 
#Note that QC flags are store bitpacks of 8bits (byte) in big endian!!!
#A data frame matching flag values is created to facilate the processing.            
#Inspiration and some code for the MODIS flag function originates from Steve Mosher:
#http://stevemosher.wordpress.com/2012/12/05/modis-qc-bits/
## MODIS WORKFLOW
# Processing of MODIS HDF files is done in 5 steps:c
# Step 1: download modis tiles for specified product and version (e.g. version 5)
# Step 2: import modis tiles for specified file format (e.g. ".tif",".rst)
# Step 3: deal with modis flags (multiple levels)
# Step 4: mosaic tiles for every time step
# Step 5: reproject and crop extent to study region
#
#AUTHOR: Benoit Parmentier                                                                       
#CREATED ON : 09/16/2013  
#MODIFIED ON : 01/28/2013
#PROJECT: NCEAS and general MODIS processing of all projects
#TODO: 
#1)Test additional Quality Flag levels for ALBEDO and other product
#2)Add function to report statistics: missing files
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

#############################
### Parameters and arguments

#in_dir<- "/Users/Parmentier/Documents/Benoit/Space_Time_Yucatan"
#out_dir<- "/Users/Parmentier/Documents/Benoit/Space_Time_Yucatan"
#in_dir<- "/Volumes/Data/Ecuador_Project" #path to parent directory where output are/will be based
#out_dir<- "/Volumes/Data/Ecuador_Project"
in_dir <- "/Volumes/Seagate Backup Plus Drive/Ecuador_Project/"
out_dir <- "/Volumes/Seagate Backup Plus Drive/Ecuador_Project/"

setwd(out_dir)

function_analyses_paper <-"MODIS_and_raster_processing_functions_01282014.R"
script_path <- "/Users/Parmentier/Dropbox/Data/NCEAS/git_space_time_lucc/scripts_queue" #path to script functions
source(file.path(script_path,function_analyses_paper)) #source all functions used in this script.

#infile_reg_outline=""  #input region outline defined by polygon: none for Venezuela
#This is the shape file of outline of the study area                                                      #It is an input/output of the covariate script
#infile_reg_outline <- "/Users/Parmentier/Documents/Benoit/Space_Time_Yucatan/GYRS_Municipios/GYRS_MX_tri-state_latlong.shp"  #input region outline defined by polygon: Oregon
infile_reg_outline <- "/Volumes/Seagate Backup Plus Drive/Ecuador_Project/region_outlines_ref_files/OR83M_state_outline.shp"  #input region outline defined by polygon: Oregon

#ref_rast_name<-""  #local raster name defining resolution, exent, local projection--. set on the fly?? 
#ref_rast_name<-"/Users/Parmentier/Documents/Benoit/Space_Time_Yucatan/GYRS_Municipios/GYRS_latlong_mask_1km.rst"  #local raster name defining resolution, exent: oregon
ref_rast_name<-"/Volumes/Seagate Backup Plus Drive/Ecuador_Project/region_outlines_ref_files/ref_reg_ecuador_250m.rst"  #local raster name defining resolution, exent: oregon

#infile_modis_grid<-"/Users/Parmentier/Documents/Benoit/Space_Time_Yucatan/modis_sinusoidal_grid_world.shp" #modis grid tiling system, global
infile_modis_grid <- "/Volumes/Seagate Backup Plus Drive/Ecuador_Project/region_outlines_ref_files/modis_sinusoidal_grid_world.shp" #modis grid tiling system, global

out_suffix <- "12282013" #output suffix for the files that are masked for quality and for 

## Other specific parameters

MODIS_product <- "MOD13Q1.005" #NDVI/EVI 250m product (monthly)
#MODIS_product <- "MOD11A1.005"
#start_date <- "2001.01.01"
#end_date <- "2001.03.05"
start_date <- "2001.01.01"
end_date <- "2012.12.31"

list_tiles_modis<- c("h09v08,h09v09,h10v09,h10v08") 
#list_tiles_modis<- NULL

file_format_download <- "hdf"
file_format <- ".rst" #output format
NA_flag_val <- -9999 #Flag used for missing values for values out of range
product_version <- 5
#temporal_granularity <- "Daily" #deal with options( 16 day, 8 day and monthly)
temporal_granularity <- "16 Day" #deal with options( 16 day, 8 day and monthly), unused at this stage...

#scaling_factors <- c(1,-273.15) #set up as slope (a) and intercept (b), if NULL, no scaling done, setting for LST 
scaling_factors <- c(0.0001,0) #set up as slope (a) and intercept (b), if NULL, no scaling done, setting for NDVI 
#scaling_factors <- NULL #set up as slope (a) and intercept (b), if NULL, no scaling done 
#modis_layer_str1 <- unlist(strsplit(modis_subdataset[1],"\""))[3] #Get day LST layer
#modis_layer_str2 <- unlist(strsplit(modis_subdataset[1],"\""))[3] #Get qc LST layer

proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
#CRS_interp <-"+proj=lcc +lat_1=43 +lat_2=45.5 +lat_0=41.75 +lon_0=-120.5 +x_0=400000 +y_0=0 +ellps=GRS80 +units=m +no_defs";
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_interp <- proj_modis_str

#run all processing steps
#steps_to_run <- list(download=TRUE,import=TRUE,apply_QC_flag=TRUE,moscaic=TRUE,reproject=TRUE)
#run specific processing steps
#steps_to_run <- list(download=FALSE,import=TRUE,apply_QC_flag=TRUE,moscaic=TRUE,reproject=TRUE)
steps_to_run <- list(download=FALSE,import=FALSE,apply_QC_flag=FALSE,moscaic=TRUE,reproject=TRUE)

######################################################
########################  BEGIN SCRIPT  #############

#####################################
#### STEP 1:  DOWNLOAD MODIS PRODUCT  ####

#Note that files are downloaded in the ouput directory in subdirectory with tile_name (e.g. h08v04)

#if(list_tiles_modis==NULL){
#  d
#}

list_tiles_modis <- unlist(strsplit(list_tiles_modis,","))  # transform string into separate element in char vector

#debug(modis_product_download)
if(steps_to_run$download==TRUE){
  download_modis_obj <- modis_product_download(MODIS_product,version,start_date,end_date,list_tiles_modis,file_format_download,out_dir,temporal_granularity)
  out_dir_tiles <- (file.path(out_dir,list_tiles_modis))
  list_files_by_tiles <-download_modis_obj$list_files_by_tiles #Use mapply to pass multiple arguments
  colnames(list_files_by_tiles) <- list_tiles_modis #note that the output of mapply is a matrix
}else{
  out_dir_tiles <- (file.path(out_dir,list_tiles_modis))
  list_files_by_tiles <-mapply(1:length(out_dir_tiles),FUN=list.files,pattern="*.hdf$",path=out_dir_tiles,full.names=T) #Use mapply to pass multiple arguments
  colnames(list_files_by_tiles) <- list_tiles_modis #note that the output of mapply is a matrix
}

####################################
##### STEP 2: IMPORT MODIS LAYERS ###

##Modify this section into a function to extract names of product and quality flag automatically!!

#infile_LST <- list_files_tiles #use only hdf!!!
#names(download_modis_obj)
infile_var <- list_files_by_tiles[,1]
#GDALinfo_hdf<-GDALinfo(infile_var[1],returnScaleOffset=FALSE)
GDALinfo_hdf<-GDALinfo(infile_var[1])
str(GDALinfo_hdf)
modis_subdataset <- attributes(GDALinfo_hdf)$subdsmdata
print(modis_subdataset)
#modis_subdataset_str1<-"SUBDATASET_1_NAME=HDF4_EOS:EOS_GRID:\"MOD11A1.A2010365.h09v04.005.2011034203707.hdf\":MODIS_Grid_Daily_1km_LST:LST_Day_1km"
#modis_subdataset_str2<-"SUBDATASET_2_NAME=HDF4_EOS:EOS_GRID:\"MOD11A1.A2010365.h09v04.005.2011034203707.hdf\":MODIS_Grid_Daily_1km_LST:QC_Day"
#modis_subdataset_str1 <- "SUBDATASET_1_NAME=HDF4_EOS:EOS_GRID:\"MOD13Q1.A2012353.h10v08.005.2013009145323.hdf\":MODIS_Grid_16DAY_250m_500m_VI:250m 16 days NDVI"
#modis_subdataset_str2<- "SUBDATASET_3_NAME=HDF4_EOS:EOS_GRID:\"MOD13Q1.A2012353.h10v08.005.2013009145323.hdf\":MODIS_Grid_16DAY_250m_500m_VI:250m 16 days VI Quality"
#modis_layer_str1 <- unlist(strsplit(modis_subdataset_str1,"\""))[3] #Get day LST layer
#modis_layer_str2 <- unlist(strsplit(modis_subdataset_str2,"\""))[3] #Get day QC layer

modis_layer_str1 <- unlist(strsplit(modis_subdataset[1],"\""))[3] #Get day LST layer
modis_layer_str2 <- unlist(strsplit(modis_subdataset[5],"\""))[3] #Get day QC layer

#modis_layer_str1 <- unlist(strsplit(modis_subdataset_str1,"\""))[3] #Get day LST layer
#modis_layer_str2 <- unlist(strsplit(modis_subdataset_str2,"\""))[3] #Get day QC layer

### import list of files before mosaicing

file_format_import <- file_format
#NA_flag_val=-9999
var_modis_name <- unlist(strsplit(modis_layer_str1,":"))[3]
qc_modis_name <- unlist(strsplit(modis_layer_str2,":"))[3]

var_modis_name <- gsub(" ","_",var_modis_name) #suffix name for product, may contain white space so replace with "_"
qc_modis_name <- gsub(" ","_",qc_modis_name)

##loop over tiles:

if(steps_to_run$import==TRUE){
  for(j in 1:length(list_tiles_modis)){
    #infile_var <- download_modis_obj$list_files_by_tiles[,j] 
    infile_var <-list_files_by_tiles[,j] #note can be any variable even thought LST presented  here
    out_dir_s <- file.path(out_dir,list_tiles_modis[j])
    out_suffix_s <- var_modis_name
    list_param_import_modis <- list(i=1,hdf_file=infile_var,subdataset=modis_layer_str1,NA_flag_val=NA_flag_val,out_dir=out_dir_s,
                                    out_suffix=out_suffix_s,file_format=file_format_import,scaling_factors=scaling_factors)
    #undebug(import_list_modis_layers_fun)
    #r1<-import_list_modis_layers_fun(1,list_param_import_modis)
    r_var_s <- lapply(1:length(infile_var),FUN=import_list_modis_layers_fun,list_param=list_param_import_modis)
    #r_var_s <- mclapply(1:11,FUN=import_list_modis_layers_fun,list_param=list_param_import_modis,mc.preschedule=FALSE,mc.cores = 11) #This is the end bracket from mclapply(...) statement
    
    #r_var_s <- mclapply(1:length(infile_var),FUN=import_list_modis_layers_fun,list_param=list_param_import_modis,mc.preschedule=FALSE,mc.cores = 1) #This is the end bracket from mclapply(...) statement
    
    out_suffix_s <- qc_modis_name
    list_param_import_modis <- list(i=1,hdf_file=infile_var,subdataset=modis_layer_str2,NA_flag_val=NA_flag_val,out_dir=out_dir_s,
                                    out_suffix=out_suffix_s,file_format=file_format_import,scaling_factors=NULL)
    #r1<-import_list_modis_layers_fun(1,list_param_import_modis)
    r_qc_s <- lapply(1:length(infile_var),FUN=import_list_modis_layers_fun,list_param=list_param_import_modis)  
    #r_qc_s <-mclapply(1:11,FUN=import_list_modis_layers_fun,list_param=list_param_import_modis,mc.preschedule=FALSE,mc.cores = 11) #This is the end bracket from mclapply(...) statement    
    #r_qc_s <-mclapply(1:length(infile_var),FUN=import_list_modis_layers_fun,list_param=list_param_import_modis,mc.preschedule=FALSE,mc.cores = 1) #This is the end bracket from mclapply(...) statement
    
  }
}

r_var_s <- unlist(r_var_s) #list of files as character vector
r_qc_s <- unlist(r_qc_s) #list of files as character vector
plot(raster(r_qc_s[1]))
plot(raster(r_var_s[2]))
print(r_var_s)
print(r_qc_s)

#################################
##### STEP 3: APPLY/DEAL WITH QC FLAG AND SCREEN VALUE FOR VALID RANGE ###

## Get QC information for lST and mask values
QC_obj <- create_MODIS_QC_table(LST=TRUE, NDVI=TRUE) #Get table corresponding to QC for LST
names(QC_obj)
#QC_data_lst <- QC_obj$LST
QC_data_ndvi <- QC_obj$NDVI

#For LST
#Select level 1:
#qc_lst_valid <- subset(x=QC_data_lst,QA_word1 == "LST Good Quality" | QA_word1 =="LST Produced,Check QA")
#Select level 2:
#qc_lst_valid <- subset(x=qc_lst_valid,QA_word2 == "Good Data" | QA_word2 =="Other Quality")
#Select level 3:
#...

#For NDVI: use this section to process
#Select level 1:
qc_lst_valid <- subset(x=QC_data_ndvi,QA_word1 == "VI Good Quality" | QA_word1 =="VI Produced,check QA")
#Select level 2:
qc_lst_valid <- subset(x=qc_lst_valid,QA_word2 %in% unique(QC_data_ndvi$QA_word2)[1:8]) #"Highest quality, 1","Lower quality, 2","Decreasing quality, 3",...,"Decreasing quality, 8" 
#Select level 3:
#...
names(qc_lst_valid)
qc_valid<-qc_lst_valid$Integer_Value #value integer values
#NA_flag_val <- -9999

#r_lst_LST <- download_modis_obj$list_files_by_tiles[,j]
list_r_lst <- vector("list",length(list_tiles_modis))
list_r_qc <- vector("list",length(list_tiles_modis))
r_stack <- vector("list",length(list_tiles_modis))

#var_modis_name <- unlist(strsplit(modis_layer_str1,":"))[3]
#qc_modis_name <- unlist(strsplit(modis_layer_str2,":"))[3]

#apply_QC_flag=TRUE,moscaic=TRUE,reproject=TRUE)
if(steps_to_run$apply_QC_flag==TRUE){
  for(j in 1:length(list_tiles_modis)){
    out_dir_s <- file.path(out_dir,list_tiles_modis)[j]
    out_suffix_s <- paste(var_modis_name,sep="") #for MODIS product (var)
    file_format_s <-file_format
    #debug(create_raster_list_from_file_pat)
    list_r_lst[[j]] <- create_raster_list_from_file_pat(out_suffix_s,file_pat="",in_dir=out_dir_s,out_prefix="",file_format=file_format_s)
    
    out_dir_s <- file.path(out_dir,list_tiles_modis)[j]
    out_suffix_s <- paste(qc_modis_name,sep="") #for qc flags
    list_r_qc[[j]] <- create_raster_list_from_file_pat(out_suffix_s,file_pat="",in_dir=out_dir_s,out_prefix="",file_format=file_format_s)
    
    ###
    
    list_param_screen_qc <- list(qc_valid,list_r_qc[[j]], list_r_lst[[j]],rast_mask=TRUE,NA_flag_val,out_dir_s,out_suffix) 
    names(list_param_screen_qc) <- c("qc_valid","rast_qc", "rast_var","rast_mask","NA_flag_val","out_dir","out_suffix") 
    #debug(screen_for_qc_valid_fun)
    #r_stack <- screen_for_qc_valid_fun(1,list_param=list_param_screen_qc)
    r_stack[[j]] <- lapply(1:length(list_r_qc[[j]]),FUN=screen_for_qc_valid_fun,list_param=list_param_screen_qc)
    #r_stack[[j]] <-mclapply(1:11,FUN=screen_for_qc_valid_fun,list_param=list_param_screen_qc,mc.preschedule=FALSE,mc.cores = 11) #This is the end bracket from mclapply(...) statement
    
    #r_stack[[j]] <-mclapply(1:length(list_r_qc[[j]]),FUN=screen_for_qc_valid_fun,list_param=list_param_screen_qc,mc.preschedule=FALSE,mc.cores = 2) #This is the end bracket from mclapply(...) statement
    
  }
}

#r_lst_by_tiles <-mapply(1:length(list_tiles_modis),FUN=list.files,pattern=paste".*.day_LST.*.rst$",path=out_dir_s,full.names=T) #Use mapply to pass multiple arguments
#r_lst <- mapply(1:length(out_suffix_s),FUN=create_raster_list_from_file_pat,
 #               file_pat="",in_dir=out_dir_s,out_prefix="",file_format=file_format_s)

### Now loop over a series of files...
#extract_list_from_list_obj(unlist(r_stack),"var")
  
#################################
##### STEP 4: MOSAIC TILES  ###

#debug(create_raster_list_from_file_pat)

list_m_var <- vector("list",length(list_tiles_modis))  
names(list_m_var)<- list_tiles_modis
list_m_qc <- vector("list",length(list_tiles_modis))  
names(list_m_qc)<- list_tiles_modis

for (j in 1:length(list_tiles_modis)){
  #out_suffix_s <- paste(list_tiles_modis[j],"_",sprintf( "%03d", product_version),"_",var_modis_name,"_",out_suffix,sep="")
  file_format_s <-file_format
  out_suffix_s <- paste(list_tiles_modis[j],"_",sprintf( "%03d", product_version),"_",var_modis_name,"_",out_suffix,file_format_s,sep="")
  out_dir_s <- file.path(out_dir,list_tiles_modis)[j]
  list_m_var[[j]]<-list.files(pattern=paste(out_suffix_s,"$",sep=""),path=out_dir_s,full.names=TRUE)
  #list_m_var[[j]] <-create_raster_list_from_file_pat(out_suffix_s,file_pat="",in_dir=out_dir_s,out_prefix="",file_format=file_format_s)
}

#Prepare list of modis tiles to mosaic
#mosaic_list_var <-mapply(FUN="c",list_m_var,SIMPLIFY=T)
x <-mapply(FUN="c",list_m_var,SIMPLIFY=T)
mosaic_list_var<-lapply(seq_len(nrow(x)), function(i) x[i,]) #list of tiles by batch to mosaic
#Prepare list of output names without extension
out_rastnames_var <- (basename(gsub(list_tiles_modis[1],"",list_m_var[[1]])))
out_rastnames_var <- gsub(extension(out_rastnames_var),"",out_rastnames_var)
j<-1
list_param_mosaic<-list(j,mosaic_list_var,out_rastnames_var,out_dir,file_format,NA_flag_val)
names(list_param_mosaic)<-c("j","mosaic_list","out_rastnames","out_path","file_format","NA_flag_val")
#debug(mosaic_m_raster_list)
#list_var_mosaiced <- mosaic_m_raster_list(1,list_param_mosaic)
#Parallelization,this works on MAC laptop too
#list_var_mosaiced <-mclapply(1:length(mosaic_list_var), list_param=list_param_mosaic, mosaic_m_raster_list,mc.preschedule=FALSE,mc.cores = 2) #This is the end bracket from mclapply(...) statement

list_var_mosaiced <- lapply(1:length(mosaic_list_var), list_param=list_param_mosaic, mosaic_m_raster_list) #This is the end bracket from mclapply(...) statement

#test_rast <- stack(list_var_mosaiced)

#################################
##### STEP 5: REPROJECT AND CROP TO STUDY REGION  ###

# FIRST SET UP STUDY AREA ####

#if (infile_reg_outline!=""){
#  filename<-sub(".shp","",basename(infile_reg_outline))   #Removing path and the extension from file name.
#  reg_outline<-readOGR(dsn=dirname(infile_reg_outline), filename) # Read in the region outline
#}
#if no shapefile defining the study/processing area then create one using modis grid tiles
#if (infile_reg_outline==""){
#  filename<-sub(".shp","",basename(infile_modis_grid))       #Removing path and the extension from file name.
#  modis_grid<-readOGR(dsn=dirname(infile_modis_grid), filename)     #Reading shape file using rgdal library
#  reg_outline_modis <-create_modis_tiles_region(modis_grid,list_tiles_modis) #problem...this does not 
# #align with extent of modis LST!!!
#now add projection on the fly
#  infile_reg_outline <-paste("modis_outline",out_region_name,"_",out_suffix,".shp",sep="")
#  writeOGR(reg_outline_modis,dsn= out_path,layer= sub(".shp","",infile_reg_outline), 
#           driver="ESRI Shapefile",overwrite_layer="TRUE")
#  reg_outline_obj <- define_crs_from_extent_fun(reg_outline_modis,buffer_dist)
#  reg_outline <-reg_outline_obj$reg_outline
#  CRS_interp <-reg_outline_obj$CRS_interp
#  infile_reg_outline <-paste("outline",out_region_name,"_",out_suffix,".shp",sep="")
#  writeOGR(reg_outline,dsn= out_path,layer= sub(".shp","",infile_reg_outline), 
#          driver="ESRI Shapefile",overwrite_layer="TRUE")
#}

# NOW PROJECT AND CROP WIHT REF REGION ####

if (ref_rast_name==""){
  #Use one mosaiced modis tile as reference image...We will need to add a function 
  ref_rast_temp <-raster(list_var_mosaiced[[1]]) 
  ref_rast <-projectRaster(from=ref_rast_temp,crs=CRS_interp,method="ngb")
  #to define a local reference system and reproject later!!
  #Assign new projection system here in the argument CRS_interp (!it is used later)
}else{
  ref_rast<-raster(ref_rast_name) #This is the reference image used to define the study/processing area
  projection(ref_rast) <- CRS_interp #Assign given reference system from master script...
}

##Create output names for region
out_suffix_var <-paste(out_suffix,file_format,sep="")          
var_list_outnames <- change_names_file_list(list_var_mosaiced,out_suffix_var,"reg_",file_format,out_path=out_dir)     

#list_param_create_region<-list(j,raster_name=list_var_mosaiced,reg_ref_rast=ref_rast,out_rast_name=var_list_outnames)
list_param_create_region<-list(j,list_var_mosaiced,ref_rast,var_list_outnames,NA_flag_val)
names(list_param_create_region) <-c("j","raster_name","reg_ref_rast","out_rast_name","NA_flag_val")

#debug(create__m_raster_region)
#test<-create__m_raster_region(1,list_param_create_region)

#reg_var_list <-mclapply(1:length(var_list_outnames), list_param=list_param_create_region, create__m_raster_region,mc.preschedule=FALSE,mc.cores = 1) #This is the end bracket from mclapply(...) statement
reg_var_list <-lapply(1:length(var_list_outnames), list_param=list_param_create_region, create__m_raster_region) 

#Still need to deal with rescaling !!!

test<-stack(reg_var_list[1:4])
plot(test)

########### END OF SCRIPT ##############