########################################  MODIS QC FOR MOD12Q1 and MOD11A2 #######################################
########################################### Read QC flags in R for MODIS #####################################
#This script provides an example of Quality Flag processing for twor MODIS product in R. 
#MODIS currently stores information in HDF4 format. Layers to be extracted must be listed first
#using for example gdalinfo. Note that QC flags are store bitpacks of 8bits (byte) in big endian!!!
#A data frame matching flag values is created to facilate the processing.            
#Much of the inspiration and code originates from Steve Mosher:
#http://stevemosher.wordpress.com/2012/12/05/modis-qc-bits/
#AUTHOR: Benoit Parmentier                                                                       
#DATE: 09/16/2013                                                                                
#PROJECT: Space Time project and NCEAS                                 
###################################################################################################

###Loading R library and packages                                                      
#library(gtools)    # loading some useful tools 

library(sp)
library(raster)
library(rgdal)
library(BMS) #contains hex2bin and bin2hex
library(bitops)

### Parameters and arguments

in_dir<- "/Users/benoitparmentier/Dropbox/Data/NCEAS/MODIS_processing"
out_dir<- "/Users/benoitparmentier/Dropbox/Data/NCEAS/MODIS_processing"
setwd(out_dir)

infile_LC<-"MCD12Q1A2001001.h10v09.051.2012157220925.hdf"
infile_LST <- list.files(path=in_dir,pattern=".*.hdf$")
infile_IDRISI <- list.files(path=in_dir,pattern=".*.rst$")  

function_analyses_paper <-"MODIS_and_raster_processing_functions_09162013.R"
script_path<-in_dir #path to script functions
source(file.path(script_path,function_analyses_paper)) #source all functions used in this script.

### BEGIN ###

rg <-GDAL.open(infile_LST[1])
getDriver(rg)
getDriverLongName(rg) #this does not work
#GDALinfo(rg)
#test<-readGDAL(rg,band=1)
GDALinfo_hdf<-GDALinfo(infile_LST[1])
str(GDALinfo_hdf)
#modis_subdataset <-attr(GDALinfo_hdf,"subdsmdata") #get modis subdataset
#GDALinfo_hdf["columns"]
#GDALinfo_hdf["rows"]

modis_subdataset <- attributes(GDALinfo_hdf)$subdsmdata
print(modis_subdataset)

###### PART I: EXPLORATION OF MODIS HDF PRODUCT: IMPORT, SCREENING OF QC FLAGS

## First: Reading LST and Land cover layers (subset) ######

#modis_subset_layer_LST_Day <- paste("HDF4_EOS:EOS_GRID:",hdf_file,":MODIS_Grid_Daily_1km_LST:LST_Day_1km",sep="")
#HDF4_EOS:EOS_GRID:"MOD11A1.A2001001.h09v04.005.2006343034412.hdf":MODIS_Grid_Daily_1km_LST:LST_Day_1km
#HDF4_EOS:EOS_GRID:"MOD11A1.A2001001.h09v04.005.2006343034412.hdf":MODIS_Grid_Daily_1km_LST:QC_Day
#modis_subset_layer <- paste("HDF4_EOS:EOS_GRID:",f20,":MODIS_Grid_Daily_1km_LST:LST_Day_1km",sep='')
modis_subset_layer_LST_Day <- paste("HDF4_EOS:EOS_GRID:",infile_LST[1],":MODIS_Grid_Daily_1km_LST:LST_Day_1km",sep="")
modis_subset_layer_LST_QC <- paste("HDF4_EOS:EOS_GRID:",infile_LST[1],":MODIS_Grid_Daily_1km_LST:QC_Day",sep="")

#modis_subset_layer <- file.path(in_dir,modis_subset_layer)
r  <- readGDAL(modis_subset_layer_LST_Day)
r  <- raster(r)

r_qc <-readGDAL(modis_subset_layer_LST_QC)
r_qc  <-raster(r_qc)

#system("gdalinfo MOD11A1.A2001020.h09v04.005.2006347194212.hdf")
#system(paste("gdalinfo"," MCD12Q1A2001001.h10v09.051.2012157220925.hdf",sep=""))
quartz(13,28)
plot(r)
plot(r_qc)
#note the storage in FLT4S, i.e. float but in effect empty for 
#for anything greater than 255!!!
freq(r_qc)

#0 is good quality
#Here are the most frequent categories found in the QC lqyers...

r_qc2 <- r_qc==2      # QC not produced clouds
r_qc0 <- r_qc==0      # QC good quality
r_qc3 <- r_qc==3      # QC not produced
r_qc65 <- r_qc==65    # LST produced

r_qcl1 <- stack(r_qc0,r_qc2,r_qc3,r_qc65)
layerNames(r_qcl1) <- c("r_qc0 Good quality","r_qc2 Not produced-Cloud","r_qc3 not produced","r_qc65 LST produced")
plot(r_qcl1)

### Second : Reading and Handling QC flags ######

## CONVERSION TO RAW FORMAT: rawToBits,inToBits,packBits: this is in {base} package
intToBits(65) #integer INT four bytes, not that the notation include 01 for 1
#rawToBits(65) for vector
length(intToBits(65))
intToBits(65)[1:8] #integer INT four bytes

as.integer(intToBits(65)[1:8])
#[1] 1 0 0 0 0 0 1 0  #this is little endian binary notation 

fg<-as.integer(intToBits(65)[1:8]) #flag reversed into big endian to match MODIS
fg[8:1] #This is number 65 in big endian format!!!    BITS 0-1 : LST produced check other QA
#[1] 0 1 0 0 0 0 0 1

rev(as.integer(intToBits(65)[1:8]))

##### Quick test

unique_val<-unique(r_qc) #unique values

f_values <- as.data.frame(freq(r_qc)) # frequency values in the raster...
head(f_values)
f_values
rev(as.integer(intToBits(65)[1:8]))
val<-(as.integer(intToBits(65)[1:8]))
val<-val[8:1]

#sprintf("%x",65) #convert decimal to hexadecimal using C-style string formatting
#sprintf("%x",123)

#r_qc2 <- r_qc==2      # LST not produced due to clouds (1-0)
#r_qc0 <- r_qc==0      # LST good quality (0-0)
#r_qc3 <- r_qc==3      # LST not produced (1-1) (in hex decimal 0x03)
#r_qc65 <- r_qc==65    # LST produced check QA (0-1)

#[1] 0 1 0 0 0 0 0 1 : this is 65
val[1:2] # LST produced check QA (0-1)
val[3:4] # good quality (0-0) #this is in the second level info, data quality flag
val[5:6] # average emissivity error <= 0.01 , Emis Error Flag
val[7:8] # average LST error <= 2K , LST Error Flag

rev(as.integer(intToBits(65)[1:8]))

rev(hex2bin("0x03")) #if 110000000000 then it is not produced and should be removed...?
#if((qc_this_day & 0x03)==0, ${lst}, null())'   # & is bitwise-and, 
#e.g qc_this day: 00000010
#e.g.    0x03   : 00000011
# result "AND"  : 00000010 hence FALSE and value is set to null in GRASS...
# Check in GRASS...

# In R use bitAnd(a,b) for the bitwise operator &

#there are 13 unique values, the most frequent one is value 2 with 3,0 and 65 following...

## REWRITE INTO A FUNCTION with options for LST,LC and NDVI/EVI

## Step 1: list values in raster

f_values <- as.data.frame(freq(r_qc)) # frequency values in the raster...
head(f_values)

## Step 2: convert integer values into relevant binary
t44 <- (sapply(f_values$value,function(x){rev(as.integer(intToBits(x)[1:8]))}))
t44 <- (lapply(f_values$value,function(x){rev(as.integer(intToBits(x)[1:8]))}))

#f_values$bin_val <- unlist(lapply(f_values$value,function(x){rev(as.integer(intToBits(x)[1:8]))}))
#This is currently created for LST (see S. Mosher blog)

########### PART II: PROCESSING MODIS HDF FILES IN WORKFLOW ##########
########START  HERE EXAMPLE OF USE OF FUNCTIONS FOR MODIS AND RASTER PROCESSING 

#Step 1: download modis tiles for specified product
#Step 2: import modis tiles for specified product
#Step 3: deal with modis flags
#Step 4: mosaic tiles
#Step 5: reproject and crop extent to study region

#### Set parameters and arguements

in_dir<- "/Users/benoitparmentier/Dropbox/Data/NCEAS/MODIS_processing"
out_dir<- "/Users/benoitparmentier/Dropbox/Data/NCEAS/MODIS_processing"
setwd(out_dir)

#infile_reg_outline=""  #input region outline defined by polygon: none for Venezuela
#This is the shape file of outline of the study area                                                      #It is an input/output of the covariate script
infile_reg_outline <- "/Users/benoitparmentier/Dropbox/Data/NCEAS/MODIS_processing/OR83M_state_outline.shp"  #input region outline defined by polygon: Oregon
#ref_rast_name<-""  #local raster name defining resolution, exent, local projection--. set on the fly?? 
ref_rast_name<-"/Users/benoitparmentier/Dropbox/Data/NCEAS/MODIS_processing/mean_day244_rescaled.rst"  #local raster name defining resolution, exent: oregon
out_suffix <-"_09202013"
#MODIS_product <- "MOD11A1.005"
#start_date <- "2001.01.01"
#end_date <- "2001.01.05"
#list_tiles<- c("h08v04","h09v04")
#out_dir<- "/Users/benoitparmentier/Dropbox/Data/NCEAS/MODIS_processing"
#file_format <- "hdf"
#product_version <- "5"
#NA_flag_val <- -9999
#import_file_format <- ".rst"
#modis_layer_str1 <- unlist(strsplit(modis_subdataset[1],"\""))[3] #Get day LST layer

##########################################
#### STEP 1:  DOWNLOAD MODIS PRODUCT  ####

MODIS_product <- "MOD11A1.005"
start_date <- "2001.01.01"
end_date <- "2001.01.05"
list_tiles<- c("h08v04","h09v04")
out_dir<- "/Users/benoitparmentier/Dropbox/Data/NCEAS/MODIS_processing"
file_format <- "hdf"
product_version <- "5"
#out_suffix <-
#Note that files are downloaded in the ouput directory in subdirectory with tile_name (e.g. h08v04)
#debug(modis_product_download)

download_modis_obj <- modis_product_download(MODIS_product,version,start_date,end_date,list_tiles,file_format,out_dir)

####################################
##### STEP 2: IMPORT MODIS LAYERS ###

#infile_LST <- list_files_tiles #use only hdf!!!
names(download_modis_obj)
infile_LST <- download_modis_obj$list_files_by_tiles[,1]
GDALinfo_hdf<-GDALinfo(infile_LST[1])
str(GDALinfo_hdf)
modis_subdataset <- attributes(GDALinfo_hdf)$subdsmdata
print(modis_subdataset)

modis_layer_str1 <- unlist(strsplit(modis_subdataset[1],"\""))[3] #Get day LST layer
modis_layer_str2 <- unlist(strsplit(modis_subdataset[3],"\""))[3] #Get day QC layer

#import one file for exploration and checking values
r_LST <-import_modis_layer_fun(infile_LST[1],subdataset=modis_layer_str1,NA_flag=-9999,out_rast_name="test.rst",memory=T)
r_qc <-import_modis_layer_fun(infile_LST[1],subdataset=modis_layer_str2,NA_flag=-9999,out_rast_name="test_qc.rst",memory=T)
plot(r_LST)

### import list of files before mosaicing

file_format_import <- ".rst"
NA_flag_val=-9999
out_suffix_s <- "day_LST"
list_param_import_modis <- list(i=1,hdf_file=infile_LST,subdataset=modis_layer_str1,NA_flag_val=NA_flag_val,out_dir=out_dir,
                                out_suffix=out_suffix_s,file_format=file_format_import)
#debug(import_list_modis_layers_fun)
import_list_modis_layers_fun(1,list_param_import_modis)
r_LST_s <- lapply(1:length(infile_LST),FUN=import_list_modis_layers_fun,list_param=list_param_import_modis)
out_suffix_s <- "day_qc"
list_param_import_modis <- list(i=1,hdf_file=infile_LST,subdataset=modis_layer_str2,NA_flag_val=NA_flag_val,out_dir=out_dir,
                                out_suffix=out_suffix_s,file_format=file_format_import)
r_LSTqc_s <- lapply(1:length(infile_LST),FUN=import_list_modis_layers_fun,list_param=list_param_import_modis)

r_LST_s <- unlist(r_LST_s) #list of files as character vector
r_LSTqc_s <- unlist(r_LSTqc_s) #list of files as character vector
r_LST_s
r_LSTqc_s

#################################
##### STEP 3: APPLY/DEAL WITH QC FLAG  ###

## Get QC information for lST and mask values
QC_obj <- create_MODIS_QC_table(LST=TRUE, NDVI=FALSE) #Get table corresponding to QC for LST
names(QC_obj)
QC_data_lst <- QC_obj$LST
QC_data_ndvi <- QC_obj$NDVI

#qc_lst_valid <- subset(x=QC_data_lst,Bit1 == 0 & Bit0 ==1 & Bit3 !=1) #only keep Good Quality (QA_word1)
### Now that QC table has been generated, screen table for desired values
#Select level 1:
qc_lst_valid <- subset(x=QC_data_lst,QA_word1 == "LST Good Quality" | QA_word1 =="LST Produced,Check QA")
#Select level 2:
qc_lst_valid <- subset(x=qc_lst_valid,QA_word2 == "Good Data" | QA_word2 =="Other Quality")
#Select level 3:
#...

names(qc_lst_valid)
qc_valid<-qc_lst_valid$Integer_Value #value integer values
NA_flag_val <- -9999
out_rast_name <- "test.tif"
#debug(qc_valid_modis_fun)
#Quick look at one image
#r_stack <- qc_valid_modis_fun(qc_valid,rast_qc=r_qc,rast_var=r_LST,rast_mask=TRUE,NA_flag_val,out_dir=".",out_rast_name)

list_param_screen_qc <- list(qc_valid,r_LSTqc_s, r_LST_s,rast_mask=TRUE,NA_flag_val,out_dir,out_suffix) 
names(list_param_screen_qc) <- c("qc_valid","rast_qc", "rast_var","rast_mask","NA_flag_val","out_dir","out_suffix") 

#debug(screen_for_qc_valid_fun)
r_stack <- screen_for_qc_valid_fun(1,list_param=list_param_screen_qc)
  
r_stack <- lapply(1:length(r_LST_s),FUN=screen_for_qc_valid_fun,list_param=list_param_screen_qc)

### Now loop over a series of files...

#################################
##### STEP 4: MOSAIC TILES  ###


#################################
##### STEP 5: REPROJECT AND CROP TO STUDY REGION  ###


