########################################  MODIS QC FOR MOD12Q1 and MOD11A2 #######################################
########################################### Read QC flags in R for MODIS #####################################
#This script provides an example of Quality Flag processing for twor MODIS product in R. 
#MODIS currently stores information in HDF4 format. Layers to be extracted must be listed first
#using for example gdalinfo. Note that QC flags are store bitpacks of 8bits (byte) in big endian!!!
#A data frame matching flag values is created to facilate the processing.            
#Much of the inspiration and code originates from Steve Mosher:
#http://stevemosher.wordpress.com/2012/12/05/modis-qc-bits/
#AUTHOR: Benoit Parmentier                                                                       
#DATE: 09/11/2013                                                                                
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

function_analyses_paper <-"MODIS_and_raster_processing_functions_09112013.R"
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

###### PART I: Reading LST and Land cover layers (subset) ######

#HDF4_EOS:EOS_GRID:"MOD11A1.A2001001.h09v04.005.2006343034412.hdf":MODIS_Grid_Daily_1km_LST:LST_Day_1km
#HDF4_EOS:EOS_GRID:"MOD11A1.A2001001.h09v04.005.2006343034412.hdf":MODIS_Grid_Daily_1km_LST:QC_Day
#modis_subset_layer <- paste("HDF4_EOS:EOS_GRID:",f20,":MODIS_Grid_Daily_1km_LST:LST_Day_1km",sep='')
modis_subset_layer_LST_Day <- paste("HDF4_EOS:EOS_GRID:",infile_LST[1],":MODIS_Grid_Daily_1km_LST:LST_Day_1km",sep="")
modis_subset_layer_LST_QC <- paste("HDF4_EOS:EOS_GRID:",infile_LST[1],":MODIS_Grid_Daily_1km_LST:QC_Day",sep="")

#modis_subset_layer <- file.path(in_dir,modis_subset_layer)
r <-readGDAL(modis_subset_layer_LST_Day)
r  <-raster(r)

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

###### PART II: Reading and Handling QC flags ######

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

QC_obj <- create_MODIS_QC_table(LST=TRUE)
names(QC_obj)
QC_data_lst <- QC_obj$LST
qc_lst_valid <- subset(x=QC_data_lst,Bit1 == 0 & Bit0 ==1 & Bit3 !=1)
names(qc_lst_valid)
qc_lst_valid$Integer_Value

qc_valid<-qc_lst_valid$Integer_Value

qc_valid_modis_fun(qc_valid,rast_qc,rast_var){
  
## function to download modis product??

#add here...





