########################################  MODIS QC FOR MOD12Q1 and MOD11A2 #######################################
########################################### Read QC flags in R for MODIS #####################################
#This script provides an example of Quality Flag processing for twor MODIS product in R. 
#MODIS currently stores information in HDF4 format. Layers to be extracted must be listed first
#using for example gdalinfo. Note that QC flags are store bitpacks of 8bits (byte) in big endian!!!
#A data frame matching flag values is created to facilate the processing.            
#Much of the inspiration and code originates from Steve Mosher:
#http://stevemosher.wordpress.com/2012/12/05/modis-qc-bits/
#AUTHOR: Benoit Parmentier                                                                       
#DATE: 08/29/2013                                                                                
#PROJECT: Space Time project and NCEAS                                 
###################################################################################################

###Loading R library and packages                                                      
#library(gtools)    # loading some useful tools 

library(sp)
library(raster)
library(rgdal)

### Parameters and arguments

in_dir<- "/Users/benoitparmentier/Dropbox/Data/NCEAS/MODIS_processing"
out_dir<- "/Users/benoitparmentier/Dropbox/Data/NCEAS/MODIS_processing"
setwd(out_dir)

infile_LC<-"MCD12Q1A2001001.h10v09.051.2012157220925.hdf"
infile_LST <- list.files(path=in_dir,pattern=".*.hdf$")
  
### BEGIN ###

rg <-GDAL.open(infile_LST[1])
getDriver(rg)
getDriverLongName(rg) #this does not work
GDALinfo(rg)

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

r_qc2 <- r_qc==2      # QC clouds
r_qc0 <- r_qc==0      # QC not produced
r_qc3 <- r_qc==3      
r_qc65 <- r_qc==65

r_qcl1 <- stack(r_qc0,r_qc2,r_qc3,r_qc65)
layerNames(r_qcl1) <- c("r_qc0","r_qc2","r_qc3","r_qc65")
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
fg[8:1] #This is number 65 in big endian format!!!

rc2<-r_qc
rc2[rc2==112]<-NA #water?
plot(rc2,colNA=c("red"))

rev(as.integer(intToBits(112)[1:8]))

##### Quick test

unique_val<-unique(rc) #unique values
f_values <-freq(rc) # frequency values in the raster...
head(f_values)
plot(rc)

rev(as.integer(intToBits(16)[1:8]))
val<-(as.integer(intToBits(16)[1:8]))
val<-val[8:1]

val[1:2] # good quality because 00
val[3:4] # quarter 1 beccause 01
val[5:8] # 0000 nothing else that land should be incorrect...need to check documentation...

rev(as.integer(intToBits(112)[1:8]))
#01110000

#0:shallow water
#1:

#there are 32 unique values, the most frequent one is value 16 with 3,295,253

## REWRITE INTO A FUNCTION with options for LST,LC and NDVI/EVI

#This is currently created for LST (see S. Mosher blog)
QC_Data <- data.frame(Integer_Value = 0:255,
                      Bit7 = NA,
                      Bit6 = NA,
                      Bit5 = NA,
                      Bit4 = NA,
                      Bit3 = NA,
                      Bit2 = NA,
                      Bit1 = NA,
                      Bit0 = NA,
                      QA_word1 = NA,
                      QA_word2 = NA,
                      QA_word3 = NA,
                      QA_word4 = NA)

for(i in QC_Data$Integer_Value){
  AsInt <- as.integer(intToBits(i)[1:8])
  QC_Data[i+1,2:9]<- AsInt[8:1]
}

QC_Data$QA_word1[QC_Data$Bit1 == 0 & QC_Data$Bit0==0] <- "LST GOOD"
QC_Data$QA_word1[QC_Data$Bit1 == 0 & QC_Data$Bit0==1] <- "LST Produced,Other Quality"
QC_Data$QA_word1[QC_Data$Bit1 == 1 & QC_Data$Bit0==0] <- "No Pixel,clouds"
QC_Data$QA_word1[QC_Data$Bit1 == 1 & QC_Data$Bit0==1] <- "No Pixel, Other QA"

QC_Data$QA_word2[QC_Data$Bit3 == 0 & QC_Data$Bit2==0] <- "Good Data"
QC_Data$QA_word2[QC_Data$Bit3 == 0 & QC_Data$Bit2==1] <- "Other Quality"
QC_Data$QA_word2[QC_Data$Bit3 == 1 & QC_Data$Bit2==0] <- "TBD"
QC_Data$QA_word2[QC_Data$Bit3 == 1 & QC_Data$Bit2==1] <- "TBD"

QC_Data$QA_word3[QC_Data$Bit5 == 0 & QC_Data$Bit4==0] <- "Emiss Error <= .01"
QC_Data$QA_word3[QC_Data$Bit5 == 0 & QC_Data$Bit4==1] <- "Emiss Err >.01 <=.02"
QC_Data$QA_word3[QC_Data$Bit5 == 1 & QC_Data$Bit4==0] <- "Emiss Err >.02 <=.04"
QC_Data$QA_word3[QC_Data$Bit5 == 1 & QC_Data$Bit4==1] <- "Emiss Err > .04"

QC_Data$QA_word4[QC_Data$Bit7 == 0 & QC_Data$Bit6==0] <- "LST Err <= 1"
QC_Data$QA_word4[QC_Data$Bit7 == 0 & QC_Data$Bit6==1] <- "LST Err > 2 LST Err <= 3"
QC_Data$QA_word4[QC_Data$Bit7 == 1 & QC_Data$Bit6==0] <- "LST Err > 1 LST Err <= 2"
QC_Data$QA_word4[QC_Data$Bit7 == 1 & QC_Data$Bit6==1] <- "LST Err > 4"

#Screen data
FINAL <- QC_Data[QC_Data$Bit1 == 0 & QC_Data$Bit0 ==1 & QC_Data$Bit3 !=1,]
