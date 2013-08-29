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
library(BMS) #contains hex2bin and bin2hex
library(bitops)

### Parameters and arguments

in_dir<- "/Users/benoitparmentier/Dropbox/Data/NCEAS/MODIS_processing"
out_dir<- "/Users/benoitparmentier/Dropbox/Data/NCEAS/MODIS_processing"
setwd(out_dir)

infile_LC<-"MCD12Q1A2001001.h10v09.051.2012157220925.hdf"
infile_LST <- list.files(path=in_dir,pattern=".*.hdf$")
infile_IDRISI <- list.files(path=in_dir,pattern=".*.rst$")  

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

## function to download modis product??



## Function to mosaic modis

mosaic_m_raster_list<-function(j,list_param){
  #This functions returns a subset of tiles from the modis grid.
  #Arguments: modies grid tile,list of tiles
  #Output: spatial grid data frame of the subset of tiles
  #Note that rasters are assumed to be in the same projection system!!
  
  #rast_list<-vector("list",length(mosaic_list))
  #for (i in 1:length(mosaic_list)){  
  # read the individual rasters into a list of RasterLayer objects
  # this may be changed so that it is not read in the memory!!!
  
  #parse output...
  
  #j<-list_param$j
  mosaic_list<-list_param$mosaic_list
  out_path<-list_param$out_path
  out_names<-list_param$out_rastnames
  ## Start
  
  input.rasters <- lapply(as.character(mosaic_list[[j]]), raster)
  mosaiced_rast<-input.rasters[[1]]
  
  for (k in 2:length(input.rasters)){
    mosaiced_rast<-mosaic(mosaiced_rast,input.rasters[[k]], fun=mean)
    #mosaiced_rast<-mosaic(mosaiced_rast,raster(input.rasters[[k]]), fun=mean)
  }
  
  data_name<-paste("mosaiced_",sep="") #can add more later...
  raster_name<-paste(data_name,out_names[j],".tif", sep="")
  writeRaster(mosaiced_rast, filename=file.path(out_path,raster_name),overwrite=TRUE)  
  #Writing the data in a raster file format...  
  rast_list<-file.path(out_path,raster_name)
  
  return(rast_list)
}

## Function to reproject and crop modis tile

create__m_raster_region <-function(j,list_param){
  #This functions returns a subset of tiles from the modis grdid.
  #Arguments: raster name of the file,reference file with
  #Output: spatial grid data frame of the subset of tiles
  
  ## Parse input arguments
  raster_name <- list_param$raster_name[[j]] #list of raster ot project and crop, this is a list!!
  reg_ref_rast <- list_param$reg_ref_rast #This must have a coordinate system defined!!
  out_rast_name <- list_param$out_rast_name[j]
  
  ## Start #
  layer_rast<-raster(raster_name)
  new_proj<-proj4string(layer_rast)                  #Extract current coordinates reference system in PROJ4 format
  region_temp_projected<-projectExtent(reg_ref_rast,CRS(new_proj))     #Project from ref to current region coord. system
  layer_crop_rast<-crop(layer_rast, region_temp_projected) #crop using the extent from the region tile
  #layer_projected_rast<-projectRaster(from=layer_crop_rast,crs=proj4string(reg_outline),method="ngb")
  layer_projected_rast<-projectRaster(from=layer_crop_rast,to=reg_ref_rast,method="ngb")
  
  writeRaster(layer_projected_rast, filename=out_rast_name,overwrite=TRUE)  
  
  return(out_rast_name)
}

## function to import modis in tif or other format...

import_modis_layer_fun <-function(hdf_file,subdataset,file_format,memory=TRUE){
  #modis_subset_layer_LST_Day <- paste("HDF4_EOS:EOS_GRID:",hdf_file,":MODIS_Grid_Daily_1km_LST:LST_Day_1km",sep="")
  modis_subset_layer_LST_Day <- paste("HDF4_EOS:EOS_GRID:",hdf_file,subdataset,sep="")
  r <-readGDAL(modis_subset_layer_LST_Day)
  r  <-raster(r)
  
  #Finish this part...write out
  #raster_name<- paste("raster","_",row, "r_", col,"c_","r_spat_s","_",out_suffix,file_format, sep="")
  #writeRaster(r_spat, NAflag=NA_flag_val,filename=raster_name,bylayer=TRUE,bandorder="BSQ",overwrite=TRUE)   
  #writeRaster(r_spat, NAflag=NA_flag_val,filename=raster_name,bylayer=TRUE,bandorder="BSQ",overwrite=TRUE)     
  
  if(memory==TRUE}{
    return(r)
  }else{
    return(raster_name)
  }  
}


## Function to  reclass value in 

qc_valid_modis_fun <-function(qc_valid,rast_qc,rast_var){
  f_values <- as.data.frame(freq(rast_qc)) # frequency values in the raster...
  if(f_values$value %in% qc_valid){ #not working...change here...
    f_values$qc_mask <- 1
  }else{
    f_values$qc_mask <- NA
  }
  
  r_qc_m <- subs(x=rast_qc,y=f_values,by=1,which=3)
  rast_lst <-mask(rast_var,r_qc_m)
  return(rast_var)
}

### FUNCTIONS

create_MODIS_QC_table <-function(LST=TRUE){
  #Function to generate MODIS QC  flag table
  #created by Benoit Parmentier with most of the lines from S. Mosher!!
  
  list_QC_Data <- vector("list", length=2)
  names(list_QC_Data) <- c("LST","NDVI")
    
  ## Generate generic product table
  QC_Data <- data.frame(Integer_Value = 0:255,
                        Bit7 = NA,Bit6 = NA,Bit5 = NA,Bit4 = NA,Bit3 = NA,Bit2 = NA,Bit1 = NA,Bit0 = NA,
                        QA_word1 = NA,QA_word2 = NA,QA_word3 = NA,QA_word4 = NA)

  for(i in QC_Data$Integer_Value){
    AsInt <- as.integer(intToBits(i)[1:8])
    QC_Data[i+1,2:9]<- AsInt[8:1]
  } 
  QC_table <- QC_Data
  
  ## PRODUCT 1: LST
  #This can be seen from table defined at LPDAAC: https://lpdaac.usgs.gov/products/modis_products_table/mod11a2
  
  QC_Data$QA_word1[QC_Data$Bit1 == 0 & QC_Data$Bit0==0] <- "LST Good Quality"    #(0-0)
  QC_Data$QA_word1[QC_Data$Bit1 == 0 & QC_Data$Bit0==1] <- "LST Produced,Check QA"
  QC_Data$QA_word1[QC_Data$Bit1 == 1 & QC_Data$Bit0==0] <- "Not Produced,clouds"
  QC_Data$QA_word1[QC_Data$Bit1 == 1 & QC_Data$Bit0==1] <- "No Produced, check Other QA"
  
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
  
  list_QC_Data[[1]] <- QC_Data
  
  ## PRODUCT 2: NDVI
  #This can be seen from table defined at LPDAAC: https://lpdaac.usgs.gov/products/modis_products_table/mod11a2
  
#   QC_Data <- data.frame(Integer_Value = 0:255,
#                         Bit7 = NA,Bit6 = NA,Bit5 = NA,Bit4 = NA,Bit3 = NA,Bit2 = NA,Bit1 = NA,Bit0 = NA,
#                         Bit15 = NA,Bit14 = NA,Bit13 = NA,Bit12 = NA,Bit11 = NA,Bit10 = NA,Bit9 = NA,Bit8 = NA,
#                         QA_word1 = NA,QA_word2 = NA,QA_word3 = NA,QA_word4 = NA,
#                         QA_word8 = NA,QA_word7 = NA,QA_word6 = NA,QA_word5 = NA
#                         QA_word9 = NA)
#   
#   QC_Data <- QC_table
#   
#   QC_Data$QA_word1[QC_Data$Bit1 == 0 & QC_Data$Bit0==0] <- "VI Good Quality"    #(0-0)
#   QC_Data$QA_word1[QC_Data$Bit1 == 0 & QC_Data$Bit0==1] <- "VI Produced,check QA"
#   QC_Data$QA_word1[QC_Data$Bit1 == 1 & QC_Data$Bit0==0] <- "Not Produced,because of clouds"
#   QC_Data$QA_word1[QC_Data$Bit1 == 1 & QC_Data$Bit0==1] <- "Not Produced, other reasons"
#   
#   QC_Data$QA_word2[QC_Data$Bit3 == 0 & QC_Data$Bit2==0] <- "Good Data"
#   QC_Data$QA_word2[QC_Data$Bit3 == 0 & QC_Data$Bit2==1] <- "Other Quality"
#   QC_Data$QA_word2[QC_Data$Bit3 == 1 & QC_Data$Bit2==0] <- "TBD"
#   QC_Data$QA_word2[QC_Data$Bit3 == 1 & QC_Data$Bit2==1] <- "TBD"
#   
#   QC_Data$QA_word3[QC_Data$Bit5 == 0 & QC_Data$Bit4==0] <- "Emiss Error <= .01"
#   QC_Data$QA_word3[QC_Data$Bit5 == 0 & QC_Data$Bit4==1] <- "Emiss Err >.01 <=.02"
#   QC_Data$QA_word3[QC_Data$Bit5 == 1 & QC_Data$Bit4==0] <- "Emiss Err >.02 <=.04"
#   QC_Data$QA_word3[QC_Data$Bit5 == 1 & QC_Data$Bit4==1] <- "Emiss Err > .04"
#   
#   QC_Data$QA_word4[QC_Data$Bit7 == 0 & QC_Data$Bit6==0] <- "LST Err <= 1"
#   QC_Data$QA_word4[QC_Data$Bit7 == 0 & QC_Data$Bit6==1] <- "LST Err > 2 LST Err <= 3"
#   QC_Data$QA_word4[QC_Data$Bit7 == 1 & QC_Data$Bit6==0] <- "LST Err > 1 LST Err <= 2"
#   QC_Data$QA_word4[QC_Data$Bit7 == 1 & QC_Data$Bit6==1] <- "LST Err > 4"
  
#  list_QC_Data[[2]]<- QC_Data
  
  return(list_QC_Data)
}

#Screen data: use only : # level 1: LST Produced good quality, LST Produced other Quality Check QA, 
                         # level 2: good data , Other quality 





