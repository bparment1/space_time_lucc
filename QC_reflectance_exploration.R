#######################################    MODIS processing   #######################################
############################ Exploration of MOD09 reflectance QC Flags #######################################
#This script explores the processing of QC flags for MODIS reflectance (MOD09).       
#It relies on functions to help the processing of QC flags for MODIS.
#The current available products are:
# - reflectance (MOD09).       
# - land surface temperature
# - NDVI 
#
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 02/20/2018 
#DATE MODIFIED: 02/23/2018
#Version 1
#PROJECT: General use for MOD09A1 quality flag processing
#TO DO:
#
#COMMIT: testing the application of reflectance flags function
#

#################################################################################################

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
library(miscFuncs) #contains binary/bit processing helper functions
library(data.table)

###### Functions used in this script sourced from other files

script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts" #path on bpy50 #PARAM 2
functions_qc_modis_processing <- "QC_modis_processing_functions_02232018f.R"
source(file.path(script_path,functions_qc_modis_processing)) #source all functions used in this script 1.

##### Functions used in this script 

create_dir_fun <- function(out_dir,out_suffix){
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    out_dir <- file.path(out_dir,out_name)
  }
  #create if does not exists
  if(!file.exists(out_dir)){
    dir.create(out_dir)
  }
  return(out_dir)
}

#####  Parameters and argument set up ###########

#ARG1: input dir
in_dir <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_RITA_reflectance/import_h09v06/"
out_dir <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_RITA_reflectance/mask_qc_h09v06/"

## PRODUCT 3: Reflectance
# This is for MOD09
#
#https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod09a1
#This is 32 bits
r_qc_s1 <- raster(file.path(in_dir,"MOD09A1_A2005001_h09v06_006_sur_refl_qc_500m.tif"))
r_var_s1 <- raster(file.path(in_dir,"MOD09A1_A2005009_h09v06_006_reflectance.tif"))

list_r_qc_s <- list.files(path=in_dir,
                          pattern="MOD09A1_A200.*._h09v06_006_sur_refl_qc_500m.tif",
                          full.names=T)
list_r_var_s <- list.files(path=in_dir,
                           pattern="MOD09A1_A200.*._h09v06_006_reflectance.tif",
                           full.names=T)

dataType(r_qc_s1) #is FLT8S, 64 bits real numbers
out_dir_s <- getwd()
#INT4U
#https://www.rdocumentation.org/packages/raster/versions/2.6-7/topics/dataType
#INT4U is available but they are best avoided as R does not support 32-bit unsigned integers.
data_type_str <- "INT4U"
raster_name_tmp <- "MOD09A1_A2005001_h09v06_006_sur_refl_qc_500m_INTU4.tif"
writeRaster(r_qc_s1,
            filename=file.path(out_dir_s,raster_name_tmp),
            #bylayer=multiband,
            #suffix=paste(names(r),"_",out_suffix,sep=""),
            #format=format_raster,
            #suffix=paste(names(r)),
            overwrite=TRUE,
            #NAflag=NA_flag_val,
            datatype=data_type_str)
r <- raster(file.path(out_dir_s,raster_name_tmp))
unique_vals_test <- unique(r)
#The flags are 32 bits int
#library(binaryLogic)
unique_vals <- unique(r_qc_s1) #first get unique values
unique_vals==unique_vals_test
#convert decimal to intu32
bin(255)
1*2^7 + 1*2^6 + 1*2^5 + 1*2^4 + 1*2^3 + 1*2^2 + 1*2^1 + 1*2^0 
length(bin(255))

##hdf are in big endian format
bin_val

(2^30)
val

### Here is Big Endian format: least significant bit is on the right
### x* 2^31 + x* 2^30 + x* 2^29 + ... + x* 2^1 + x* 2^0

bin(61)
#[1] 1 1 1 1 0 1

length(bin(61))

0*2^7 + 0*2^6 + 1*2^5 + 1*2^4 + 1*2^3 + 1*2^2 + 0*2^1 + + 1*2^0 

#So you need to pad on the left!!!

intToBits(as.integer(65))
val <- (2^32)-1
intToBits(val)

test<- (intToBits(61))    
length(intToBits(61))
#rawToBits()
class(test)  
test[1]
#debug(generate_qc_MODIS_reflectance_MOD09_table)
qc_table_modis <- generate_qc_MODIS_reflectance_MOD09_table()

#View(qc_table_modis)

#### Now read in the values and process to match the selected flags!!!

i <- 1
val <- unique_vals[i]

#undebug(convert_decimal_to_uint32)
bin_val <- convert_decimal_to_uint32(val)
convert_to_decimal(bin_val)

### Do this for unique value 1:
i <- 1
val <- unique_vals[i]

bin_val <- convert_decimal_to_uint32(val)
convert_to_decimal(bin_val)

#### This is where you set up the desired qc flags:
desired_qc_rows <- c(1,2,5,14,23,32,41,50,59,69,71)

qc_table_modis_selected <- qc_table_modis[desired_qc_rows,]
#View(qc_table_modis_selected)

#debug(apply_mask_from_qc_layer)
#### Generate function to create mask from qc
#unique_bit_range <- unique(qc_table_modis$bitNo)

#apply_mask_from_qc_layer(r_qc= r_qc_s1, qc_table_modis=qc_table_modis,out_dir=out_dir_s,out_suffix="")
#
i <- 1
out_suffix_s <- "masked"

list_obj_mask <- apply_mask_from_qc_layer(i,
                         rast_qc=r_qc_s1,
                         rast_var=r_var_s1,
                         qc_table_modis_selected,
                         NA_flag_val= NULL,
                         rast_mask=TRUE,
                         qc_info=F,
                         out_dir=".",
                         out_suffix=out_suffix_s)

list_obj_mask <- mclapply(1:length(list_r_var_s)[1:2],
                          FUN= apply_mask_from_qc_layer,
                          rast_qc=list_r_qc_s,
                          rast_var=list_r_var_s,
                          qc_table_modis_selected,
                          NA_flag_val= NULL,
                          rast_mask=TRUE,
                          qc_info=F,
                          out_dir=".",
                          out_suffix=out_suffix_s,
                          mc.cores=num_cores,
                          mc.preschedule=T)

r_var_m <- raster(list_obj_mask$var)
r_mask <- raster(list_obj_mask$mask)

plot(stack(r_var_m,r_mask))

################################## End of script  #################################