#######################################    MODIS processing   #######################################
############################ Exploration of MOD09 reflectance QC Flags #######################################
#This script explores the processing of QC flags for MODIS reflectance (MOD09).       
#
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/07/2017 
#DATE MODIFIED: 02/21/2018
#Version: 1

#PROJECT: Space beats time Framework
#TO DO:
#
#COMMIT: splitting function and main script for the assessment
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

#function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_11242015_functions.R" #PARAM 1
function_paper_figures_analyses <- "space_beats_time_sbt_paper_figures_functions_02102018.R" #PARAM 1
function_space_and_time_assessment <- "space_and_time_assessment_functions_02202018.R" #PARAM 1

script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts" #path on bpy50 #PARAM 2
#script_path <- "/home/parmentier/Data/Space_beats_time/sbt_scripts" #path on Atlas
#source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.
source(file.path(script_path,function_paper_figures_analyses)) #source all functions used in this script 1.
source(file.path(script_path,function_space_and_time_assessment)) #source all functions used in this script 1.

#zonal stratum for NLU
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

convert_decimal_to_uint32 <- function(x){
  bin_val <- bin(x)
  bin_val <- as.numeric(as.logical(bin_val))
  
  if(length(bin_val)< 32){
    add_zero <- 32-length(bin_val) 
  }
  
  bin_val <- c(rep(0,add_zero),bin_val)
  
  return(bin_val)
}

### This is little endian?
convert_to_decimal <- function(bin_val){
  n <- length(bin_val)
  vals<- lapply(1:n,function(i,coef){coef[i]*2^(n-i)},coef=bin_val)
  sum(unlist(vals))
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

dataType(r_qc_s1) #is FLT8S, 64 bits real numbers
writeRaster(r_qc_s1,data) #is FLT8S, 64 bits real numbers
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
class(as.integer(65))
#line1<-c(readBin(to.read,"int",5), 
#         readBin(to.read,"double",1,size=4),
#         readBin(to.read,"int",2))
test <- (unique_vals[1]) #no need to do this
test <- 10
bin_val<- bin(test)
bin_val
bin_val <- as.numeric(as.logical(bin_val))

debug(convert_decimal_to_uint32)
convert_decimal_to_uint32(61)
bin(61)

#add_zero <- 2
#rep(0,add_zero)

#debug(convert_to_decimal)
test
bin(test)
convert_to_decimal(bin_val)
test

rep(NA,64)
#qc_table_modis <- data.table(bitNo=NA,param_name=NA,BitComb=NA,Sur_refl_qc_500m=NA)
qc_table_modis <- data.frame(bitNo=NA,param_name=NA,BitComb=NA,Sur_refl_qc_500m=NA)
### Make this a function
### now add band quality
band_no <- 1
list_band_no <- 1:7
list_qc_table_modis <- vector("list",length=length(list_band_no))
start_bit <- 0
i <- 1
for(i in 1:length(list_band_no)){
  band_no <- list_band_no[i]
  if(i==1){
    start_bit <- 2
  }else{
    start_bit <- end_bit + 1
  }
  #start_bit <- 2
  end_bit <- start_bit + 3 
  bit_range <- paste(start_bit,end_bit,sep="-")
  qc_table_modis[1:9,c("bitNo")] <- rep(bit_range,9) 
  
  qc_table_modis[1:9,c("param_name")] <- rep(paste("band",band_no,"quality four bit range"),9) 
  qc_table_modis[1:9,c("BitComb")] <- c("0000","1000","1001","1010","1011","1100","1101","1110","1111") 
  qc_table_modis[1:9,c("Sur_refl_qc_500m")] <- c("highest quality",
                                                  "dead detector; data interpolated in L1B",
                                                  "solar zenith >=86 degrees",
                                                  "solar zenith >=85 degrees and < 86 degrees",
                                                  "missing input",
                                                  "internal constant used in place of climatological data for at least one atmospheric constant",
                                                  "correction out of bounds pixel constrained to extreme allowable value",
                                                  "L1B data faulty",
                                                  "not processed due to deep ocean or clouds") 
  list_qc_table_modis[[i]] <- qc_table_modis
}

qc_modis_table_bands <- do.call(rbind,list_qc_table_modis)
View(qc_modis_table_bands)

qc_table_modis_modland <- data.frame(bitNo=NA,param_name=NA,BitComb=NA,Sur_refl_qc_500m=NA)

qc_table_modis_modland[1:4,c("bitNo")] <- rep("0-1",4) 
qc_table_modis_modland[1:4,c("param_name")] <- rep("MODLAND QA bits",4) 
qc_table_modis_modland[1:4,c("BitComb")] <- c("00","01","10","11") 
qc_table_modis_modland[1:4,c("Sur_refl_qc_500m")] <- c("corrected product produced at ideal quality all bands",
                                               "corrected product produced at less than ideal quality some or all bands",
                                               "corrected product not produced due to cloud effects all bands",
                                               "corrected producted not produced due to other reasons some or all bands may be fill value [Note that a value of (11) overrides a value of (01)]") 

qc_table_modis_other <- data.frame(bitNo=NA,param_name=NA,BitComb=NA,Sur_refl_qc_500m=NA)
qc_table_modis_other[1:4,c("bitNo")] <- c(rep("30-30",2),
                                          rep("31-31",2))
qc_table_modis_other[1:4,c("param_name")] <- c(rep("atmospheric correction performed",2),
                                               rep("adjacency correction performed",2))
qc_table_modis_other[1:4,c("BitComb")] <- c("1","0","1","0") 
qc_table_modis_other[1:4,c("Sur_refl_qc_500m")] <- c("yes","no","yes","no") 

#qc_table_modis <- rbind(qc_table_modis_modland,qc_modis_table_bands)
qc_table_modis <- rbind(qc_table_modis_modland,qc_modis_table_bands,qc_table_modis_other)

View(qc_table_modis)

write.table(qc_table_modis,file=file.path(out_dir,"qc_table_modis.txt"),sep=",")

### 

#selected_flags <- list(QA_word1 ="VI Good Quality",QA_word1 ="VI Produced,check QA") #if NULL use default
selected_flags <- list(Sur_refl_qc_500m ="highest quality",
                       Sur_refl_qc_500m = "corrected product produced at ideal quality all bands",
                       Sur_refl_qc_500m = "corrected product produced at less than ideal quality some or all bands")

#sbuset(qc_table_modis
#lapply(selected_flags,function(x){subset(qc_table_modis,})

#### Now read in the values and process to match the selected flags!!!

i <- 1
val <- unique_vals[i]

bin_val <- convert_decimal_to_uint32(val)
convert_to_decimal(bin_val)

##hdf are in big endian format
bin_val

(2^30)
val

### Here is Big Endian format: least significant bit is on the right
### x* 2^31 + x* 2^30 + x* 2^29 + ... + x* 2^1 + x* 2^0

?bitShiftR()
bitShiftR(val,1)
(2^30)/2

bitShiftR(val,1)
(2^30)/2
?bitwShiftL()
unique_bit_range <- unique(qc_table_modis$bitNo)

i <- 1 #modland

bit_range_processed <- unique_bit_range[i]
start_bit <- as.integer(unlist(strsplit(bit_range_processed,"-"))[1]) +1 #start at 1 in R
end_bit <- as.integer(unlist(strsplit(bit_range_processed,"-"))[2]) +1
bin_val[start_bit:end_bit] # corrected produced at less than ...

i <- 2 #band 1

bit_range_processed <- unique_bit_range[i]
start_bit <- as.integer(unlist(strsplit(bit_range_processed,"-"))[1]) +1 #start at 1 in R
end_bit <- as.integer(unlist(strsplit(bit_range_processed,"-"))[2]) +1
bin_val[start_bit:end_bit] # corrected produced at less than ...

i <- 3 #band 2

bit_range_processed <- unique_bit_range[i]
start_bit <- as.integer(unlist(strsplit(bit_range_processed,"-"))[1]) +1 #start at 1 in R
end_bit <- as.integer(unlist(strsplit(bit_range_processed,"-"))[2]) +1
bin_val[start_bit:end_bit] # corrected produced at less than ...

extract_qc_bit_info <- function(bit_range,bin_val){
  bit_range_processed <- unique_bit_range[i]
  start_bit <- as.integer(unlist(strsplit(bit_range_processed,"-"))[1]) +1 #start at 1 in R
  end_bit <- as.integer(unlist(strsplit(bit_range_processed,"-"))[2]) +1
  bin_val[start_bit:end_bit] # corrected produced at less than ...
  
}

#####################

#Now generate reference table
#(create_qa_mask)

#' A function to translate QA information into a datamask
#' 
#' @param qa_bits numeric. Number of bits of info in the QA
#' @param wordinbits logical. Should the bit word parsing stay in bit format, or translate to numeric within the word
#' @param bit_interpreter data.table. data.table with columns word_number, word_start, word_end
#'
create_qa_mask = function(qa_bits, wordinbits = T, bit_interpreter){
  
  #create a table that translates integer values to bits
  cw = data.table::data.table(intval = 0:(2^qa_bits-1))
  
  #GEE converts the bits into their integer value per word. That conversion is better if least significant digit is on the left
  #modis appears to be the opposite
  if(wordinbits){
    bitorder = qa_bits:1
  }else{
    bitorder = 1:qa_bits
  }
  
  bits = lapply(cw[,intval], function(x) as.integer(intToBits(x))[bitorder]) #starting val is on the lhs.
  bits = do.call(rbind, bits)
  cw = cbind(cw, bits)
  setnames(cw, c('intval',paste0('V',bitorder - 1)))
  
  #create words
  for(qqq in seq(nrow(bit_interpreter))){
    
    
    start = bit_interpreter[qqq,word_start]
    end = bit_interpreter[qqq,word_end]
    i = bit_interpreter[qqq,word_number]
    
    if(wordinbits){
      cw[, paste0('word',i) := apply(.SD, 1, paste, collapse = ""), .SDcols = paste0('V',end:start)]
    } else{
      #convert the word to its integer counterpart
      cw[, paste0('word',i) := rowSums(.SD[, lapply(1:ncol(.SD), function(x) 2^(x-1) * .SD[[x]])]), .SDcols = paste0('V',start:end)] 
    }
    
  }
  
  #return the grid
  return(cw)
}

build_lst_qa = function(){
  bi = data.table(word_number = 1:4 ,
                  word_start = c(0,2,4,6),
                  word_end = c(1,3,5,7))
  qa = create_qa_mask(8, F, bi)
  
  #return rows where the data is good at word 1 or word 2
  qa = qa[word1 ==0 | (word1 == 1 & word2 == 0), intval]
  
  return(qa)
  
}

build_vi_qa = function(){
  bi = data.table(word_number = 1:9 ,
                  word_start = c(0,2,6,8,9,10,11,14,15),
                  word_end =   c(1,5,7,8,9,10,13,14,15))
  qa = create_qa_mask(16, F, bi)
  
  #select for data quality
  qa = qa[word1 ==0 | (word1 == 1 & word2 %in% c(0,1,2,4,8,9,10)), ]
  
  #select for land
  qa = qa[word7 %in% c(1,2),intval]
  
  return(qa)
  
}




###########

#length(intToBits(test[1])) 
length(intToBits(as.integer(unique_vals[1]))) #that is 32
str_bit <- as.character(intToBits(as.integer(unique_vals[1])))
str_bit

length(str_bit)
str_bit[1]

#range for 32 bits is 0 to 2^32 -1
(2^32)-1
unique_vals > (2^32)-1

test <- data.table(INTU4_val=0:2^32-1)
library(data.table)
#https://stackoverflow.com/questions/43274706/apply-function-bitwise-and-on-each-cell-of-a-raster-in-r/43274922
#as.binary(test)
#class(intToBits(test))
#potential logic:
# Write table of bits as described in the documentation
# Extract unique values from raster
# Convert unique values into bit sequence
# Match bit sequence to table
# screen for only sequence that we want
# 
# Other logic:
# 
  

create_MODIS_QC_table <-function(LST=TRUE, NDVI=TRUE,reflectance=TRUE){
  #Function to generate MODIS QC  flag table
  #Author: Benoit Parmentier (with some lines from S.Mosher)
  #Date CREATED: 09/16/2013
  #Date MODIFIED: 02/15/2018
  #
  #Some of the inspiration and code originates from Steve Mosher' s blog:
  #http://stevemosher.wordpress.com/2012/12/05/modis-qc-bits/
  
  #r_qc_s1 <- raster("/home/bparmentier/Google Drive/Space_beats_time/Data/data_RITA_reflectance/import_h09v06/MOD09A1_A2005001_h09v06_006_sur_refl_qc_500m.tif")
  
  list_QC_Data <- vector("list", length=3)
  names(list_QC_Data) <- c("LST","NDVI","reflectance")
  
  ## PRODUCT 1: LST
  #This can be seen from table defined at LPDAAC: https://lpdaac.usgs.gov/products/modis_products_table/mod11a2
  #LST MOD11A2 has 4 levels/indicators of QA:
  
  ## Generate product table
  if (LST==TRUE){
    QC_Data <- data.frame(Integer_Value = 0:255,
                          Bit7 = NA,Bit6 = NA,Bit5 = NA,Bit4 = NA,Bit3 = NA,Bit2 = NA,Bit1 = NA,Bit0 = NA,
                          QA_word1 = NA,QA_word2 = NA,QA_word3 = NA,QA_word4 = NA)
    #Populate table/data frame
    for(i in QC_Data$Integer_Value){
      AsInt <- as.integer(intToBits(i)[1:8])
      QC_Data[i+1,2:9]<- AsInt[8:1]
    } 
    #Level 1: Overal MODIS Quality which is common to all MODIS product
    QC_Data$QA_word1[QC_Data$Bit1 == 0 & QC_Data$Bit0==0] <- "LST Good Quality"    #(0-0)
    QC_Data$QA_word1[QC_Data$Bit1 == 0 & QC_Data$Bit0==1] <- "LST Produced,Check QA"
    QC_Data$QA_word1[QC_Data$Bit1 == 1 & QC_Data$Bit0==0] <- "Not Produced,clouds"
    QC_Data$QA_word1[QC_Data$Bit1 == 1 & QC_Data$Bit0==1] <- "No Produced, check Other QA"
    
    #Level 2: Information on quality of product (i.e. LST produced, Check QA) for LST
    QC_Data$QA_word2[QC_Data$Bit3 == 0 & QC_Data$Bit2==0] <- "Good Data"
    QC_Data$QA_word2[QC_Data$Bit3 == 0 & QC_Data$Bit2==1] <- "Other Quality"
    QC_Data$QA_word2[QC_Data$Bit3 == 1 & QC_Data$Bit2==0] <- "TBD"
    QC_Data$QA_word2[QC_Data$Bit3 == 1 & QC_Data$Bit2==1] <- "TBD"
    
    #Level 3: Information on quality of of emissitivity 
    QC_Data$QA_word3[QC_Data$Bit5 == 0 & QC_Data$Bit4==0] <- "Emiss Error <= .01"
    QC_Data$QA_word3[QC_Data$Bit5 == 0 & QC_Data$Bit4==1] <- "Emiss Err >.01 <=.02"
    QC_Data$QA_word3[QC_Data$Bit5 == 1 & QC_Data$Bit4==0] <- "Emiss Err >.02 <=.04"
    QC_Data$QA_word3[QC_Data$Bit5 == 1 & QC_Data$Bit4==1] <- "Emiss Err > .04"
    
    #Level 4: Uncertaing for LST error
    QC_Data$QA_word4[QC_Data$Bit7 == 0 & QC_Data$Bit6==0] <- "LST Err <= 1"
    QC_Data$QA_word4[QC_Data$Bit7 == 0 & QC_Data$Bit6==1] <- "LST Err > 2 LST Err <= 3"
    QC_Data$QA_word4[QC_Data$Bit7 == 1 & QC_Data$Bit6==0] <- "LST Err > 1 LST Err <= 2"
    QC_Data$QA_word4[QC_Data$Bit7 == 1 & QC_Data$Bit6==1] <- "LST Err > 4"
    
    list_QC_Data[[1]] <- QC_Data
  }
  
  ## PRODUCT 2: NDVI
  #This can be seen from table defined at LPDAAC: https://lpdaac.usgs.gov/products/modis_products_table/mod11a2
  
  if(NDVI==TRUE){
    QC_Data <- data.frame(Integer_Value = 0:65535,
                          Bit15 = NA,Bit14 = NA,Bit13 = NA,Bit12 = NA,Bit11 = NA,Bit10 = NA,Bit9 = NA,Bit8 = NA,
                          Bit7 = NA,Bit6 = NA,Bit5 = NA,Bit4 = NA,Bit3 = NA,Bit2 = NA,Bit1 = NA,Bit0 = NA,
                          QA_word1 = NA,QA_word2 = NA,QA_word3 = NA,QA_word4 = NA,
                          QA_word5 = NA,QA_word6 = NA,QA_word7 = NA,QA_word8 = NA,
                          QA_word9 = NA)
    #Populate table...this is extremely slow...change???
    for(i in QC_Data$Integer_Value){
      AsInt <- as.integer(intToBits(i)[1:16]) #16bit unsigned integer
      QC_Data[i+1,2:17]<- AsInt[16:1]
    } 
    
    #Level 1: Overal MODIS Quality which is common to all MODIS product
    QC_Data$QA_word1[QC_Data$Bit1 == 0 & QC_Data$Bit0==0] <- "VI Good Quality"    #(0-0)
    QC_Data$QA_word1[QC_Data$Bit1 == 0 & QC_Data$Bit0==1] <- "VI Produced,check QA"
    QC_Data$QA_word1[QC_Data$Bit1 == 1 & QC_Data$Bit0==0] <- "Not Produced,because of clouds"
    QC_Data$QA_word1[QC_Data$Bit1 == 1 & QC_Data$Bit0==1] <- "Not Produced, other reasons"
    
    #Level 2: VI usefulness (read from right to left)
    QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==0 & QC_Data$Bit3 == 0 & QC_Data$Bit2==0] <- "Highest quality, 1"
    QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==0 & QC_Data$Bit3 == 0 & QC_Data$Bit2==1] <- "Lower quality, 2"
    QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==0 & QC_Data$Bit3 == 1 & QC_Data$Bit2==0] <- "Decreasing quality, 3 "
    QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==0 & QC_Data$Bit3 == 1 & QC_Data$Bit2==1] <- "Decreasing quality, 4"
    QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==1 & QC_Data$Bit3 == 0 & QC_Data$Bit2==0] <- "Decreasing quality, 5"
    QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==1 & QC_Data$Bit3 == 0 & QC_Data$Bit2==1] <- "Decreasing quality, 6"
    QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==1 & QC_Data$Bit3 == 1 & QC_Data$Bit2==0] <- "Decreasing quality, 7"
    QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==1 & QC_Data$Bit3 == 1 & QC_Data$Bit2==1] <- "Decreasing quality, 8"
    QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==0 & QC_Data$Bit3 == 0 & QC_Data$Bit2==0] <- "Decreasing quality, 9"
    QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==0 & QC_Data$Bit3 == 0 & QC_Data$Bit2==1] <- "Decreasing quality, 10"
    QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==0 & QC_Data$Bit3 == 1 & QC_Data$Bit2==0] <- "Decreasing quality, 11"
    QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==0 & QC_Data$Bit3 == 1 & QC_Data$Bit2==1] <- "Decreasing quality, 12"
    QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==1 & QC_Data$Bit3 == 0 & QC_Data$Bit2==0] <- "Lowest quality, 13"
    QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==1 & QC_Data$Bit3 == 0 & QC_Data$Bit2==1] <- "Quality so low that not useful, 14"
    QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==1 & QC_Data$Bit3 == 1 & QC_Data$Bit2==0] <- "L1B data faulty, 15"
    QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==1 & QC_Data$Bit3 == 1 & QC_Data$Bit2==1] <- "Not useful/not processed, 16"
    
    # Level 3: Aerosol quantity 
    QC_Data$QA_word3[QC_Data$Bit7 == 0 & QC_Data$Bit6==0] <- "Climatology"
    QC_Data$QA_word3[QC_Data$Bit7 == 0 & QC_Data$Bit6==1] <- "Low"
    QC_Data$QA_word3[QC_Data$Bit7 == 1 & QC_Data$Bit6==0] <- "Average"
    QC_Data$QA_word3[QC_Data$Bit7 == 1 & QC_Data$Bit6==1] <- "High"
    
    # Level 4: Adjacent cloud detected
    QC_Data$QA_word4[QC_Data$Bit8==0] <- "No"
    QC_Data$QA_word4[QC_Data$Bit8==1] <- "Yes"
    
    # Level 5: Atmosphere BRDF correction performed
    QC_Data$QA_word5[QC_Data$Bit9 == 0] <- "No"
    QC_Data$QA_word5[QC_Data$Bit9 == 1] <- "Yes"
    
    # Level 6: Mixed Clouds
    QC_Data$QA_word6[QC_Data$Bit10 == 0] <- "No"
    QC_Data$QA_word6[QC_Data$Bit10 == 1] <- "Yes"
    
    #Level 7: Land/Water Flag (read from right to left)
    QC_Data$QA_word7[QC_Data$Bit13==0 & QC_Data$Bit12 == 0 & QC_Data$Bit11==0] <- "Shallow Ocean"
    QC_Data$QA_word7[QC_Data$Bit13==0 & QC_Data$Bit12 == 0 & QC_Data$Bit11==1] <- "Land"
    QC_Data$QA_word7[QC_Data$Bit13==0 & QC_Data$Bit12 == 1 & QC_Data$Bit11==0] <- "Ocean coastlines and lake shorelines"
    QC_Data$QA_word7[QC_Data$Bit13==0 & QC_Data$Bit12 == 1 & QC_Data$Bit11==1] <- "Shallow inland water"
    QC_Data$QA_word7[QC_Data$Bit13==1 & QC_Data$Bit12 == 0 & QC_Data$Bit11==0] <- "Ephemeral water"
    QC_Data$QA_word7[QC_Data$Bit13==1 & QC_Data$Bit12 == 0 & QC_Data$Bit11==1] <- "Deep inland water"
    QC_Data$QA_word7[QC_Data$Bit13==1 & QC_Data$Bit12 == 1 & QC_Data$Bit11==0] <- "Moderate or continental water"
    QC_Data$QA_word7[QC_Data$Bit13==1 & QC_Data$Bit12 == 1 & QC_Data$Bit11==1] <- "Deep ocean"
    
    # Level 8: Possible snow/ice
    QC_Data$QA_word8[QC_Data$Bit14 == 0] <- "No"
    QC_Data$QA_word8[QC_Data$Bit14 == 1] <- "Yes"
    
    # Level 9: Possible shadow
    QC_Data$QA_word9[QC_Data$Bit15 == 0] <- "No"
    QC_Data$QA_word9[QC_Data$Bit15 == 1] <- "Yes"
    
    list_QC_Data[[2]]<- QC_Data
  }
  
  ## PRODUCT 3: Reflectance
  # This is for MOD09
  #
  if(reflectance==TRUE){
    #https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod09a1
    #This is 32 bits
    #r_qc_s1 <- raster("/home/bparmentier/Google Drive/Space_beats_time/Data/data_RITA_reflectance/import_h09v06/MOD09A1_A2005001_h09v06_006_sur_refl_qc_500m.tif")
    #dataType(r_qc_s1) is FLT8S, 64 bits real numbers
    #The flags are 32 bits int
    #library(binaryLogic)
    test <- unique(r_qc_s1) #first get unique values
    intToBits(65)    
    length(intToBits(65))
    #rawToBits()
    
    line1<-c(readBin(to.read,"int",5), 
             readBin(to.read,"double",1,size=4),
             readBin(to.read,"int",2))
    as.integer(test[1]) #no need to do this
    #length(intToBits(test[1])) 
    length(intToBits(test[1])) #that is 32
    str_bit <- as.character(intToBits(test[1]))
    
    
    #https://stackoverflow.com/questions/43274706/apply-function-bitwise-and-on-each-cell-of-a-raster-in-r/43274922
    #as.binary(test)
    #class(intToBits(test))
    #potential logic:
    # Write table of bits as described in the documentation
    # Extract unique values from raster
    # Convert unique values into bit sequence
    # Match bit sequence to table
    # screen for only sequence that we want
    # 
    # Other logic:
    # 
    
    QC_Data <- data.frame(Integer_Value = 0:65535,
                          Bit15 = NA,Bit14 = NA,Bit13 = NA,Bit12 = NA,Bit11 = NA,Bit10 = NA,Bit9 = NA,Bit8 = NA,
                          Bit7 = NA,Bit6 = NA,Bit5 = NA,Bit4 = NA,Bit3 = NA,Bit2 = NA,Bit1 = NA,Bit0 = NA,
                          QA_word1 = NA,QA_word2 = NA,QA_word3 = NA,QA_word4 = NA,
                          QA_word5 = NA,QA_word6 = NA,QA_word7 = NA,QA_word8 = NA,
                          QA_word9 = NA)
    
    
    #QC_Data <- data.frame(,
    #                      0,1, NA,NA,..,NA #
    #                      0,0
    #                      Bit15 = NA,Bit14 = NA,Bit13 = NA,Bit12 = NA,Bit11 = NA,Bit10 = NA,Bit9 = NA,Bit8 = NA,
    #                      Bit7 = NA,Bit6 = NA,Bit5 = NA,Bit4 = NA,Bit3 = NA,Bit2 = NA,Bit1 = NA,Bit0 = NA,
    #                      QA_word1 = NA,QA_word2 = NA,QA_word3 = NA,QA_word4 = NA,
    #                      QA_word5 = NA,QA_word6 = NA,QA_word7 = NA,QA_word8 = NA,
    #                      QA_word9 = NA)
    
    #Populate table...this is extremely slow...change???
    for(i in QC_Data$Integer_Value){
      AsInt <- as.integer(intToBits(i)[1:16]) #16bit unsigned integer
      QC_Data[i+1,2:17]<- AsInt[16:1]
    } 
    
    #Level 1: Overal MODIS Quality which is common to all MODIS product
    QC_Data$QA_word1[QC_Data$Bit1 == 0 & QC_Data$Bit0==0] <- "VI Good Quality"    #(0-0)
    QC_Data$QA_word1[QC_Data$Bit1 == 0 & QC_Data$Bit0==1] <- "VI Produced,check QA"
    QC_Data$QA_word1[QC_Data$Bit1 == 1 & QC_Data$Bit0==0] <- "Not Produced,because of clouds"
    QC_Data$QA_word1[QC_Data$Bit1 == 1 & QC_Data$Bit0==1] <- "Not Produced, other reasons"
    
    #Level 2: VI usefulness (read from right to left)
    QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==0 & QC_Data$Bit3 == 0 & QC_Data$Bit2==0] <- "Highest quality, 1"
    QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==0 & QC_Data$Bit3 == 0 & QC_Data$Bit2==1] <- "Lower quality, 2"
    QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==0 & QC_Data$Bit3 == 1 & QC_Data$Bit2==0] <- "Decreasing quality, 3 "
    QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==0 & QC_Data$Bit3 == 1 & QC_Data$Bit2==1] <- "Decreasing quality, 4"
    QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==1 & QC_Data$Bit3 == 0 & QC_Data$Bit2==0] <- "Decreasing quality, 5"
    QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==1 & QC_Data$Bit3 == 0 & QC_Data$Bit2==1] <- "Decreasing quality, 6"
    QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==1 & QC_Data$Bit3 == 1 & QC_Data$Bit2==0] <- "Decreasing quality, 7"
    QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==1 & QC_Data$Bit3 == 1 & QC_Data$Bit2==1] <- "Decreasing quality, 8"
    QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==0 & QC_Data$Bit3 == 0 & QC_Data$Bit2==0] <- "Decreasing quality, 9"
    QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==0 & QC_Data$Bit3 == 0 & QC_Data$Bit2==1] <- "Decreasing quality, 10"
    QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==0 & QC_Data$Bit3 == 1 & QC_Data$Bit2==0] <- "Decreasing quality, 11"
    QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==0 & QC_Data$Bit3 == 1 & QC_Data$Bit2==1] <- "Decreasing quality, 12"
    QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==1 & QC_Data$Bit3 == 0 & QC_Data$Bit2==0] <- "Lowest quality, 13"
    QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==1 & QC_Data$Bit3 == 0 & QC_Data$Bit2==1] <- "Quality so low that not useful, 14"
    QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==1 & QC_Data$Bit3 == 1 & QC_Data$Bit2==0] <- "L1B data faulty, 15"
    QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==1 & QC_Data$Bit3 == 1 & QC_Data$Bit2==1] <- "Not useful/not processed, 16"
    
    # Level 3: Aerosol quantity 
    QC_Data$QA_word3[QC_Data$Bit7 == 0 & QC_Data$Bit6==0] <- "Climatology"
    QC_Data$QA_word3[QC_Data$Bit7 == 0 & QC_Data$Bit6==1] <- "Low"
    QC_Data$QA_word3[QC_Data$Bit7 == 1 & QC_Data$Bit6==0] <- "Average"
    QC_Data$QA_word3[QC_Data$Bit7 == 1 & QC_Data$Bit6==1] <- "High"
    
    # Level 4: Adjacent cloud detected
    QC_Data$QA_word4[QC_Data$Bit8==0] <- "No"
    QC_Data$QA_word4[QC_Data$Bit8==1] <- "Yes"
    
    # Level 5: Atmosphere BRDF correction performed
    QC_Data$QA_word5[QC_Data$Bit9 == 0] <- "No"
    QC_Data$QA_word5[QC_Data$Bit9 == 1] <- "Yes"
    
    # Level 6: Mixed Clouds
    QC_Data$QA_word6[QC_Data$Bit10 == 0] <- "No"
    QC_Data$QA_word6[QC_Data$Bit10 == 1] <- "Yes"
    
    #Level 7: Land/Water Flag (read from right to left)
    QC_Data$QA_word7[QC_Data$Bit13==0 & QC_Data$Bit12 == 0 & QC_Data$Bit11==0] <- "Shallow Ocean"
    QC_Data$QA_word7[QC_Data$Bit13==0 & QC_Data$Bit12 == 0 & QC_Data$Bit11==1] <- "Land"
    QC_Data$QA_word7[QC_Data$Bit13==0 & QC_Data$Bit12 == 1 & QC_Data$Bit11==0] <- "Ocean coastlines and lake shorelines"
    QC_Data$QA_word7[QC_Data$Bit13==0 & QC_Data$Bit12 == 1 & QC_Data$Bit11==1] <- "Shallow inland water"
    QC_Data$QA_word7[QC_Data$Bit13==1 & QC_Data$Bit12 == 0 & QC_Data$Bit11==0] <- "Ephemeral water"
    QC_Data$QA_word7[QC_Data$Bit13==1 & QC_Data$Bit12 == 0 & QC_Data$Bit11==1] <- "Deep inland water"
    QC_Data$QA_word7[QC_Data$Bit13==1 & QC_Data$Bit12 == 1 & QC_Data$Bit11==0] <- "Moderate or continental water"
    QC_Data$QA_word7[QC_Data$Bit13==1 & QC_Data$Bit12 == 1 & QC_Data$Bit11==1] <- "Deep ocean"
    
    # Level 8: Possible snow/ice
    QC_Data$QA_word8[QC_Data$Bit14 == 0] <- "No"
    QC_Data$QA_word8[QC_Data$Bit14 == 1] <- "Yes"
    
    # Level 9: Possible shadow
    QC_Data$QA_word9[QC_Data$Bit15 == 0] <- "No"
    QC_Data$QA_word9[QC_Data$Bit15 == 1] <- "Yes"
    
    list_QC_Data[[2]]<- QC_Data
  }
  
  
  
  
  ## PRODUCT 4: Albedo
  #This can be seen from table defined at LPDAAC: https://lpdaac.usgs.gov/products/modis_products_table/mod11a2
  #To be added...
  
  ###Now return and save object:
  #Prepare object to return
  
  save(list_QC_Data,file= file.path(".",paste("list_QC_Data",".RData",sep="")))
  
  return(list_QC_Data)
}


################################## End of script  #################################