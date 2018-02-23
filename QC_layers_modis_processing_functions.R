#######################################    MODIS processing   #######################################
############################  QC Flags processing related functions #######################################
#This script contains functions to help the processing of QC flags for MODIS.
#The current available products are:
# - reflectance (MOD09).       
# - land surface temperature
# - NDVI 
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 02/23/2017 
#DATE MODIFIED: 02/23/2018
#Version: 1

#PROJECT: General use
#TO DO:
#
#COMMIT: first commit for general modis flag processing functions
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

extract_qc_bit_info <- function(bit_range_processed,bin_val){
  #
  start_bit <- as.integer(unlist(strsplit(bit_range_processed,"-"))[1]) +1 #start at 1 in R
  end_bit <- as.integer(unlist(strsplit(bit_range_processed,"-"))[2]) +1
  bin_val_extracted <- bin_val[start_bit:end_bit] # corrected produced at less than ...
  bitNo <- paste(bit_range_processed,collapse = "-")
  bin_val_extracted <- paste(as.character(bin_val_extracted),collapse = "")
  #test <- paste(test,collapse="")
  #bin_val_extracted <- as.character(bin_val_extracted)
  #data.frame(bitNo=bit_range_processed,bin_val_extracted)
  df_bit_qc_info <-data.frame(bitNo=bit_range_processed,bin_val_extracted)
  return(df_bit_qc_info)
}

### This is little endian?
convert_to_decimal <- function(bin_val){
  n <- length(bin_val)
  vals<- lapply(1:n,function(i,coef){coef[i]*2^(n-i)},coef=bin_val)
  sum(unlist(vals))
}

screen_qc_bitNo <- function(qc_val,qc_table_modis_selected){
  
  unique_bitNo  <- unique(qc_table_modis_selected$bitNo)
  
  list_matching_val <- vector("list",length=nrow(qc_val))
  #for(i in 1:length(unique_bitNo)){
  for(i in 1:nrow(qc_val)){
    #bitNo_val <- qc_table_modis_selected$bitNo[1]
    #bitNo_val <- unique_bitNo[i]
    bitNo_val <- qc_val$bitNo[i]
    list_match_val[[i]] <- qc_val[qc_val$bitNo==bitNo_val,2]%in%qc_table_modis_selected[qc_table_modis_selected$bitNo==bitNo_val,c("BitComb")]
  }
  qc_val$screen_qc <- as.numeric(list_match_val)
  sum_val <- sum(qc_val$screen_qc)
  
  if(sum_val==nrow(qc_val)){
    mask_val <- 0 #do not mask, keep value
  }else{
    mask_val <- 0 #mask
  }
  
  obj_screen_qc_bitNo <- list(qc_val,mask_val)
  names(obj_screen_qc_bitNo) <- c("qc_val","mask_val")
  
  return(obj_screen_qc_bitNo)
}

generate_qc_MODIS_reflectance_MOD09_table <- function(){
  #This function generates a MODIS quality flag table based on the information
  #from LPDAAC:
  #https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod09a1
  #AUTHORS: Benoit Parmentier
  #CREATED: 02/23/2018
  #MODIFIED: 02/23/2018
  #
  
  ######## Start ######
  
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
  #View(qc_modis_table_bands)
  
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
  
  #View(qc_table_modis)
  
  write.table(qc_table_modis,file=file.path(out_dir,"qc_table_modis.txt"),sep=",")
  return(qc_table_modis)
}

generate_qc_val_from_int32_reflectance <- function(val){
  #
  #
  
  #i <- 1
  #val <- unique_vals[i]
  
  bin_val <- convert_decimal_to_uint32(val)
  convert_to_decimal(bin_val)
  
  test_bin_val_qc <- lapply(unique_bit_range,FUN=extract_qc_bit_info,bin_val=bin_val)
  test_bin_val_qc[[1]]
  
  qc_val <- do.call(rbind,test_bin_val_qc)
  #View(qc_val)
  return(qc_val)
}


generate_mask_from_qc_layer <- function(r_qc,qc_table_modis_selected,out_dir=".",out_suffix=""){
  #
  #
  
  ### Begin function ###
  
  if(class(r_qc)=="character"){
    r_qc <- raster(r_qc)
  }
  
  unique_vals <- unique(r_qc) #first get unique values
  
  ### Do this for unique value 1:
  
  list_qc_val <- lapply(unique_vals,FUN=generate_qc_val_from_int32_reflectance)
  #length(list_qc_val)
  length(unique_vals)
  
  ### Find if need to mask the values:
  #### Now compare:
  list_qc_val_screened <- lapply(list_qc_val,
                                 FUN=screen_qc_bitNo,
                                 qc_table_modis = qc_table_modis)
  names(list_qc_val_screened[[1]])
  rm(list_qc_val)
  
  mask_val <- unlist(lapply(list_qc_val_screened,FUN=function(x){x$mask_val}))
  
  df_qc_val <- data.frame(val=unique_vals,mask_val=mask_val)
  #View(df_qc_val)
  
  ###### Reclassify QC image if necessary:
  
  if(sum(df_qc_val$mask_val)==0){
    qc_mask <- init(r_qc,0)
  }else{
    qc_mask <- subs(r_qc, df_qc_val,by="val",which="mask_val")
  }
  
  #####
  writeRaster()
  
  return()
}

################################## End of script  #################################