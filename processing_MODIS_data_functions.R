################################################  PROCESSING MODIS DATA #######################################
########################################### Function to process MODIS datasets #####################################
#The current version processes MODIS dataset given a set of inputs.
#This script download and processes MODIS tiles using Quality Flag. 
#Tiles are mosaiced and reprojected for a specific study region.
#MODIS currently stores information in HDF4 format. Layers must be extracted and must be listed first
#using for example gdalinfo to identify the relevant subdatasets and QC flag options. 
#Note that QC flags are store bitpacks of 8bits (byte) in big endian!!!
#A data frame matching flag values is created to facilitate the processing.            
#Inspiration and some code for the MODIS flag function originates from Steve Mosher:
#http://stevemosher.wordpress.com/2012/12/05/modis-qc-bits/
#Note that the downloading step requires a .netrc file and and login to EARTHDATA.
#See for more information on the LPDAAC data pool: https://lpdaac.usgs.gov/data_access/data_pool
#
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
#MODIFIED ON : 07/10/2018
#PROJECT: General MODIS processing of all projects
#COMMIT: dealing with multibands outputs in import
#
#TODO: 
#1)Test additional Quality Flag levels for ALBEDO and other products (MOD09)
#2) Add function to report statistics: missing files
#3) Currently 20 input arguments (param), reduce to 15 or less
#4) Make this script a function callable from shell!!
#5) This script can be transformed to process other datasets using the https://lpdaac.usgs.gov/data_access/data_pool
#   e.g."https://e4ftl01.cr.usgs.gov/WELD/" for WELD datasets.
#6) adding multiband option for import of MOD09


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


processing_modis_data <- function(in_dir,
                                  out_dir,
                                  CRS_reg,
                                  method_proj_val,
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
                                  multiband,
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
                                  qc_info,
                                  out_dir_processing_steps){
  

  #AUTHOR: Benoit Parmentier                                                                       
  #CREATED ON : 02/08/2018  
  #MODIFIED ON : 07/10/2018
  
  
  ##################################
  
  #Create output directory
  
  if(create_out_dir_param==TRUE){
    out_dir <- create_dir_fun(out_dir,out_suffix)
    setwd(out_dir)
  }else{
    setwd(out_dir) #use previoulsy defined directory
  }
  
  #####################################
  #### STEP 1:  DOWNLOAD MODIS PRODUCT  ####
  
  date_param <- unlist(strsplit(date_param,";"))
  start_date <- date_param[1] 
  end_date <- date_param[2]
  time_step <- as.integer(date_param[3])
  
  ### If tiles to download, determine using region outline
  if(is.null(list_tiles_modis)){
    # generate figures of selected tiles later on?
    obj <- get_modis_tiles_list(infile_modis_grid,
                                reg_outline=infile_reg_outline,
                                CRS_reg)
    list_tiles_modis <- obj$tiles_modis
    
  }
  list_tiles_modis <- unlist(strsplit(list_tiles_modis,","))  # transform string into separate element in char vector
  
  ### Plot area with kml before donwloading:
  ### This can be displayed using Google Earth
  if(!is.null(ref_rast_name)){
    r <- raster(ref_rast_name)
    filename_kml <- file.path(out_dir,"ref_rast.kml")
    projection(r) <- CRS_reg
    r_WGS84 <- projectRaster(r,crs=CRS_WGS84)
    KML(r_WGS84,filename_kml , col=rev(terrain.colors(255)), 
        colNA=NA, zip='', overwrite=TRUE)
  }
  
  if(!is.null(infile_reg_outline)){
    # S4 method for Spatial
    reg_sf <- st_read(infile_reg_outline)
    reg_spdf <- as(reg_sf,"Spatial")
    filename_kml <- file.path(out_dir,"reg_outline.kml")
    KML(reg_spdf, filename_kml, zip='', overwrite=TRUE)
  }
  
  #debug(modis_product_download)
  #browser()
  if(steps_to_run$download==TRUE){
    #debug(modis_product_download)
    #9:36 to 10: for MOD09 4 tiles
    
    download_modis_obj <- modis_product_download(MODIS_product,
                                                 product_version,
                                                 start_date,
                                                 end_date,
                                                 list_tiles_modis,
                                                 file_format_download,
                                                 out_dir,
                                                 temporal_granularity)
    #browser()
    out_dir_tiles <- (file.path(in_dir,list_tiles_modis))
    list_files_by_tiles <- download_modis_obj$list_files_by_tiles #Use mapply to pass multiple arguments
    colnames(list_files_by_tiles) <- list_tiles_modis #note that the output of mapply is a matrix
  }else{
    
    ### Add download_dir here
    out_dir_tiles <- (file.path(in_dir,list_tiles_modis))
    #list_files_by_tiles <- mapply(1:length(out_dir_tiles),FUN=list.files,MoreArgs=list(pattern="*.hdf$",path=out_dir_tiles,full.names=T),SIMPLIFY=T) #Use mapply to pass multiple arguments 
    list_files_by_tiles <-mapply(1:length(out_dir_tiles),
                                 FUN=function(i,x){list.files(path=x[[i]],pattern="*.hdf$",full.names=T)},MoreArgs=(list(x=out_dir_tiles)),SIMPLIFY=T) #Use mapply to pass multiple arguments
    colnames(list_files_by_tiles) <- list_tiles_modis #note that the output of mapply is a matrix
  }
  
  ##### Add a check for missing here??: clean up and make a function:
  #undebug(extract_dates_from_raster_name)
  ##This assumes MODIS data!!!
  #browser()
  
  df_m_list_hdf <- extract_dates_from_raster_name(1,list_files_by_tiles,split_char=".")
  df_m_list_hdf <- lapply(1:nrow(list_files_by_tiles),
                          FUN=extract_dates_from_raster_name,
                          list_files=list_files_by_tiles,
                          split_char=".")
  df_m_list_hdf <- do.call(rbind,df_m_list_hdf)
  df_m_list_hdf <- data.frame(lapply(df_m_list_hdf , as.character), stringsAsFactors=FALSE)
  #View(df_m_list_hdf)
  #df_m_list_hdf <- as.data.frame(df_m_list_hdf,StringAsFactor=F)
  str(df_m_list_hdf)
  #names(df_m_list_hdf) <- c("raster_name","doy")
  names(df_m_list_hdf)
  df_m_list_hdf$doy
  #df_m_list_hdf$dates <- as.Date(strptime(as.character(df_m_list_hdf$doy), format="%Y %j"))
  
  df_dates <- as.data.frame(generate_dates_by_step(start_date,end_date,time_step))
  df_dates$doy <- as.character(df_dates$doy)
  df_dates$dates <- as.character(df_dates$dates)
  df_dates$missing <- 0
  #missing_dates <- setdiff(as.character(df_dates$doy),as.character(df_m_list_hdf$doy))
  missing_dates <- setdiff(as.character(df_dates$dates),as.character(df_m_list_hdf$dates))
  ### Assign 1 if a date is missing
  df_dates$missing[df_dates$date %in% missing_dates] <- 1
  class(df_dates$doy)
  class(df_m_list_hdf$doy)
  
  str(df_dates)
  df_hdf_products <- merge(df_dates,df_m_list_hdf,by="doy",all=T)
  #df_time_series <- merge(df_time_series,df_files,by="date",all=T) #outer join to keep missing dates
  df_hdf_products <- merge(df_dates,df_m_list_hdf,by="dates",all=T)
  
  number_missing_dates <- sum(df_hdf_products$missing)
  message(paste("Number of missing files is :",number_missing_dates))
  ## Still need to modify this for each tile!!!
  #browser()
  
  ####################################
  ##### STEP 2: IMPORT MODIS LAYERS ###
  
  ##Modify this section into a function to extract names of product and quality flag automatically!!
  infile_var <- list_files_by_tiles[[1]]
  GDALinfo_hdf<-GDALinfo(infile_var[1],returnScaleOffset=FALSE)
  str(GDALinfo_hdf)
  modis_subdataset <- attributes(GDALinfo_hdf)$subdsmdata
  print(modis_subdataset)
  
  hdf_df <- (strsplit(modis_subdataset,":"))
  hdf_df <- as.data.frame(do.call(rbind,hdf_df),stringsAsFactors=F)
  
  names(hdf_df) <- c("subdataset_name","description","dir","product","var_name")
  #Select automatically QC flag!!
  #View(hdf_df)
  
  write.table(hdf_df,"hdf_subdataset.txt",sep=",")
  
  browser()
  
  if(product_type=="NDVI"){
    #Will Need to change if EVI!!!
    #1 km 16 days EVI
    #modis_layer_str1 <- unlist(strsplit(modis_subdataset[3],"\""))[3] #Get day EVI layer
    #var_name_index <- which(hdf_df$var_name==var_name) #find matching variable in hdf file
    
    modis_layer_str1 <- unlist(strsplit(modis_subdataset[1],"\""))[3] #Get day NDVI layer
    modis_layer_str2 <- unlist(strsplit(modis_subdataset[5],"\""))[3] #Get day VI QC layer
  }
  
  #If product type EVI...
  #Implement here
  
  if(product_type=="LST"){
    #var_name, make a data.frame and pick the correct line...
    var_name_index <- which(hdf_df$var_name==var_name) #find matching variable in hdf file
    qc_name_index <- which(hdf_df$var_name==qc_name)
    #modis_layer_str1 <- unlist(strsplit(modis_subdataset[1],"\""))[3] #Get day LST layer
    #modis_layer_str2 <- unlist(strsplit(modis_subdataset[3],"\""))[3] #Get day QC layer
    modis_layer_str1 <- unlist(strsplit(modis_subdataset[var_name_index],"\""))[3] #Get day LST layer
    modis_layer_str2 <- unlist(strsplit(modis_subdataset[qc_name_index],"\""))[3] #Get day QC layer
    
  }

  if(product_type=="reflectance"){
    
    #nb_subdatasets <- length(modis_subdataset)/2
    #seq(1,nb_subdatasets,by=2)
    index_var_layers<- seq(1,length(modis_subdataset),by=2)
    index_qc_layers<- seq(2,length(modis_subdataset),by=2)
    
    modis_subdataset[index_var_layers] #13 bands for MOD09 reflectance
    
    #var_name_index <- index_var_layers # if get all...this inlcudes other things
    #qc_name_index <- index_qc_layers
    var_name_index <- index_var_layers[1:7] #get reflectance bands 1 to 7 
    qc_name_index <- 15 #qc band
    
    #modis_layer_str1 <- unlist(strsplit(modis_subdataset[1],"\""))[3] #Get day sur_refl_b01
    #modis_layer_str2 <- unlist(strsplit(modis_subdataset[1],"\""))[3] #Get day sur_refl_b01
    
    #modis_layer_str2 <- unlist(strsplit(modis_subdataset[5],"\""))[3] #Get day VI QC layer
    #var_name, make a data.frame and pick the correct line...
    #var_name_index <- which(hdf_df$var_name==var_name) #find matching variable in hdf file
    #qc_name_index <- which(hdf_df$var_name==qc_name)
    #modis_layer_str1 <- unlist(strsplit(modis_subdataset[1],"\""))[3] #Get day LST layer
    #modis_layer_str2 <- unlist(strsplit(modis_subdataset[3],"\""))[3] #Get day QC layer
    
    modis_layer_str1 <- unlist(lapply(var_name_index,function(i){unlist(strsplit(modis_subdataset[i],"\""))[3]}))
    modis_layer_str2 <- unlist(lapply(qc_name_index,function(i){unlist(strsplit(modis_subdataset[i],"\""))[3]}))
    
    #modis_layer_str1 <- unlist(strsplit(modis_subdataset[var_name_index],"\""))[3] #Get day LST layer
    #modis_layer_str2 <- unlist(strsplit(modis_subdataset[qc_name_index],"\""))[3] #Get day QC layer
    
  }
  ### import list of files before mosaicing
  
  file_format_import <- file_format
  #var_modis_name <- unlist(strsplit(modis_layer_str1,":"))[3]
  #qc_modis_name <- unlist(strsplit(modis_layer_str2,":"))[3]
  var_modis_name <- unlist(lapply(modis_layer_str1, function(x){unlist(strsplit(x,":"))[3]}))
  qc_modis_name <-  unlist(lapply(modis_layer_str2, function(x){unlist(strsplit(x,":"))[3]}))
  
  browser()
  
  #### Remote white space if present
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
      #browser()
      out_suffix_s <- var_modis_name
      list_param_import_modis <- list(i=1,
                                      hdf_file=infile_var,
                                      subdataset=modis_layer_str1,
                                      NA_flag_val=NA_flag_val,
                                      out_dir=out_dir_s,
                                      out_suffix=out_suffix_s,
                                      file_format=file_format_import,
                                      scaling_factors=scaling_factors,
                                      product_type=product_type,
                                      multiband=multiband)
      #debug(import_list_modis_layers_fun)
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
                                      hdf_file=infile_var,
                                      subdataset=modis_layer_str2,
                                      NA_flag_val=NA_flag_val,
                                      out_dir=out_dir_s,
                                      out_suffix=out_suffix_s,
                                      file_format=file_format_import,
                                      scaling_factors=NULL,
                                      product_type=product_type,
                                      multiband=multiband)
      #r1<-import_list_modis_layers_fun(1,list_param_import_modis)
      #r_qc_s <-mclapply(1:12,
      #                  FUN=import_list_modis_layers_fun,
      #                  list_param=list_param_import_modis,
      #                 mc.preschedule=FALSE,
      #                  mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
      
      list_r_qc_s <- mclapply(1:length(infile_var),
                              FUN=import_list_modis_layers_fun,
                              list_param=list_param_import_modis,
                              mc.preschedule=FALSE,
                              mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
      
      #### Report on errors when importing data
      
      df_import <- data.frame(var=unlist(list_r_var_s),qc=unlist(list_r_qc_s))
      list_error <- unlist(lapply(list_r_var_s, FUN=function(x){class(x)=="try-error"}))
      index_error <- which(list_error==TRUE)
      df_import$var[index_error] <- NA
      list_error <- unlist(lapply(list_r_qc_s, FUN=function(x){class(x)=="try-error"}))
      index_error <- which(list_error==TRUE)
      df_import$qc[index_error] <- NA
      #View(df_import)
      out_file <- paste0("df_import",".txt")
      write.table(df_import,file.path(out_dir_s,out_file))
      
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
  
  browser()
  ### Should report on errors:
  #import_error_list <- unlist(lapply(1:length(r_qc_s),FUN=function(x){class(x)=="try-error"}))
  #import_error_list[as.numeric(import_error_list)==1]
  #sum(as.numeric(import_error_list))
  
  #list_r_var_s <- unlist(list_r_var_s) #list of files as character vector
  #list_r_qc_s <- unlist(list_r_qc_s) #list of files as character vector
  list_r_var_s1 <- unlist(list_imported_files[[1]]$var[1]) #first tile, variable
  list_r_qc_s1 <- unlist(list_imported_files[[1]]$qc[1]) #first tile, quality flags
  
  if(length(var_modis_name) > 1){
    r_var_s1 <- brick(list_r_var_s1)
  }else{
    r_var_s1 <- raster(list_r_var_s1)
  }
  
  if(length(qc_modis_name) > 1){
    r_qc_s1 <- brick(list_r_qc_s1)
  }else{
    r_qc_s1 <- raster(list_r_qc_s1)
  }
  
  plot(r_var_s1)
  plot(r_qc_s1)
  
  #################################
  ##### STEP 3: APPLY/DEAL WITH QC FLAG AND SCREEN VALUE FOR VALID RANGE ###
  
  ## Get QC information for lST/NDVI and mask values: imporove and automate this later
  if(product_type=="NDVI"){
    #debug(create_MODIS_QC_table)
    QC_obj <- create_MODIS_QC_table(product_type) #Get table corresponding to QC for LST
    names(QC_obj)
    QC_data_ndvi <- QC_obj$NDVI
    #For NDVI: use this section to process. This is the default processing quality, this should go top in the parameters!!!
    #Select level 1:
    qc_lst_valid <- subset(x=QC_data_ndvi,QA_word1 == "VI Good Quality" | QA_word1 =="VI Produced,check QA")
    #Select level 2:
    qc_lst_valid <- subset(x=qc_lst_valid,QA_word2 %in% unique(QC_data_ndvi$QA_word2)[1:8]) #"Highest quality, 1","Lower quality, 2","Decreasing quality, 3",...,"Decreasing quality, 8" 
    
    #Select level 3:
    #...Not implemented at this stage
    ## Now select valid integer values
    if(selected_flags=="None"){
      #qc_lst_valid
      qc_valid <- QC_data_ndvi$Integer_Value
    }else{
      ##Now select the valid integer values:
      qc_valid <- qc_lst_valid$Integer_Value #valid integer values
      ## change the method above using similar to reflectance qc screening later on
    }
    
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
    
    qc_valid <- qc_lst_valid$Integer_Value #valid integer values
    
  }
  #QC_obj <- create_MODIS_QC_table(LST=TRUE, NDVI=TRUE) #Get table corresponding to QC for LST

  if(product_type=="reflectance"){
    
    #### This is where you set up the desired qc flags:
    desired_qc_rows <- c(1,2,5,14,23,32,41,50,59,69,71)
    
    ### Note this is using a different function than for LST and NDVI
    #debug(generate_qc_MODIS_reflectance_MOD09_table)
    qc_table_modis <- generate_qc_MODIS_reflectance_MOD09_table()
    
    #View(qc_table_modis)
    
    #### Now read in the values and process to match the selected flags!!!
    
    qc_table_modis_selected <- qc_table_modis[desired_qc_rows,]
    #View(qc_table_modis_selected)
    
    #debug(apply_mask_from_qc_layer)
    #### Generate function to create mask from qc
    #unique_bit_range <- unique(qc_table_modis$bitNo)
    
    #apply_mask_from_qc_layer(r_qc= r_qc_s1, qc_table_modis=qc_table_modis,out_dir=out_dir_s,out_suffix="")
    #
  }
  
  #r_lst_LST <- download_modis_obj$list_files_by_tiles[,j]
  list_r_var <- vector("list",length(list_tiles_modis)) #to contain image
  list_r_qc <- vector("list",length(list_tiles_modis)) #to contain qc mask image
  list_r_stack <- vector("list",length(list_tiles_modis)) #to contain results
  
  #26 minutes for 230 files to apply NDVI mask
  #12 minutes per tile for MOD09
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
      
      if(product_type=="reflectance"){
        out_suffix_s <- paste(product_type,sep="") #for MODIS product (var)
      }else{
        out_suffix_s <- paste(var_modis_name,sep="") #for MODIS product (var)
      }
      file_format_s <-file_format
      #undebug(create_raster_list_from_file_pat)
      list_r_var[[j]] <- create_raster_list_from_file_pat(out_suffix_s,file_pat="",
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
      
      ##Set output dir for QC mask operation
      out_dir_tmp <- paste0("mask_qc_",list_tiles_modis[j])
      #out_dir_s <- file.path(out_dir,list_tiles_modis[j])
      out_dir_s <- file.path(out_dir,out_dir_tmp) #input dir is import out dir
      
      if(product_type%in%c("NDVI","LST)")){
        list_param_screen_qc <- list(qc_valid,
                                     list_r_qc[[j]], 
                                     list_r_var[[j]],
                                     rast_mask=TRUE,
                                     NA_flag_val,out_dir_s,out_suffix) 
        names(list_param_screen_qc) <- c("qc_valid",
                                         "rast_qc", 
                                         "rast_var","rast_mask",
                                         "NA_flag_val","out_dir","out_suffix") 
        #undebug(screen_for_qc_valid_fun)
        test <- screen_for_qc_valid_fun(5,list_param=list_param_screen_qc)
        #r_stack[[j]] <- lapply(1:length(list_r_qc[[j]]),FUN=screen_for_qc_valid_fun,list_param=list_param_screen_qc)
        #r_test <-mclapply(1:11,FUN=screen_for_qc_valid_fun,list_param=list_param_screen_qc,mc.preschedule=FALSE,mc.cores = 11) #This is the end bracket from mclapply(...) statement
        
        list_r_stack[[j]] <-mclapply(1:length(list_r_qc[[j]]),
                                     FUN=screen_for_qc_valid_fun,
                                     list_param=list_param_screen_qc,
                                     mc.preschedule=FALSE,
                                     mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
        #r_stack <- stack(unlist(lapply(list_r_stack[[j]],function(x){x$var})))
        #plot(r_stack,y=5)
      }
      
      if(product_type=="reflectance"){
        i <- 1
        #out_suffix_s <- "masked"
        out_suffix_s <- out_suffix
        
        ## error in 37 for h09v06
        #undebug(apply_mask_from_qc_layer)
        list_r_stack[[j]] <- apply_mask_from_qc_layer(i,
                                                  rast_qc=list_r_qc[[j]],
                                                  rast_var=list_r_var[[j]],
                                                  qc_table_modis_selected,
                                                  NA_flag_val= NULL,
                                                  rast_mask=TRUE,
                                                  qc_info=F,
                                                  multiband=multiband,
                                                  out_dir=out_dir_s,
                                                  out_suffix=out_suffix_s)
        
        list_r_stack[[j]] <- mclapply(1:length(list_r_qc[[j]]),
                                  FUN= apply_mask_from_qc_layer,
                                  rast_qc=list_r_qc[[j]],
                                  rast_var=list_r_var[[j]],
                                  qc_table_modis_selected,
                                  NA_flag_val= NULL,
                                  rast_mask=TRUE,
                                  qc_info=F,
                                  multiband=multiband,
                                  out_dir=out_dir_s,
                                  out_suffix=out_suffix_s,
                                  mc.cores=num_cores,
                                  mc.preschedule=FALSE)
        
      }
      
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
  #browser()
  
  #################################
  ##### STEP 4: MOSAIC TILES  ###
  
  ##No mosaicing if only one tile!
  if(length(list_tiles_modis)==1){
    steps_to_run$mosaic <- FALSE
  }
  
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
      file_pattern <- paste0(".*.",product_type,".*.",
                             out_suffix,file_format,"$")
      
      ### Assume that the data is in the mask_qc folder
      
      in_dir_tmp <- paste0("mask_qc_",list_tiles_modis[j])
      in_dir_s <- file.path(out_dir,in_dir_tmp) #input dir is import out dir
      
      list_r_var_s <-list.files(pattern=file_pattern,
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
      #names(df_m_var) <- c(paste("raster_name",j,sep="_"),"date")
      names(df_m_var) <- c(paste("raster_name",j,sep="_"),"doy","date")
      
      df_m_var[,1] <- as.character(df_m_var[,1]) #make sure it is not factor
      #drop "doy column"?
      df_m_var <- df_m_var[,-2]
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
                            out_suffix="",
                            file_format,
                            NA_flag_val,
                            multiband)
    
    names(list_param_mosaic)<-c("j",
                                "mosaic_list",
                                "out_rastnames",
                                "out_dir",
                                "out_suffix",
                                "file_format",
                                "NA_flag_val",
                                "multiband")
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
    
    if(length(list_tiles_modis)==1){
      #dlook into import
      #out_dir_s
      #file_pattern <- paste0(".*.",product_type,".*.",
      #                       out_suffix,file_format,"$")
      
      ### Assume that the data is in the mask_qc folder
      
      in_dir_tmp <- paste0("mask_qc_",list_tiles_modis[j])
      
      in_dir_s <- file.path(out_dir,in_dir_tmp) #input dir is import out dir
      out_dir_s <- in_dir_s
      #list_r_var_s <-list.files(pattern=file_pattern,
      #                          path=in_dir_s,
      #                          full.names=TRUE) #inputs for moasics
      #list_m_var[[j]] <- list_r_var_s
      #list_m_var[[j]] <-create_raster_list_from_file_pat(out_suffix_s,file_pat="",in_dir=out_dir_s,out_prefix="",file_format=file_format_s)
      
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
      #method_proj_val <- "bilinear" 
      #method_proj_val <- "ngb" 
      
      ref_rast_prj <-projectRaster(from=ref_rast_tmp,
                                   res=res(ref_rast_tmp), #set resolution to the same as input
                                   crs=CRS_reg,
                                   method=method_proj_val)
      #to define a local reference system and reproject later!!
      #Assign new projection system here in the argument CRS_reg (!it is used later)
      if(!is.null(infile_reg_outline)){
        reg_sf <- st_read(infile_reg_outline)
        reg_sf <- st_transform(reg_sf,crs=CRS_reg)
        reg_sp <-as(reg_sf, "Spatial") 
        ref_rast <- crop(ref_rast_prj,reg_sp)  
        ref_rast_name_generated <- paste("ref_rast_",out_suffix,file_format,sep="")
        writeRaster(ref_rast,file.path(out_dir,ref_rast_name_generated))
      }
      
      ##This is the mosaiced and reproject tiles matching:
      ref_rast_prj_name_generated <- paste("ref_mosaiced_input_rast_",out_suffix,file_format,sep="")
      writeRaster( ref_rast_prj,file.path(out_dir,ref_rast_prj_name_generated))
      
    }  
      #Use the reference raster
    if(!is.null(ref_rast_name)){
        ref_rast<-raster(ref_rast_name) #This is the reference image used to define the study/processing area
        projection(ref_rast) <- CRS_reg #Assign given reference system from master script...
    }
      
    ##Create output names for region
    list_var_mosaiced_tmp <- remove_from_list_fun(list_var_mosaiced,condition_class ="try-error")$list
      
    out_suffix_var <-paste(out_suffix,file_format,sep="")          
    var_list_outnames <- change_names_file_list(list_var_mosaiced_tmp,
                                                  out_suffix_var,
                                                  "reg_",
                                                  file_format,
                                                  out_path=out_dir)     
      
    #list_param_create_region<-list(j,raster_name=list_var_mosaiced,reg_ref_rast=ref_rast,out_rast_name=var_list_outnames)
    #j<-1
    #list_param_create_region<-list(j,list_var_mosaiced,ref_rast,var_list_outnames,NA_flag_val)
    #names(list_param_create_region) <-c("j","raster_name","reg_ref_rast","out_rast_name","NA_flag_val")
      
    #undebug(create__m_raster_region)
      
    #### pasted from other file:
    out_rast_name <- NULL
    lf_r <- list_var_mosaiced
    list_param_create_region <- list(as.list(lf_r),
                                       ref_rast, 
                                       out_rast_name,
                                       agg_param,
                                       file_format,
                                       NA_flag_val,
                                       input_proj_str=NULL,
                                       multiband,
                                       method_proj_val,
                                       out_suffix="",
                                       out_dir_s)
    names(list_param_create_region) <- c("raster_name",
                                           "reg_ref_rast", 
                                           "out_rast_name",
                                           "agg_param",
                                           "file_format",
                                           "NA_flag_val",
                                           "input_proj_str",
                                           "multiband",
                                           "method_proj_val",
                                           "out_suffix",
                                           "out_dir")
    #undebug(create__m_raster_region)
    #r_filename <- create__m_raster_region(1,list_param=list_param_create_region)
    #r <- raster(r_filename)
    reg_var_list <- mclapply(1:length(lf_r),
                               FUN=create__m_raster_region,
                               list_param=list_param_create_region,
                               mc.preschedule=FALSE,
                               mc.cores = num_cores)
      
    ###Note: add x,y raster and mask defining the study area for the stack below!!
      
    lf <- list.files(path=out_dir_s,
                       pattern=paste0("*.",file_format,"$"),
                       full.names = T)
    write.table(lf,
                  "raster_list_data.txt",
                  sep=",",
                  row.names=F)
      
    if(save_textfile==TRUE){
        r_reg_var<- stack(reg_var_list)
        r_x <-init(r_reg_var,v="x")
        r_y <-init(r_reg_var,v="y")
        r_stack <- stack(r_x,r_y,r_reg_var)
        
        dat_reg_var_spdf <- as(r_stack,"SpatialPointsDataFrame")
        dat_reg_var <- as.data.frame(dat_reg_var_spdf) 
        #dat_out <- as.data.frame(r_reg_var)
        #dat_out <- na.omit(dat_out)
      
        out_filename <- paste0("dat_reg_var_list_",product_type,"_",out_suffix,".txt")
        out_filename <- file.path(out_dir_s,out_filename)
        write.table(dat_reg_var,
                    out_filename,
                    row.names=F,sep=",",col.names=T)
        
        #write.table(dat_reg_var)
    }else{
      out_filename <- NULL
    }
      
  }### End of step 5: reproject  
  
  ####### NOW prepare object to Return here:
  
  processing_modis_obj <- list(lf,out_filename)
  names(processing_modis_obj) <- c("project_files","out_filename")
  
  return(processing_modis_obj)
}


############################## END OF SCRIPT #################################

