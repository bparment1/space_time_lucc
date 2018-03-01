########################################  MODIS AND RASTER PROCESSING #######################################
########################################### Read, project, crop and process rasters #####################################
#This script contains general functions to processs raster images, raster time series as well MODIS specific functions.
#This script will form the basis of a library of functions for raster processing of for GIS and Remote Sensing applications.
#AUTHOR: Benoit Parmentier                                                                       
#CREATED ON: 09/16/2013
#MODIFIED ON: 03/01/2018
#PROJECT: None, general utility functions for raster (GIS) processing. 
#COMMIT: multiband option changes for mosaic of MOD09A1
#
#TODO:
#1)Modify generation of CRS for additional projected system (only LCC, Lambert Conformal at this stage)
#2)Add plotting function for raster stack
#3)Add procedures to fill in missing values using temporal and spatial interpolation for raster time series.
#4)Add proper documentation for general use.
#4)Test additional Quality Flag levels for Reflectand, Fire, LST,ALBEDO and other products
#7)Add function to report statistics: missing files in modis time series downloaded (make it general)
#
###################################################################################################

### List of 22 functions currently available:
# Add documentation later

#[1] "assign_projection_crs"            
#[2] "change_names_file_list"           
#[3] "create_MODIS_QC_table"            
#[4] "create__m_raster_region"         
#[5] "create_dir_fun"                   
#[6] "create_idrisi_rgf"                
#[7] "create_modis_tiles_region"        
#[8] "create_polygon_from_extent"      
#[9] "create_raster_list_from_file_pat" 
#[10] "define_crs_from_extent_fun"       
#[11] "extract_list_from_list_obj"
#[12] "extract_dates_from_raster_name"
#[13] "generate_dates_by_step"
#[14] "get_modis_tiles_list"             
#[15] "import_list_modis_layers_fun"     
#[16] "import_modis_layer_fun"           
#[17] "load_obj" 
#[18] "merge_multiple_df"
#[19] "modis_product_download"           
#[20] "mosaic_m_raster_list"             
#[21] "remove_from_list_fun"             
#[22] "screen_for_qc_valid_fun"         
#[23] "screening_val_r_stack_fun" 

   
###Loading R library and packages                                                      
#library(gtools)    # loading some useful tools 

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

## Function to mosaic modis or other raster images

mosaic_m_raster_list<-function(j,list_param){
  #This functions returns a subset of tiles from the modis grid.
  #Arguments: modies grid tile,list of tiles
  #Output: spatial grid data frame of the subset of tiles
  #Note that rasters are assumed to be in the same projection system!!
  
  #rast_list<-vector("list",length(mosaic_list))
  #for (i in 1:length(mosaic_list)){  
  # read the individual rasters into a list of RasterLayer objects
  # this may be changed so that it is not read in the memory!!!
  
  ## CREATED ON: 09/16/2013
  ## MODIFIED ON: 03/01/2018
  ## AUTHOR: Benoit Parmentier
  ##
  ## INPUTS
  # 1)
  # 2)
  # 3)
  ## OUTPUTS
  # 1)
  #
  
  #### Start function ####
  
  ### Step 1: Parse inputs
  
  #j<-list_param$j
  mosaic_list<-list_param$mosaic_list
  out_dir<-list_param$out_dir
  out_suffix <- list_param$out_suffix
  out_names<-list_param$out_rastnames
  file_format <- list_param$file_format
  NA_flag_val <- list_param$NA_flag_val
  multiband <- list_param$multiband
  
  
  ### Step 2: Mosaic lists of files
  
  ## Check input to see if this is a multiband file
  r_in <- brick(as.character(mosaic_list[[j]])[1])
  n_layer <- nlayers(r_in)
  
  if(n_layer>1){
    input.rasters <- lapply(as.character(mosaic_list[[j]]), brick)
  }else{
    input.rasters <- lapply(as.character(mosaic_list[[j]]), raster)
  }
  
  mosaiced_rast <- input.rasters[[1]]
  
  for (k in 2:length(input.rasters)){
    mosaiced_rast <- mosaic(mosaiced_rast,input.rasters[[k]], fun=mean)
    #mosaiced_rast<-mosaic(mosaiced_rast,raster(input.rasters[[k]]), fun=mean)
  }
  
  #### Step 3: Write out images and return values
  
  data_name<-paste("mosaiced_",sep="") #can add more later...
  #raster_name<-paste(data_name,out_names[j],".tif", sep="")
  #raster_name<-paste(data_name,out_names[j],file_format, sep="")
  raster_name<-paste(data_name,out_names[j], sep="") #don't add file format here
  
  #raster_name <- basename(sub(extension(rast_name_var),"",rast_name_var))
  #raster_name <- paste(raster_name,"_",out_suffix,extension(rast_name_var),sep="")
  
  #### change to compress if format is tif!! add this later...
  #file_format <- extension(mosaic_list[[j]])[1] #assumes that all inputs have the same file type, take first
  
  #Write out as brick
  data_type_str <- dataType(mosaiced_rast) #find the dataType, this should be a future input param
  
  if(is.null(NA_flag_val)){
    NA_flag_val <- NAvalue(mosaiced_rast)
  }
  
  #browser()

  if(n_layer>1){
    suffix_str <- 1:n_layer
    if(out_suffix!=""){
      suffix_str <- paste(out_suffix,suffix_str,sep="_")
    }
    #if not, don't add out_sufffix
    
    if(multiband==TRUE){
      #raster_name_tmp <- basename(rast_name_var)
      #raster_name <- basename(sub(file_format,"",raster_name))
      if(out_suffix!=""){
        raster_name_tmp <- paste(raster_name,"_",out_suffix,file_format,sep="")
      }else{
        raster_name_tmp <- paste(raster_name,file_format,sep="")
      }
      bylayer_val <- FALSE #don't write out separate layer files for each "band"
      rast_list <- file.path(out_dir,raster_name_tmp) #as return from function
    }
    if(multiband==FALSE){
      raster_name_tmp <- paste(raster_name,file_format,sep="") #don't add output suffix because in suffix_str
      bylayer_val <- TRUE #write out separate layer files for each "band"
      rast_list <- file.path(out_dir,(paste(raster_name,"_",suffix_str,file_format,sep=""))) 
    }
    
    if(file_format==".tif"){
      #Use compression option for tif
      writeRaster(mosaiced_rast,
                  filename=file.path(out_dir,raster_name_tmp),
                  bylayer=bylayer_val,
                  suffix=suffix_str,
                  overwrite=TRUE,
                  NAflag=NA_flag_val,
                  datatype=data_type_str,
                  options=c("COMPRESS=LZW"))
    }else{
      #Don't use compression option if not tif
      writeRaster(mosaiced_rast,
                  filename=file.path(out_dir,raster_name_tmp),
                  bylayer=multiband,
                  suffix=suffix_str,
                  overwrite=TRUE,
                  NAflag=NA_flag_val,
                  datatype=data_type_str)
    }
    
  }
  
  if(n_layer==1){
    #raster_name_tmp <- basename(rast_name_var)
    #raster_name <- basename(sub(extension(rast_name_var),"",rast_name_var))
    raster_name_tmp <- paste(raster_name,"_",out_suffix,file_format,sep="")
    writeRaster(mosaiced_rast, 
                NAflag=NA_flag_val,
                filename=file.path(out_dir,raster_name_tmp),
                bylayer=FALSE,
                bandorder="BSQ",
                datatype=data_type_str,
                overwrite=TRUE)
    rast_list <- file.path(out_dir,raster_name_tmp)
  }
  
  ## The Raster and rgdal packages write temporary files on the disk when memory is an issue. This can potential build up
  ## in long  loops and can fill up hard drives resulting in errors. The following  sections removes these files 
  ## as they are created in the loop. This code section  can be transformed into a "clean-up function later on
  ## Start remove
  tempfiles <- list.files(tempdir(),full.names=T) #GDAL transient files are not removed
  files_to_remove<-grep(out_suffix,tempfiles,value=T) #list files to remove
  if(length(files_to_remove)>0){
    file.remove(files_to_remove)
  }
  #now remove temp files from raster package located in rasterTmpDir
  removeTmpFiles(h=0) #did not work if h is not set to 0
  ## end of remove section
  
  return(rast_list)
}

## Function to reproject and crop modis tile or other raster images

### This is a very general function to process raster to match a raster of reference,
create__m_raster_region <-function(j,list_param){
  #This processes a list of raster to match to a region of interest defined by a reference raster
  #INPUT Arguments: raster name of the file,reference file with
  # j: file to be processed with input parameters
  # raster_name: list of raster to process i.e. match the region of interest
  # reg_ref_rast: reference raster used to defined the region of interest and spatial parameters
  # out_rast_name: output raster name, if NULL then use out_suffix to the input name to generate output name
  # agg_param: aggregation parameters: this is a vector used in in the aggregate function. It has three items:
  #                                     -TRUE/FALSE: if true then aggregate
  #                                     -agg_fact: aggregation factor, if NULL compute on the fly
  #                                     -agg_fun: aggregation function to use, the default is mean
  # file_format: output format used in the raster e.g. .tif, .rst
  # NA_flag_val: flag value used for no data
  # input_proj_str: defined projection,default null in which case it is extract from the input raster
  # out_suffix : output suffix added to output names if no output raster name is given
  # out_dir:  <- list_param$out_dir
  # Output: spatial grid data frame of the subset of tiles
  #
  # Authors: Benoit Parmentier
  # Created: 10/01/2015
  # Modified: 01/20/2016
  #TODO:
  # - Add option to disaggregate...
  # - Modify agg param to be able to use different ones by file j for the mcapply function
  #
  ################################################
  ## Parse input arguments
  raster_name <- list_param$raster_name[[j]] #list of raster ot project and crop, this is a list!!
  reg_ref_rast <- list_param$reg_ref_rast #This must have a coordinate system defined!!
  out_rast_name <- list_param$out_rast_name[j] #if NULL then use out_suffix to add to output name
  agg_param <- list_param$agg_param #TRUE,agg_fact,agg_fun
  file_format <- list_param$file_format #.tif, .rst
  NA_flag_val <- list_param$NA_flag_val #flag value used for no data
  input_proj_str <- list_param$input_proj_str #default null?
  out_suffix <- list_param$out_suffix
  out_dir <- list_param$out_dir
  
  ## Start #
  
  ## Create raster object if not already present
  if(class(raster_name)!="RasterLayer"){
    layer_rast<-raster(raster_name)
  }else{
    layer_rast <- raster_name
    raster_name <- filename(layer_rast)
  }
  
  ## Create output raster name if out_rast_name is null
  if(is.null(out_rast_name)){
    extension_str <- extension(raster_name)
    raster_name_tmp <- gsub(extension_str,"",basename(raster_name))
    if(out_suffix!=""){
      out_rast_name <- file.path(out_dir,paste(raster_name_tmp,"_crop_proj_reg_",out_suffix,file_format,sep="")) #for use in function later...
    }else{
      out_rast_name <- file.path(out_dir,paste(raster_name_tmp,"_crop_proj_reg",out_suffix,file_format,sep="")) #for use in function later...
    }
  }
  
  ## Get the input raster projection information if needed
  if(is.null(input_proj_str)){
    input_proj_str <-projection(layer_rast)   #Extract current coordinates reference system in PROJ4 format
  }else{
    projection(layer_rast) <- input_proj_str #assign projection info
  }
  region_temp_projected <- projectExtent(reg_ref_rast,CRS(input_proj_str))     #Project from ref to current region coord. system
  
  layer_crop_rast <- crop(layer_rast, region_temp_projected) #crop using the extent from the region tile
  #layer_projected_rast<-projectRaster(from=layer_crop_rast,crs=proj4string(reg_outline),method="ngb")
  if(agg_param[1]==TRUE){
    agg_fact <- as.numeric(agg_param[2]) #in case we have a string/char type
    agg_fun <- agg_param[3]
    #debug(aggregate_raster)
    r_agg_raster_name <- aggregate_raster(reg_ref_rast, #reference raster with the desired resolution
                                          agg_fact=agg_fact, #given aggregation factor
                                          r_in=layer_crop_rast, #raster to be aggregated
                                          agg_fun="mean", #aggregation function
                                          out_suffix=out_suffix,
                                          file_format=".tif",
                                          out_dir=out_dir)
    layer_crop_rast <- raster(r_agg_raster_name)
  }
  #Should check if different projection!!!
  layer_projected_rast <- projectRaster(from=layer_crop_rast,
                                        to=reg_ref_rast,
                                        method="ngb",
                                        NAflag=NA_flag_val,
                                        filename=out_rast_name,
                                        overwrite=TRUE)
  
  #NAvalue() #set above
  #Need cleanup of tmp files here!!! building up to 19gb!
  removeTmpFiles(h=0)
  
  return(out_rast_name)
}


#####

change_names_file_list<-function(list_name,out_suffix,out_prefix,extension,out_path=""){
  #Function to add suffix and prefix to list of file names
  lf_new_names_list<-vector("list",length(list_name)) #this will contain new names for files
  for (i in 1:length(list_name)){
    
    lf_name<-basename(list_name[[i]])
    lf_out_path<-dirname(list_name[[i]])
    data_name<-paste(out_prefix,sub(extension,"",lf_name),"_",sep="") #can add more later...
    raster_name<-paste(data_name,out_suffix, sep="") #out_suffix must include extension!!!
    if((lf_out_path!="") & (out_path=="")){
      lf_new_names_list[[i]]<-file.path(lf_out_path,raster_name)
    }else{
      lf_new_names_list[[i]]<-file.path(out_path,raster_name)
    }
    
  }
  return(unlist(lf_new_names_list))
}

screening_val_r_stack_fun<-function(list_val_range,r_stack){
  #Screening values for a raster stack by providing a valid range. Values outside the valid
  #range are assigned NA. Layers in the stack/brick are only screened if a name valid range is provided.
  #input: list_val_range: list of character strings comma separated
  #        e.g.: "mm_12,-15,50","mm_12,-15,50"
  #               variable name, min value, max value
  #The user must include the name of the variable matching the names in the raster brick/stack.
  #Values are assigned NA if they are less than the mini value or greater than the maximum value.
  #Output: stack with screened values. Note that the original order of layer names is not preserved!!!
  
  ## Parameters: parsing
  
  tab_range_list<-do.call(rbind,as.list(list_val_range))
  
  #tab_range <- strsplit(tab_range_list[[j]],",")
  
  tab_range <- strsplit(tab_range_list,",")
  tab_range <-as.data.frame(do.call(rbind, tab_range))
  names(tab_range)<-c("varname","vmin","vmax")
  tab_range$vmin <- as.numeric(as.character(tab_range$vmin)) #transform to character first to avoid values being considered as factor
  tab_range$vmax <- as.numeric(as.character(tab_range$vmax))
  tab_range$varname <- as.character(tab_range$varname)
  val_rst<-vector("list",nrow(tab_range)) #list of one row data.frame
  
  for (k in 1:nrow(tab_range)){
    #avl<-c(-Inf,tab_range$vmin[k],NA, tab_range$vmax[k],+Inf,NA)   #This creates a input vector...val 1 are -9999, 2 neg, 3 positive
    #avl<-c(tab_range$vmin[k],tab_range$vmax[k],NA)   #This creates a input vector...val 1 are -9999, 2 neg, 3 positive
    
    #rclmat<-matrix(avl,ncol=3,byrow=TRUE)
    #s_raster_r<-raster(r_stack,match(tab_range$varterm[k],names(r_stack))) #select relevant layer from stack
    s_raster_r<-raster(r_stack,match(tab_range$varname[k],names(r_stack)))
    #s_raster_r<-reclassify(s_raster_r,rclmat)  #now reclass values 
    #s_raster_r<-reclassify(s_raster_r,rclmat,include.lowest=TRUE,right=FALSE)  #now reclass values 
    #s_raster_r<-reclassify(s_raster_r,rclmat,include.lowest=FALSE,right=FALSE)  #now reclass values 
    #s_raster_r<-reclassify(s_raster_r,rclmat,include.lowest=TRUE,right=TRUE)  #now reclass values
    #s_raster_r<-reclassify(s_raster_r,rclmat,include.lowest=FALSE,right=TRUE)  #now reclass values
    #r_stack<-dropLayer(r_stack,match(tab_range$varname[k],names(r_stack)))
    s_raster_r[s_raster_r < tab_range$vmin[k]] <- NA #Assign NA if less than the minimum value in the valid range
    s_raster_r[s_raster_r > tab_range$vmax[k]] <- NA #Assign NA if greater than the maxim value in the valid range
    
    names(s_raster_r)<-tab_range$varname[k] #Loss of layer names when using reclass
    val_rst[[k]]<-s_raster_r
  }
  #could be taken out of function for parallelization
  s_rst_m<-stack(val_rst) #This a raster stack with valid range of values
  retained_names<-setdiff(names(r_stack),tab_range$varname)
  r_stack <- dropLayer(r_stack,match(tab_range$varname,names(r_stack)))
  names(r_stack) <-retained_names
  r_stack <- addLayer(r_stack,s_rst_m) #add back layers that were screened out
  
  return(r_stack)
}

define_crs_from_extent_fun<-function(reg_outline,buffer_dist){
  #Screening values for raster stack
  #input: list_val_range: list of character strings comma separated
  #        e.g.: "mm_12,-15,50","mm_12,-15,50"
  #               variable name, min value, max value
  library(rgeos)
  
  #Buffer function not in use yet!! need query for specific matching MODIS tile !!! use gIntersection
  if (buffer_dist!=0){
    reg_outline_dissolved <- gUnionCascaded(reg_outline)  #dissolve polygons
    reg_outline <- gBuffer(reg_outline_dissolved,width=buffer_dist*1000)
  }
  
  #CRS_interp <-"+proj=lcc +lat_1=43 +lat_2=45.5 +lat_0=41.75 +lon_0=-120.5 +x_0=400000 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs";
  reg_centroid <- gCentroid(reg_outline)
  reg_centroid_WGS84 <- spTransform(reg_centroid,CRS_locs_WGS84) #get cooddinates of center of region in lat, lon
  reg_outline_WGS84 <- spTransform(reg_outline,CRS_locs_WGS84) #get cooddinates of center of region in lat, lon
  reg_extent <-extent( reg_outline_WGS84) #get boudning box of extent
  #  xy_latlon<-project(xy, CRS_interp, inv=TRUE) # find lat long for projected coordinats (or pixels...)
  
  #Calculate projection parameters
  reg_lat_1 <- ymin(reg_extent)+((ymax(reg_extent)- ymin(reg_extent))/4)
  reg_lat_2 <- ymax(reg_extent)-((ymax(reg_extent)- ymin(reg_extent))/4)
  
  reg_lon_0 <- coordinates(reg_centroid_WGS84)[1]
  reg_lat_0 <- coordinates(reg_centroid_WGS84)[2]
  reg_x_0 <- 0
  reg_y_0 <- 0
  
  #Add false northing and false easting calucation for y_0,x_0
  #CRS_interp <- paste("+proj=lcc +lat_1=",43," +lat_2=",45.5," +lat_0=",41.75," +lon_0=",-120.5,
  #                    " +x_0=",0,"+y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  
  CRS_interp <- paste("+proj=lcc +lat_1=",reg_lat_1," +lat_2=",reg_lat_2," +lat_0=",reg_lat_0," +lon_0=",reg_lon_0,
                      " +x_0=",reg_x_0," +y_0=",reg_y_0," +ellps=WGS84 +datum=WGS84 +units=m +no_defs",sep="")
  
  reg_outline_interp <- spTransform(reg_outline,CRS(CRS_interp)) #get cooddinates of center of region in lat, lon
  
  #add part to save projected file??
  #return reg_outline!!!
  reg_outline_obj <-list(reg_outline_interp,CRS_interp)
  names(reg_outline_obj) <-c("reg_outline","CRS_interp")
  return(reg_outline_obj)
} 

### Assing projection system to raster layer
assign_projection_crs <-function(i,list_param){
  #assign projection to list of raster
  #proj_str: proj4 information
  #filename: raster file 
  proj_str<-list_param$proj_str
  list_filename<-list_param$list_filename
  
  filename <-list_filename[[i]]
  r<-raster(readGDAL(filename))
  projection(r)<-proj_str
  writeRaster(r,filename=filename,overwrite=TRUE)
}

## Function to  reclass value in 

#qc_valid_modis_fun <-function(qc_valid,rast_qc,rast_var,rast_mask=FALSE,NA_flag_val,out_dir=".",out_rast_name){
#  f_values <- as.data.frame(freq(rast_qc)) # frequency values in the raster...
#  f_values$qc_mask <- as.integer(f_values$value %in% qc_valid)
#  f_values$qc_mask[f_values$qc_mask==0] <- NA
#  
#  r_qc_m <- subs(x=rast_qc,y=f_values,by=1,which=3)
#  rast_var_m <-mask(rast_var,r_qc_m)
#  
#  if(rast_mask==FALSE){
#    raster_name<- out_rast_name
#   writeRaster(rast_var_m, NAflag=NA_flag_val,filename=file.path(out_dir,raster_name)
#                ,bylayer=FALSE,bandorder="BSQ",overwrite=TRUE)
#    return(rast_var_m)
#  }else{
#    r_stack <-stack(rast_var_m,r_qc_m)
#    raster_name<- out_rast_name
#    writeRaster(r_stack, NAflag=NA_flag_val,filename=file.path(out_dir,raster_name)
#                ,bylayer=FALSE,bandorder="BSQ",overwrite=TRUE)
#    return(r_stack)
#  }
#}

extract_list_from_list_obj<-function(obj_list,list_name){
  #Create a list of an object from a given list of object using a name prodived as input
  list_tmp<-vector("list",length(obj_list))
  for (i in 1:length(obj_list)){
    tmp<-obj_list[[i]][[list_name]] #double bracket to return data.frame
    list_tmp[[i]]<-tmp
  }
  return(list_tmp) #this is  a data.frame
}

remove_from_list_fun <- function(l_x,condition_class ="try-error"){
  index <- vector("list",length(l_x))
  for (i in 1:length(l_x)){
    if (inherits(l_x[[i]],condition_class)){
      index[[i]] <- FALSE #remove from list
    }else{
      index[[i]] <- TRUE
    }
  }
  l_x<-l_x[unlist(index)] #remove from list all elements using subset
  
  obj <- list(l_x,index)
  names(obj) <- c("list","valid")
  return(obj)
}

screen_for_qc_valid_fun <-function(i,list_param){
  ##Function to assign NA given qc flag values from MODIS or other raster
  #Author: Benoit Parmentier
  #Created On: 09/20/2013
  #Modified On: 10/09/2017
  
  #Parse arguments:
  
  qc_valid <- list_param$qc_valid # valid pixel values as dataframe
  rast_qc <- list_param$rast_qc[i] #raster with integer values reflecting quality flag layer e.g. day qc
  rast_var <- list_param$rast_var[i] #raster with measured/derived variable e.g. day LST
  rast_mask <- list_param$rast_mask #return raster mask as separate layer, if TRUE then raster is written and returned?
  NA_flag_val <- list_param$NA_flag_val #value for NA
  out_dir <- list_param$out_dir #output dir
  out_suffix <- out_suffix # suffix used for raster files written as outputs
  
  #########################
  #### Start script:
  
  if(!file.exists(out_dir)){
    dir.create(out_dir)
  }
  
  if(is.character(rast_qc)==TRUE){
    rast_name_qc <- rast_qc
    rast_qc <-raster(rast_qc)
  }
  if(is.character(rast_var)==TRUE){
    rast_name_var <- rast_var
    rast_var<-raster(rast_var)
  }
  
  f_values <- as.data.frame(freq(rast_qc)) # frequency values in the raster...as a dataframe
  f_values$qc_mask <- as.integer(f_values$value %in% qc_valid) # values that should be masked out
  f_values$qc_mask[f_values$qc_mask==0] <- NA #NA for masked out values
  
  #Use "subs" function to assign NA to values that are masked, column 1 contains the identifiers i.e. values in raster
  r_qc_m <- subs(x=rast_qc,y=f_values,by=1,which=3) #Use column labeled as qc_mask (number 3) to assign value
  rast_var_m <-mask(rast_var,r_qc_m)
  
  if(rast_mask==FALSE){  #then only write out variable that is masked out
    raster_name <-basename(sub(extension(rast_name_var),"",rast_name_var))
    raster_name<- paste(raster_name,"_",out_suffix,extension(rast_name_var),sep="")
    writeRaster(rast_var_m, NAflag=NA_flag_val,filename=file.path(out_dir,raster_name)
                ,bylayer=FALSE,bandorder="BSQ",overwrite=TRUE)
    
    ## The Raster and rgdal packages write temporary files on the disk when memory is an issue. This can potential build up
    ## in long  loops and can fill up hard drives resulting in errors. The following  sections removes these files 
    ## as they are created in the loop. This code section  can be transformed into a "clean-up function later on
    ## Start remove
    rm(rast_var)
    tempfiles<-list.files(tempdir(),full.names=T) #GDAL transient files are not removed, only RASTER, tempdir() from R basic
    files_to_remove<-grep(out_suffix,tempfiles,value=T) #list files to remove
    if(length(files_to_remove)>0){
      file.remove(files_to_remove)
    }
    #now remove temp files from raster package located in rasterTmpDir
    removeTmpFiles(h=0) #did not work if h is not set to 0
    ## end of remove section
    return(raster_name)
  }else{ #if keep mask true
    raster_name <-basename(sub(extension(rast_name_var),"",rast_name_var))
    raster_name<- paste(raster_name,"_",out_suffix,extension(rast_name_var),sep="")
    writeRaster(rast_var_m, NAflag=NA_flag_val,filename=file.path(out_dir,raster_name)
                ,bylayer=FALSE,bandorder="BSQ",overwrite=TRUE)
    raster_name_qc <-basename(sub(extension(rast_name_qc),"",rast_name_qc))
    raster_name_qc <- paste(raster_name_qc,"_","mask","_",out_suffix,extension(rast_name_qc),sep="")
    writeRaster(r_qc_m, NAflag=NA_flag_val,filename=file.path(out_dir,raster_name_qc)
                ,bylayer=FALSE,bandorder="BSQ",overwrite=TRUE)
    r_stack_name <- list(file.path(out_dir,raster_name),file.path(out_dir,raster_name_qc))
    names(r_stack_name) <- c("var","mask")
    
    ## The Raster and rgdal packages write temporary files on the disk when memory is an issue. This can potential build up
    ## in long  loops and can fill up hard drives resulting in errors. The following  sections removes these files 
    ## as they are created in the loop. This code section  can be transformed into a "clean-up function later on
    ## Start remove
    rm(rast_var_m)
    rm(r_qc_m)
    tempfiles<-list.files(tempdir(),full.names=T) #GDAL transient files are not removed
    files_to_remove<-grep(out_suffix,tempfiles,value=T) #list files to remove
    if(length(files_to_remove)>0){
      file.remove(files_to_remove)
    }
    #now remove temp files from raster package located in rasterTmpDir
    removeTmpFiles(h=0) #did not work if h is not set to 0
    ## end of remove section
    return(r_stack_name)
  }
}

create_raster_list_from_file_pat <- function(out_suffix_s,file_pat="",in_dir=".",out_prefix="",file_format=".rst"){
  #create a list of raster files to creater R raster stacks
  if(file_pat==""){
    file_pat <- paste(".*.",out_suffix_s,file_format,"$",sep="")
    #list_raster_name <- list.files(path=in_dir,pattern=,full.names=T)
  }
  
  list_raster_name <- list.files(path=in_dir,
                                 pattern=file_pat,
                                 full.names=T)
  #list_raster_name <- list.files(path=in_dir,
  #                              pattern=".*.LST_Day_1km.rst$",
  #                               full.names=T)
  
  dat_list<-c(mixedsort(unlist(list_raster_name)))
  #dat_list <- sub("[.][^.]*$", "", dat_list, perl=TRUE) 
  #writeLines(dat_list,con=paste(out_prefix,out_suffix_s,".rgf",sep=""))
  return(dat_list)
}

create_idrisi_rgf <- function(out_suffix_s,file_pat="",in_dir=".",out_prefix="",ending=FALSE){
  #create a list of raster idrisi  files...
  #out_suffix_s: ending pattern to wich .rst is attached
  #file_pat: string found in the list of files in what ever place
  if(file_pat==""){
    list_raster_name <- list.files(path=in_dir,pattern=paste(out_suffix_s,".*.rst$",sep=""),full.names=F)
  }else{
    list_raster_name <- list.files(path=in_dir,pattern=file_pat,full.names=F)
  }
  if(ending==TRUE){
    list_raster_name <- grep(paste(out_suffix_s,".rst$",sep=""),list_raster_name,value=TRUE)
  }
  dat_list<-c(as.integer(length(list_raster_name)),mixedsort(unlist(list_raster_name)))
  dat_list <- sub("[.][^.]*$", "", dat_list, perl=TRUE) #remove extension using regular expression
  writeLines(dat_list,con=file.path(in_dir,paste(out_prefix,out_suffix_s,".rgf",sep="")))
  
  return(file.path(in_dir,list_raster_name))
}

#This function is very very slow not be used most likely
create_polygon_from_extent<-function(reg_ref_rast,outDir=NULL,outSuffix=NULL){
  #This functions returns polygon sp from input rast
  #Arguments: input ref rast
  #Output: spatial polygon
  if(is.null(outDir)){
    outDir=getwd()
  }
  if(is.null(outSuffix)){
    outSuffix=""
  }
  set1f <- function(x){rep(1, x)}
  
  tmp_rast <- init(reg_ref_rast, fun=set1f, overwrite=TRUE)
  reg_outline_poly<-rasterToPolygons(tmp_rast,dissolve=T)
  infile_reg_outline <- paste("reg_out_line_",out_suffix,".shp",sep="")
  writeOGR(reg_outline_poly,dsn= outDir,layer= sub(".shp","",infile_reg_outline), 
           driver="ESRI Shapefile",overwrite_layer="TRUE")

  return(reg_outline_poly)
}

#merge a list of data.frame 
#function created since merge_all requires ordering from most to least 
merge_multiple_df<-function(df_list,by_name){
  for (i in 1:(length(df_list)-1)){
    if (i==1){
      df1=df_list[[i]]
    }
    if (i!=1){
      df1=df_m
    }
    df2<-df_list[[i+1]]
    df_m<-merge(df1,df2,by=by_name,all=T)
  }
  return(df_m)
}

### MODIS SPECIFIC FUNCTIONS

create_modis_tiles_region<-function(modis_grid,tiles){
  #This functions returns a subset of tiles from the modis grdi.
  #Arguments: modies grid tile,list of tiles
  #Output: spatial grid data frame of the subset of tiles
  
  h_list<-lapply(tiles,substr,start=2,stop=3) #passing multiple arguments
  v_list<-lapply(tiles,substr,start=5,stop=6) #passing multiple arguments
  
  selected_tiles<-subset(subset(modis_grid,subset = h %in% as.numeric (h_list) ),
                         subset = v %in% as.numeric(v_list)) 
  return(selected_tiles)
}

get_modis_tiles_list <-function(modis_grid,reg_outline,CRS_interp){
  #
  #Usage:This function finds the matching modis tiles given a vector polygon in shapefile format.
  #Inputs:
  #modis_grid: shapefiles of modis grid.
  #reg_outline: file name or sf object of processing region with extent
  #CRS_interp: projection system
  #Outputs:
  #list of modis tiles in the hxxvxx format eg h09v06
  
  ## SCRIPT BEGIN ##
  
  if((class(reg_outline)[1]!="sf")){
    reg_outline <- st_read(reg_outline)
  }
  
  modis_grid<-st_read(infile_modis_grid)    #Reading shape file using rgdal library
  reg_outline_sin <- st_transform(reg_outline,st_crs(modis_grid)$proj4string)

  l_poly <- st_intersects(reg_outline_sin,modis_grid) #intersected poly
  l_poly <- unique(unlist(l_poly))
  modis_grid_selected <- modis_grid[l_poly,]
  plot(modis_grid_selected$geometry)
  plot(reg_outline_sin,col="red",add=T)
  df_tmp <- as.data.frame(modis_grid_selected)
  tiles_modis <- paste(sprintf("h%02d", df_tmp$h),sprintf("v%02d", df_tmp$v),sep="")
  tiles_modis <- paste(tiles_modis,collapse=",")
  
  ##### Prepare object to return
  
  obj <- list(tiles_modis,modis_grid_selected)
  names(obj) <- c("tiles_modis","modis_grid_selected")
  return(obj)
}
## function to download modis product??

## For some time the ftp access does not work for MOLT!! now use curl and list from http.

modis_product_download <- function(MODIS_product,version,start_date,end_date,list_tiles,file_format,out_dir,temporal_granularity){
  
  ##Functions used in the script
  
  extractFolders=function(urlString) {
    htmlString=getURL(urlString)
    ret=gsub("]", "", str_replace_all(str_extract_all(htmlString, paste('DIR',".([^]]+).", '/\">',sep=""))[[1]], "[a-zA-Z\"= <>/]", ""))
    return(ret[which(nchar(ret)>0)])
  }
  
  #list_folders_files[[i]] <- extractFiles(url_folders_str[i], list_tiles)[file_format]
  extractFiles=function(urlString, list_tiles_str) {
    #Slight modifications by Benoit
    #list_tiles: modis tiles as character vectors eg c("h10v06","h09v07")
    #urlString: character vector with url folder to specific dates for product
    
    # get filename strings
    htmlString=getURL(urlString)
    #htmlString=getURL(urlString[2])
    allVec=gsub('\">', '', gsub('<a href=\"', "", str_extract_all(htmlString, paste('<a href=\"',"([^]]+)", '\">',sep=""))[[1]]))
    #allVec: this contains list of all files! need to select the correct tiles...
    #ret=c()
    #for (currSel in list_tiles_str) {
    #  ret=c(ret, grep(currSel, allVec, value=TRUE))
    #}
    #list_tiles_str <- c("h10v06","h11v07")
    ret <- lapply(list_tiles_str,function(x){grep(x,allVec,value=T)})
    
    # select specific files
    #ret <- paste(urlString,ret,sep="") #append the url of folder
    ret <-file.path(urlString,unlist(ret))
    jpg=sapply(ret, FUN=endswith, char=".jpg")
    xml=sapply(ret, FUN=endswith, char=".xml")
    hdf=sapply(ret, FUN=endswith, char=".hdf")
    
    retList=list(jpg=ret[which(jpg)], xml=ret[which(xml)], hdf=ret[which(hdf)])
    return(retList)
  }
  
  endswith=function(x, char) {
    currSub = substr(x, as.numeric(nchar(x)-nchar(char))+1,nchar(x))
    if (currSub==char) {return(TRUE)}
    return(FALSE)
  }
  
  ########## BEGIN SCRIPT #######
  
  #step 1: parse input elements
  
  st <- as.Date(start_date,format="%Y.%m.%d") #start date
  en <- as.Date(end_date,format="%Y.%m.%d") #end date
  ll <- seq.Date(st, en, by="1 day") #sequence of dates
  dates_queried <- format(ll,"%Y.%m.%d") #formatting queried dates
  
  #This is where it can be changed for other product
  url_product <-paste("https://e4ftl01.cr.usgs.gov/MOLT/",MODIS_product,"/",sep="") #URL is a constant...
  #url_product <- file.path("http://e4ftl01.cr.usgs.gov/MOLT/",MODIS_product)
  #debug(extractFolders)
  dates_available <- extractFolders(url_product)  #Get the list of available driectory dates for the product, from 2000 to now
  
  list_folder_dates <- intersect(as.character(dates_queried), as.character(dates_available)) #list of remote folders to access
  #list_folder_dates <-setdiff(as.character(dates_available), as.character(dates_queried))
  
  #step 2: list content of specific day folder to obtain specific file...  #parse by tile name!!!

  url_folders_str <-paste(url_product,list_folder_dates,"/",sep="") #url for the folders matching dates to download
  
  ## loop over
  #debug(extractFiles)
  #generate outdir for each tile!!!
  
  list_folders_files <- vector("list",length(url_folders_str))
  d_files <- vector("list",length(list_folders_files))
  file_format<-c("hdf","xml")
  #can make this faster using parallelization...
  for (i in 1:length(url_folders_str)){
    #debug(extractFiles)
    list_folders_files[[i]] <- try(extractFiles(url_folders_str[i], list_tiles)[file_format]) 
    #d_files[[i]] <- list_folders_files[[i]][[file_format]]                      
  }
  #list_folders_files <- lapply(url_folders_str,extractFiles,list_tiles_str=list_tiles)
  #Now remove error objects...
  d_files_tmp <-remove_from_list_fun(l_x=list_folders_files,condition_class ="try-error")
  d_files <- as.character(unlist(d_files_tmp$list)) #all the files to download...
  n_dir_available <- length(d_files)
  n_dir_requested <- length(list_folders_files)

  #Step 3: download file to the directory 
  #browser()
  #prepare files and directories for download
  out_dir_tiles <- file.path(out_dir,list_tiles)
  list_files_tiles <- vector("list",length(list_tiles))
  for(j in 1:length(out_dir_tiles)){
    if (!file.exists(out_dir_tiles[j])){
      dir.create(out_dir_tiles[j])
    }
    list_files_tiles[[j]] <- grep(pattern=list_tiles[j],x=d_files,value=TRUE) 
  }
    
  #Now download per tiles: can be parallelized
  for (j in 1:length(list_files_tiles)){ #loop around tils
    file_items <- list_files_tiles[[j]]
    for (i in 1:length(file_items)){
      file_item <- file_items[i]
      #download.file(file_item,destfile=file.path(out_dir_tiles[j],basename(file_item)))
      #download.file(file_item,destfile="test.hdf")
      ## need a .netrc file,
      ##set the file to at least Read (400) or Read/Write (600)
      ##chmod 0600 ~/.netrc
      #curl -n -L -c cookiefile -b cookiefile http://e4ftl01.cr.usgs.gov/MOLT/MOD09A1.006/2001.01.09/MOD09A1.A2001009.h13v01.006.2015140120258.hdf.xml 
      #system("curl -n -L -c cookiefile -b cookiefile http://e4ftl01.cr.usgs.gov/MOLT/MOD09A1.006/2001.01.09/MOD09A1.A2001009.h13v01.006.2015140120258.hdf.xml") 
      
      #system("curl -n -L -c cookiefile -b cookiefile https://e4ftl01.cr.usgs.gov/MOLT/MOD09A1.006/2001.01.09/MOD09A1.A2001009.h13v01.006.2015140120258.hdf --output MOD09A1.A2001009.h13v01.006.2015140120258.hdf")
      cmd_curl_str <- paste("curl -n -L -c cookiefile -b cookiefile",
                            file_item,
                            "--output",
                            paste0("'",file.path(out_dir_tiles[j],basename(file_item)),"'")
                            ) 
      system(cmd_curl_str)
      #system("curl -n -L -c cookiefile -b cookiefile https://e4ftl01.cr.usgs.gov/MOLT/MOD11A1.006/2001.01.01/MOD11A1.A2001001.h08v05.006.2015111170727.hdf")
      #myopts <- RCurl::curlOptions(netrc=TRUE, netrc.file=path.expand("~/.netrc"), 
      #                             cookiefile=path.expand("~/.urs_cookies"), 
      #                             followlocation=TRUE) 
      
      #your_url <- file_item
      #getBinaryURL(your_url, .opts=myopts) 
      
      #where the .netrc file should look like this: 
        
      #machine urs.earthdata.nasa.gov login your_login password your_pass 
      #machine e4ftl01.cr.usgs.gov login your_login password your_pass
    }
  }
  
  #Prepare return object: list of files downloaded with http and list downloaded of files in tiles directories
  
  list_files_by_tiles <-mapply(1:length(out_dir_tiles),
                               FUN=function(i,x){list.files(path=x[[i]],pattern="*.hdf$",full.names=T)},MoreArgs=(list(x=out_dir_tiles))) #Use mapply to pass multiple arguments
  #list_files_by_tiles <-mapply(1:length(out_dir_tiles),FUN=list.files,MoreArgs=list(pattern="*.hdf$",path=out_dir_tiles,full.names=T)) #Use mapply to pass multiple arguments
  
  colnames(list_files_by_tiles) <- list_tiles #note that the output of mapply is a matrix
  download_modis_obj <- list(list_files_tiles,list_files_by_tiles)
  names(download_modis_obj) <- c("downloaded_files","list_files_by_tiles")
  return(download_modis_obj)
}

#######
## function to import modis in tif or other format...
import_modis_layer_fun <-function(hdf_file,subdataset,NA_flag,out_rast_name="test.tif",memory=TRUE){
  
  #PARSE input arguments/parameters
  
  modis_subset_layer_Day <- paste("HDF4_EOS:EOS_GRID:",hdf_file,subdataset,sep="")
  r <-readGDAL(modis_subset_layer_Day)
  r  <-raster(r)
  
  if(memory==TRUE){
    return(r)
  }else{
    #Finish this part...write out
    raster_name<- out_rast_name
    writeRaster(r_spat, NAflag=NA_flag_val,filename=raster_name,bylayer=TRUE,bandorder="BSQ",overwrite=TRUE)       
    return(raster_name)
  }  
}

## function to import modis in tif or other format...
import_list_modis_layers_fun <-function(i,list_param){
  
  #PARSE input arguments/parameters
  
  hdf_file <- list_param$hdf_file
  subdataset <- list_param$subdataset
  NA_flag_val <- list_param$NA_flag_val
  out_dir_s <- list_param$out_dir
  out_suffix <- list_param$out_suffix
  file_format <- list_param$file_format
  scaling_factors <- list_param$scaling_factors
  product_type <- list_param$product_type #this is the variable name: LST, NDVI, reflectance
  multiband <- list_param$multiband # 
    
  ######## Begin script #####
  
  if(!file.exists(out_dir_s)){
    dir.create(out_dir_s)
  }
  
  #setwd(out_dir_s)

  #Now get file to import
  hdf_filename <-hdf_file[i] # must include input path!!
  
  modis_subset_layer_Day <- paste("HDF4_EOS:EOS_GRID:",
                                  hdf_filename,
                                  subdataset,sep="")
  
  n_layer <- length(modis_subset_layer_Day)

  if(n_layer==1){
    r <- readGDAL(modis_subset_layer_Day) 
    r  <-raster(r)
  }
  if(n_layer>1){
    list_r <- lapply(modis_subset_layer_Day,function(x){raster(readGDAL(x))})
    #r <- readGDAL(modis_subset_layer_Day) 
    r <- stack(list_r)
    get_names_layers <- function(x){val_extracted <- unlist(strsplit(x,":")); val_extracted[length(val_extracted)]}
    names_layers <- unlist(lapply(modis_subset_layer_Day,get_names_layers))
    names(r) <- names_layers
  }
  
  if(!is.null(scaling_factors)){ #if scaling factor exists, scale values...(not applied for QC flags!!!)
    r <- scaling_factors[1]*r + scaling_factors[2]
  }
  #Finish this part...write out
  names_hdf <- as.character(unlist(strsplit(x=basename(hdf_filename), split="[.]")))
  
  char_nb<-length(names_hdf)-2
  names_hdf <- names_hdf[1:char_nb]
  names_hdf <- paste(names_hdf,collapse="_") #this is the name of the hdf file with "." replaced by "_"
  raster_name <- paste(names_hdf,"_",out_suffix,file_format,sep="")
  #out_dir_str <-  dirname(hdf)
  #set output dir from input above
  if(n_layer==1){
    raster_name_tmp <- raster_name
    if(file_format==".tif"){
      
      writeRaster(r, 
                  NAflag=NA_flag_val,
                  filename=file.path(out_dir_s,raster_name_tmp),
                  bylayer=TRUE,
                  bandorder="BSQ",
                  overwrite=TRUE,
                  #datatype=data_type_str, #this should be a future option for reduced size!!!
                  options=c("COMPRESS=LZW")) #compress by default
    }else{
      
      writeRaster(r, 
                  NAflag=NA_flag_val,
                  filename=file.path(out_dir_s,raster_name_tmp),
                  bylayer=TRUE,
                  bandorder="BSQ",
                  overwrite=TRUE)   
    }
  }
  
  #### Now deal with multiband e.g. MOD09 reflectance
  if(n_layer>1){
    
    #Write out as brick
    data_type_str <- dataType(r) #find the dataType, this should be a future input param
    if(is.null(NA_flag_val)){
      NA_flag_val <- NAvalue(r)
    }
    
    if(multiband==TRUE){
      raster_name_tmp <- paste(names_hdf,"_",product_type,file_format,sep="")
      bylayer_val <- FALSE #don't write out separate layer files for each "band"
    }
    if(multiband==FALSE){
      raster_name_tmp <- raster_name
      bylayer_val <- TRUE #write out separate layer files for each "band"
    }
    
    if(file_format==".tif"){
      writeRaster(r,
                  filename=file.path(out_dir_s,raster_name_tmp),
                  bylayer=bylayer_val,
                  #suffix=paste(names(r),"_",out_suffix,sep=""),
                  #format=format_raster,
                  suffix=paste(names(r)),
                  overwrite=TRUE,
                  NAflag=NA_flag_val,
                  datatype=data_type_str,
                  options=c("COMPRESS=LZW"))

    }else{
    #Don't use compression option if not tif
    writeRaster(r,
                filename=file.path(out_dir_s,raster_name_tmp),
                bylayer=multiband,
                #suffix=paste(names(r),"_",out_suffix,sep=""),
                #format=format_raster,
                suffix=paste(names(r)),
                overwrite=TRUE,
                NAflag=NA_flag_val,
                datatype=data_type_str)
      
    }
    
  }
  
  ## The Raster and rgdal packages write temporary files on the disk when memory is an issue. This can potential build up
  ## in long  loops and can fill up hard drives resulting in errors. The following  sections removes these files 
  ## as they are created in the loop. This code section  can be transformed into a "clean-up function later on
  ## Start remove
  rm(r)
  product <- names_hdf[1]
  tempfiles<-list.files(tempdir(),full.names=T) #GDAL transient files are not removed
  files_to_remove<-grep(product,tempfiles,value=T) #list files to remove
  if(length(files_to_remove)>0){
    file.remove(files_to_remove)
  }
  #now remove temp files from raster package located in rasterTmpDir
  removeTmpFiles(h=0) #did not work if h is not set to 0
  #tempfiles<-list.files(tempdir(),full.names=T) #GDAL transient files are not removed
  ## end of remove section
  
  return(file.path(out_dir_s,raster_name_tmp)) 
}

create_MODIS_QC_table <-function(LST=TRUE, NDVI=TRUE,reflectance=TRUE){
  #Function to generate MODIS QC  flag table
  #Author: Benoit Parmentier (with some lines from S.Mosher)
  #Date CREATED: 09/16/2013
  #Date MODIFIED: 02/15/2018
  #
  #Some of the inspiration and code originates from Steve Mosher' s blog:
  #http://stevemosher.wordpress.com/2012/12/05/modis-qc-bits/
  
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
    test <- unique(r_qc_s1)
    intToBits(65)    
    length(intToBits(65))
    rawToBits()
    line1<-c(readBin(to.read,"int",5), 
             readBin(to.read,"double",1,size=4),
             readBin(to.read,"int",2))
    as.integer(test[1]) #no need to do this
    #length(intToBits(test[1])) 
    length(intToBits(test[1])) #that is 32
    #intToBits()
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

    list_QC_Data[[3]]<- QC_Data
  }
  
  ## PRODUCT 4: Albedo
  #This can be seen from table defined at LPDAAC: https://lpdaac.usgs.gov/products/modis_products_table/mod11a2
  #To be added...
  
  ###Now return and save object:
  #Prepare object to return
  
  save(list_QC_Data,file= file.path(".",paste("list_QC_Data",".RData",sep="")))
  
  return(list_QC_Data)
}

#This function extract dates from MODIS exported raster files, used in the mosaics stage... 
extract_dates_from_raster_name <- function(i,list_files,split_char="_"){
  #Prepare list of modis tiles to mosaic
  raster_name <- list_files[i]
  #raster_name <- list_m_var[[1]][1]
  split_char_val <- paste0("[",split_char,"]")
  
  #names_part_raster <- as.character(unlist(strsplit(x=basename(raster_name), split="[_]")))
  #names_part_raster <- as.character(unlist(strsplit(x=basename(raster_name), split="[.]")))
  names_part_raster <- as.character(unlist(strsplit(x=basename(raster_name), split=split_char_val)))
  
  char_nb<-length(names_part_raster)-2
  #names_hdf <- names_hdf[1:char_nb]

  doy_raster <- names_part_raster[2] #this is MODIS DOY (day of year)
  doy_raster <- strsplit(doy_raster,"A")[[1]][[2]]
  dates_val <- as.Date(strptime(as.character(doy_raster), format="%Y %j"))
  df_raster_name <- data.frame(raster_name=raster_name,doy=doy_raster,dates=dates_val)
  return(df_raster_name)
}


#Screen data: use only : # level 1: LST Produced good quality, LST Produced other Quality Check QA, 
                         # level 2: good data , Other quality 

#classify_raster_fun <- function(list_rast){
#  library(raster)
#  #raster_data <- list.files(path=getwd())    #promt user for dir containing raster files
#  rast_s <- stack(list_rast)
#  f <- function(x) { rowSums(x >= 4 & x <= 9) }
#  x <- calc(rast_s, f, progress='text', filename='output.tif')
#}

generate_dates_by_step <-function(start_date,end_date,step_date){
  #library(xts) declare out of this function
  #library(zoo)
  #library(lubridate)
  
  st <- as.Date(start_date,format="%Y.%m.%d")
  en <- as.Date(end_date,format="%Y.%m.%d")
  #year_list <-seq(format(st,"%Y"),format(en,"%Y")) #extract year
  year_list <- seq(as.numeric(strftime(st,"%Y")),as.numeric(strftime(en,"%Y"))) #extract year
  
  ll_list <- vector("list",length=length(year_list))
  for (i in 1:length(year_list)){
    if(i==1){
      first_date <-st
    }else{
      first_date<-paste(year_list[[i]],"-01","-01",sep="")
    }
    if(i==length(year_list)){
      last_date <-en
    }else{
      last_date<-paste(year_list[[i]],"-12","-31",sep="")
    }
    #ll <- seq.Date(st, en, by=step)
    ll <- seq.Date(as.Date(first_date), as.Date(last_date), by=step_date)
    ll_list[[i]]<-as.character(ll)
    #paste(yday(ll,)
  }
  
  #
  dates_modis <-as.Date(unlist((ll_list))) 
  
  dates_DOY_modis <- as.character(paste(year(dates_modis),sprintf("%03d", yday(dates_modis)),sep=""))
  dates_obj <- list(dates_modis,dates_DOY_modis)
  names(dates_obj) <- c("dates","doy")  
  return(dates_obj)
}
 ################################## END OF SCRIPT #######################################

