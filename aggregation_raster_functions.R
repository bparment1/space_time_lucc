####################################  Aggregation functions for general use  #######################################
###########################################  Processing and Aggregating raster  #######################################
#This script contains genreal functions for aggregation of categorical or continuous raster layers.
#Functions are to be called from other scripts. 
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/03/2015 
#DATE MODIFIED: 06/01/2018
#Version: 1
#PROJECT: General utility functions for raster aggregation.            

#
#COMMENTS: -  
#          - 
#TO DO:
# - Add function to reclassify inputs from aggregation
# - 
# - 
#
#################################################################################################

###Loading R library and packages                                                      

library(raster)                 # loading the raster package
library(gtools)                 # loading R helper programming tools/functions
library(sp)                     # spatial objects in R
library(gplots)                 # plotting functions such as plotCI
library(rgdal)                  # gdal driver for R
library(RColorBrewer)           # color scheme, palettes used for plotting
library(gdata)                  # read different format (including .xlsx)
library(plotrix)                # plot options and functions 
library(rasterVis)              # raster visualization
library(colorRamps)             # contains matlab.like palette
library(zoo)                    # time series objects and methods
library(maptools)               #
library(rgeos)                  # spatial analysis, topological and geometric operations e.g. interesect, union, contain etc.
library(parallel)               # mclapply with cores...

##### Functions used in this script 

aggregate_raster_fun <- function(l_rast,cat_names,agg_method_cat="majority",agg_fact=2,agg_fun=mean,file_format=".tif",rast_ref=NULL,num_cores=1,out_suffix=NULL, out_dir=NULL){
  #
  #Aggregate raster from raster input and reference file
  #INPUT arguments:
  #1) l_rast : set of input raster layers as list, can be filenames or raster layers objects
  #2) cat_names: within the list, give names of raster that contain categorical variables  
  #3) agg_method_cat: aggregation rule for categorical variable, default is majority
  #4) agg_fact: factor to aggregate
  #5) agg_fun: default is mean
  #6) file_format: e.g. ".tif"
  #7) rast_ref: reference raster to match in resolution, if NULL then send a message
  #6) file_Format: raster format used e.g. .tif
  #5) num_cores: number of cores to use:
  #5) out_suffix: output suffix
  #7) out_dir: output directory
  #8) out_rast_name: output raster name if null it is created from the input file name
  #OUTPUT:
  # out_raster_name: name of the file containing the aggregated raster
  #
  # Authors: Benoit Parmentier
  # Created: 03/02/2017
  # Modified: 06/01/2018
  # To Do: 
  # - Add option to disaggregate
  # - add additional options to do aggregation
  ################################
  
  
  #
  #Function to aggregate input raster stack
  #if use majority then the zonal layer is aggregated and then reclassfied based by the majority rule
  
  if(!is.null(cat_names)==TRUE){
    selected_continuous_layers <- names(l_rast)!=cat_names #use set it will be cleaner?
    #continuous_layers <- l_rast[selected_continuous_layers]
    selected_cat_layers <- names(l_rast)==cat_names
  }else{
    selected_continuous_layers <- rep(1,length=length(l_rast))
  }
  
  
  #if(!is.null(l_rast)){
  if(sum(selected_continuous_layers)>0){
    
    #debug(aggregate_raster)
    #lf_agg_test <- aggregate_raster(l_rast[[1]],
    #                   #r_in=raster(lf_layerized_bool[1]),
    #                   agg_fact=agg_fact,
    #                   reg_ref_rast=NULL,
    #                   #agg_fun="mean",
    #                   agg_fun=agg_fun,
    #                   out_suffix=NULL,
    #                   file_format=file_format,
    #                   out_dir=out_dir,
    #                   out_rast_name = NULL) 
    
    lf_agg <- mclapply(l_rast[selected_continuous_layers],
                       FUN=aggregate_raster,
                       #r_in=raster(lf_layerized_bool[1]),
                       agg_fact=agg_fact,
                       reg_ref_rast=rast_ref,
                       #agg_fun="mean",
                       agg_fun=agg_fun,
                       out_suffix=NULL,
                       file_format=file_format,
                       out_dir=out_dir,
                       out_rast_name = NULL,
                       mc.preschedule=FALSE,
                       mc.cores = num_cores) 
    
    #l_rast_original <- l_rast
    #l_rast <- unlist(lf_agg) 
    l_rast_continuous <- lf_agg
  }else{
    l_rast_continuous <- NULL
  }
  
  ###Break out and get mean per class and do majority rule!
  
  #if(use_majority==TRUE){
  #if(!is.null(cat_names)==TRUE){
  if(sum(selected_cat_layers)>0){
    
    ## Use loop because we already have a num_cores
    l_rast_cat <- vector("list",length=length(selected_cat_layers))
    l_rast_cat <- l_rast[selected_cat_layers]
    for(j in 1:length(l_rast_cat)){
      
      #debug(generate_soft_cat_aggregated_raster_fun)
      #out_suffix_str <- out_suffix #may change this to included "majority" in the name
      out_suffix_str <- paste(cat_names[j],"_",out_suffix,sep="")
      raster_name_cat <- l_rast[selected_cat_layers][[j]]
      
      lf_agg_soft <- generate_soft_cat_aggregated_raster_fun(raster_name_cat,
                                                             reg_ref_rast=rast_ref,
                                                             agg_fact,
                                                             agg_fun,
                                                             num_cores,
                                                             NA_flag_val=NA_flag_val,
                                                             file_format,
                                                             out_dir,
                                                             out_suffix_str)
      
      if(class(raster_name_cat)=="character"){
        raster_name_cat <- raster(raster_name_cat)
      }
      
      reclass_val <- unique(raster_name_cat) #unique zonal values to reassign
      
      r_stack <- stack(lf_agg_soft)
      
      if(agg_method_cat=="majority"){
        r_reclass_obj <- reclass_in_majority(r_stack= r_stack,
                                             threshold_val=NULL,
                                             max_aggregation = TRUE,
                                             reclass_val = reclass_val)
        plot(r_reclass_obj$r_rec)
        rast_zonal <- r_reclass_obj$r_rec
      }
      
      #Other methods should bee added here here
      
      ### Compute agg factor
      if(is.null(agg_fact)){
        if(is.character(rast_ref)){
          rast_ref <- raster(rast_ref)
        }
        res_ref <- res(rast_ref) #assumes square cells, and decimal degrees from WGS84 for now...
        res_in <- res(raster_name_cat) #input resolution, assumes decimal degrees
        agg_fact_val <- unique(round(res_ref/res_in)) #find the factor needed..
        
        #fix this to add other otpions e.g. aggregating down
      }
      
      #output aggregated categorical layer:
      raster_name <- paste0("agg_",agg_fact_val,"_","r_",out_suffix_str,file_format)
      
      #### Need to add multiband support
      writeRaster(rast_zonal,
                  filename=file.path(out_dir,raster_name),
                  overwrite=TRUE)  
      
      #l_rast_cat[[i]] <- rast_zonal 
      l_rast_cat[[j]] <- raster_name 
      
    }
    
  }else{
    l_rast_cat <- NULL
  }
  
  ###
  #zonal_colnames <- gsub(extension(raster_name),"",raster_name)
  ##
  
  ##########################
  #### prepare return object
  
  obj_agg <- list(cat_names,l_rast_cat,l_rast_continuous)
  names(obj_agg) <- c("cat_names","l_rast_cat","l_rast_continuous")
  
  return(obj_agg)
}


reclass_in_majority <- function(r_stack,threshold_val=0.5,max_aggregation=FALSE,reclass_val){
  ##
  #This function reclassify a set of soft layers using the majority or maximum value rule.
  #When max_aggregation is TRUE, the max value rule is used in the aggregation.
  #
  #INPUTS
  #1) r_stack
  #2) threshold_val
  #3) max_aggregation
  #4) reclass_val
  #
  
  ## Reclass
  if(!is.null(threshold_val) & (max_aggregation==FALSE)){
    r_rec_threshold <- r_stack > threshold_val
    #use calc later to save directly the file
    #
    
    r_rec_val_s <- lapply(1:nlayers(r_rec_threshold),
                          function(i,r_stack){df_subs <- data.frame(id=c(0,1),v=c(0,reclass_val[i]));
                          x <- subs(subset(r_stack,i), df_subs)},r_stack=r_rec_threshold)
    r_rec_val_s <- stack(r_rec_val_s) #this contains pixel above 0.5 with re-assigned values
    r_rec <- calc(r_test,function(x){sum(x)})
    
    ### prepare return object
    reclass_obj <- list(r_rec,r_rec_val_s)
    names(reclass_obj) <- c("r_rec","r_rec_val_s")
    
  }
  
  if(max_aggregation==TRUE){
    #r_zonal_agg_soft <- stack(lf_agg_soft)
    #Find the max, in stack of pixels (can be used for maximum compositing)
    r_max_s <- calc(r_stack, function(x) max(x, na.rm = TRUE))
    #maxStack <- stackApply(r_zonal_agg_soft, indices=rep(1,nlayers(r_zonal_agg_soft)), fun = max, na.rm=TRUE)
    r_max_rec_s <- overlay(r_stack,r_max_s, fun=function(x,y){as.numeric(x==y)})
    r_ties <- sum(r_max_rec_s) #find out ties
    #this may be long
    #freq_r_rec_df <- freq(r_rec_max,merge=T)
    r_ties_mask <- r_ties > 1
    r_max_rec_s <- mask(r_max_rec_s,r_ties_mask,maskvalue=1)
    r_rec_val_s <- lapply(1:nlayers(r_max_rec_s),
                          function(i,r_stack){df_subs <- data.frame(id=c(0,1),v=c(0,reclass_val[i]));
                          x <- subs(subset(r_stack,i), df_subs)},r_stack=r_max_rec_s)
    r_rec_val_s <- stack(r_rec_val_s)
    r_rec <- calc(r_rec_val_s,function(x){sum(x)})#overlays the layer with sum, 
    #x2 <- subs(r, df, subsWithNA=FALSE)
    
    ### prepare return object
    reclass_obj <- list(r_rec,r_rec_val_s,r_max_rec_s,r_ties)
    names(reclass_obj) <- c("r_rec","r_rec_val_s","r_max_rec_s","r_ties")
  }
  
  ###
  return(reclass_obj)
}

generate_soft_cat_aggregated_raster_fun <- function(r,reg_ref_rast,agg_fact,agg_fun,num_cores,NA_flag_val,file_format,out_dir,out_suffix){
  ## Function to aggregate categories
  
  ##INPUTS
  #1) r: raster to aggregate and crop/project if necessary
  #2) reg_ref_rast: reference raster, it must have a coordinate system defined
  #3) agg_fact:
  #4) agg_fun:
  #5) num_cores: number of cores used in the parallel processsing
  #6) NA_flag_val: flag value used for no data
  #7) file_format: raster file format e.g. .tif, .rst
  #8) out_dir: output directory
  #9) out_suffix: output suffix added to files
  #
  #OUTPUTS
  #
  #
  #
  
  
  ##### STEP 1: Check input
  
  if(class(r)!="RasterLayer"){
    r <- raster(r)
  }
  
  
  NAvalue(r) <- NA_flag_val #make sure we have a flag value assigned
  
  ###### STEP 1: BREAK OUT
  ## Breakout layers: separate categories in individual layers
  
  freq_tb <- as.data.frame(freq(r)) #zero is NA?
  out_filename_freq <- paste("freq_tb_",out_suffix,".txt",sep="")
  write.table(freq_tb,file= file.path(out_dir,out_filename_freq))
  
  ## get the names
  names_val <- freq_tb$value
  names_val <- names_val[!is.na(names_val)] #remove NA
  
  ## Make a brick composed of multiple layers: one layer per category (in one unique multiband file)
  out_raster_name <- file.path(out_dir,paste("r_layerized_bool_",out_suffix,file_format,sep=""))
  r_layerized <- layerize(r,
                          classes=names_val,
                          filename= out_raster_name,
                          overwrite=T)
  list_out_raster_name <- file.path(out_dir,paste("r_layerized_bool_",names_val,"_",out_suffix,file_format,sep=""))
  
  ## Now write out separate layers in one unique file (one per categories)
  #notice some issue here, looping through
  #for(i in 1:nlayers){
  #  writeRaster(r_layerized,
  #              bylayer=T,
  #              suffix=paste0(names_val,"_",out_suffix),
  #              filename=paste("r_layerized_bool",
  #                             file_format,sep="")
  #             ,overwrite=T)
  #}
  
  writeRaster(r_layerized,
              bylayer=T,
              suffix=paste0(names_val,"_",out_suffix),
              filename=paste("r_layerized_bool",
                             file_format,sep="")
              ,overwrite=T)
  
  #browser()
  
  lf_layerized_bool <- paste("r_layerized_bool_",names_val,"_",out_suffix,file_format,sep="")
  #names(r_layerized) <- 
  #inMemory(r_date_layerized)
  #filename(r_date_layerized)
  #r_test <- raster(lf_layerized_bool[1])
  
  ### STEP 2: aggregate
  ## Aggregate
  
  ##To do
  ##Set correct names based on categories of input
  ##Set output name
  
  #r_agg_fname <-aggregate_raster(agg_fact=agg_fact,
  #                               r_in=r_date_layerized,reg_ref_rast=NULL,agg_fun="mean",out_suffix=NULL,file_format=".tif",out_dir=NULL)
  #debug(aggregate_raster)
  #r_agg_fname <-aggregate_raster(r_in=raster(lf_layerized_bool[1]),
  #                               agg_fact=agg_fact,
  #                               reg_ref_rast=NULL,
  #                               #agg_fun="mean",
  #                               agg_fun=agg_fun,
  #                               out_suffix=NULL,
  #                               file_format=file_format,
  #                               out_dir=out_dir,
  #                               out_rast_name = NULL)
  #r_var_s <- mclapply(1:length(infile_var),FUN=import_list_modis_layers_fun,list_param=list_param_import_modis,mc.preschedule=FALSE,mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
  
  #debug(aggregate_raster)
  lf_agg <- aggregate_raster(lf_layerized_bool[1],
                             #r_in=raster(lf_layerized_bool[1]),
                             agg_fact=agg_fact,
                             reg_ref_rast=reg_ref_rast,
                             #agg_fun="mean",
                             agg_fun=agg_fun,
                             out_suffix=NULL,
                             file_format=file_format,
                             out_dir=out_dir,
                             out_rast_name = NULL)
  
  lf_agg <- mclapply(lf_layerized_bool,
                     FUN=aggregate_raster,
                     #r_in=raster(lf_layerized_bool[1]),
                     agg_fact=agg_fact,
                     reg_ref_rast=reg_ref_rast,
                     #agg_fun="mean",
                     agg_fun=agg_fun,
                     out_suffix=NULL,
                     file_format=file_format,
                     out_dir=out_dir,
                     out_rast_name = NULL,
                     mc.preschedule=FALSE,
                     mc.cores = num_cores) 
  #r_agg <- stack(unlist(lf_agg))
  raster_outname <- unlist(lf_agg)
  #apply function by layer using lapply.
  
  ## Reclassify by labels
  return(raster_outname)
}



#Function to aggregate from fine to coarse resolution, this will change accordingly once the input raster ref is given..
#
aggregate_raster <- function(r_in, agg_fact, reg_ref_rast=NULL,agg_fun="mean",out_suffix=NULL,file_format=".tif",out_dir=NULL,out_rast_name=NULL){
  #Aggregate raster from raster input and reference file
  #INPUT arguments:
  #1) r_in: input raster layer
  #2) agg_fact: factor to aggregate
  #3) reg_ref_rast: reference raster to match in resolution, if NULL then send a message
  #4) agg_fun: default is mean
  #5) out_suffix: output suffix
  #6) file_Format: raster format used e.g. .tif
  #7) out_dir: output directory
  #8) out_rast_name: output raster name if null it is created from the input file name
  #OUTPUT:
  # out_raster_name: name of the file containing the aggregated raster
  #
  # Authors: Benoit Parmentier
  # Created: 10/15/2015
  # Modified: 03/15/2017
  # To Do: 
  # - Add option to disaggregate
  #
  ################################
  
  ##If file is provided rather than RasterLayer
  if(class(r_in)!="RasterLayer"){
    r_in <- raster(r_in)
  }
  
  if(!is.null(reg_ref_rast)){
    ##If file is provided rather than RasterLayer
    if(class(reg_ref_rast)!="RasterLayer"){
      reg_ref_rast <- raster(reg_ref_rast)
    }
  }
  
  ### Compute aggregation factor if no ref image...
  if(is.null(agg_fact)){
    res_ref <- res(reg_ref_rast)[1] #assumes square cells, and decimal degrees from WGS84 for now...
    res_in <- res(r_in)[1] #input resolution, assumes decimal degrees
    agg_fact <-round(res_ref/res_in) #find the factor needed..
    #fix this to add other otpions e.g. aggregating down
  }
  
  #Default values...
  if(is.null(out_suffix)){
    out_suffix <- ""
  }
  
  if(is.null(out_dir)){
    out_dir <- "."
  }
  
  ## Create output raster name if out_rast_name is null
  if(is.null(out_rast_name)){
    raster_name <- filename(r_in)
    extension_str <- extension(raster_name)
    raster_name_tmp <- gsub(extension_str,"",basename(raster_name))
    out_rast_name <- file.path(out_dir,paste("agg_",agg_fact,"_",raster_name_tmp,out_suffix,file_format,sep="")) #output name for aster file
    #out_rast_name <- raster_name #for use in function later...
  }
  
  r_agg <- aggregate(r_in, 
                     fact=agg_fact,
                     FUN=agg_fun,
                     filename=out_rast_name,
                     overwrite=TRUE)
  
  return(out_rast_name)
  
}
