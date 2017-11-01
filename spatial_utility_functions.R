

r_test <- generate_outline_and_grid(r_val)

generate_outline_and_grid <- function(reg_layer,n_grid=NULL,out_suffix="",out_dir="."){
  #This function generates a grid from an input layer
  if(class(reg_layer)=="RasterLayer"){
    r <- reg_layer
    r[] <- 1 #change this later
    r <- mask(r,reg_layer)
    reg_sp <- rasterToPolygons(r,dissolve=T)
    reg_layer <- as(reg_sp, "sf")
  }
  if(class(reg_layer)=="SpatialPolygons"){
    reg_layer <- as(reg_layer, "sf")
  }
  ## generate grid: default is 2x2
  if(is.null(n_grid)){
    n_grid <- c(2,2)
  }
  reg_grid <- st_make_grid(reg_layer,n=n_grid,offset = n_distance)
  
  ## save grid figures
  plot(reg_sf)
  plot(reg_grid,border="red",add=T)
  
  ### buffer: make overlapping area later!!!
  #n_distance <- 5000
  #reg_buffer <- st_buffer(reg_layer,dist=n_distance)
  
  out_filename <- file.path(out_dir,
                            paste0("reg_grid_",out_suffix,".shp"))
  st_write(reg_grid,out_filename)
  
  return(reg_grid)
}


    