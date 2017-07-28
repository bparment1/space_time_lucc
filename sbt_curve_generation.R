library(raster)
library(rasterVis)
library(colorRamps)

n_row <- 3
n_col <- 3

r <- raster(nrows=n_row,ncols=n_col)

vtn1 <- c(0.8,0.8,0.2,
         0.8,0.2,0.2,
         0.2,0.8,0.8)

vt0 <- c(0.8,0.8,0.2,
         0.8,0.2,0.2,
         0.2,0.8,0.8)

vt1 <- c(0.3,0.3,0.2,
         0.3,0.2,0.3,
         0.2,0.3,0.3)

vt2 <- c(0.4,0.4,0.2,
         0.4,0.2,0.4,
         0.2,0.4,0.4)

vt3 <- c(0.6,0.6,0.2,
         0.6,0.2,0.6,
         0.2,0.6,0.6)

vt4 <- c(0.6,0.6,0.2,
         0.6,0.2,0.6,
         0.2,0.6,0.6)

list_vt <- list(vtn1,vt0,vt1,vt2,vt3,vt4)

list_r <-(lapply(list_vt, function(x,rast){rast[]<-x;return(rast)},rast=r))
r_var <- stack(list_r)

## generating additional variables
r_zonal <- subset(r_var,1)
r_zonal[] <- rep(1,ncell(r_zonal))
r_x <-init(r_var,v="x")
r_y <-init(r_var,v="y")

col_palette <- matlab.like(100)
#col_palette <- matlab.like(100)


levelplot(r_var,col.regions=col_palette,main="Time series")



#?Moran
f <- matrix(c(1,1,1,
              1,0,1,
              1,1,1), nrow=3)

t_corr_fun <- function(i,list_rast){
  r_subset <- stack(list_rast[[i-1]],list_rast[[i]]) 
  cor_val<-layerStats(r_subset,'pearson')
  cor_val <- cor_val[[1]][1,2]
  return(cor_val)
}

#Call function


################# END OF SCRIPT ##################

