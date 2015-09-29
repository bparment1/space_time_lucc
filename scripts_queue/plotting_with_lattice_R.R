install.packages("gridExtra")
require(gridExtra) # also loads grid
require(lattice)
x <- seq(pi/4, 5 * pi, length.out = 100)
y <- seq(pi/4, 5 * pi, length.out = 100)
r <- as.vector(sqrt(outer(x^2, y^2, "+")))

grid <- expand.grid(x=x, y=y)
grid$z <- cos(r^2) * exp(-r/(pi^3))
plot1 <- levelplot(z~x*y, grid, cuts = 50, scales=list(log="e"), xlab="",
                   ylab="", main="Weird Function", sub="with log scales",
                   colorkey = FALSE, region = TRUE)

plot2 <- levelplot(z~x*y, grid, cuts = 50, scales=list(log="e"), xlab="",
                   ylab="", main="Weird Function", sub="with log scales",
                   colorkey = FALSE, region = TRUE)
grid.arrange(plot1,plot2, ncol=2)

### Wthie line issues
library(lattice)
library(Cairo) #install it if you don’t have it.

# makes some interesting data to plot
x <- y <- seq(-1,1,.005)
dataGrid <- data.frame(expand.grid(x=x,y=y))

# Base graphics image() is similar to lattice levelplot().
# But it takes x, y and z where z is a matrix
z <- with(dataGrid, x * exp(-x^2 - y^2))
dim(z) <- c(length(y), length(y))

# Plot the image. It has the ugly white lines (at least on my win7/64 PC).
image(x,y,z, xlim=c(-1,1), ylim=c(-1,1), asp=1)

# Plot the image with an option. Voila, lines are gone!
# But there is no color key on this plot :(
image(x,y,z, xlim=c(-1,1), ylim=c(-1,1), asp=1, useRaster=TRUE)


# levelplot() from lattice will take either the matrix or a dataframe.
# Let’s try with the matrix first, since we have it already.
# This also has ugly white lines in the plot :(
lp1 <- levelplot(z, aspect = "iso", xlim=c(-1,1), ylim=c(-1,1),
                 row.values = x, column.values = y)
plot(lp1)

# Here is the levelplot method for data frames. What you are using. 
# Syntax is bit easier, but still has the white lines.
lp2 <- levelplot(z ~ x * y, aspect = "iso", data = dataGrid)
plot(lp2)

# Using the non-default panel function for levelplot makes the white
# The lines disappear. I highly recommend the book on Lattice if you 
# plan to do statistical plotting frequently.
lp3 <- levelplot(z ~ x * y, aspect = "iso", data = dataGrid,
                 raster = TRUE, # this only rasters the key! (not working?? Key still has white lines)
                 panel = function (x,y,z, ...){
                   panel.levelplot.raster(x, y, z,...)
                 })
plot(lp3)

# Note. Other graphics devices will automatically raster, but they may 
# struggle plotting really big data. 
CairoWin()
plot(lp2) # the "bad" lattice levelplot is now good!

