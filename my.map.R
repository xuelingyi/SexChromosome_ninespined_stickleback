library(mapdata)
library(maptools)
library(dplyr)
library(ggsn)
library(PBSmapping)
library(data.table)
library(ggrepel)
library(ggplot2)
library(ggpubr)
my.map = function(sif, ...) {
  xmin = min(sif$long_ling) - 1
  xmax = max(sif$long_ling) + 1
  ymin = min(sif$lat_ling) - 1
  ymax = max(sif$lat_ling) + 1
  
  world = map_data("world")
  setnames(world, c("X","Y","PID","POS","region","subregion"))
  map = clipPolys(world, xlim=c(xmin, xmax), ylim=c(ymin, ymax), keepExtra=TRUE)
  
  mercator = ggplot() + coord_map(xlim=c(xmin, xmax), ylim=c(ymin, ymax)) + 
    geom_polygon(data=map, aes(X,Y,group=PID), fill="grey97", color="grey50") + 
    theme_void() 
  return(mercator)
}