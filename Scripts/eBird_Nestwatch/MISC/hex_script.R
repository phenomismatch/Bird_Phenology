###############
#Script to create hex cells, filter them, and plot them
#
###############


# load packages -----------------------------------------------------------

library(dplyr)
library(dggridR)
library(maps)
library(rgdal)
library(raster)
library(sp)
library(ggplot2)


# create grid -------------------------------------------------------------

#create hexgrid
hexgrid6 <- dggridR::dgconstruct(res = 6)

#save hexgrid as shp file
#dggridR::dgearthgrid(hexgrid6, savegrid = 'global_hex.shp')

#read in grid
hge <- rgdal::readOGR('global_hex.shp', verbose = FALSE)



# read in bird range map --------------------------------------------------

#get shp filename and read in
sp_rng <- rgdal::readOGR('Vireo_olivaceus_22705243.shp', verbose = FALSE)

#crop to area of interest
sp_rng2 <- raster::crop(sp_rng, raster::extent(-95, -50, 24, 90))


#filter by breeding (2) and migration (4) range - need to convert spdf to sp
nrng <- sp_rng2[which(sp_rng2$SEASONAL == 2 | sp_rng2$SEASONAL == 4),]

#filter by resident (1) and non-breeding (3) to exclude hex cells that contain 2/4 and 1/3
nrng_rm <- sp_rng2[which(sp_rng2$SEASONAL == 1 | sp_rng2$SEASONAL == 3),]

#add CRS
nrng_sp <- sp::SpatialPolygons(nrng@polygons)
sp::proj4string(nrng_sp) <- sp::CRS(sp::proj4string(nrng))

#find intersections with code from here: https://gis.stackexchange.com/questions/140504/extracting-intersection-areas-in-r
poly_int <- rgeos::gIntersects(hge, nrng_sp, byid=TRUE)

tpoly <- which(poly_int == TRUE, arr.ind = TRUE)[,2]
cell_mrg <- as.numeric(tpoly[!duplicated(tpoly)])



# cells and cell centers --------------------------------------------------

#make data frame with cells and lat/lons of cell centers
cell_centers <- dggridR::dgSEQNUM_to_GEO(hexgrid6, cell_mrg)
cc_df <- data.frame(cell = cell_mrg, 
                    lon = cell_centers$lon_deg, 
                    lat = cell_centers$lat_deg)


#transform cells to grid
cell_grid <- dggridR::dgcellstogrid(hexgrid6, cc_df$cell)
cell_grid$cell <- as.numeric(cell_grid$cell)

#select particular cell of interest
cell_676 <- dplyr::filter(cell_grid, cell == 676)


# plot grid ---------------------------------------------------------------

usamap <- data.frame(maps::map("world", "USA", plot = FALSE)[c("x", "y")])
canadamap <- data.frame(maps::map("world", "Canada", plot = FALSE)[c("x", "y")])
mexicomap <- data.frame(maps::map("world", "Mexico", plot = FALSE)[c("x", "y")])


#can comment out lines to highlight cell of interest/add cell numbers below

p <- ggplot() +
  geom_path(data = usamap,
            aes(x = x, y = y), color = 'black', size = 0.5) +
  geom_path(data = canadamap,
            aes(x = x, y = y), color = 'black', size = 0.5) +
  geom_path(data = mexicomap,
            aes(x = x, y = y), color = 'black', size = 0.5) +
  coord_map("ortho", orientation = c(35, -80, 0),
            xlim = c(-105, -55), ylim = c(20, 66)) +
  geom_polygon(data = cell_grid, aes(x = long, y = lat, group = group),
               alpha = 0.2) +
  #highlight cell of interest
  geom_polygon(data = cell_676, aes(x = long, y = lat, group = group),
               alpha = 0.5, fill = 'red') +
  #add cell number on plot
  annotate('text', x = cc_df$lon, y = cc_df$lat,
           label = cc_df$cell, col = 'blue', alpha = 0.8,
           size = 3) +
  geom_path(data = cell_grid, aes(x = long, y = lat, group = group),
            alpha = 0.4, color = 'black') +
  theme_bw() +
  xlab('Longitude') +
  ylab('Latitude')

p
