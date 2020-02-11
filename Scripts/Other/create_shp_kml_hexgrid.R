###############
# Create hex grid and write to shp and kml files
#
###############


# load packages -----------------------------------------------------------

library(dggridR)
library(sp)
library(raster)
library(maps)


# create grid -------------------------------------------------------------

#create grid
hexgrid6 <- dggridR::dgconstruct(res = 6)
global <- dggridR::dgearthgrid(hexgrid6, frame = FALSE)
#add cell names
global.cell <- data.frame(cell = sp::getSpPPolygonsIDSlots(global), 
                          row.names = sp::getSpPPolygonsIDSlots(global))
#create SPDF
global <- sp::SpatialPolygonsDataFrame(global, global.cell)



# Fix line issues and crop ----------------------------------------------------

#sort out lines issues - from here: https://github.com/r-barnes/dggridR/issues/35
for (i in 1:length(global@polygons)) 
{
  if(max(global@polygons[[i]]@Polygons[[1]]@coords[,1]) - 
     min(global@polygons[[i]]@Polygons[[1]]@coords[,1]) > 270) 
  {
      global@polygons[[i]]@Polygons[[1]]@coords[,1] <- (global@polygons[[i]]@Polygons[[1]]@coords[,1] + 360) %% 360
  }
}


# crop --------------------------------------------------------------------

#get cell ids
cells <- as.numeric(as.vector.factor(global@data$cell))
#get lat/lon of cell centers
cell_centers <- dggridR::dgSEQNUM_to_GEO(hexgrid6, cells)

#only cell centers than fall within specified range
cell_idx <- which(cell_centers$lon_deg < -50 & 
                    cell_centers$lon_deg > -95 &
                    cell_centers$lat_deg > 24 &
                    cell_centers$lat_deg < 90)

#filter SPDF by idx
out <- global[cell_idx,]


# visualize on map --------------------------------------------------------

world <- data.frame(maps::map("world", plot = FALSE)[c("x", "y")])
plot(world, type = 'l')
lines(out, col = 'red')



# save out files ----------------------------------------------------------

setwd('~/Desktop/')

#save as shp file
writeOGR(out, "hex_grid_crop.shp", "", "ESRI Shapefile")
#save as kml file
writeOGR(out, "hex_grid_crop.kml", "", "KML")

