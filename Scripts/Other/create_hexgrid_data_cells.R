cells <- unique(ucell_grid$cell)

N <- length(cells)
tt <- vector("list", N)
for (i in 1:N)
{
  #i <- 1
  uc_t <- dplyr::filter(ucell_grid, cell == cells[i])
  tpol <- sp::Polygon(uc_t[,1:2])
  ps <- sp::Polygons(list(tpol), ID = cells[i])
  sps <- sp::SpatialPolygons(list(ps))

  tt[[i]] <- sps
}

joined <- sp::SpatialPolygons(lapply(tt, function(x){x@polygons[[1]]}))
jdata <- sp::SpatialPolygonsDataFrame(Sr = joined, 
                                        data = data.frame(cell = cells[1:N]),FALSE)

proj4string(jdata) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#proj4string(jdata) = CRS("+proj=longlat +datum=WGS84 +no_defs")


setwd("~/Google_Drive/R/Bird_Phenology/Data/hex_grid_crop")
rgdal::writeOGR(jdata, "hex_grid_crop_data.shp", "", "ESRI Shapefile")



# check -------------------------------------------------------------------

setwd("~/Google_Drive/R/Bird_Phenology/Data/hex_grid_crop")
jdata2 <- readOGR('hex_grid_crop_data.shp')
hge <- rgdal::readOGR('hex_grid_crop.shp', verbose = FALSE)

plot(hge)
plot(jdata2, col = 'red', add = TRUE)

jdata2
hge

hge_cells <- as.numeric(as.character(hge@data[,1]))
jdata2_cells <- as.numeric(as.character(jdata2@data[,1]))

