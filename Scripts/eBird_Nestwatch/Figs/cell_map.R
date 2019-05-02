##############
#create cell maps
#
##############



# load maps ---------------------------------------------------------------

usa_m <- maps::map('usa', plot = FALSE, fill = TRUE)
IDs <- sapply(strsplit(usa_m$names, ":"), function(x) x[1])
usa <- maptools::map2SpatialPolygons(usa_m, IDs = IDs, proj4string = sp::CRS('+init=epsg:4326'))

canada_m <- maps::map('world', regions = 'Canada', plot = FALSE, fill = TRUE)
IDs <- sapply(strsplit(canada_m$names, ":"), function(x) x[1])
canada <- maptools::map2SpatialPolygons(canada_m, IDs = IDs, proj4string = sp::CRS('+init=epsg:4326'))



# distribute points uniformly on sphere -----------------------------------

#from: http://mathworld.wolfram.com/SpherePointPicking.html
set.seed(1)

N <- 10000000    #How many cells to sample
u     <- runif(N)
v     <- runif(N)
theta <- 2*pi*u      * 180/pi
phi   <- acos(2*v-1) * 180/pi
lon   <- theta-180
lat   <- phi-90
df    <- data.frame(lat = lat, lon = lon)



# convert to spatial points - WGS84 ---------------------------------------

pts <- sp::SpatialPoints(cbind(df$lon, df$lat), 
                         proj4string = sp::CRS('+init=epsg:4326'))




# particle filter to select points inside US ------------------------------

nn <- which(!is.na(sp::over(pts, usa)) | !is.na(sp::over(pts, canada)))
npts <- pts[nn]
ndf <- df[nn,]

ndf2 <- dplyr::filter(ndf, lat > 26 & lon < -50 & lon > -95)

# construct hex cells using points in US ----------------------------------

hexgrid6 <- dggridR::dgconstruct(res = 6)
cells <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                  in_lon_deg = ndf2$lon, in_lat_deg = ndf2$lat)[[1]]
cell_grid <- dggridR::dgcellstogrid(hexgrid6, cells)

cell_3673 <- dplyr::filter(cell_grid, cell == 3673)
cell_676 <- dplyr::filter(cell_grid, cell == 676)

# ggplot ------------------------------------------------------------------

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
  # geom_polygon(data = cell_676, aes(x = long, y = lat, group = group),
  #              alpha = 0.5, fill = 'red') +
  geom_path(data = cell_grid, aes(x = long, y = lat, group = group),
            alpha = 0.4, color = 'black') +
  theme_bw() +
  xlab('Longitude') +
  ylab('Latitude')


setwd('~/Desktop/')
ggsave(plot = p,
       filename = 'cell_map_676.pdf')
  