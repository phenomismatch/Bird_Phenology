##############
#create cell maps
#
##############


# Load packages -----------------------------------------------------------

library(dplyr)
library(ggplot2)
library(maps)


# load maps ---------------------------------------------------------------

usa_m <- maps::map('usa', plot = FALSE, fill = TRUE)
IDs <- sapply(strsplit(usa_m$names, ":"), function(x) x[1])
usa <- maptools::map2SpatialPolygons(usa_m, IDs = IDs, proj4string = sp::CRS('+init=epsg:4326'))

canada_m <- maps::map('world', regions = 'Canada', plot = FALSE, fill = TRUE)
IDs <- sapply(strsplit(canada_m$names, ":"), function(x) x[1])
canada <- maptools::map2SpatialPolygons(canada_m, IDs = IDs, proj4string = sp::CRS('+init=epsg:4326'))

#load maps
usamap <- data.frame(maps::map("world", "USA", plot = FALSE)[c("x", "y")])
canadamap <- data.frame(maps::map("world", "Canada", plot = FALSE)[c("x", "y")])
mexicomap <- data.frame(maps::map("world", "Mexico", plot = FALSE)[c("x", "y")])

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

ndf2 <- dplyr::filter(ndf, lat > 24 & lon < -50 & lon > -95)
#ndf2 <- df[which(df$lat > 24 & df$lat < 90 & df$lon < -50 & df$lon > -95),]

# construct hex cells using points in US ----------------------------------

hexgrid6 <- dggridR::dgconstruct(res = 6)
cells <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                  in_lon_deg = ndf2$lon, in_lat_deg = ndf2$lat)[[1]]
#get cell centers
cell_centers <- dggridR::dgSEQNUM_to_GEO(hexgrid6, cells)
cc_df <- data.frame(cell = cells, lon = cell_centers$lon_deg, 
                    lat = cell_centers$lat_deg)

#cells only within the range that ebird surveys were filtered to
n_cc_df <- unique(cc_df[which(cc_df$lon > -95 & cc_df$lon < -50 & cc_df$lat > 24),])
cells2 <- n_cc_df$cell

#from near top of 2-halfmax script
#cells2 <- cells

#cells to grid
cell_grid <- dggridR::dgcellstogrid(hexgrid6, cells2)

# cell_3673 <- dplyr::filter(cell_grid, cell == 3673)
# cell_676 <- dplyr::filter(cell_grid, cell == 676)




# ggplot ------------------------------------------------------------------

# args <- 'XXXX'
# #reference key for species synonyms
# setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/metadata/'))
# sp_key <- read.csv('species_filenames_key.csv')
# 
# #change dir to shp files
# setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/shapefiles/'))
# 
# #filter by breeding/migration cells
# #match species name to shp file name
# g_ind <- grep(args, sp_key$file_names_2016)
# 
# #check for synonyms if there are no matches
# if (length(g_ind) == 0)
# {
#   g_ind2 <- grep(args, sp_key$BL_Checklist_name)
# } else {
#   g_ind2 <- g_ind
# }
# 
# #get filename and read in
# fname <- as.character(sp_key[g_ind2,]$filenames[grep('.shp', sp_key[g_ind2, 'filenames'])])
# sp_rng <- rgdal::readOGR(fname[1], verbose = FALSE)
# 
# #filter by breeding (2) and migration (4) range - need to convert spdf to sp
# nrng <- sp_rng[which(sp_rng$SEASONAL == 2 | sp_rng$SEASONAL == 4),]
# 
# #filter by resident (1) and non-breeding (3) to exclude hex cells that contain 2/4 and 1/3
# nrng_rm <- sp_rng[which(sp_rng$SEASONAL == 1 | sp_rng$SEASONAL == 3),]
# 
# nrng_sp <- sp::SpatialPolygons(nrng@polygons)
# sp::proj4string(nrng_sp) <- sp::CRS(sp::proj4string(nrng))


p <- ggplot() +
  geom_path(data = usamap,
            aes(x = x, y = y), color = 'black', size = 0.5) +
  geom_path(data = canadamap,
            aes(x = x, y = y), color = 'black', size = 0.5) +
  geom_path(data = mexicomap,
            aes(x = x, y = y), color = 'black', size = 0.5) +
  coord_map("ortho", orientation = c(35, -80, 0),
            xlim = c(-105, -55), ylim = c(20, 66)) +
  # geom_polygon(data = nrng_sp,
  #              aes(x = long, y = lat, group=group), fill = 'green', alpha = 0.4) +
  geom_polygon(data = cell_grid, aes(x = long, y = lat, group = group),
               alpha = 0.2) +
  # #highlight cell of interest
  # geom_polygon(data = cell_676, aes(x = long, y = lat, group = group),
  #              alpha = 0.5, fill = 'red') +
  #add cell number on plot
  annotate('text', x = n_cc_df$lon, y = n_cc_df$lat,
           label = n_cc_df$cell, col = 'blue', alpha = 0.8,
           size = 3) +
  geom_path(data = cell_grid, aes(x = long, y = lat, group = group),
            alpha = 0.1, color = 'black', size = 0.1) +
  theme_bw() +
  xlab('Longitude') +
  ylab('Latitude')

p
setwd('~/Desktop/')
ggsave(plot = p,
       filename = 'cell_map.pdf')
  