
# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/home/CAM/cyoungflesh/phenomismatch/'




# load packages -----------------------------------------------------------

library(ggplot2)
library(rgdal)
library(sp)


# db/hm query dir ------------------------------------------------------------

db_dir <- 'db_query_2018-10-15'
hm_dir <- 'halfmax_species_2018-10-16'
IAR_dir <- 'IAR_2018-10-26'



setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_dir))

IAR_date <- substr(IAR_dir, start = 5, stop = 15)

diagnostics_frame_p <- readRDS(paste0('diagnostics_frame-', IAR_date,'.rds'))
cells_frame <- readRDS(paste0('cells_frame-', IAR_date, '.rds'))
yrs_frame <- readRDS(paste0('yrs_frame-', IAR_date, '.rds'))

diagnostics_frame <- diagnostics_frame_p

#which species/year has the most cells - to model 'data rich species'
DR_sp <- as.character(yrs_frame[which.max(yrs_frame[,1:3]$n_cells),1])
DR_filt <- dplyr::filter(yrs_frame, species == DR_sp)[,1:3]
DR_yr <- c(2002, 2010, 2017)



# plot cell numbers on map ------------------------------------------------

#create hex cell grid
hexgrid6 <- dggridR::dgconstruct(res = 6)

#get boundaries of all cells over earth
hge <- dggridR::dgearthgrid(hexgrid6)

#load species range map
setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/shapefiles/'))
sp_rng <- rgdal::readOGR('Vireo_olivaceus_22705243.shp', verbose = FALSE)

#filter by breeding range
nrng <- subset(sp_rng, SEASONAL == 2)

#convert cell vetices to spatial points - need to convert spdf to sp
nrng_sp <- sp::SpatialPolygons(nrng@polygons)
sp::proj4string(nrng_sp) <- sp::CRS(sp::proj4string(nrng))
ecells_sp <- sp::SpatialPoints(cbind(hge$long, hge$lat), 
                               proj4string = sp::CRS(sp::proj4string(nrng)))

#which vertices overlap species breeding range
ovrlp <- which(!is.na(sp::over(ecells_sp, nrng_sp)))
overlap_cells <- unique(as.numeric(hge[ovrlp,]$cell))

#get cell centers
cell_centers <- dggridR::dgSEQNUM_to_GEO(hexgrid6, overlap_cells)
cc_df <- data.frame(cell = overlap_cells, lon = cell_centers$lon_deg, 
                    lat = cell_centers$lat_deg)

#cells only within the range that ebird surveys were filtered to
n_cc_df <- cc_df[which(cc_df$lon > -100 & cc_df$lon < -50 & cc_df$lat > 26),]

#create cell grid
cell_grid <- dggridR::dgcellstogrid(hexgrid6, n_cc_df$cell)



p <- ggplot() + 
  coord_map("ortho", orientation = c(35, -80, 0), 
            xlim = c(-100, -55), ylim = c(20, 90)) + 
  geom_polygon(data = cell_grid, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  #geom_polygon(data = nrng, aes(x = long, y = lat), fill = 'red', alpha = 0.5) + 
  theme_bw()
p



cells <- unique(diagnostics_frame$cell)
ncel <- length(cells)



#transform cells to grid
cell_grid <- dggridR::dgcellstogrid(hexgrid6, cells)
cell_grid$cell <- as.numeric(cell_grid$cell)
cell_centers <- dggridR::dgSEQNUM_to_GEO(hexgrid6, cells)
to_plt <- data.frame(cell = cells, lon_deg = cell_centers$lon_deg, lat_deg = cell_centers$lat_deg)

to_plt2 <- dplyr::inner_join(to_plt, cell_grid, by = 'cell')


#load maps
usamap <- data.frame(maps::map("world", "USA", plot = FALSE)[c("x", "y")])
canadamap <- data.frame(maps::map("world", "Canada", plot = FALSE)[c("x", "y")])
mexicomap <- data.frame(maps::map("world", "Mexico", plot = FALSE)[c("x", "y")])

#plot
p <- ggplot() +
  geom_path(data = usamap, 
            aes(x = x, y = y), color = 'black', alpha = 0.5) + 
  geom_path(data = canadamap, 
            aes(x = x, y = y), color = 'black', alpha = 0.5) + 
  geom_path(data = mexicomap, 
            aes(x = x, y = y), color = 'black', alpha = 0.5) + 
  coord_map("ortho", orientation = c(35, -80, 0), 
            xlim = c(-100, -55), ylim = c(20, 90)) + 
  geom_polygon(data = to_plt2, aes(x = long, y = lat, group = group), fill = 'red', alpha = 0.3) +
  geom_path(data = to_plt2, aes(x = long, y = lat, group = group), col = 'black', alpha = 0.3) + 
  annotate('text', x = to_plt2$lon_deg, y = to_plt2$lat_deg, 
           label = to_plt2$cell, col = 'black', alpha = 0.3,
           size = 4) +
  ggtitle('Cell Number') +
  theme_bw() +
  xlab('Longitude') +
  ylab('Latitude')

p





# filter cells ------------------------------------------------------------

#cells identified manually from plot (plot_code.R) to filter
cells_to_keep <- c(622, 595, 567, 594, 621, 648, 649, 676, 675, 647, 619, 
                   620, 621, 593, 592, 619, 618, 591, 564, 565, 566, 646, 
                   674, 702, 703)


to_plt3 <- to_plt2[which(to_plt2$cell %in% cells_to_keep),]


#plot
p2 <- ggplot() +
  geom_path(data = usamap, 
            aes(x = x, y = y), color = 'black') + 
  geom_path(data = canadamap, 
            aes(x = x, y = y), color = 'black') + 
  geom_path(data = mexicomap, 
            aes(x = x, y = y), color = 'black') + 
  coord_map("ortho", orientation = c(35, -80, 0), 
            xlim = c(-100, -55), ylim = c(20, 90)) + 
  geom_polygon(data = to_plt3, aes(x = long, y = lat, group = group), fill = 'red', alpha = 0.3) +
  geom_path(data = to_plt3, aes(x = long, y = lat, group = group), col = 'black', alpha = 0.3) + 
  annotate('text', x = to_plt3$lon_deg, y = to_plt3$lat_deg, 
           label = to_plt3$cell, col = 'black', alpha = 0.3,
           size = 4) +
  ggtitle('Cell Number') +
  theme_bw() +
  xlab('Longitude') +
  ylab('Latitude')

p2

