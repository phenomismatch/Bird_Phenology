
# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/home/CAM/cyoungflesh/phenomismatch/'




# load packages -----------------------------------------------------------

library(ggplot2)
library(rgdal)
library(sp)

# plot cell numbers on map ------------------------------------------------

#create hex cell grid
hexgrid6 <- dggridR::dgconstruct(res = 6)

#get boundaries of all cells over earth
hge <- dggridR::dgearthgrid(hexgrid6)

all_cells <- unique(as.numeric(hge$cell))

#get cell centers
cell_centers <- dggridR::dgSEQNUM_to_GEO(hexgrid6, all_cells)
cc_df <- data.frame(cell = all_cells, lon = cell_centers$lon_deg, lat = cell_centers$lat_deg)

#cells only within the range that ebird surveys were filtered to
n_cc_df <- cc_df[which(cc_df$lon > -100 & cc_df$lon < -50 & cc_df$lat > 26 & cc_df$lat < 65),]

#create cell grid
cell_grid <- dggridR::dgcellstogrid(hexgrid6, n_cc_df$cell)


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
  geom_polygon(data = cell_grid, aes(x = long, y = lat, group = group), fill = 'red', alpha = 0.3) +
  geom_path(data = cell_grid, aes(x = long, y = lat, group = group), col = 'black', alpha = 0.3) + 
  annotate('text', x = n_cc_df$lon, y = n_cc_df$lat, 
           label = n_cc_df$cell, col = 'black', alpha = 0.6,
           size = 3) +
  ggtitle('Cell Number') +
  theme_bw() +
  xlab('Longitude') +
  ylab('Latitude')

p

