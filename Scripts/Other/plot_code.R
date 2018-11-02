

# Plot pre-IAR halfmax estimates ------------------------------------------


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
            aes(x = x, y = y), color = 'black') + 
  geom_path(data = canadamap, 
            aes(x = x, y = y), color = 'black') + 
  geom_path(data = mexicomap, 
            aes(x = x, y = y), color = 'black') + 
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


#cells that should probbaly be removed (can be filtered more formally at the query stage in the future)
cells_to_rm <- c(208, 845, 289, 316, 343, 345, 346, 319, 320, 373, 
                 401, 402, 375, 376, 402, 403, 404, 370, 770, 796, 
                 766, 765, 763, 735, 734, 731, 758, 813, 3669, 397,
                 398)





'%ni%' <- Negate('%in%')
to_plt3 <- to_plt2[which(to_plt2$cell %ni% cells_to_rm),]


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

