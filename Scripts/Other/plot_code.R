
tt <- data.frame(lng = rep(NA, 3), lat = rep(NA, 3))
tt$lng <- c(-73, -83, -107)
tt$lat <- c(40, 30, 37)

#get cell #s
cell_ns <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                    in_lon_deg = tt$lng, 
                                    in_lat_deg = tt$lat)[[1]]
#get lat/lon cell centers
cell_centers <- dggridR::dgSEQNUM_to_GEO(hexgrid6, cell_ns)

#transform cells to grid
cell_grid <- dggridR::dgcellstogrid(hexgrid6, cells)




library(ggplot2)

all_states <- map_data("state")
w2hr <- map_data("world")

usamap <- data.frame(map("world", "USA", plot = FALSE)[c("x", "y")])
canadamap <- data.frame(map("world", "Canada", plot = FALSE)[c("x", "y")])
mexicomap <- data.frame(map("world", "Mexico", plot = FALSE)[c("x", "y")])

p <- ggplot() +
  geom_path(data = usamap, 
            aes(x = x, y = y), color = 'black') + 
  geom_path(data = canadamap, 
            aes(x = x, y = y), color = 'black') + 
  geom_path(data = mexicomap, 
            aes(x = x, y = y), color = 'black') + 
  coord_map("ortho", orientation = c(35, -95, 0)) + 
  geom_polygon(data = cell_grid, aes(x = long, y = lat, group = group, fill = 'red'), 
               alpha = 0.4) 
p
