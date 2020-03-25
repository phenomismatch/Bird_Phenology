
#dir <- '~/Google_Drive/R/'
dir <- '/labs/Tingley/phenomismatch/'

library(dplyr)
library(dggridR)
library(sp)
library(raster)
library(rgeos)
library(rgdal)

# create grid -------------------------------------------------------------

#make hexgrid
hexgrid6 <- dggridR::dgconstruct(res = 6)

setwd(paste0(dir, 'Bird_Phenology/Data/hex_grid_crop/'))

#read in grid
hge <- rgdal::readOGR('hex_grid_crop.shp', verbose = FALSE)
hge_cells <- as.numeric(as.character(hge@data[,1]))


# filter cells by range  ---------------------------------------------------

#reference key for species synonyms
setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/metadata/'))
sp_key <- read.csv('species_filenames_key.csv')




setwd(paste0(dir, 'Bird_Phenology/Data/'))

species_list_i <- read.table('eBird_species_list.txt', stringsAsFactors = FALSE)
species_list <- species_list_i[,1]
nsp <- length(species_list)



out <- list()
for (i in 1:nsp)
{
  #i <- 20
  #i <- 44
  args <- species_list[i]
  
  print(args)
  #change dir to shp files
  setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/shapefiles/'))
  
  #filter by breeding/migration cells
  #match species name to shp file name
  g_ind <- grep(args, sp_key$file_names_2016)
  
  #check for synonyms if there are no matches
  if (length(g_ind) == 0)
  {
    g_ind2 <- grep(args, sp_key$BL_Checklist_name)
  } else {
    g_ind2 <- g_ind
  }
  
  #get filename and read in
  fname <- as.character(sp_key[g_ind2,]$filenames[grep('.shp', sp_key[g_ind2, 'filenames'])])
  sp_rng <- rgdal::readOGR(fname[1], verbose = FALSE)
  
  #filter by breeding (2) and migration (4) range - need to convert spdf to sp
  nrng <- sp_rng[which(sp_rng$SEASONAL == 2 | sp_rng$SEASONAL == 4),]
  
  #filter by resident (1) and non-breeding (3) to exclude hex cells that contain 2/4   and 1/3
  nrng_rm <- sp_rng[which(sp_rng$SEASONAL == 1 | sp_rng$SEASONAL == 3),]
  
  #remove unneeded objects
  rm(sp_rng)
  rm(fname)
  
  #if there is a legitimate range
  if (NROW(nrng@data) > 0 & extent(nrng)@xmax > -95)
  {
    #good cells
    nrng_sp <- sp::SpatialPolygons(nrng@polygons)
    sp::proj4string(nrng_sp) <- sp::CRS(sp::proj4string(nrng))
    #find intersections with code from here: https://gis.stackexchange.com/questions  /140504/extracting-intersection-areas-in-r
    poly_int <- rgeos::gIntersects(hge, nrng_sp, byid=TRUE)
    tpoly <- which(poly_int == TRUE, arr.ind = TRUE)[,2]
    br_mig_cells <- hge_cells[as.numeric(tpoly[!duplicated(tpoly)])]
    
    #bad cells
    if (length(nrng_rm) > 0)
    {
      nrng_rm_sp <- sp::SpatialPolygons(nrng_rm@polygons)
      sp::proj4string(nrng_rm_sp) <- sp::CRS(sp::proj4string(nrng_rm))
      poly_int_rm <- rgeos::gIntersects(hge, nrng_rm_sp, byid=TRUE)
      tpoly_rm <- which(poly_int_rm == TRUE, arr.ind = TRUE)[,2]
      res_ovr_cells <- hge_cells[as.numeric(tpoly_rm[!duplicated(tpoly_rm)])]
      
      #remove cells that appear in resident and overwinter range that also appear in   breeding range
      cell_mrg <- c(br_mig_cells, res_ovr_cells)
      to_rm <- c(cell_mrg[duplicated(cell_mrg)], 812, 813, 841)
      
      rm(nrng_rm)
      rm(nrng_rm_sp)
      rm(res_ovr_cells)
      
    } else {
      cell_mrg <- br_mig_cells
      to_rm <- c(812, 813, 841)
    }
    
    c_rm <- which(br_mig_cells %in% to_rm)
    if (length(c_rm) > 0)
    {
      overlap_cells <- br_mig_cells[-c_rm]  
    } else {
      overlap_cells <- br_mig_cells
    }

    #get cell centers
    cell_centers <- dggridR::dgSEQNUM_to_GEO(hexgrid6, overlap_cells)
    cc_df <- data.frame(cell = overlap_cells, lon = cell_centers$lon_deg, 
                        lat = cell_centers$lat_deg)
    
    
    cells <- cc_df$cell
    if (length(cells) > 0)
    out <- c(out, list(cells))
  } else {
    print('not for this species')
  }
}

setwd(paste0(dir, 'Bird_Phenology/Data'))
saveRDS(out_old, 'new-process.rds')

