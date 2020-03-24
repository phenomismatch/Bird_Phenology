
#dir <- '~/Google_Drive/R/'
dir <- '/labs/Tingley/phenomismatch/'

# create grid -------------------------------------------------------------

hexgrid6 <- dggridR::dgconstruct(res = 6)

#get boundaries of all cells over earth
setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/shapefiles/'))

#dggridR::dgearthgrid(hexgrid6, savegrid = 'global_hex.shp')
#read in grid
hge <- rgdal::readOGR('global_hex.shp', verbose = FALSE)


# filter cells by range  ---------------------------------------------------

#reference key for species synonyms
setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/metadata/'))
sp_key <- read.csv('species_filenames_key.csv')



setwd(paste0(dir, 'Bird_Phenology/Data/'))

species_list_i <- read.table('eBird_species_list.txt', stringsAsFactors = FALSE)
species_list <- species_list_i[,1]
nsp <- length(species_list)



out_old <- list()
for (i in 1:nsp)
{
  #i <- 1
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
sp_rng <- rgdal::readOGR(fname, verbose = FALSE)
#crop to area of interest
sp_rng2 <- raster::crop(sp_rng, raster::extent(-95, -50, 24, 90))

#filter by breeding (2) and migration (4) range - need to convert spdf to sp
nrng <- sp_rng2[which(sp_rng2$SEASONAL == 2 | sp_rng2$SEASONAL == 4),]

#filter by resident (1) and non-breeding (3) to exclude hex cells that contain 2/4 and 1/3
nrng_rm <- sp_rng2[which(sp_rng2$SEASONAL == 1 | sp_rng2$SEASONAL == 3),]

#remove unneeded objects
rm(sp_rng)
rm(sp_rng2)
rm(fname)


#if there is a legitimate range
if (NROW(nrng@data) > 0)
{
  #good cells
  nrng_sp <- sp::SpatialPolygons(nrng@polygons)
  sp::proj4string(nrng_sp) <- sp::CRS(sp::proj4string(nrng))
  #find intersections with code from here: https://gis.stackexchange.com/questions/140504/extracting-intersection-areas-in-r
  poly_int <- rgeos::gIntersects(hge, nrng_sp, byid=TRUE)
  tpoly <- which(poly_int == TRUE, arr.ind = TRUE)[,2]
  br_mig_cells <- as.numeric(tpoly[!duplicated(tpoly)])
  
  #bad cells - also exclude cells 812, 813, and 841 (Bahamas)
  if (length(nrng_rm) > 0)
  {
    nrng_rm_sp <- sp::SpatialPolygons(nrng_rm@polygons)
    sp::proj4string(nrng_rm_sp) <- sp::CRS(sp::proj4string(nrng_rm))
    poly_int_rm <- rgeos::gIntersects(hge, nrng_rm_sp, byid=TRUE)
    tpoly_rm <- which(poly_int_rm == TRUE, arr.ind = TRUE)[,2]
    res_ovr_cells <- as.numeric(tpoly_rm[!duplicated(tpoly_rm)])
    
    #remove cells that appear in resident and overwinter range that also appear in breeding range
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
  
  
  #cells only within the range that ebird surveys were filtered to
  n_cc_df <- cc_df[which(cc_df$lon > -95 & cc_df$lon < -50 & cc_df$lat > 24),]
  cells <- n_cc_df$cell
  out_old <- c(out_old, list(cells))
  
} else {
  print('not for this species')
}
}

setwd(paste0(dir, 'Data'))
saveRDS(out, 'old-process.rds')
