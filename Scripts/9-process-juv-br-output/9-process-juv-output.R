######################
# 9 - process juv GAM output
#
# Aggregate posterior info and diagnostic info from 7-halfmax-juvs.R
#
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/labs/Tingley/phenomismatch/'


# db/juv query dir ------------------------------------------------------------

#input dir
juv_date <- '2020-03-09'
juv_dir <- paste0(dir, 'Bird_Phenology/Data/Processed/halfmax_juvs_', juv_date)

#output dir
juv_master_dir <- paste0(dir, 'Bird_Phenology/Data/Processed/juv_master_', juv_date)


# runtime -----------------------------------------------------------------

tt <- proc.time()


# Load packages -----------------------------------------------------------

library(dplyr)
library(dggridR)
library(sp)
library(maptools)
library(rgdal)
library(rgeos)


# calculate grid/land overlap ---------------------------------------------

#make hexgrid
hexgrid6 <- dggridR::dgconstruct(res = 6)

setwd(paste0(dir, 'Bird_Phenology/Data/hex_grid_crop/'))

#read in grid
hge <- rgdal::readOGR('hex_grid_crop.shp', verbose = FALSE)
hge_cells <- as.numeric(as.character(hge@data[,1]))

#load maps
usamap <- data.frame(maps::map("world", "USA", plot = FALSE)[c("x", "y")])
canadamap <- data.frame(maps::map("world", "Canada", plot = FALSE)[c("x", "y")])
mexicomap <- data.frame(maps::map("world", "Mexico", plot = FALSE)[c("x", "y")])

#combine maps
combmap <- maps::map("world", c("USA", "Canada", "Mexico"), fill = TRUE, plot = FALSE)

#make maps into spatial polygon object
IDs <- sapply(strsplit(combmap$names, ":"), function(x) x[1])
cmap <- maptools::map2SpatialPolygons(combmap, IDs = IDs, 
                                      proj4string = CRS(proj4string(hge)))

#convert map to equal area projection
cmap2 <- sp::spTransform(cmap, sp::CRS("+proj=laea +lat_0=0 +lon_0=-70 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

#calculate overlap between each hex cell and map
ovr_df <- data.frame(cell = hge_cells, per_ovr = NA)
for (i in 2:length(hge))
{
  #just one hex cell
  hge_spoly <- sp::SpatialPolygons(list(hge@polygons[[i]]))
  proj4string(hge_spoly) <- sp::CRS(proj4string(hge))
  
  #tranform to equal area projection
  hge_spoly2 <- sp::spTransform(hge_spoly, sp::CRS(proj4string(cmap2)))
  
  #area of hex cell (should all be the same)
  hex_area <- rgeos::gArea(hge_spoly2)
  #get intersection of hex cell and map
  int <- rgeos::gIntersection(hge_spoly2, cmap2)
  if (!is.null(int))
  {
    #get area of intersection
    ovr_area <- rgeos::gArea(int)
  } else {
    ovr_area <- 0
  }
  
  #percent of hex cell that is covered by land
  ovr_df$per_ovr[i] <- round(ovr_area / hex_area, 2)
}


# import eBird species list -----------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/'))

species_list_i <- read.table('eBird_species_list.txt', stringsAsFactors = FALSE)
species_list <- species_list_i[,1]
nsp <- length(species_list)
#years for MAPS data
years <- 1989:2018
nyr <- length(years)


# Proces halfmax results -----------------------------------------------------------------

#get number of cell/years
tu_cells <- rep(NA, nsp)
for (i in 1:nsp)
{
  #i <- 20
  
  #import halfmax estimates and diagnostics from logit cubic model
  setwd(juv_dir)
  temp_halfmax <- readRDS(paste0('halfmax_juvs_', species_list[i], '.rds'))
   
  tu_cells[i] <- length(unique(temp_halfmax$cell))
}


#reference key for species synonyms for ranges
setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/metadata/'))
sp_key <- read.csv('species_filenames_key.csv')


counter <- 1
for (i in 1:nsp)
{
  #i <- 23
  
  #import halfmax estimates and diagnostics from GAM
  setwd(juv_dir)
  temp_halfmax <- readRDS(paste0('halfmax_juvs_', species_list[i], '.rds'))
  
  if (i == 1)
  {
    na_reps <- rep(NA, max(tu_cells) * nyr * nsp)
    
    diagnostics_frame <- data.frame(species = na_reps,
                                    year = na_reps,
                                    cell = na_reps,
                                    juv_mean = na_reps,
                                    juv_sd = na_reps,
                                    max_Rhat = na_reps,
                                    min_neff = na_reps,
                                    mlmax = na_reps,
                                    plmax = na_reps,
                                    num_diverge = na_reps,
                                    num_tree = na_reps,
                                    num_BFMI = na_reps,
                                    delta = na_reps,
                                    tree_depth = na_reps,
                                    t_iter = na_reps,
                                    n1 = na_reps,
                                    #n1W = na_reps,
                                    n0 = na_reps,
                                    n0i = na_reps,
                                    njd = na_reps,
                                    njd1 = na_reps,
                                    njd0 = na_reps,
                                    njd0i = na_reps)
  }
  
  #filter by breeding/migration cells
  #match species name to shp file name
  g_ind <- grep(species_list[i], sp_key$file_names_2016)
  
  #check for synonyms if there are no matches
  if (length(g_ind) == 0)
  {
    g_ind2 <- grep(species_list[i], sp_key$BL_Checklist_name)
  } else {
    g_ind2 <- g_ind
  }
  
  #change dir to shp files
  setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/shapefiles/'))
  #get filename and read in
  fname <- as.character(sp_key[g_ind2,]$filenames[grep('.shp', sp_key[g_ind2, 'filenames'])])
  #read in shp file
  sp_rng <- rgdal::readOGR(fname[1], verbose = FALSE)
  
  #filter by breeding (2) range - need to convert spdf to sp
  nrng <- sp_rng[which(sp_rng$SEASONAL == 2),]
  
  #filter by resident (1), non-breeding (3), and migratory (4) ranges
  nrng_rm <- sp_rng[which(sp_rng$SEASONAL == 1 | sp_rng$SEASONAL == 3 | sp_rng$SEASONAL == 4),]
  
  #remove unneeded objects
  rm(sp_rng)
  
  #if there is a legitimate range
  if (NROW(nrng@data) > 0)
  {
    #good cells
    nrng_sp <- sp::SpatialPolygons(nrng@polygons)
    sp::proj4string(nrng_sp) <- sp::CRS(sp::proj4string(nrng))
    #find intersections with code from here: https://gis.stackexchange.com/questions/140504/extracting-intersection-areas-in-r
    poly_int <- rgeos::gIntersects(hge, nrng_sp, byid=TRUE)
    tpoly <- which(poly_int == TRUE, arr.ind = TRUE)[,2]
    br_cells <- hge_cells[as.numeric(tpoly[!duplicated(tpoly)])]
    
    #bad cells
    if (length(nrng_rm) > 0)
    {
      nrng_rm_sp <- sp::SpatialPolygons(nrng_rm@polygons)
      sp::proj4string(nrng_rm_sp) <- sp::CRS(sp::proj4string(nrng_rm))
      poly_int_rm <- rgeos::gIntersects(hge, nrng_rm_sp, byid=TRUE)
      tpoly_rm <- which(poly_int_rm == TRUE, arr.ind = TRUE)[,2]
      res_ovr_cells <- hge_cells[as.numeric(tpoly_rm[!duplicated(tpoly_rm)])]
      
      #remove cells that appear in resident and overwinter range that also appear in breeding range
      cell_mrg <- c(br_cells, res_ovr_cells)
      to_rm <- c(cell_mrg[duplicated(cell_mrg)])
      
      rm(nrng_rm)
      rm(nrng_rm_sp)
      rm(res_ovr_cells)
      
    } else {
      cell_mrg <- br_cells
      to_rm <- NA
    }
    
    #remove unneeded objects
    rm(nrng)
    rm(nrng_sp)
    
    c_rm <- which(br_cells %in% to_rm)
    if (length(c_rm) > 0)
    {
      overlap_cells <- br_cells[-c_rm]  
    } else {
      overlap_cells <- br_cells
    }
    
    #remove unneeded objects
    rm(br_cells)
    rm(cell_mrg)
    rm(to_rm)
    
    #get cell centers
    cell_centers <- dggridR::dgSEQNUM_to_GEO(hexgrid6, overlap_cells)
    cc_df <- data.frame(cell = overlap_cells, lon = cell_centers$lon_deg, 
                        lat = cell_centers$lat_deg)
    
    #remove unneeded objects
    rm(cell_centers)
    rm(overlap_cells)
    
    cells <- cc_df$cell
    ncell <- length(cells)
    
    #remove unneeded objects
    rm(cc_df)
  } else {
    ncell <- NA
  } 
  
  #loop through years
  for (j in 1:length(years))
  {
    #j <- 1
    tt_halfmax1 <- dplyr::filter(temp_halfmax, 
                                 year == years[j])
    
    if (!is.na(ncell))
    {
      for (k in 1:ncell)
      {
        #k <- 1
        
        diagnostics_frame$species[counter] <- species_list[i]
        diagnostics_frame$year[counter] <- years[j]
        diagnostics_frame$cell[counter] <- cells[k]
        diagnostics_frame$shp_fname[counter] <- fname
        
        #get model fits
        tt_halfmax2 <- dplyr::filter(tt_halfmax1, 
                             cell == cells[k])
        
        print(paste0('species: ', species_list[i], ', ',
                     'year: ', years[j], ', ', 
                     'cell: ', cells[k]))

        if (NROW(tt_halfmax2) > 0)
        {
          diagnostics_frame$min_neff[counter] <- tt_halfmax2$min_neff
          diagnostics_frame$max_Rhat[counter] <- tt_halfmax2$max_Rhat
          #did mean have local max
          diagnostics_frame$mlmax[counter] <- tt_halfmax2$mlmax
          #percent of iter with local max
          diagnostics_frame$plmax[counter] <- tt_halfmax2$plmax
          diagnostics_frame$num_diverge[counter] <- tt_halfmax2$num_diverge
          diagnostics_frame$num_tree[counter] <- tt_halfmax2$num_tree
          diagnostics_frame$num_BFMI[counter] <- tt_halfmax2$num_BFMI
          diagnostics_frame$delta[counter] <- tt_halfmax2$delta
          diagnostics_frame$tree_depth[counter] <- tt_halfmax2$tree_depth
          #total iterations
          diagnostics_frame$t_iter[counter] <- tt_halfmax2$t_iter
          #number of surveys where species was detected
          diagnostics_frame$n1[counter] <- tt_halfmax2$n1
          #number of detections that came before jday 60
          #diagnostics_frame$n1W[counter] <- tt_halfmax2$n1W
          #number of surveys where species was not detected
          diagnostics_frame$n0[counter] <- tt_halfmax2$n0
          #number of non-detections before first detection
          diagnostics_frame$n0i[counter] <- tt_halfmax2$n0i
          #number of julian days with obs
          diagnostics_frame$njd[counter] <- tt_halfmax2$njd
          #number of unique days with detections
          diagnostics_frame$njd1[counter] <- tt_halfmax2$njd1
          #number of unique days with non-detection
          diagnostics_frame$njd0[counter] <- tt_halfmax2$njd0
          #number of unique days of non-detections before first detection
          diagnostics_frame$njd0i[counter] <- tt_halfmax2$njd0i
          
          #posterior for halfmax
          cndf <- colnames(tt_halfmax2)
          iter_ind <- grep('iter', cndf)
          to.rm <- which(cndf == 't_iter')
          halfmax_posterior <- as.numeric(tt_halfmax2[,iter_ind[-which(iter_ind == to.rm)]])
          
          #calculate posterior mean and sd
          if (sum(!is.na(halfmax_posterior)) > 0)
          {
            diagnostics_frame$juv_mean[counter] <- mean(halfmax_posterior)
            diagnostics_frame$juv_sd[counter] <- sd(halfmax_posterior)
          }
        }
        counter <- counter + 1
      } # if loop for at least one cell - species without sufficient ranges have 0 cells
    } # k -cell
  } # j - year
} # i - species


# strip excess NAs --------------------------------------------------------

to.rm <- min(which(is.na(diagnostics_frame$species))):NROW(diagnostics_frame)
diagnostics_frame2 <- diagnostics_frame[-to.rm,]


# filter bad results ------------------------------------------------------

### add NA for both juv_mean and juv_sd if any of the following conditions are met

to.NA <- which(diagnostics_frame2$num_diverge > 0 | 
                 diagnostics_frame2$max_Rhat >= 1.05 |
                 diagnostics_frame2$min_neff < 350 |
                 diagnostics_frame2$num_BFMI > 0 |
                 diagnostics_frame2$juv_sd > 15 | 
                 diagnostics_frame2$mlmax == FALSE) #must have hump

# #53% of cells are bad
# length(to.NA)/sum(!is.na(diagnostics_frame2$juv_mean))
# diagnostics_frame2[to.NA,c('species', 'cell', 'year',
#                           'juv_mean', 'juv_sd', 'min_neff',
#                           'num_diverge', 'max_Rhat', 'mlmax')]

if (length(to.NA) > 0)
{
  diagnostics_frame2[to.NA,'juv_mean'] <- NA
  diagnostics_frame2[to.NA,'juv_sd'] <- NA
}


# combine data with overlap df --------------------------------------------

diagnostics_frame3 <- dplyr::left_join(diagnostics_frame2, ovr_df, by = 'cell')


# order -------------------------------------------------------------------

#order diagnostics frame by species, year, and cell #
df_master <- diagnostics_frame3[with(diagnostics_frame3, order(species, year, cell)),]


# add cell lat/lon --------------------------------------------------------

#get hexgrid cell centers
cellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid6, df_master$cell)

df_master$cell_lat <- round(cellcenters$lat_deg, digits = 2)
df_master$cell_lng <- round(cellcenters$lon_deg, digits = 2)


# write to RDS --------------------------------------------------

dir.create(juv_master_dir)
setwd(juv_master_dir)

saveRDS(df_master, paste0('juv_master_', juv_date, '2.rds'))
