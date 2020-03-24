######################
# 3 - proces arrival GAM output
#
# Aggeregate posterior info and diagnostic info from 2-halfmax-arr.R to be used in IAR model
# Determine which cells should be used in IAR model (cells that overlap breeding and migratory ranges, but do not overlap resident or non-breeding ranges) - 18 hour run time?
#
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/labs/Tingley/phenomismatch/'



# db/hm query dir ------------------------------------------------------------

hm_dir <- 'halfmax_arrival_2020-02-26'
hm_dir <- 'halfmax_species_2019-05-03'
hm_date <- substr(hm_dir, start = 17, stop = 26)


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

setwd(paste0(dir, 'Bird_Phenology/Data/hex_grid_crop'))

#read in cropped hex shp file
hex_shp <- rgdal::readOGR('hex_grid_crop.shp')

#load maps
usamap <- data.frame(maps::map("world", "USA", plot = FALSE)[c("x", "y")])
canadamap <- data.frame(maps::map("world", "Canada", plot = FALSE)[c("x", "y")])
mexicomap <- data.frame(maps::map("world", "Mexico", plot = FALSE)[c("x", "y")])

#combine maps
combmap <- maps::map("world", c("USA", "Canada", "Mexico"), fill = TRUE, plot = FALSE)

#make maps into spatial polygon object
IDs <- sapply(strsplit(combmap$names, ":"), function(x) x[1])
cmap <- maptools::map2SpatialPolygons(combmap, IDs = IDs, 
                                      proj4string = CRS(proj4string(hex_shp)))

#convert map to equal area projection
cmap2 <- sp::spTransform(cmap, sp::CRS("+proj=laea +lat_0=0 +lon_0=-70 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

#calculate overlap between each hex cell and map
hex_shp_cells <- as.numeric(as.character(hex_shp@data$cell))
ovr_df <- data.frame(cell = hex_shp_cells, per_ovr = NA)
for (i in 1:length(hex_shp))
{
  #just one hex cell
  hex_spoly <- sp::SpatialPolygons(list(hex_shp@polygons[[i]]))
  proj4string(hex_spoly) <- sp::CRS(proj4string(hex_shp))
  
  #tranform to equal area projection
  hex_spoly2 <- sp::spTransform(hex_spoly, sp::CRS(proj4string(cmap2)))
  
  #plot(hex_spoly2, add = TRUE, col = 'red')
  
  #area of hex cell (should all be the same)
  hex_area <- rgeos::gArea(hex_spoly2)
  #get intersection of hex cell and map
  int <- rgeos::gIntersection(hex_spoly2, cmap2)
  if (!is.null(int))
  {
    #get area of intersection
    ovr_area <- rgeos::gArea(int)
  } else {
    ovr_area <- 0
  }
  
  #percent of hex cell that is covered by land
  ovr_df$per_ovr[i] <- ovr_area / hex_area
}



# import eBird species list -----------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/'))

species_list_i <- read.table('eBird_species_list.txt', stringsAsFactors = FALSE)
species_list <- species_list_i[,1]
nsp <- length(species_list)
years <- 2002:2018
nyr <- length(years)


# Proces halfmax results -----------------------------------------------------------------

#get number of cell/years
tu_cells <- rep(NA, nsp)
for (i in 1:nsp)
{
  #i <- 20
  
  #import halfmax estimates and diagnostics from logit cubic model
  setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', hm_dir))
  temp_halfmax <- readRDS(paste0('halfmax_arrival_', species_list[i], '.rds'))

  tu_cells[i] <- length(unique(temp_halfmax$cell))
}

'%ni%' <- Negate('%in%')

#reference key for species synonyms for ranges
setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/metadata/'))
sp_key <- read.csv('species_filenames_key.csv')


counter <- 1
for (i in 1:nsp)
{
  #i <- 20
  
  #import halfmax estimates and diagnostics from GAM
  setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', hm_dir))
  temp_halfmax <- readRDS(paste0('halfmax_arrival_', species_list[i], '.rds'))
  
  if (i == 1)
  {
    na_reps <- rep(NA, max(tu_cells) * nyr * nsp)
    
    diagnostics_frame <- data.frame(species = na_reps,
                                    year = na_reps,
                                    cell = na_reps,
                                    HM_mean = na_reps,
                                    HM_sd = na_reps,
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
                                    n1W = na_reps,
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
  #to avoid self-intersection problem
  sp_rng2 <- rgeos::gBuffer(sp_rng, byid = TRUE, width = 0)
  #crop to area of interest
  sp_rng3 <- raster::crop(sp_rng2, raster::extent(-95, -50, 24, 90))
  
  #filter by breeding (2) and migration (4) range - need to convert spdf to sp
  nrng <- sp_rng3[which(sp_rng3$SEASONAL == 2 | sp_rng3$SEASONAL == 4),]
  
  #filter by resident (1) and non-breeding (3) to exclude hex cells that contain 2/4 and 1/3
  nrng_rm <- sp_rng3[which(sp_rng3$SEASONAL == 1 | sp_rng3$SEASONAL == 3),]

  #remove unneeded objects
  rm(sp_rng)
  rm(sp_rng2)
  rm(sp_rng3)

  #if there is a legitimate range
  if (NROW(nrng@data) > 0)
  {
    #good cells
    nrng_sp <- sp::SpatialPolygons(nrng@polygons)
    sp::proj4string(nrng_sp) <- sp::CRS(sp::proj4string(nrng))
    #find intersections with code from here: https://gis.stackexchange.com/questions/140504/extracting-intersection-areas-in-r
    poly_int <- rgeos::gIntersects(hex_shp, nrng_sp, byid = TRUE)
    tpoly <- which(poly_int == TRUE, arr.ind = TRUE)[,2]
    br_mig_cells <- hex_shp_cells[as.numeric(tpoly[!duplicated(tpoly)])]
  
    #bad cells - also exclude cells 812, 813, and 841 (Bahamas)
    if (length(nrng_rm) > 0)
    {
      nrng_rm_sp <- sp::SpatialPolygons(nrng_rm@polygons)
      sp::proj4string(nrng_rm_sp) <- sp::CRS(sp::proj4string(nrng_rm))
      poly_int_rm <- rgeos::gIntersects(hex_shp, nrng_rm_sp, byid = TRUE)
      tpoly_rm <- which(poly_int_rm == TRUE, arr.ind = TRUE)[,2]
      res_ovr_cells <- hex_shp_cells[as.numeric(tpoly_rm[!duplicated(tpoly_rm)])]
      
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
    
    #remove unneeded objects
    rm(nrng)
    rm(nrng_sp)
    
    c_rm <- which(br_mig_cells %in% to_rm)
    if (length(c_rm) > 0)
    {
      overlap_cells <- br_mig_cells[-c_rm]  
    } else {
      overlap_cells <- br_mig_cells
    }
    
    #remove unneeded objects
    rm(br_mig_cells)
    rm(cell_mrg)
    rm(to_rm)
    
    #get cell centers
    cell_centers <- dggridR::dgSEQNUM_to_GEO(hexgrid6, overlap_cells)
    cc_df <- data.frame(cell = overlap_cells, lon = cell_centers$lon_deg, 
                        lat = cell_centers$lat_deg)
    
    #remove unneeded objects
    rm(cell_centers)
    rm(overlap_cells)
    
    #cells only within the range that ebird surveys were filtered to
    n_cc_df <- cc_df[which(cc_df$lon > -95 & cc_df$lon < -50 & cc_df$lat > 24),]
    cells <- n_cc_df$cell
    ncell <- length(cells)
    
    #remove unneeded objects
    rm(cc_df)
    rm(n_cc_df)
  } else {
    ncell <- NA
  } 
  
  #loop through years
  for (j in 1:nyr)
  {
    #j <- 1
    tt_halfmax1 <- dplyr::filter(temp_halfmax, 
                                 year == years[j])
    if (ncell > 0)
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
          diagnostics_frame$n1W[counter] <- tt_halfmax2$n1W
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
            diagnostics_frame$HM_mean[counter] <- mean(halfmax_posterior)
            diagnostics_frame$HM_sd[counter] <- sd(halfmax_posterior)
          }
        }
        counter <- counter + 1
      } # if loop for at least one cell - species without sufficient ranges have 0 cells
    } # k -cell
  } # j - year
} # i - species


#strip excess NAs


#ADD HERE



#add 'meets criteria' column
diagnostics_frame$MODEL <- NA
diagnostics_frame$shp_fname <- NA

diagnostics_frame2 <- diagnostics_frame


### add NA for both HM_mean and HM_sd if any of the following conditions are met
#these params used for most recent IAR input run (as of 2020-02-XX)

#filter by mlmax and plmax?
to.NA <- which(diagnostics_frame2$num_diverge > 0 | 
                 diagnostics_frame2$max_Rhat >= 1.05 |
                 diagnostics_frame2$min_neff < 350 |
                 diagnostics_frame2$num_BFMI > 0 |
                 diagnostics_frame2$HM_sd > 15)

# #XX% of cells are bad
# length(to.NA)/sum(!is.na(diagnostics_frame2$HM_mean))
# diagnostics_frame2[to.NA,c('species', 'cell', 'year',
#                            'HM_mean', 'HM_sd', 'min_neff', 'num_diverge')]

if (length(to.NA) > 0)
{
  diagnostics_frame2[to.NA,'HM_mean'] <- NA
  diagnostics_frame2[to.NA,'HM_sd'] <- NA
}




# combine data with overlap df --------------------------------------------

diagnostics_frame3 <- dplyr::left_join(diagnostics_frame2, ovr_df, by = 'cell')



# Filter data based on criteria -----------------------------------------------------------

# Which species-years are worth modeling? At a bare minimum:
#     Species with at least 'NC' cells in all three years from 2016-2018
#     Species-years with at least 'NC' cells for those species

NC <- 3

df_out <- data.frame()
#which species/years meet criteria for model
for (i in 1:length(species_list))
{
    #i <- 1
    
    print(i)
    #filter by species
    t_sp <- dplyr::filter(diagnostics_frame3, species == species_list[i])
    
    if (NROW(t_sp) > 0)
    {
      #number of cells with good data in each year from 2016-2018
      nobs_yr <- c()
      for (j in 2016:2018)
      {
        #j <- 2018
        ty_sp3 <- dplyr::filter(t_sp, year == j)
        ind <- which(!is.na(ty_sp3$HM_mean))
        nobs_yr <- c(nobs_yr, length(ind))
        #ty_sp[ind,]
      }

      #if all three years have greater than or = to 'NC' cells of data, figure 
      #...out which years have at least 'NC' cells
      yrs_kp <- c()
      if (sum(nobs_yr >= NC) == 3)
      {
        #see which years have more than 3 cells of data
        nobs_yr2 <- c()
        for (j in min(years):max(years))
        {
          #j <- 2012
          ty2_sp3 <- dplyr::filter(t_sp, year == j)
          ind2 <- which(!is.na(ty2_sp3$HM_mean))
          nobs_yr2 <- c(nobs_yr2, length(ind2))
        }
      
        #years to keep (more than three cells of data)
        yrs_kp <- years[which(nobs_yr2 >= NC)]
      }
    
      if (length(yrs_kp) > 0)
      {
        t_sp[which(t_sp$year %in% yrs_kp),]$MODEL <- TRUE
      }
    
      df_out <- rbind(df_out, t_sp)
    }
}
run_time <- (proc.time()[3] - tt[3]) / 60



# order -------------------------------------------------------------------

#order diagnostics frame by species, year, and cell #
df_master <- df_out[with(df_out, order(species, year, cell)),]


# write to RDS --------------------------------------------------

IAR_dir_path <- paste0(dir, 'Bird_phenology/Data/Processed/IAR_input_', hm_date)

dir.create(IAR_dir_path)
setwd(IAR_dir_path)

saveRDS(df_master, paste0('IAR_input-', hm_date, '.rds'))



# create list of species to run through IAR model -------------------------

species_tm <- aggregate(MODEL ~ species, data = df_master, FUN = function(x) sum(x, na.rm = TRUE))$species

setwd(paste0(dir, 'Bird_Phenology/Data/'))
write.table(species_tm, file = paste0('IAR_species_list.txt'), row.names = FALSE, col.names = FALSE)




# # create dfs that show # cells with data in each year/species, and # years with data in each cell/species -----------------
# 
# #create df with species/cells/n_yrs per sp_cell
# cells_frame <- data.frame(species = rep(NA, cell_years), 
#                           cell = rep(NA, cell_years), 
#                           n_yrs = rep(NA, cell_years))
# 
# yrs_frame <- data.frame(species = rep(NA, cell_years), 
#                         year = rep(NA, cell_years), 
#                         n_cells = rep(NA, cell_years))
# 
# counter_cell <- 1
# counter_year <- 1
# #fill cells_frame and yrs_frame
# for (i in 1:nsp)
# {
#   #i <- 101
#   print(i)
#   
#   tspf <- dplyr::filter(df_master, species == species_list[i])
#   tcells <- unique(tspf$cell)
#   tyears <- unique(tspf$year)
#   
#   if (NROW(tspf) > 0)
#   {
#     for (k in 1:length(tcells))
#     {
#       #k <- 1
#       t_cell <- dplyr::filter(tspf, cell == tcells[k])
#     
#       cells_frame[counter_cell, 'species'] <- species_list[i]
#       cells_frame[counter_cell, 'cell'] <- tcells[k]
# 
#       #insert number of yrs with data
#       yrs_d <- t_cell$year[which(!is.na(t_cell$HM_mean))]
#       cells_frame[counter_cell,'n_yrs'] <- length(yrs_d)
#       counter_cell <- counter_cell + 1
#     }
#   
#     for (j in 1:length(tyears))
#     {
#       #j <- 1
#       t_year <- dplyr::filter(tspf, year == tyears[j])
#     
#       yrs_frame[counter_year, 'species'] <- species_list[i]
#       yrs_frame[counter_year, 'year'] <- tyears[j]
#     
#       #insert number of yrs with data
#       yrs_d <- t_cell$year[which(!is.na(t_year$HM_mean))]
#       yrs_frame[counter_year,'n_cells'] <- length(yrs_d)
#       counter_year <- counter_year + 1
#     }
#   }
# }
# 
# 
# #remove NA from end of df
# to.rm.cell <- min(which(is.na(cells_frame$species))):NROW(cells_frame)
# to.rm.year <- min(which(is.na(yrs_frame$species))):NROW(yrs_frame)
# 
# cells_frame2 <- cells_frame[-to.rm.cell,]
# yrs_frame2 <- yrs_frame[-to.rm.year,]
# 
# #write to RDS
# setwd(IAR_dir_path)
# saveRDS(cells_frame2, paste0('cells_frame-', hm_date, '.rds'))
# saveRDS(yrs_frame2, paste0('yrs_frame-', hm_date, '.rds'))



# explore data ------------------------------------------------------------

# aggregate(n_cells ~ species, data = yrs_frame, FUN = max)
# aggregate(n_cells ~ species, data = yrs_frame, FUN = mean)


