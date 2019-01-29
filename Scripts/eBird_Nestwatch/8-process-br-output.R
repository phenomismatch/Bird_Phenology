######################
# 8 - Process breeding code logit cubic output
#
# Compiles dataframe with breeding data from Nestwatch (XXX), eBird breeding codes (7-logit-cubic-breeding.R), and MAPS observations (6-query-nesting-data.R) for every year/cell
######################

#filter cells by only those in breeding range (those modeled in 4-IAR [regardles whether there were data in those cells])



# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/home/CAM/cyoungflesh/phenomismatch/'



# db/hm query dir ------------------------------------------------------------

ebird_date <- '2019-01-16'
bc_query_date <- '2019-01-15'
NW_date <- '2019-01-28'
MAPS_date <- '2019-01-28'
db_query_dir <- 'db_query_2018-10-15'


# Load packages -----------------------------------------------------------

library(dplyr)
library(dggridR)


# Set wd ------------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/'))



# import eBird species list -----------------------------------------------------

species_list_i <- read.table('IAR_species_list.txt', stringsAsFactors = FALSE)
species_list <- species_list_i[,1]
nsp <- length(species_list)


#DATA ONLY VALID THROUGH 2017 (2018 data only goes to ~ jday 60 as of 2018-10-15 query)
years <- 2002:2017
nyr <- length(years)


# create grid -------------------------------------------------------------

#construct grid
hexgrid6 <- dggridR::dgconstruct(res = 6)



# read in NW and MAPS data ------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))

NW_data <- readRDS(paste0('breeding_NW_', NW_date, '.rds'))
MAPS_data <- readRDS(paste0('breeding_MAPS_obs_', MAPS_date, '.rds'))


# create and fill df ------------------------------------------------------

#read in ebird surveys to get full list of cells
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', db_query_dir))
t_data <- readRDS(paste0('ebird_NA_phen_proc_', species_list[1], '.rds'))

na_reps <- rep(NA, (nsp*nyr*length(unique(t_data$cell6))))
rm(t_data)


m_breeding_df <- data.frame(SPECIES = na_reps,
                            CELL = na_reps,
                            YEAR = na_reps,
                            EB_HM_mean = na_reps,
                            EB_HM_sd = na_reps,
                            EB_HM_LCI = na_reps,
                            EB_HM_UCI = na_reps,
                            EB_n1 = na_reps,
                            EB_n1W = na_reps,
                            EB_n0 = na_reps,
                            EB_n0i = na_reps,
                            EB_njd1 = na_reps,
                            EB_njd0 = na_reps,
                            EB_njd0i = na_reps,
                            EB_max_Rhat = na_reps,
                            EB_min_neff = na_reps,
                            EB_sh_pv = na_reps,
                            EB_nphen_bad = na_reps,
                            NW_mean_cid = na_reps,
                            NW_sd_cid = na_reps,
                            NW_num_obs = na_reps,
                            MAPS_midpoint = na_reps,
                            MAPS_l_bounds = na_reps,
                            MAPS_u_bounds = na_reps,
                            MAPS_n_stations = na_reps)

counter <- 1
for (i in 1:nsp)
{
  #i <- 1
  
  #readin halfmax data
  setwd(paste0(dir, 'Bird_Phenology/Data/Processed/halfmax_breeding_', ebird_date))
  temp_halfmax <- readRDS(paste0('halfmax_df_breeding_', species_list[i], '.rds'))
  
  setwd(paste0(dir, 'Bird_Phenology/Data/Processed/breeding_cat_query_', bc_query_date))
  temp_bc <- readRDS(paste0('ebird_NA_breeding_cat_', species_list[i], '.rds'))
  
  cells <- unique(temp_halfmax$cell)
  years <- unique(temp_halfmax$year)
  nyr <- length(years)
  ncell <- length(cells)
  
  #loop through years
  for (j in 1:nyr)
  {
    #j <- 1
    
    for (k in 1:ncell)
    {
      #k <- 1
      print(paste0('species: ', species_list[i], ', ',
                   'year: ', years[j], ', ',
                   'cell: ', cells[k]))
      
      #####################
      #ebird breeding codes
      #####################
      
      cysdata <- dplyr::filter(temp_bc, 
                               year == years[j] & 
                               cell == cells[k])
      
      tt_halfmax <- filter(temp_halfmax, 
                           year == years[j], 
                           cell == cells[k])
      
      #write species/cell/year data
      m_breeding_df$SPECIES[counter] <- species_list[i]
      m_breeding_df$CELL[counter] <- cells[k]
      m_breeding_df$YEAR[counter] <- years[j]
        
      #add column with 1/0 breeding or not
      cysdata$br <- as.numeric(cysdata$bba_category == 'C3' | cysdata$bba_category == 'C4')
      
      #bird not seen - fill with zeros
      na.ind <- which(is.na(cysdata$br))
      cysdata$br[na.ind] <- 0
      
      #number of surveys where breeding was detected (confirmed or probable)
      m_breeding_df$EB_n1[counter] <- sum(cysdata$br)
      #number of surveys where breeding was not detected (bird not seen breeding or not seen)
      m_breeding_df$EB_n0[counter] <- sum(cysdata$br == 0)
      #number of detections that came before jday 60
      m_breeding_df$EB_n1W[counter] <- sum(cysdata$br * as.numeric(cysdata$day < 60))
      #number of unique days with detections
      m_breeding_df$EB_njd1[counter] <- length(unique(cysdata$day[which(cysdata$br == 1)]))
      #number of unique days with non-detection
      m_breeding_df$EB_njd0[counter] <- length(unique(cysdata$day[which(cysdata$br == 0)]))
      
      
      if (m_breeding_df$EB_n1[counter] > 0)
      {
        #number of unique days of non-detections before first detection
        m_breeding_df$EB_njd0i[counter] <- length(unique(cysdata$day[which(cysdata$br == 0 & cysdata$day < 
                                                                        min(cysdata$day[which(cysdata$br == 1)]))]))
        #number of non-detections before first detection
        m_breeding_df$EB_n0i[counter] <- length(which(cysdata$br == 0 & 
                                                         cysdata$day < min(cysdata$day[which(cysdata$br == 1)])))
      } else {
        m_breeding_df$EB_njd0i[counter] <- 0
        m_breeding_df$EB_n0i[counter] <- 0
      }
      
      m_breeding_df$EB_min_neff[counter] <- tt_halfmax$min_neff
      m_breeding_df$EB_max_Rhat[counter] <- tt_halfmax$max_Rhat
      m_breeding_df$EB_sh_pv[counter] <- tt_halfmax$sh
        
      iter_ind <- grep('iter', colnames(tt_halfmax))
      halfmax_posterior <- as.vector(tt_halfmax[,iter_ind])
      
      #determine how many estimates are 1 and not 1 (estimates of 1 are bogus)
      m_breeding_df$EB_nphen_bad[counter] <- sum(halfmax_posterior == 1)

      #calculate posterior mean and sd
      if (sum(!is.na(halfmax_posterior)) > 0)
      {
        m_breeding_df$EB_HM_mean[counter] <- mean(halfmax_posterior)
        m_breeding_df$EB_HM_sd[counter] <- sd(halfmax_posterior)
        m_breeding_df$EB_HM_LCI[counter] <- quantile(halfmax_posterior, probs = 0.025)
        m_breeding_df$EB_HM_UCI[counter] <- quantile(halfmax_posterior, probs = 0.975)
      }
      
      
      ########
      #NW data
      ########
      
      t_NW <- dplyr::filter(NW_data, 
                            CELL == cells[k],
                            YEAR == years[j])
      
      if (NROW(t_NW) > 0)
      {
        m_breeding_df$NW_mean_cid[counter] <- t_NW$MEAN_FIRST_LAY
        m_breeding_df$NW_sd_cid[counter] <- t_NW$SD_FIRST_LAY
        m_breeding_df$NW_num_obs[counter] <- t_NW$NUM_OBS
      }
      
      #########
      #MAPS obs
      #########

      t_MAPS <- dplyr::filter(MAPS_data, 
                            cell == cells[k],
                            YR == years[j])
      
      if (NROW(t_MAPS) > 0)
      {
        m_breeding_df$MAPS_midpoint[counter] <- t_MAPS$midpoint
        m_breeding_df$MAPS_l_bounds[counter] <- t_MAPS$l_bounds
        m_breeding_df$MAPS_u_bounds[counter] <- t_MAPS$u_bounds
        m_breeding_df$MAPS_n_stations[counter] <- t_MAPS$n_stations
      }
      
      counter <- counter + 1
    } # k -cell
  } # j - year
} # i - species








#vvv OLD vvv

#add 'meets criteria' column
diagnostics_frame$MODEL <- NA
diagnostics_frame$shp_fname <- NA


# Filter data based on criteria -----------------------------------------------------------

# Which species-years are worth modeling? At a bare minimum:
#     Species with at least 'NC' cells in all three years from 2015-2017
#     Species-years with at least 'NC' cells for those species

NC <- 3

'%ni%' <- Negate('%in%')

#reference key for species synonyms
setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/metadata/'))
sp_key <- read.csv('species_filenames_key.csv')

#change dir to shp files
setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/shapefiles/'))

df_out <- data.frame()
#which species/years meet criteria for model
for (i in 1:length(species_list))
{
  #i <- 16
  
  print(i)
  #filter by species
  t_sp <- dplyr::filter(diagnostics_frame, species == species_list[i])
  
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
  
  #get filename and read in
  fname <- as.character(sp_key[g_ind2,]$filenames[grep('.shp', sp_key[g_ind2, 'filenames'])])
  t_sp$shp_fname <- fname
  sp_rng <- rgdal::readOGR(fname, verbose = FALSE)
  
  #filter by breeding (2) and migration (4) range - need to convert spdf to sp
  nrng <- sp_rng[which(sp_rng$SEASONAL == 2 | sp_rng$SEASONAL == 4),]
  
  #filter by resident (1) and non-breeding (3) to exclude hex cells that contain 2/4 and 1/3
  nrng_rm <- sp_rng[which(sp_rng$SEASONAL == 1 | sp_rng$SEASONAL == 3),]
  
  #only process if there is a seasonal range
  if (NROW(nrng@data) > 0)
  {
    #good cells
    nrng_sp <- sp::SpatialPolygons(nrng@polygons)
    sp::proj4string(nrng_sp) <- sp::CRS(sp::proj4string(nrng))
    ptsreg <- sp::spsample(nrng, 50000, type = "regular")
    br_mig_cells <- as.numeric(which(!is.na(sp::over(hge, ptsreg))))
    
    #bad cells
    nrng_rm_sp <- sp::SpatialPolygons(nrng_rm@polygons)
    sp::proj4string(nrng_rm_sp) <- sp::CRS(sp::proj4string(nrng_rm))
    ptsreg_rm <- sp::spsample(nrng_rm_sp, 50000, type = "regular")
    res_ovr_cells <- as.numeric(which(!is.na(sp::over(hge, ptsreg_rm))))
    
    #remove cells that appear in resident and overwinter range that also appear in breeding range
    cell_mrg <- c(br_mig_cells, res_ovr_cells)
    to_rm <- cell_mrg[duplicated(cell_mrg)]
    
    if (length(to_rm) > 0)
    {
      overlap_cells <- br_mig_cells[-which(br_mig_cells %in% to_rm)]
    } else {
      overlap_cells <- br_mig_cells
    }
    
    #get cell centers
    cell_centers <- dggridR::dgSEQNUM_to_GEO(hexgrid6, overlap_cells)
    cc_df <- data.frame(cell = overlap_cells, lon = cell_centers$lon_deg, 
                        lat = cell_centers$lat_deg)
    
    #cells only within the range that ebird surveys were filtered to
    n_cc_df <- cc_df[which(cc_df$lon > -100 & cc_df$lon < -50 & cc_df$lat > 26),]
    cells <- n_cc_df$cell
    
    #retain rows that match selected cells
    t_sp2 <- t_sp[which(t_sp$cell %in% cells),]
    
    #create rows for cells that were missing in ebird data
    missing_cells <- cells[which(cells %ni% t_sp2$cell)]
    
    temp_dff <- t_sp2[1,]
    temp_dff[,2:20] <- NA
    
    nmc <- length(missing_cells)
    nreps <- nmc * nyr
    
    temp_dff2 <- temp_dff[rep(row.names(temp_dff), nreps),]
    rownames(temp_dff2) <- NULL
    
    temp_dff2$year <- rep(years, nmc)
    temp_dff2$cell <- rep(missing_cells, each = nyr)
    
    #combine filtered data with missing cells
    t_sp3 <- rbind(t_sp2, temp_dff2)
    
    
    #number of cells with good data in each year from 2015-2017
    nobs_yr <- c()
    for (j in 2015:2017)
    {
      #j <- 2017
      ty_sp3 <- dplyr::filter(t_sp3, year == j)
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
        ty2_sp3 <- dplyr::filter(t_sp3, year == j)
        ind2 <- which(!is.na(ty2_sp3$HM_mean))
        nobs_yr2 <- c(nobs_yr2, length(ind2))
      }
      
      #years to keep (more than three cells of data)
      yrs_kp <- years[which(nobs_yr2 >= NC)]
    }
    
    if (length(yrs_kp) > 0)
    {
      t_sp3[which(t_sp3$year %in% yrs_kp),]$MODEL <- TRUE
    }
    
    df_out <- rbind(df_out, t_sp3)
  } else {
    #merge unchanged data if there isn't a seasonal range
    df_out <- rbind(df_out, t_sp)
  }
}


# order -------------------------------------------------------------------

#order diagnostics frame by species, year, and cell #
df_master <- df_out[with(df_out, order(species, year, cell)),]

#create scaling factor column
df_master$scaling_factor <- NA


# add scaling factor to df ------------------------------------------------

for (k in 1:length(species_list))
{
  #filter by species
  f_in <- dplyr::filter(df_master, species == species_list[k])
  
  #filter by year
  f_out <- f_in[which(f_in$MODEL == TRUE),]
  
  #define cells and years to be modeled
  cells <- unique(f_out$cell)
  years <- unique(f_out$year)
  nyr <- length(years)
  ncel <- length(cells)
  
  # create adjacency matrix
  #make hexgrid
  hexgrid6 <- dggridR::dgconstruct(res = 6)
  
  #get hexgrid cell centers
  cellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid6, cells)
  
  #create adjacency matrix - 1 if adjacent to cell, 0 if not
  adjacency_matrix <- matrix(data = NA, nrow = length(cells), ncol = length(cells))
  
  for (i in 1:length(cells))
  {
    #i <- 1
    for (j in i:length(cells))
    {
      #j <- 4
      dists <- geosphere::distm(c(cellcenters$lon_deg[i], cellcenters$lat_deg[i]),
                                c(cellcenters$lon_deg[j], cellcenters$lat_deg[j]))
      adjacency_matrix[i,j] <- as.numeric((dists/1000) > 0 & (dists/1000) < 311)
    }
  }
  
  #indices for 1s
  ninds <- which(adjacency_matrix == 1, arr.ind = TRUE)
  
  # Estimate scaling factor for BYM2 model with INLA
  #Build the adjacency matrix using INLA library functions
  adj.matrix <- Matrix::sparseMatrix(i = ninds[,1], j = ninds[,2], x = 1, symmetric = TRUE)
  
  #The IAR precision matrix (note! This is singular)
  Q <- Matrix::Diagonal(ncel, Matrix::rowSums(adj.matrix)) - adj.matrix
  #Add a small jitter to the diagonal for numerical stability (optional but recommended)
  Q_pert <- Q + Matrix::Diagonal(ncel) * max(diag(Q)) * sqrt(.Machine$double.eps)
  
  # Compute the diagonal elements of the covariance matrix subject to the 
  # constraint that the entries of the ICAR sum to zero.
  # See the inla.qinv function help for further details.
  Q_inv <- INLA::inla.qinv(Q_pert, 
                           constr = list(A = matrix(1, 1, ncel), e = 0))
  
  #Compute the geometric mean of the variances, which are on the diagonal of Q.inv
  scaling_factor <- exp(mean(log(diag(Q_inv))))
  
  t_sc_ind <- which(df_master$species == species_list[k])
  df_master$scaling_factor[t_sc_ind] <- scaling_factor
}



# write to RDS --------------------------------------------------

IAR_dir_path <- paste0(dir, 'Bird_phenology/Data/Processed/IAR_', Sys.Date())

dir.create(IAR_dir_path)
setwd(IAR_dir_path)

saveRDS(df_master, paste0('IAR_input-', Sys.Date(), '.rds'))



# create list of species to run through IAR model -------------------------

species_tm <- aggregate(MODEL ~ species, data = df_master, FUN = function(x) sum(x, na.rm = TRUE))$species

setwd(paste0(dir, 'Bird_Phenology/Data/'))
write.table(species_tm, file = 'IAR_species_list.txt', row.names = FALSE, col.names = FALSE)










# create dfs that show # cells with data in each year/species, and # years with data in each cell/species -----------------

#create vector of year names
yrs_vec <- c()
for (i in min(years):max(years))
{
  yrs_vec <- c(yrs_vec, paste0('yr_', i))
}

#create vector of cell names
cell_vec <- c()
for (i in 1:ncel)
{
  cell_vec <- c(cell_vec, paste0('cell_', cells[i]))
}


#create df with species/cells/n_yrs per sp_cell
cells_frame <- data.frame(species = rep(species_list, each = length(cells)), 
                          cell = rep(cells, length(species_list)), 
                          n_yrs = NA)

yrs_frame <- data.frame(species = rep(species_list, each = length(years)), 
                        year = rep(years, length(species_list)), 
                        n_cells = NA)


#add years and cell cols
cells_frame[yrs_vec] <- NA
yrs_frame[cell_vec] <- NA

#fill cells_frame and yrs_frame
for (i in 1:nsp)
{
  #i <- 80
  
  for (k in 1:ncel)
  {
    #k <- 30
    t_cell <- dplyr::filter(diagnostics_frame, species == species_list[i], cell == cells[k])
    
    #get index for species/cell row
    idx <- which(cells_frame$species == species_list[i] & 
                   cells_frame$cell == cells[k])
    
    #insert number of yrs with data
    yrs_d <- t_cell$year[which(!is.na(t_cell$HM_mean))]
    cells_frame[idx,'n_yrs'] <- length(yrs_d)
    
    #find out which years have data and insert TRUE into df - leave NA otherwise
    temp_yrs <- paste0('yr_', yrs_d)
    cells_frame[idx, which(colnames(cells_frame) %in% temp_yrs)] <- TRUE
  }
  
  for (j in min(years):max(years))
  {
    #j <- 2015
    t_yr <- dplyr::filter(diagnostics_frame, species == species_list[i], year == j)
    
    #get index for species/year row
    idx_yr <- which(yrs_frame$species == species_list[i] & 
                      yrs_frame$year == j)
    
    #insert number of cells with data
    cells_d <- t_yr$cell[which(!is.na(t_yr$HM_mean))]
    yrs_frame[idx_yr,'n_cells'] <- length(cells_d)
    
    #find out which years have data and insert TRUE into df - leave NA otherwise
    temp_cells <- paste0('cell_', cells_d)
    yrs_frame[idx_yr, which(colnames(yrs_frame) %in% temp_cells)] <- TRUE
  }
}

#write to RDS
setwd(IAR_dir_path)
saveRDS(cells_frame, paste0('cells_frame-', Sys.Date(), '.rds'))
saveRDS(yrs_frame, paste0('yrs_frame-', Sys.Date(), '.rds'))



# explore data ------------------------------------------------------------

aggregate(n_cells ~ species, data = yrs_frame, FUN = max)
aggregate(n_cells ~ species, data = yrs_frame, FUN = mean)


