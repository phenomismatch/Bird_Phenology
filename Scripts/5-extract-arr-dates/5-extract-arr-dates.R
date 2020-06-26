######################
# 5 - Extract arrival dates from IAR output
#
# Compiles dataframe with pre (2-halfmax-arr.R) and post (4-IAR-arr.R) IAR data for every year/cell
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'


# other dir ---------------------------------------------------------------

IAR_in_dir <- 'IAR_input_2020-05-07'
IAR_out_dir <- 'IAR_output_2020-06-08'
master_out_dir <- 'arrival_master_2020-06-08'


# Load packages -----------------------------------------------------------

library(MCMCvis)
library(rstan)
library(dplyr)
library(dggridR)
library(rgdal)
library(sp)
library(rgeos)


# setwd -------------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_in_dir))

IAR_in_date <- substr(IAR_in_dir, start = 11, stop = 20)
IAR_out_date <- substr(IAR_out_dir, start = 12, stop = 21)


# Filter data by species/years ------------------------------------------------------

#read in master df
df_master <- readRDS(paste0('IAR_input-', IAR_in_date, '.rds'))

species <- as.character(read.table(paste0(dir, 'Bird_Phenology/Data/IAR_species_list.txt'))[,1])


# create empty dataframes to fill -----------------------------------------

#combine pre and post IAR data for every year/cell that was modeled (including cell lat/lon)

out <- data.frame(species = rep(NA, NROW(df_master)), cell = NA, 
                  mig_cell = NA, breed_cell = NA,
                  cell_lat = NA, cell_lng = NA, per_ovr = NA,
                  year = NA, arr_GAM_mean = NA, arr_GAM_sd = NA, 
                  VALID_GAM = NA, arr_IAR_mean = NA, arr_IAR_sd = NA, plmax = NA,
                  gamma_mean = NA, gamma_sd = NA,
                  beta0_mean = NA, beta0_sd = NA,
                  sigma_beta0_mean = NA, sigma_beta0_sd = NA,
                  beta_gamma_mean = NA, beta_gamma_sd = NA,
                  num_diverge = NA, max_Rhat = NA, min_neff = NA,
                  PPC_mn_dis = NA, PPC_mn_pval = NA)


# run loop to fill empty dfs ----------------------------------------------------------------

#make hexgrid
hexgrid6 <- dggridR::dgconstruct(res = 6)

#read in cropped grid .shp file
setwd(paste0(dir, 'Bird_Phenology/Data/hex_grid_crop/'))
hge <- rgdal::readOGR('hex_grid_crop.shp', verbose = FALSE)
hge_cells <- as.numeric(as.character(hge@data[,1]))

counter <- 1
for (i in 1:length(species))
{
  #i <- 94 #Vireo olivaceus
  #i <- 1
  
  #filter by species
  sp <- species[i]
  
  #switch to out dir
  setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_out_dir))
  
  #if that species RDS object exists in dir
  if (length(grep(paste0(sp, '-iar-stan_output-', IAR_out_date, '.rds'), list.files())) > 0)
  {
    
    # filter and read in data -------------------------------------------------
    f_in_p <- dplyr::filter(df_master, species == sp & MODEL == TRUE)
    
    #read in IAR model output and input
    t_fit <- readRDS(paste0(sp, '-iar-stan_output-', IAR_out_date, '.rds'))
    t_data <- readRDS(paste0(sp, '-iar-stan_input-', IAR_out_date, '.rds'))
    
    #only cells and years that were modeled (to account for any lone cells that were dropped in 4-IAR-arr.R)
    f_in <- dplyr::filter(f_in_p, cell %in% t_data$cells)
    
    t_cells <- unique(f_in$cell)
    t_years <- unique(f_in$year)
  
    #get hexgrid cell centers
    cellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid6, t_cells)
  
    # filter cells ------------------------------------------------------------
    
    #from 2-halfmax-arr.R
    #reference key for species synonyms
    setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/metadata/'))
    sp_key <- read.csv('species_filenames_key.csv')
    
    #change dir to shp files
    setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/shapefiles/'))
    
    #filter by breeding/migration cells
    #match species name to shp file name
    g_ind <- grep(sp, sp_key$file_names_2016)
    
    #check for synonyms if there are no matches
    if (length(g_ind) == 0)
    {
      g_ind2 <- grep(sp, sp_key$BL_Checklist_name)
    } else {
      g_ind2 <- g_ind
    }
    
    #get filename and read in
    fname <- as.character(sp_key[g_ind2,]$filenames[grep('.shp', sp_key[g_ind2, 'filenames'])])
    sp_rng <- rgdal::readOGR(fname[1], verbose = FALSE)
    
    #filter by breeding (2) range
    br_rng <- sp_rng[which(sp_rng$SEASONAL == 2),]
    #filter by migration (4) range
    mig_rng <- sp_rng[which(sp_rng$SEASONAL == 4),]
    
    if (length(br_rng) > 0)
    {
      br_rng_sp <- sp::SpatialPolygons(br_rng@polygons)
      sp::proj4string(br_rng_sp) <- sp::CRS(sp::proj4string(br_rng))
      #find intersections with code from here: https://gis.stackexchange.com/questions/1504/extracting-intersection-areas-in-r
      hge2 <- sp::spTransform(hge, sp::proj4string(br_rng_sp))
      poly_int_br <- rgeos::gIntersects(hge2, br_rng_sp, byid = TRUE)
      tpoly_br <- which(poly_int_br == TRUE, arr.ind = TRUE)[,2]
      br_cells <- hge_cells[as.numeric(tpoly_br[!duplicated(tpoly_br)])]
      #see which of these cells were actually used in modeling (don't overlap non-migratory/breeding ranges)
      br_cells_f <- br_cells[which(br_cells %in% t_cells)]
      #add to master
      f_in$breed_cell <- f_in$cell %in% br_cells_f
    } else {
      f_in$breed_cell <- rep(FALSE, NROW(f_in))
    }
    
    if (length(mig_rng) > 0)
    {
      mig_rng_sp <- sp::SpatialPolygons(mig_rng@polygons)
      sp::proj4string(mig_rng_sp) <- sp::CRS(sp::proj4string(mig_rng))
      hge2 <- sp::spTransform(hge, sp::proj4string(mig_rng_sp))
      poly_int_mig <- rgeos::gIntersects(hge2, mig_rng_sp, byid = TRUE)
      tpoly_mig <- which(poly_int_mig == TRUE, arr.ind = TRUE)[,2]
      mig_cells <- hge_cells[as.numeric(tpoly_mig[!duplicated(tpoly_mig)])]
      mig_cells_f <- mig_cells[which(mig_cells %in% t_cells)]
      f_in$mig_cell <- f_in$cell %in% mig_cells_f
    } else {
      f_in$mig_cell <- rep(FALSE, NROW(f_in))
    }
    
    
    # extract posteriors ------------------------------------------------------
    
    #extract median and sd for IAR arrival dates
    fit_mean <- MCMCvis::MCMCpstr(t_fit, params = 'y_true', func = mean)[[1]]
    fit_med <- MCMCvis::MCMCpstr(t_fit, params = 'y_true', func = median)[[1]]
    fit_sd <- MCMCvis::MCMCpstr(t_fit, params = 'y_true', func = sd)[[1]]
    
    #extract cell random effect (gamma)
    gamma_mean <- MCMCvis::MCMCpstr(t_fit, params = 'gamma', func = mean)[[1]]
    gamma_sd <- MCMCvis::MCMCpstr(t_fit, params = 'gamma', func = sd)[[1]]
    
    #extract year effect (beta0)
    beta0_mean <- MCMCvis::MCMCpstr(t_fit, params = 'beta0', func = mean)[[1]]
    beta0_sd <- MCMCvis::MCMCpstr(t_fit, params = 'beta0', func = sd)[[1]]
    
    #extract variance year effect (sigma_beta0)
    sigma_beta0_mean <- MCMCvis::MCMCpstr(t_fit, params = 'sigma_beta0', func = mean)[[1]]
    sigma_beta0_sd <- MCMCvis::MCMCpstr(t_fit, params = 'sigma_beta0', func = sd)[[1]]

    #extract migration speed (beta_gamma)
    beta_gamma_mean <- MCMCvis::MCMCpstr(t_fit, params = 'beta_gamma', 
                                          func = mean)[[1]]
    beta_gamma_sd <- MCMCvis::MCMCpstr(t_fit, params = 'beta_gamma', 
                                        func = sd)[[1]]
    # beta_gamma_ch <- MCMCvis::MCMCchains(t_fit, params = 'beta_gamma')
    # colnames(beta_gamma_ch) <- sp
    
    #diagnostics
    num_diverge <- rstan::get_num_divergent(t_fit)
    model_summary <- MCMCvis::MCMCsummary(t_fit, excl = 'y_rep', round = 3)
    max_Rhat <- max(as.vector(model_summary[, grep('Rhat', colnames(model_summary))]))
    min_neff <- min(as.vector(model_summary[, grep('n.eff', colnames(model_summary))]))
    
    #PPC
    y_rep_ch <- MCMCvis::MCMCpstr(t_fit, params = 'y_rep', type = 'chains')[[1]]
    t_y_rep <- t(y_rep_ch)
    
    #remove NA vals
    na.y.rm <- which(is.na(t_data$y_PPC))
    n_y_PPC <- t_data$y_PPC[-na.y.rm]
    n_y_rep <- t_y_rep[, -na.y.rm]

    #PPC
    #p-val
    PPC_fun <- function(FUN, YR = n_y_rep, D = n_y_PPC)
    {
      fout <- sum(apply(YR, 1, FUN) > FUN(D)) / NROW(YR)
      return(fout)
    }
    PPC_mn_pval <- PPC_fun(mean)
    
    #discrepancy
    PPC_dis_fun <- function(FUN, YR = n_y_rep, D = n_y_PPC)
    {
      fout <- apply(YR, 1, FUN) - FUN(D)
      mn_out <- mean(fout)
      return(mn_out)
    }
    PPC_mn_dis <- PPC_dis_fun(mean)

    #loop through years
    for (j in 1:length(t_years))
    {
      #j <- 1
      print(paste0('species: ', sp, ', ', 
                   'year: ', t_years[j]))
      
      t_f_in <- dplyr::filter(f_in, year == t_years[j])
      
      t_full <- data.frame(t_f_in[,c('species','cell', 'mig_cell', 'breed_cell')], 
                         cell_lat = round(cellcenters$lat_deg, digits = 2), 
                         cell_lng = round(cellcenters$lon_deg, digits = 2),
                         t_f_in[,c('per_ovr', 'year', 'arr_GAM_mean', 'arr_GAM_sd')],
                         VALID_GAM = t_f_in$VALID,
                         arr_IAR_mean = fit_mean[j,], 
                         arr_IAR_sd = fit_sd[j,],
                         plmax = t_f_in$plmax,
                         gamma_mean,
                         gamma_sd,
                         beta0_mean = beta0_mean[j],
                         beta0_sd = beta0_sd[j],
                         sigma_beta0_mean,
                         sigma_beta0_sd,
                         beta_gamma_mean,
                         beta_gamma_sd,
                         num_diverge,
                         max_Rhat,
                         min_neff,
                         PPC_mn_dis,
                         PPC_mn_pval)
     
      #fill empty df
      out[counter:(counter + NROW(t_full) - 1),] <- t_full
      
      #advance counter
      counter <- counter + NROW(t_full)
      
    } #end year loop
  } else {
    print(paste0('.rds file for ', sp, ' not found in directory'))
  }
} #end species loop


#remove zeros
out2 <- out[-c(min(which(is.na(out$species))):NROW(out)),]


# write to file -----------------------------------------------------------

#create dir if it does not exist

ifelse(!dir.exists(paste0(dir, 'Bird_Phenology/Data/Processed/', master_out_dir)),
       dir.create(paste0(dir, 'Bird_Phenology/Data/Processed/', master_out_dir)),
       FALSE)

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', master_out_dir))

saveRDS(out2, file = paste0('arrival_master_', IAR_out_date, '.rds'))

print('I completed!')
