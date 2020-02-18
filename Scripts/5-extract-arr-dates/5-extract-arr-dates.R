######################
# 5 - Extract arrival dates from IAR output
#
# Compiles dataframe with pre (2-halfmax-arr.R) and post (4-IAR-arr.R) IAR data for every year/cell
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/labs/Tingley/phenomismatch/'



# other dir ---------------------------------------------------------------

IAR_in_dir <- 'IAR_input_2019-05-03'
IAR_out_dir <- '~/Google_Drive/R/Bird_Phenology/Data/Processed/IAR_output_2020-02-10'
master_out_dir <- '~/Google_Drive/R/Bird_Phenology/Data/Processed/arrival_master_2020-02-10'


# Load packages -----------------------------------------------------------

library(MCMCvis)
library(dplyr)
library(dggridR)
library(reshape2)


# setwd -------------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_in_dir))

IAR_in_date <- substr(IAR_in_dir, start = 11, stop = 20)
IAR_out_date <- substr(IAR_out_dir, start = 12, stop = 21)

nc_out_dir <- nchar(IAR_out_dir)
IAR_out_date <- substr(IAR_out_dir, start = (nc_out_dir - 9), stop = nc_out_dir)


# Filter data by species/years ------------------------------------------------------

#read in master df
df_master <- readRDS(paste0('IAR_input-', IAR_in_date, '.rds'))

species <- as.character(read.table(paste0(dir, 'Bird_Phenology/Data/IAR_species_list.txt'))[,1])




# create empty dataframes to fill -----------------------------------------

#combine pre and post IAR data for every year/cell that was modeled (including cell lat/lon)
out <- data.frame(species = rep(NA, NROW(df_master)), cell = NA,
                                mig_cell = NA, breed_cell = NA,
                                cell_lat = NA, cell_lon = NA,
                                year = NA, mean_pre_IAR = NA, sd_pre_IAR = NA,
                                mean_post_IAR = NA, sd_post_IAR = NA,
                                mean_gamma = NA, sd_gamma = NA,
                                mean_beta0 = NA, sd_beta0 = NA,
                                mean_alpha_gamma = NA, sd_alpha_gamma = NA,
                                mean_beta_gamma = NA, sd_beta_gamma = NA,
                                num_diverge = NA, max_rhat = NA, min_neff = NA)



# #create NA matrix for posterior iter
# iter_mat <- matrix(NA, nrow = NROW(df_master), ncol = 24000)
# colnames(iter_mat) <- paste0('iter_', 1:24000)
# 
# #create empty dataframe for posterior
# arr_post <- cbind(data.frame(species = rep(NA, NROW(df_master)), cell = NA, 
#                              mig_cell = NA, breed_cell = NA, 
#                              year = NA, cell_lat = NA, cell_lon = NA), iter_mat)
# 
# rm(iter_mat)
# gc()


# run loop to fill empty dfs ----------------------------------------------------------------

counter <- 1
for (i in 1:length(species))
{
  #i <- 94 #Vireo olivaceus
  #i <- 1
  
  #filter by species
  sp <- species[i]
  
  #switch to out dir
  setwd(IAR_out_dir)
  
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
  
    #make hexgrid
    hexgrid6 <- dggridR::dgconstruct(res = 6)
  
    #get hexgrid cell centers
    cellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid6, t_cells)
  
  
    # filter cells ------------------------------------------------------------
    
    #from 2-halfmax-arr.R
    '%ni%' <- Negate('%in%')
    
    #reference key for species synonyms
    setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/metadata/'))
    sp_key <- read.csv('species_filenames_key.csv')
    
    #change dir to shp files
    setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/shapefiles/'))
    
    hge <- rgdal::readOGR('global_hex.shp', verbose = FALSE)
    
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
      poly_int_br <- rgeos::gIntersects(hge, br_rng_sp, byid = TRUE)
      tpoly_br <- which(poly_int_br == TRUE, arr.ind = TRUE)[,2]
      br_cells <- as.numeric(tpoly_br[!duplicated(tpoly_br)])
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
      poly_int_mig <- rgeos::gIntersects(hge, mig_rng_sp, byid = TRUE)
      tpoly_mig <- which(poly_int_mig == TRUE, arr.ind = TRUE)[,2]
      mig_cells <- as.numeric(tpoly_mig[!duplicated(tpoly_mig)])
      mig_cells_f <- mig_cells[which(mig_cells %in% t_cells)]
      f_in$mig_cell <- f_in$cell %in% mig_cells_f
    } else {
      f_in$mig_cell <- rep(FALSE, NROW(f_in))
    }
    
    
    # extract posteriors ------------------------------------------------------
    
    #extract median and sd for IAR arrival dates
    mean_fit <- MCMCvis::MCMCpstr(t_fit, params = 'y_true', func = mean)[[1]]
    med_fit <- MCMCvis::MCMCpstr(t_fit, params = 'y_true', func = median)[[1]]
    sd_fit <- MCMCvis::MCMCpstr(t_fit, params = 'y_true', func = sd)[[1]]
    
    #extract cell random effect (gamma)
    mean_gamma <- MCMCvis::MCMCpstr(t_fit, params = 'gamma', func = mean)[[1]]
    sd_gamma <- MCMCvis::MCMCpstr(t_fit, params = 'gamma', func = sd)[[1]]
    
    #extract year effect (beta0)
    mean_beta0 <- MCMCvis::MCMCpstr(t_fit, params = 'beta0', func = mean)[[1]]
    sd_beta0 <- MCMCvis::MCMCpstr(t_fit, params = 'beta0', func = sd)[[1]]
    
    #extract arrival date of species at lat 0 (alpha_gamma)
    mean_alpha_gamma <- MCMCvis::MCMCpstr(t_fit, params = 'alpha_gamma', 
                                          func = mean)[[1]]
    sd_alpha_gamma <- MCMCvis::MCMCpstr(t_fit, params = 'alpha_gamma', 
                                        func = sd)[[1]]
    alpha_gamma_ch <- MCMCvis::MCMCchains(t_fit, params = 'alpha_gamma')
    colnames(alpha_gamma_ch) <- sp
    
    #extract migration speed (beta_gamma)
    mean_beta_gamma <- MCMCvis::MCMCpstr(t_fit, params = 'beta_gamma', 
                                          func = mean)[[1]]
    sd_beta_gamma <- MCMCvis::MCMCpstr(t_fit, params = 'beta_gamma', 
                                        func = sd)[[1]]
    beta_gamma_ch <- MCMCvis::MCMCchains(t_fit, params = 'beta_gamma')
    colnames(beta_gamma_ch) <- sp
    
    #diagnostics
    num_diverge <- rstan::get_num_divergent(t_fit)
    model_summary <- MCMCvis::MCMCsummary(t_fit, excl = 'y_rep', round = 2)
    max_rhat <- max(as.vector(model_summary[, grep('Rhat', colnames(model_summary))]))
    min_neff <- min(as.vector(model_summary[, grep('n.eff', colnames(model_summary))]))
    
    # #extract posteriors for arrival dates
    # yt_ch <- MCMCvis::MCMCpstr(t_fit, params = 'y_true', type = 'chains')[[1]]
    # #colnames for posterior df
    # iter_lab <- paste0('iter_', 1:dim(yt_ch)[3])
    
    # #matrix objects with alpha_gamma and beta_gamma posteriors
    # if (i == 1)
    # {
    #   alpha_gamma_post <- alpha_gamma_ch
    #   beta_gamma_post <- beta_gamma_ch
    # } else {
    #   alpha_gamma_post <- cbind(alpha_gamma_post, alpha_gamma_ch)
    #   beta_gamma_post <- cbind(beta_gamma_post, beta_gamma_ch)
    # }
    # 
    #loop through years
    for (j in 1:length(t_years))
    {
      #j <- 1
      print(paste0('species: ', sp, ', ', 
                   'year: ', t_years[j]))
      
      t_f_in <- dplyr::filter(f_in, year == t_years[j])
      
      t_full <- data.frame(t_f_in[,c('species','cell', 'mig_cell', 'breed_cell')], 
                         cell_lat = round(cellcenters$lat_deg, digits = 2), 
                         cell_lon = round(cellcenters$lon_deg, digits = 2),
                         t_f_in[,c('year', 'HM_mean', 'HM_sd')],
                         mean_post_IAR = mean_fit[j,], 
                         sd_post_IAR = sd_fit[j,],
                         mean_gamma,
                         sd_gamma,
                         mean_beta0[j],
                         sd_beta0[j],
                         mean_alpha_gamma,
                         sd_alpha_gamma,
                         mean_beta_gamma,
                         sd_beta_gamma,
                         num_diverge,
                         max_rhat,
                         min_neff)
      
      colnames(t_full)[c(8,9,14,15)] <- c('mean_pre_IAR', 'sd_pre_IAR', 'mean_beta0', 'sd_beta0')
     
      #fill empty df
      out[counter:(counter + NROW(t_full) - 1),] <- t_full
      
      #run time is excessive due to mem constraints on desktop
      # t_post <- data.frame(t_f_in[,c('species','cell', 'mig_cell', 'breed_cell', 'year')], 
      #                      cell_lat = round(cellcenters$lat_deg, digits = 2), 
      #                      cell_lon = round(cellcenters$lon_deg, digits = 2),
      #                      yt_ch[,j,])
      # colnames(t_post)[-c(1:7)] <- iter_lab
      # 
      # #fill empty df
      # arr_post[counter:(counter + NROW(t_full) - 1),] <- t_post
      
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
ifelse(!dir.exists(master_out_dir),
       dir.create(master_out_dir),
       FALSE)

setwd(master_out_dir)

saveRDS(out2, file = paste0('arrival_master_', IAR_out_date, '.rds'))
# saveRDS(alpha_gamma_post, file = paste0('arrival_alpha_gamma_post_', IAR_out_date, '.rds'))
# saveRDS(beta_gamma_post, file = paste0('arrival_beta_gamma_post_', IAR_out_date, '.rds'))
# saveRDS(arr_post, file = paste0('arrival_post_', IAR_out_date, '.rds'))

print('I completed!')
