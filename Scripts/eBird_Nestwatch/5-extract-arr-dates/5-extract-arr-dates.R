######################
# 5 - Extract arrival dates from IAR output
#
# Compiles dataframe with pre (2-halfmax-arr.R) and post (4-IAR-arr.R) IAR data for every year/cell
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/UCHC/LABS/Tingley/phenomismatch/'



# other dir ---------------------------------------------------------------

IAR_in_dir <- 'IAR_input_2019-05-03'
IAR_out_dir <- '/Users/caseyyoungflesh/Desktop/Bird_Phenology_Offline/Data/Processed/IAR_output_2019-05-26'
master_out_dir <- '~/Google_Drive/R/Bird_Phenology/Data/Processed/arrival_master_2019-05-26'


# Load packages -----------------------------------------------------------

library(MCMCvis)
library(dplyr)
library(dggridR)


# setwd -------------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_in_dir))

IAR_in_date <- substr(IAR_in_dir, start = 11, stop = 20)
IAR_out_date <- substr(IAR_out_dir, start = 12, stop = 21)

nc_out_dir <- nchar(IAR_out_dir)
IAR_out_date <- substr(IAR_out_dir, start = (nc_out_dir - 9), stop = nc_out_dir)


# Filter data by species/years ------------------------------------------------------

#read in master df
df_master <- readRDS(paste0('IAR_input-', IAR_in_date, '.rds'))

species <- as.character(read.table('../../IAR_species_list.txt')[,1])

#switch to out dir
setwd(IAR_out_dir)
#setwd(paste0('~/Desktop/Bird_Phenology_Offline/Data/Processed/', IAR_out_dir))


#combine pre and post IAR data for every year/cell that was modeled (including cell lat/lon)
out <- data.frame()
for (i in 1:length(species))
{
  #i <- 94 #Vireo olivaceus
  #i <- 13
  
  #filter by species
  sp <- species[i]
  
  #if that species RDS object exists in dir
  if (length(grep(paste0(sp, '-', IAR_out_date, '-iar-stan_output.rds'), list.files())) > 0)
  {
    f_in_p <- dplyr::filter(df_master, species == sp & MODEL == TRUE)
    
    #read in IAR model output and input
    t_fit <- readRDS(paste0(sp, '-', IAR_out_date, '-iar-stan_output.rds'))
    t_data <- readRDS(paste0(sp, '-', IAR_out_date, '-iar-stan_input.rds'))
    
    #only cells and years that were modeled (to account for any lone cells that were dropped in 4-IAR-arr.R)
    f_in <- dplyr::filter(f_in_p, cell %in% t_data$cells)
    
    t_cells <- unique(f_in$cell)
    t_years <- unique(f_in$year)
  
    #make hexgrid
    hexgrid6 <- dggridR::dgconstruct(res = 6)
  
    #get hexgrid cell centers
    cellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid6, t_cells)
  
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
    
    #matrix objects with alpha_gamma and beta_gamma posteriors
    if (i == 1)
    {
      alpha_gamma_post <- alpha_gamma_ch
      beta_gamma_post <- beta_gamma_ch
    } else {
      alpha_gamma_post <- cbind(alpha_gamma_post, alpha_gamma_ch)
      beta_gamma_post <- cbind(beta_gamma_post, beta_gamma_ch)
    }
    
    #loop through years
    for (j in 1:length(t_years))
    {
      #j <- 1
      print(paste0('species: ', sp, ', ', 
                   'year: ', t_years[j]))
      
      t_f_in <- dplyr::filter(f_in, year == t_years[j])
      
      t_full <- data.frame(t_f_in[,c('species','cell')], 
                         cell_lat = round(cellcenters$lat_deg, digits = 2), 
                         cell_lon = round(cellcenters$lon_deg, digits = 2),
                         t_f_in[,c('year', 'HM_mean', 'HM_sd')],
                         mean_post_IAR = mean_fit[,j], 
                         sd_post_IAR = sd_fit[,j],
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
      
      colnames(t_full)[c(6,7,12,13)] <- c('mean_pre_IAR', 'sd_pre_IAR', 'mean_beta0', 'sd_beta0')
     
      out <- rbind(out, t_full)
    } #end year loop
  } else {
    print(paste0('.rds file for ', sp, ' not found in directory'))
  }
} #end species loop


#create dir if it does not exist
ifelse(!dir.exists(master_out_dir),
       dir.create(master_out_dir),
       FALSE)

setwd(master_out_dir)

saveRDS(out, file = paste0('arrival_master_', IAR_out_date, '.rds'))
saveRDS(alpha_gamma_post, file = paste0('arrival_alpha_gamma_post_', IAR_out_date, '.rds'))
saveRDS(beta_gamma_post, file = paste0('arrival_beta_gamma_post_', IAR_out_date, '.rds'))

print('I completed!')
