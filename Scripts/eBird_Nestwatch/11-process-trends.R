######################
# Extract trends for arr and juv
#
# sigma - interannual variability in arrival date after accounting for trend
# gamma - intercept for phenological event
# theta - effect of lat on phenological event
# pi - magnitude of phenological change over time
# nu - effect of lat on magnitude of phenological change over time
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/labs/Tingley/phenomismatch/'



# other dir ---------------------------------------------------------------

arr_run_date <- '2019-10-15'
juv_run_date <- '2019-10-15'
IAR_date <- '2019-05-26'


# arr_dir <- paste0(dir, 'Bird_Phenology/Data/Processed/arr_trends_output',
#                   arr_run_date)
# juv_dir <- paste0(dir, 'Bird_Phenology/Data/Processed/juv_trends_output_',
#                   juv_run_date)

arr_dir <- paste0('/Users/caseyyoungflesh/Desktop/Bird_Phenology_Offline/Data/Processed/arr_trends_output_',
                  arr_run_date)
juv_dir <- paste0('/Users/caseyyoungflesh/Desktop/Bird_Phenology_Offline/Data/Processed/juv_trends_output_',
                  juv_run_date)

# out_dir <- paste0(dir, 'Bird_Phenology/Data/Processed/trends_summary_', arr_run_date)
out_dir <- paste0('/Users/caseyyoungflesh/Desktop/Bird_Phenology_Offline/Data/Processed/trends_summary_', arr_run_date)

am_dir <- paste0('~/Google_Drive/R/Bird_Phenology/Data/Processed/arrival_master_', IAR_date)



# read in arrival master ---------------------------------------------------

setwd(am_dir)

am_out <- readRDS(paste0('arrival_master_', IAR_date, '.rds'))
am_out_f <- unique(am_out[,c('species', 'mean_alpha_gamma', 'sd_alpha_gamma', 
                                'mean_beta_gamma' ,'sd_beta_gamma')])


# Load packages -----------------------------------------------------------

library(MCMCvis)
library(dplyr)
library(dggridR)


# setwd -------------------------------------------------------------------

setwd(arr_dir)


# Filter arr results by species/years ------------------------------------------------------

files <- list.files()
species <- unique(sapply(strsplit(files, split = '-'), head , 1))

out <- data.frame()
for (i in 1:length(species))
{
  #i <- 3
  
  #filter by species
  sp <- species[i]
  print(sp)
  
  #if that species RDS object exists in dir
  if (length(grep(paste0(sp, '-', arr_run_date, '-arr_trends_stan_output.rds'), list.files())) > 0)
  {
    #read in model output
    t_fit <- readRDS(paste0(sp, '-', arr_run_date, '-arr_trends_stan_output.rds'))
    t_in <- readRDS(paste0(sp, '-', arr-run_date, '-arr_trends_stan_input.rds'))
    
    #make sure there are samples
    if (length(t_fit@par_dims) > 0)
    {
      #add year range and cell range (maybe range area calculated from range maps)?
      
      mn_gamma <- MCMCvis::MCMCpstr(t_fit, params = 'gamma', func = mean)[[1]]
      sd_gamma <- MCMCvis::MCMCpstr(t_fit, params = 'gamma', func = sd)[[1]]
      # gamma_ch <- MCMCvis::MCMCchains(t_fit, params = 'gamma')
      # colnames(gamma_ch) <- sp
      
      mn_theta <- MCMCvis::MCMCpstr(t_fit, params = 'theta', func = mean)[[1]]
      sd_theta <- MCMCvis::MCMCpstr(t_fit, params = 'theta', func = sd)[[1]]
      # theta_ch <- MCMCvis::MCMCchains(t_fit, params = 'theta')
      # colnames(theta_ch) <- sp
      
      mn_pi <- MCMCvis::MCMCpstr(t_fit, params = 'pi', func = mean)[[1]]
      sd_pi <- MCMCvis::MCMCpstr(t_fit, params = 'pi', func = sd)[[1]]
      # pi_ch <- MCMCvis::MCMCchains(t_fit, params = 'pi')
      # colnames(pi_ch) <- sp
      
      mn_nu <- MCMCvis::MCMCpstr(t_fit, params = 'nu', func = mean)[[1]]
      sd_nu <- MCMCvis::MCMCpstr(t_fit, params = 'nu', func = sd)[[1]]
      # nu_ch <- MCMCvis::MCMCchains(t_fit, params = 'nu')
      # colnames(nu_ch) <- sp
      
      mn_sigma <- MCMCvis::MCMCpstr(t_fit, params = 'sigma', func = mean)[[1]]
      sd_sigma <- MCMCvis::MCMCpstr(t_fit, params = 'sigma', func = sd)[[1]]
      # sigma_ch <- MCMCvis::MCMCchains(t_fit, params = 'sigma')
      # colnames(sigma_ch) <- sp
      
      mn_beta <- MCMCvis::MCMCpstr(t_fit, params = 'beta', func = mean)[[1]]
      sd_beta <- MCMCvis::MCMCpstr(t_fit, params = 'beta', func = sd)[[1]]
      
      n_cells <- t_in$NC
      n_years <- length(unique(t_in$year))
      rng_lat <- range(t_in$lat)[2] - range(t_in$lat)[1]
      #mn_lat <- mean(t_in$lat_usc)
      
      #which index is the most northerly
      max_idx <- which.max(t_in$lat)
      min_idx <- which.min(t_in$lat)
      
      mxl_mn_beta <- mn_beta[max_idx]
      mxl_sd_beta <- sd_beta[max_idx]
      
      mnl_mn_beta <- mn_beta[min_idx]
      mnl_sd_beta <- sd_beta[min_idx]
      
      #diagnostics
      num_diverge <- rstan::get_num_divergent(t_fit)
      model_summary <- MCMCvis::MCMCsummary(t_fit, excl = 'y_rep', round = 2)
      max_rhat <- max(as.vector(model_summary[, grep('Rhat', colnames(model_summary))]))
      min_neff <- min(as.vector(model_summary[, grep('n.eff', colnames(model_summary))]))
      
      #arrival master
      am_t <- dplyr::filter(am_out_f, species == sp)
      
      t_mean_pre_IAR <- dplyr::filter(am_out, species == sp, year == 2018)$mean_pre_IAR
      
      #percent of cells over range that data is avail for in 2018
      nc_pre_IAR <- sum(!is.na(t_mean_pre_IAR))
      per_pre_IAR <- nc_pre_IAR / length(t_mean_pre_IAR)
      
      #percent of cells over range that have 3+ years data
      per_cd <- n_cells / length(t_mean_pre_IAR)
      
      
      t_full <- data.frame(species = sp,
                           mn_alpha_gamma = am_t$mean_alpha_gamma,
                           sd_alpha_gamma = am_t$sd_alpha_gamma,
                           mn_beta_gamma = am_t$mean_beta_gamma,
                           sd_beta_gamma = am_t$sd_beta_gamma,
                           mn_gamma,
                           sd_gamma,
                           mn_theta,
                           sd_theta,
                           mn_pi,
                           sd_pi,
                           mn_nu,
                           sd_nu,
                           mn_sigma,
                           sd_sigma,
                           n_cells,
                           n_years,
                           per_pre_IAR,
                           per_cd,
                           rng_lat,
                           mxl_mn_beta,
                           mxl_sd_beta,
                           mnl_mn_beta,
                           mnl_sd_beta,
                           #mn_lat,
                           num_diverge,
                           max_rhat,
                           min_neff)
      
      out <- rbind(out, t_full)
      
      if (i == 1)
      {
        # mu_alpha_post <- mu_alpha_ch
        # mu_beta_post <- mu_beta_ch
        # sigma_post <- sigma_ch
        # alpha_beta_post <- alpha_beta_ch
        # beta_beta_post <- beta_beta_ch
      } else {
        # mu_alpha_post <- cbind(mu_alpha_post, mu_alpha_ch)
        # mu_beta_post <- cbind(mu_beta_post, mu_beta_ch)
        # sigma_post <- cbind(sigma_post, sigma_ch)
        # alpha_beta_post <- cbind(alpha_beta_post, alpha_beta_ch)
        # beta_beta_post <- cbind(beta_beta_post, beta_beta_ch)
      }
    } else {
      print('No samples for species')
    }
  } else {
    print(paste0('.rds file for ', sp, ' not found in directory'))
  }
} #end species loop


#create dir if doesn't exist
ifelse(!dir.exists(out_dir),
       dir.create(out_dir),
       FALSE)

setwd(out_dir)

saveRDS(out, file = paste0('arr_trends_master_', arr_run_date, '.rds'))
# saveRDS(mu_alpha_post, file = paste0('pheno_trends_mu_alpha_post_', arr_run_date, '.rds'))
# saveRDS(mu_beta_post, file = paste0('pheno_trends_mu_beta_post_', arr_run_date, '.rds'))
# saveRDS(sigma_post, file = paste0('pheno_trends_sigma_post_', arr_run_date, '.rds'))
# saveRDS(alpha_beta_post, file = paste0('pheno_trends_alpha_beta_post_', run_date, '.rds'))
# saveRDS(beta_beta_post, file = paste0('pheno_trends_beta_beta_post_', run_date, '.rds'))

print('I completed!')
