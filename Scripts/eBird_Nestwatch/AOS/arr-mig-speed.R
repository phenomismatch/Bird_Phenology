########################
# arrival ~ migration distance + migration speed + mean breeding lat
#
#
########################



# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/UCHC/LABS/Tingley/phenomismatch/'



# other dir ---------------------------------------------------------------

IAR_date <- '2019-05-26'


# Load packages -----------------------------------------------------------

library(rstan)
library(MCMCvis)
library(dplyr)


# setwd -------------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/arrival_master_', IAR_date))



# process data ------------------------------------------------------------

#master data
master_in <- readRDS(paste0('arrival_master_', IAR_date, '.rds'))

#species list
sp <- as.character(read.table(paste0(dir, 'Bird_Phenology/Data/IAR_species_list.txt'))[,1])

#filter master data for breeding range
master_in_f <- dplyr::filter(master_in, 
                             breed_cell == TRUE, mig_cell == FALSE)

#calculate: 
#mean arrival date (over breeding range)
#standard error of mean arrival date (over breeding range)
#mean breeding latitude
#standard error of mean breeding latitude
#number of cells of input data for IAR in 2018
#percent of cells in migratory and breeding range that had input data for IAR in 2018
#mean beta_gamma from IAR model
#sd beta_gamma from IAR model

arr_df <- data.frame(species = rep(NA, length(sp)),
                     mn_arr = NA,
                     std_err_arr = NA,
                     mn_lat = NA,
                     std_err_lat = NA,
                     rng_lat = NA,
                     nc_pre_IAR = NA,
                     per_pre_IAR = NA,
                     mn_beta_gamma = NA,
                     sd_beta_gamma = NA)

c1 <- 1
for (i in 1:length(sp))
{
  #i <- 1
  t_sp <- dplyr::filter(master_in_f, species == sp[i])
  t_cell <- unique(t_sp$cell)
  
  t_mean_pre_IAR <- dplyr::filter(master_in, species == sp[i], year == 2018)$mean_pre_IAR
  nc_pre_IAR <- sum(!is.na(t_mean_pre_IAR))
  per_pre_IAR <- nc_pre_IAR / length(t_mean_pre_IAR)
  
  #to combine posteriors:
  # out_p <- matrix(NA, nrow = NROW(t_sp), ncol = 2000)
  # counter <- 1
  # for (j in 1:length(t_cell))
  # {
  #   #j <- 1
  #   t_sp2 <- dplyr::filter(t_sp, cell == t_cell[j])
  #   t_year <- unique(t_sp2$year)
  #   
  #   for (k in 1:length(t_year))
  #   {
  #     #k <- 1
  #     t_sp3 <- dplyr::filter(t_sp2, year == t_year[k])
  #     #simulate posterior for each cell
  #     out_p[counter,] <- rnorm(2000, t_sp3$mean_post_IAR, t_sp3$sd_post_IAR)
  #     counter <- counter + 1
  #   }
  # }
  
  arr_df$species[c1] <- sp[i]
  arr_df$mn_arr[c1] <- mean(t_sp$mean_post_IAR)
  arr_df$std_err_arr[c1] <- sd(t_sp$mean_post_IAR) / sqrt(NROW(t_sp))
  
  #calculate mean and std error of mean for lat
  fl <- dplyr::filter(t_sp, year == t_sp$year[1])
  arr_df$mn_lat[c1] <- mean(fl$cell_lat)
  arr_df$std_err_lat[c1] <- sd(fl$cell_lat) / sqrt(NROW(fl))
  arr_df$rng_lat[c1] <- max(fl$cell_lat) - min(fl$cell_lat)
  
  #insert beta_gamma
  arr_df$mn_beta_gamma[c1] <- t_sp$mean_beta_gamma[1]
  arr_df$sd_beta_gamma[c1] <- t_sp$sd_beta_gamma[1]
  
  #insert info number of cells pre IAR
  arr_df$nc_pre_IAR[c1] <- nc_pre_IAR
  arr_df$per_pre_IAR[c1] <- round(per_pre_IAR, 2)
  
  c1 <- c1 + 1
}



# merge mig dis with arr -----------------------------------------------------------------

#read in traits
setwd("~/Google_Drive/R/Bird_Phenology/Data/Traits")
tdb <- read.csv('Trait_database-2019-07-07.csv')
tdb$species <- gsub(' ', '_', tdb$SCI_NAME)

arr_mrg <- dplyr::left_join(arr_df, 
                            tdb[,c('species', 'MIGRATION_DISTANCE_LASORTE')], by = 'species')

colnames(arr_mrg)[grep('MIGRATION', colnames(arr_mrg))] <- 'mig_dis'




# Stan model arr ~ speed + mig_dis + cell_lat ---------------------------------------------------------

#model missing covariate data (migratory distance) using approach in MCcElreath 2016 p. 432
#https://github.com/ssp3nc3r/rethinking/blob/master/chapter14.Rmd


#only species that have data from at least 20% of their breeding and migratory range
arr_mrg_f <- dplyr::filter(arr_mrg, per_pre_IAR > 0.2)

#|r| < 0.5 - multicolinearity not a major issue here
cor(arr_mrg_f$mig_dis, arr_mrg_f$mn_beta_gamma, use = 'pairwise.complete.obs')


#create data list for Stan
DATA <- list(N = NROW(arr_mrg_f),
             y_obs = arr_mrg_f$mn_arr,
             y_sd = arr_mrg_f$std_err_arr,
             bg_obs = arr_mrg_f$mn_beta_gamma,
             bg_sd = arr_mrg_f$sd_beta_gamma,
             lat_obs = arr_mrg_f$mn_lat,
             lat_sd = arr_mrg_f$std_err_lat,
             #fill missing slots with 0
             mig_dis_obs = ifelse(is.na(arr_mrg_f$mig_dis), 0,  arr_mrg_f$mig_dis / 1000),
             N_miss = sum(is.na(arr_mrg_f$mig_dis)),
             miss_idx = is.na(arr_mrg_f$mig_dis) * cumsum(is.na(arr_mrg_f$mig_dis)))



# Stan model --------------------------------------------------------------

#y_{obs} \sim N(y_{true}, \sigma_{y})
#y_{true} \sim N(\mu, \sigma)
#\mu = \alpha + \beta \times MD + \gamma \times IS_{true} + \pi \times LAT_{true}
#MD \sim N(\mu_{MD}, \sigma_{MD})
#IS_{obs} \sim N(IS_{true}, \sigma_{IS})
#LAT_{obs} \sim N(LAT_{true}, \sigma_{LAT})

#y_{obs} - mean arrival date over species breeding range (mean of all mean cell arrival dates over breeding range) - given as data
#y_{sd} - standard error arrival date over species breeding range (sd of all mean cell arrival dates over breeding range over the sqrt number of cell/years) - given as data
#MD - migration distance, from La Sorte et al. 2013 (Ecology). Not available for all species. When data are available, line 4 is a likelihood (estimates \mu_{MD} and \sigma_{MD}), when data are not available, it is a prior - given as data
#IS_{obs} - inverse migration speed (mean posterior for \beta_{\gamma} from IAR model) - given as data
#\sigma_{IS} - uncertainty inverse migration speed (sd posterior for \beta_{\gamma} from IAR model) - given as data
#LAT_{obs} - mean breeding latitude (mean of cell center latitudes that cover species breeding range) - given as data
#\sigma_{LAT} - standard error of LAT_{obs} (sd(lat) / sqrt(N)) - given as data


model <- "
  data {
  int<lower = 0> N;                                 // number of obs
  vector<lower = 0, upper = 200>[N] y_obs;
  vector<lower = 0>[N] y_sd;
  vector<lower = 0>[N] mig_dis_obs;                // for missing cov
  int N_miss;                                       // for missing cov
  int<lower = 0, upper = N_miss> miss_idx[N];       // for missing cov
  vector<lower = 0>[N] bg_obs;
  vector<lower = 0>[N] bg_sd;
  vector<lower = 0>[N] lat_obs;
  vector<lower = 0>[N] lat_sd;
  }
  
  parameters {
  real<lower = 0> sigma_raw;
  vector[N] bg_true_raw;
  vector[N] lat_true_raw;
  vector[N] y_true_raw;
  real alpha_raw;
  real beta_raw;
  real gamma_raw;
  real pi_raw;
  real<lower = 0> mu_mig_dis;       // for missing cov
  real<lower = 0> sigma_mig_dis;    // for missing cov
  vector<lower = 0>[N_miss] mig_dis_impute;    // for missing cov
  }
  
  transformed parameters {
  vector[N] mu;
  vector[N] bg_true;
  vector[N] lat_true;
  real<lower = 0> sigma;
  vector[N] y_true;
  real alpha;
  real beta;
  real gamma;
  real pi;
  vector<lower = 0>[N] mig_dis;     // for missing cov
  
  alpha = alpha_raw * 70 + 100;
  beta = beta_raw * 100;
  gamma = gamma_raw * 30;
  pi = pi_raw * 10;
  sigma = sigma_raw * 20;
  bg_true = bg_true_raw * 10 + 5;
  lat_true = lat_true_raw * 10 + 45;
  
  // for missing cov
  // combine observe and estimates for mising
  mig_dis = mig_dis_obs;
  for (i in 1:N)
  {
    if (miss_idx[i] > 0)
    {
      mig_dis[i] = mig_dis_impute[miss_idx[i]];
    }
  }
  
  mu = alpha + beta * mig_dis + gamma * bg_true + pi * lat_true;
  
  // implies y_true ~ N(mu, sigma)
  y_true = y_true_raw * sigma + mu;
  }
  
  model {
  alpha_raw ~ std_normal();
  beta_raw ~ std_normal();
  gamma_raw ~ std_normal();
  sigma_raw ~ std_normal();
  lat_true_raw ~ std_normal();
  bg_true_raw ~ std_normal();
  y_true_raw ~ std_normal();
  mu_mig_dis ~ normal(4, 2);
  sigma_mig_dis ~ normal(2, 2);
  
  mig_dis ~ normal(mu_mig_dis, sigma_mig_dis);   // for missing cov
  bg_obs ~ normal(bg_true, bg_sd);
  lat_obs ~ normal(lat_true, lat_sd);
  
  y_obs ~ normal(y_true, y_sd);
  }
  
  generated quantities {
  
  real y_rep[N];
  
  // PPC
  y_rep = normal_rng(y_true, y_sd);
  }
  "



# Run model ---------------------------------------------------------------

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.85
TREE_DEPTH <- 15
STEP_SIZE <- 0.0025
CHAINS <- 4
ITER <- 5000

tt <- proc.time()
fit <- rstan::stan(model_code = model,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('mu', 'sigma', 'alpha', 
                            'beta', 'gamma', 'pi', 'mu_mig_dis',
                            'sigma_mig_dis', 'mig_dis_impute',
                            'y_true', 'bg_true', 'lat_true', 
                            'y_rep'),
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60



# Model summaries ---------------------------------------------------------

MCMCvis::MCMCsummary(fit, round = 3, params = c('beta', 'gamma', 'pi'))

MCMCvis::MCMCplot(fit, params = c('beta', 'gamma', 'pi'))#, 
                  # labels = c('Migratory distance (beta)',
                  #            'Inverse migration speed (gamma)',
                  #            'Mean breeding latitude (pi)'))

#MCMCvis::MCMCplot(fit, params = 'mig_dis_impute')


# Shinystan ---------------------------------------------------------------

# library(shinystan)
# launch_shinystan(fit)



# PPC ---------------------------------------------------------------------

y_val <- DATA$y_obs
y_rep <- MCMCvis::MCMCchains(fit, params = 'y_rep')
bayesplot::ppc_dens_overlay(y_val, y_rep[1:100,])

N <- 100
plot(matrix(rep(y_val, N), nrow = N, byrow = TRUE), 
     y_rep[1:N,], pch = '.')
abline(a = 0, b = 1, col = 'red', lty = 2)



# PPO ---------------------------------------------------------------------

# alpha = alpha_raw * 70 + 100;
# beta = beta_raw * 100;
# gamma = gamma_raw * 30;
# pi = pi_raw * 10;
# sigma = sigma_raw * 20;
# bg_true = bg_true_raw * 10 + 5;
# lat_true = lat_true_raw * 10 + 45;
# mu_mig_dis ~ normal(4, 2);
# sigma_mig_dis ~ normal(2, 2);

setwd('~/Desktop')
#sigma ~ halfnormal(0, 20)
PR_p <- rnorm(10000, 0, 20)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma',
                   priors = PR,
                   pdf = TRUE,
                   filename = 'trace_sigma.pdf')

#alpha ~ normal(100, 70)
PR <- rnorm(10000, 100, 70)
MCMCvis::MCMCtrace(fit,
                   params = 'alpha',
                   priors = PR,
                   pdf = TRUE,
                   filename = 'trace_alpha.pdf')

#beta ~ normal(0, 100)
PR <- rnorm(10000, 0, 100)
MCMCvis::MCMCtrace(fit,
                   params = 'beta',
                   priors = PR,
                   pdf = TRUE,
                   filename = 'trace_beta.pdf')

#gamma ~ normal(0, 30)
PR <- rnorm(10000, 0, 30)
MCMCvis::MCMCtrace(fit,
                   params = 'gamma',
                   priors = PR,
                   pdf = TRUE,
                   filename = 'trace_gamma.pdf')

#pi ~ normal(0, 10)
PR <- rnorm(10000, 0, 10)
MCMCvis::MCMCtrace(fit,
                   params = 'pi',
                   priors = PR,
                   pdf = TRUE,
                   filename = 'trace_pi.pdf')

#bg_true ~ normal(5, 10)
PR <- rnorm(10000, 5, 10)
MCMCvis::MCMCtrace(fit,
                   params = 'bg_true',
                   priors = PR,
                   pdf = TRUE,
                   filename = 'trace_bg_true.pdf')

#lat_true ~ normal(45, 10)
PR <- rnorm(10000, 45, 10)
MCMCvis::MCMCtrace(fit,
                   params = 'lat_true',
                   priors = PR,
                   pdf = TRUE,
                   filename = 'trace_lat_true.pdf')

#mu_mig_dis ~ normal(4, 2)
PR <- rnorm(10000, 4, 2)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_mig_dis',
                   priors = PR,
                   pdf = TRUE,
                   filename = 'trace_mu_mig_dis.pdf')

#sigma_mig_dis ~ normal(2, 2)
PR_p <- rnorm(10000, 2, 2)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_mig_dis',
                   priors = PR,
                   pdf = TRUE,
                   filename = 'trace_sigma_mig_dis.pdf')
