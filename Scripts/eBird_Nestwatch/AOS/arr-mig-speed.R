########################
#calculate average arrival date for each species
#compare to migration speed across migration and breeding range
#
#
########################



# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/UCHC/LABS/Tingley/phenomismatch/'



# other dir ---------------------------------------------------------------

IAR_out_dir <- '/Users/caseyyoungflesh/Desktop/Bird_Phenology_Offline/Data/Processed/IAR_output_2019-05-26'



# Load packages -----------------------------------------------------------

library(rstan)
library(MCMCvis)
library(dplyr)


# setwd -------------------------------------------------------------------

setwd(IAR_out_dir)

nc_out_dir <- nchar(IAR_out_dir)
IAR_out_date <- substr(IAR_out_dir, start = (nc_out_dir - 9), stop = nc_out_dir)


# process data ------------------------------------------------------------

#master data
master_in <- readRDS(paste0('arrival_master_', IAR_out_date, '.rds'))

#run in loop
species <- as.character(read.table(paste0(dir, 'Bird_Phenology/Data/IAR_species_list.txt'))[,1])

out <- data.frame()
for (i in 1:length(species))
{
  #i <- 2
  d_filt <- dplyr::filter(master_in, species == species[i],
                          breed_cell == TRUE, mig_cell == FALSE)
  
  #create data list for Stan
  DATA <- list(N = NROW(d_filt),
               y_obs = d_filt$mean_post_IAR,
               y_sd = d_filt$sd_post_IAR,
               year = as.numeric(factor(d_filt$year)),
               NY = length(unique(d_filt$year)))
  
  
  # Stan model --------------------------------------------------------------
  
  model <- "
  data {
  int<lower = 0> N;                                 // number of obs
  vector<lower = 0, upper = 200>[N] y_obs;
  vector<lower = 0>[N] y_sd;
  int<lower = 1> year[N];                           // year id
  int<lower = 0> NY;                                // number of years
  }
  
  parameters {
  real mu_raw;
  real<lower = 0> sigma_raw;
  vector[NY] mu_yr_raw;
  vector<lower = 0>[NY] sigma_yr_raw;
  vector[N] y_true_raw;
  }
  
  transformed parameters {
  real mu;
  real<lower = 0> sigma;
  vector[NY] mu_yr;
  vector<lower = 0>[NY] sigma_yr;
  vector[N] y_true;
  
  mu = mu_raw * 60 + 160;
  sigma = sigma_raw * 40;
  sigma_yr = sigma_yr_raw * 40;
  
  mu_yr[year] = mu_yr_raw[year] * sigma + mu;     // implies mu_yr[j] ~ N(mu, sigma)
  
  for (i in 1:N)
  {
    y_true[i] = y_true_raw[i] * sigma_yr[year[i]] + mu_yr[year[i]];   // implies y_true[i] ~ N(mu_yr[j], sigma_yr[j])
  }
  
  }
  
  model {
  mu_raw ~ std_normal();
  sigma_raw ~ std_normal();
  mu_yr_raw ~ std_normal();
  sigma_yr ~ std_normal();
  y_true_raw ~ std_normal();
  
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
  
  DELTA <- 0.99
  TREE_DEPTH <- 15
  STEP_SIZE <- 0.01
  CHAINS <- 3
  ITER <- 2000
  
  tt <- proc.time()
  fit <- rstan::stan(model_code = model,
                     data = DATA,
                     chains = CHAINS,
                     iter = ITER,
                     cores = CHAINS,
                     pars = c('mu', 'sigma', 
                              'mu_yr', 'sigma_yr', 
                              'y_true', 'y_rep'),
                     control = list(adapt_delta = DELTA,
                                    max_treedepth = TREE_DEPTH,
                                    stepsize = STEP_SIZE))
  run_time <- (proc.time()[3] - tt[3]) / 60
  
  
  # y_rep <- MCMCvis::MCMCchains(fit, params = 'y_rep')
  # bayesplot::ppc_stat(DATA$y_obs, y_rep, stat = 'mean')
  # bayesplot::ppc_dens_overlay(DATA$y_obs, y_rep[1:500,])
  
  
  num_diverge <- rstan::get_num_divergent(fit)
  num_tree <- rstan::get_num_max_treedepth(fit)
  num_BFMI <- length(rstan::get_low_bfmi_chains(fit))
  
  mean_mu <- round(MCMCvis::MCMCpstr(fit, params = 'mu')[[1]], 4)
  sd_mu <- round(MCMCvis::MCMCpstr(fit, params = 'mu', func = sd)[[1]], 4)
  
  mn_cell_lat <- round(mean(d_filt$cell_lat), 2)
  
  summary <- MCMCvis::MCMCsummary(fit, excl = c('y_true', 'y_rep'))
  max_Rhat <- max(summary[,'Rhat'])
  min_neff <- min(summary[,'n.eff'])
  
  temp <- data.frame(species = species[i],
                     mean_mu,
                     sd_mu,
                     mn_cell_lat,
                     mn_alpha_gamma = round(d_filt$mean_alpha_gamma[1], 4),
                     sd_alpha_gamma = round(d_filt$sd_alpha_gamma[1], 4),
                     mn_beta_gamma = round(d_filt$mean_beta_gamma[1], 4),
                     sd_beta_gamma = round(d_filt$sd_beta_gamma[1], 4),
                     max_Rhat,
                     min_neff,
                     num_diverge,
                     num_tree,
                     num_BFMI)
  
  out <- rbind(out, temp)
}




# plot results ------------------------------------------------------------

#mean and +- 1 sd for data
y_true_mn <- out$mean_mu
y_true_sd <- out$sd_mu
y_true_LCI <- y_true_mn - y_true_sd
y_true_UCI <- y_true_mn + y_true_sd

x_true_mn <- out$mn_beta_gamma
x_true_sd <- out$sd_beta_gamma
x_true_LCI <- x_true_mn - x_true_sd
x_true_UCI <- x_true_mn + x_true_sd


DATA_PLOT <- data.frame(y_true_mn, y_true_LCI, y_true_UCI,
                        x_true_mn, x_true_LCI, x_true_UCI)

plt <- ggplot(data = DATA_PLOT, aes(x_true_mn, y_true_mn), 
              color = 'black', alpha = 0.6) +
  geom_point(data = DATA_PLOT,
             aes(x_true_mn, y_true_mn), color = 'black',
             inherit.aes = FALSE, size = 5, alpha = 0.4) +
  geom_errorbar(data = DATA_PLOT,
                aes(ymin = y_true_LCI, ymax = y_true_UCI), #width = 0.05,
                color = 'black', alpha = 0.2) +
  geom_errorbarh(data = DATA_PLOT,
                 aes(xmin = x_true_LCI, xmax = x_true_UCI), #height = 0.05,
                 color = 'black', alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
    axis.ticks.length= unit(0.2, 'cm')) #length of axis tick



# setwd('~/Desktop')
# ggsave(paste0('arrival_migration_speed.pdf'), plt)








