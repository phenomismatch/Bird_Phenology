####################
# Changes in weight, st. weight, and fat over time - MAPS
# Changes in weight, st. weight, and fat over age - MAPS
#
####################


# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/'


# load packages -----------------------------------------------------------

library(dplyr)
library(rstan)



# load data ---------------------------------------------------------------

setwd(paste0(dir, '/..'))
#setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))

maps_data <- readRDS('MAPS_age_filled.rds')





# data processing --------------------------------------------------------

maps_adults <- dplyr::filter(maps_data, age %in% c('1', '5', '6', '7', '8'))

#QC data - remove zeros
to.rm <- which(maps_adults$weight == 0 | maps_adults$wing_chord == 0 | is.na(maps_adults$weight) | is.na(maps_adults$weight) | is.na(maps_adults$fat_content) | is.na(maps_adults$wing_chord))
maps_adults_qc <- maps_adults[-to.rm, ]


#only species with > 200 data points
d_cnt <- plyr::count(maps_adults_qc, 'sci_name')
sp <- dplyr::filter(d_cnt, freq > 200)[,1]

maps_f <- dplyr::filter(maps_adults_qc, sci_name %in% sp)

#factor for species_id
sp_f <- data.frame(sp = maps_f$sci_name, sp_factor = as.numeric(factor(maps_f$sci_name)))
usp <- unique(sp_f)[order(unique(sp_f$sp_factor)),]

#range(unique(maps_f$year))
#range(usp[,2])


# model input -------------------------------------------------------------


DATA <- data.frame(N = NROW(maps_f),
                   L = length(unique(maps_f$year)),
                   year = maps_f$year,
                   sp = sp_f$sp_factor,
                   
                   mass_obs = mass_obs,
                   mass_mis = mass_mis,
                   smass_obs = smass_obs,
                   smass_mis = smass_mis,
                   fat_obs = fat_obs,
                   fat_mis = fat_mis,
                   N_mass_obs = length(mass_obs),
                   N_mass_mis = length(mass_mis),
                   N_smass_obs = length(smass_obs),
                   N_smass_mis = length(smass_mis),
                   N_fat_obs = length(fat_obs),
                   N_fat_mis = length(fat_mis),
                   ii_mass_obs = NA,
                   ii_mass_mis = NA,
                   ii_smass_obs = NA,
                   ii_smass_mis = NA,
                   ii_fat_obs = NA,
                   ii_fat_mis = NA,
                   
                   #FIX
                   mass = maps_f$weight,
                   smass = maps_f$weight/maps_f$wing_chord,
                   fat = maps_f$fat_content)



# model -------------------------------------------------------------------

#varying slopes, varying intercept models
#need to accomodate missing data for years

MAPS_m <- '
data {
int<lower = 0> N;                                     // number of data points
int<lower = 0> J;                                     // number of species
int<lower = 0> N_mass_obs;
int<lower = 0> N_mass_mis;
int<lower = 0> N_smass_obs;
int<lower = 0> N_smass_mis;
int<lower = 0> N_fat_obs;
int<lower = 0> N_fat_mis;
int<lower = 1992, upper = 2017> year;                 // year
int<lower = 1, upper = 88> sp;                        // species id
real<lower = 0> mass_obs;                             // mass obs
real<lower = 0> smass_obs;                            // mass standardized by wing chord obs
int<lower = 0> fat_obs;                               // fat score obs
real<lower = 0> mass_mis;                             // mass NA
real<lower = 0> smass_mis;                            // mass standardized by wing chord NA
int<lower = 0> fat_mis;                               // fat score NA
real<lower = 1> ii_mass_obs[N_mass_obs];              //indices for mass obs
real<lower = 1> ii_mass_mis[N_mass_mis];              //indices for missing obs
real<lower = 1> ii_smass_obs[N_smass_obs];            //indices for mass obs
real<lower = 1> ii_smass_mis[N_smass_mis];            //indices for missing obs
real<lower = 1> ii_fat_obs[N_fat_obs];                //indices for mass obs
real<lower = 1> ii_fat_mis[N_fat_mis];                //indices for missing obs
}

parameters {
real alpha_mass[J];                                   // intercept mass
real beta_mass[J];                                    // slope mass
real alpha_smass[J];                                  // intercept smass
real beta_smass[J];                                   // slope smass
real alpha_fat[J];                                    // intercept fat
real beta_fat[J];                                     // slope fat

real<lower = 0> sigma_mass;
real<lower = 0> sigma_smass;
real<lower = 0> sigma_fat;

real mu_alpha_mass_raw;
real mu_alpha_smass_raw;
real mu_alpha_fat_raw;
real mu_beta_mass_raw;
real mu_beta_smass_raw;
real mu_beta_fat_raw;

real<lower = 0> sigma_alpha_mass;
real<lower = 0> sigma_beta_mass;
real<lower = 0> sigma_alpha_smass;
real<lower = 0> sigma_beta_smass;
real<lower = 0> sigma_alpha_fat;
real<lower = 0> sigma_beta_fat;
}

transformed parameters {
real<lower = 0, upper = 200> y[N, J];                 // response data to be modeled
real mu_mass[N];
real mu_smass[N];
real mu_fat[N];
real mu_alpha_mass;
real mu_alpha_smass;
real mu_alpha_fat;
real mu_beta_mass;
real mu_beta_smass;
real mu_beta_fat;
real mass[N];
real smass[N];
real fat[N];

mu_alpha_mass = mu_alpha_mass_raw * 10;
mu_alpha_smass = mu_alpha_smass_raw * 10;
mu_alpha_fat = mu_alpha_fat_raw * 10;

mu_beta_mass = mu_beta_mass_raw * 5;
mu_beta_smass = mu_beta_smass_raw * 5;
mu_beta_fat = mu_beta_dat_raw * 5;

for (j in 1:J)
{
  alpha_mass[j] ~ normal(mu_alpha_mass, sigma_alpha_mass)
  beta_mass[j] ~ normal(mu_beta_mass, sigma_beta_mass)
  
  alpha_smass[j] ~ normal(mu_alpha_smass, sigma_alpha_smass)
  beta_smass[j] ~ normal(mu_beta_smass, sigma_beta_smass)

  alpha_fat[j] ~ normal(mu_alpha_fat, sigma_alpha_fat)
  beta_fat[j] ~ normal(mu_beta_fat, sigma_beta_fat)
}


for (i in 1:N)
{
  mu_mass[i] = alpha_mass[sp[i]] + beta_mass[sp[i]] * year[sp[i]];
  mu_smass[i] = alpha_smass[sp[i]] + beta_smass[sp[i]] * year[sp[i]];
  mu_fat[i] = alpha_fat[sp[i]] + beta_fat[sp[i]] * year[sp[i]];
}

// indexing to avoid NAs
mass[ii_mass_obs] = mass_obs
mass[ii_mass_mis] = mass_mis
smass[ii_smass_obs] = smass_obs
smass[ii_smass_mis] = smass_mis
fat[ii_fat_obs] = fat_obs
fat[ii_fat_mis] = fat_mis
}

model {

// priors
mu_alpha_mass_raw ~ normal(0, 1);
mu_alpha_smass_raw ~ normal(0, 1);
mu_alpha_fat_raw ~ normal(0, 1);

sigma_alpha_mass ~ normal(0, 5);
sigma_beta_mass ~ normal(0, 5);
sigma_alpha_smass ~ normal(0, 5);
sigma_beta_smass ~ normal(0, 5);
sigma_alpha_fat ~ normal(0, 5);
sigma_beta_fat ~ normal(0, 5);

sigma_mass ~ normal(0, 5);
sigma_smass ~ normal(0, 5);
sigma_fat ~ normal(0, 5);

// model
mass ~ normal(mu_mass, sigma_mass);
smass ~ normal(mu_smass, sigma_smass);
fat ~ normal(mu_fat, sigma_fat);

}

generated quantities {

vector[N] mass_rep;
vector[N] smass_rep;
vector[N] fat_rep;

mass_rep = normal_rng(mu_mass, sigma_mass);
smass_rep = normal_rng(mu_smass, sigma_smass);
fat_rep = normal_rng(mu_fat, sigma_fat);
}'



# Run model ---------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.97
TREE_DEPTH <- 19
STEP_SIZE <- 0.0005
CHAINS <- 4
ITER <- 6000

tt <- proc.time()
fit <- rstan::stan(model_code = IAR_2,
            data = DATA,
            chains = CHAINS,
            iter = ITER,
            cores = CHAINS,
            pars = c('sigma_beta0', 'beta0',
                     'alpha_gamma', 'beta_gamma', 'sigma_gamma', 'gamma',
                     'sigma_nu', 'mu_sn', 'rho', 'nu', 'theta', 'phi', 
                     'y_true', 'y_rep'),
            control = list(adapt_delta = DELTA,
                           max_treedepth = TREE_DEPTH,
                           stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60