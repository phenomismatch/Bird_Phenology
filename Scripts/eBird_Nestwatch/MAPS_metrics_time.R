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

setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))

maps_data <- readRDS('MAPS_age_filled.rds')





# data processing --------------------------------------------------------

maps_adults <- dplyr::filter(maps_data, age %in% c('1', '5', '6', '7', '8'))

#QC data - remove zeros
to.rm <- which(maps_adults$weight == 0 | maps_adults$wing_chord == 0 | is.na(maps_adults$weight) | is.na(maps_adults$fat_content))
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
                   mass = maps_f$weight,
                   smass = maps_f$weight/maps_f$wing_chord,
                   fat = maps_f$fat_content)



# model -------------------------------------------------------------------

#varying slopes, varying intercept models
#need to accomodate missing data for years

MAPS_m <- '
data {
int<lower = 0> N;                                     // number of data points
int<lower = 0> L;                                     // number of years
int<lower = 1992, upper=2017> year;                   // year
int<lower = 1, upper = 88> sp;                        // species id
real<lower = 0> mass;                                 // mass
real<lower = 0> smass;                                // mass standardized by wing chord
int<lower = 0> fat;                                   // fat score
}

parameters {
real<lower = 1, upper = 200> y_mis[N, J];             // missing response data
real alpha_mass;                                      // intercept mass
real beta_mass;                                       // slope mass
real alpha_smass;                                     // intercept smass
real beta_smass;                                      // slope smass
real alpha_fat;                                       // intercept fat
real beta_fat;                                        // slope fat

real sigma_mass;
real sigma_smass;
real sigma_fat;

real mu_alpha_mass;
real mu_beta_mass;
real mu_alpha_smass;
real mu_beta_smass;
real mu_alpha_fat;
real mu_beta_fat;

real sigma_alpha_mass;
real sigma_beta_mass;
real sigma_alpha_smass;
real sigma_beta_smass;
real sigma_alpha_fat;
real sigma_beta_fat;
}

transformed parameters {
real<lower = 0, upper = 200> y[N, J];                 // response data to be modeled
vector[N] gamma;
real alpha_gamma;
real beta_gamma;
real<lower = 0> sigma_gamma;
real mu_gamma[N];
vector<lower = 0>[J] sigma_nu;
matrix[N, J] y_true;
matrix[N, J] nu;                            // spatial and non-spatial component
real mu_sn;
real<lower = 0> sigma_beta0;
real beta0[J];

alpha_gamma = alpha_gamma_raw * 30;
beta_gamma = beta_gamma_raw * 3 + 2;
sigma_gamma = sigma_gamma_raw * 5;
mu_sn = mu_sn_raw * 1.5;
sigma_beta0 = sigma_beta0_raw * 5;


for (j in 1:sp)
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
for (j in 1:J)
{
  y[ii_obs[1:N_obs[j], j], j] = y_obs[1:N_obs[j], j];
  y[ii_mis[1:N_mis[j], j], j] = y_mis[1:N_mis[j], j];
}
}

model {

// priors

mass ~ normal(mu_mass, sigma_mass);
smass ~ normal(mu_smass, sigma_smass);
fat ~ normal(mu_fat, sigma_fat);

}

generated quantities {

// real y_rep[N, J];
vector[NJ] y_rep;
int<lower = 0> counter;

counter = 1;
for (j in 1:J)
{
  for (n in 1:N)
  {
  y_rep[counter] = normal_rng(y_true[n,j], sigma_y[n,j]);
  counter = counter + 1;
  }
}
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