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

maps_f$year_f <- as.numeric(factor(maps_f$year))

# model input -------------------------------------------------------------


DATA <- list(N = NROW(maps_f),
             J = NROW(usp),
             year = maps_f$year_f,
             sp = sp_f$sp_factor,
             mass = maps_f$weight,
             smass = maps_f$weight/maps_f$wing_chord,
             fat = maps_f$fat_content)



# model -------------------------------------------------------------------

#varying slopes, varying intercept models
#need to accomodate missing data for years

wf_time <- '
data {
int<lower = 0> N;                                     // number of data points
int<lower = 0> J;                                     // number of species
int<lower = 1, upper = 26> year[N];                   // year id
int<lower = 1, upper = 88> sp[N];                     // species id
real<lower = 0> mass[N];                              // mass
real<lower = 0> smass[N];                             // mass standardized by wing chord
real<lower = 0> fat[N];                                // fat score
}

parameters {
real alpha_mass[J];                                   // intercept mass
real beta_mass[J];                                    // slope mass
// real alpha_smass[J];                                  // intercept smass
// real beta_smass[J];                                   // slope smass
// real alpha_fat[J];                                    // intercept fat
// real beta_fat[J];                                     // slope fat

real<lower = 0> sigma_mass;
// real<lower = 0> sigma_smass;
// real<lower = 0> sigma_fat;

real mu_alpha_mass_raw;
// real mu_alpha_smass_raw;
// real mu_alpha_fat_raw;
real mu_beta_mass_raw;
// real mu_beta_smass_raw;
// real mu_beta_fat_raw;

real<lower = 0> sigma_alpha_mass;
real<lower = 0> sigma_beta_mass;
// real<lower = 0> sigma_alpha_smass;
// real<lower = 0> sigma_beta_smass;
// real<lower = 0> sigma_alpha_fat;
// real<lower = 0> sigma_beta_fat;
}

transformed parameters {
real mu_mass[N];
// real mu_smass[N];
//real mu_fat[N];
real mu_alpha_mass;
// real mu_alpha_smass;
// real mu_alpha_fat;
real mu_beta_mass;
// real mu_beta_smass;
// real mu_beta_fat;

mu_alpha_mass = mu_alpha_mass_raw * 10;
// mu_alpha_smass = mu_alpha_smass_raw * 10;
// mu_alpha_fat = mu_alpha_fat_raw * 10;

mu_beta_mass = mu_beta_mass_raw * 5;
// mu_beta_smass = mu_beta_smass_raw * 5;
// mu_beta_fat = mu_beta_fat_raw * 5;

for (i in 1:N)
{
  mu_mass[i] = alpha_mass[sp[i]] + beta_mass[sp[i]] * year[sp[i]];
  // mu_smass[i] = alpha_smass[sp[i]] + beta_smass[sp[i]] * year[sp[i]];
  // mu_fat[i] = alpha_fat[sp[i]] + beta_fat[sp[i]] * year[sp[i]];
}

}

model {

// priors
mu_alpha_mass_raw ~ normal(0, 1);
// mu_alpha_smass_raw ~ normal(0, 1);
// mu_alpha_fat_raw ~ normal(0, 1);
mu_beta_mass_raw ~ normal(0, 1);
// mu_beta_smass_raw ~ normal(0, 1);
// mu_beta_fat_raw ~ normal(0, 1);

sigma_alpha_mass ~ normal(0, 5);
sigma_beta_mass ~ normal(0, 5);
// sigma_alpha_smass ~ normal(0, 5);
// sigma_beta_smass ~ normal(0, 5);
// sigma_alpha_fat ~ normal(0, 5);
// sigma_beta_fat ~ normal(0, 5);

sigma_mass ~ normal(0, 5);
// sigma_smass ~ normal(0, 5);
// sigma_fat ~ normal(0, 5);


for (j in 1:J)
{
  alpha_mass[j] ~ normal(mu_alpha_mass, sigma_alpha_mass);
  beta_mass[j] ~ normal(mu_beta_mass, sigma_beta_mass);
  
  // alpha_smass[j] ~ normal(mu_alpha_smass, sigma_alpha_smass);
  // beta_smass[j] ~ normal(mu_beta_smass, sigma_beta_smass);
  
  // alpha_fat[j] ~ normal(mu_alpha_fat, sigma_alpha_fat);
  // beta_fat[j] ~ normal(mu_beta_fat, sigma_beta_fat);
}


// model
mass ~ normal(mu_mass, sigma_mass);
// smass ~ normal(mu_smass, sigma_smass);
// fat ~ normal(mu_fat, sigma_fat);

}

generated quantities {

// real mass_rep[N];
// real smass_rep[N];
// real fat_rep[N];

// mass_rep = normal_rng(mu_mass, sigma_mass);
// smass_rep = normal_rng(mu_smass, sigma_smass);
// fat_rep = normal_rng(mu_fat, sigma_fat);
}'



# Run model ---------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.90
TREE_DEPTH <- 15
STEP_SIZE <- 0.05
CHAINS <- 4
ITER <- 1000

tt <- proc.time()
fit <- rstan::stan(model_code = wf_time,
            data = DATA,
            chains = CHAINS,
            iter = ITER,
            cores = CHAINS,
            pars = c('alpha_mass', #'alpha_smass', 'alpha_fat',
                     'beta_mass', #'beta_smass', 'beta_fat',
                     'mu_alpha_mass', #'mu_alpha_smass', 'mu_alpha_fat',
                     'mu_beta_mass', #'mu_beta_smass', 'mu_beta_fat',
                     'sigma_mass'))#, 'sigma_smass', 'sigma_fat'),
                     #'mu_mass', 'mu_smass', 'mu_fat',
                     #'mass_rep', 'smass_rep', 'fat_rep'),
            # control = list(adapt_delta = DELTA,
            #                max_treedepth = TREE_DEPTH,
            #                stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60



#setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
#saveRDS(fit, file = 'MAPS-metrics-time-stan_output.rds'))

#save data to RDS
#saveRDS(DATA, file = 'MAPS-metrics-time-stan_input.rds'))



# Calc diagnostics and rerun if needed ------------------------------------

num_diverge <- rstan::get_num_divergent(fit)
num_tree <- rstan::get_num_max_treedepth(fit)
num_BFMI <- length(rstan::get_low_bfmi_chains(fit))



# analyze data ------------------------------------------------------------

MCMCvis::MCMCsummary(fit, round = 2)
MCMCvis::MCMCplot(fit, params = 'beta_mass')
MCMCvis::MCMCplot(fit, params = 'beta_smass')
MCMCvis::MCMCplot(fit, params = 'beta_fat')

# library(shinystan)
# launch_shinystan(fit)





# fat model -------------------------------------------------------------------

#varying slopes, varying intercept models

#random slopes, random intercepts
#fat ~ year + (year | sp)

#only species with > 9k data points
d_cnt <- plyr::count(maps_adults_qc, 'sci_name')
sp <- dplyr::filter(d_cnt, freq > 9000)[,1]

brms_data <- dplyr::filter(maps_adults_qc, sci_name %in% sp)

#add ids for species and year
brms_data$spf <- factor(brms_data$sci_name)
brms_data$yearf <- as.numeric(factor(brms_data$year))
brms_data$fatf <- ordered(brms_data$fat_content)

b_fit_fat <- brms::brm(formula = fatf ~ yearf + (yearf | spf), 
                       data = brms_data, 
                       family = cumulative('logit'),
                       cores = 3,
                       chains = 3,
                       iter = 1000)

#https://groups.google.com/forum/#!msg/brms-users/2dH6EawVrtM/X7NPi0CACAAJ
#https://discourse.mc-stan.org/t/interpretation-brms-results/6374/17
#https://discourse.mc-stan.org/t/coef-versus-fixef-and-ranef/3914/8

# weight model with brms --------------------------------------------------

b_fit_weight <- brms::brm(formula = weight ~ yearf + (yearf | spf), 
                          data = brms_data, 
                          cores = 3,
                          chains = 3,
                          iter = 1000)


