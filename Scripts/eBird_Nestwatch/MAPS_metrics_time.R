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
library(brms)


# load data ---------------------------------------------------------------

setwd(paste0(dir, '/..'))
#setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))

maps_data <- readRDS('MAPS_age_filled.rds')





# data processing --------------------------------------------------------

maps_adults <- dplyr::filter(maps_data, age %in% c('1', '5', '6', '7', '8'))

#QC data - remove zeros
to.rm <- which(maps_adults$weight == 0 | maps_adults$wing_chord == 0 | is.na(maps_adults$weight) | is.na(maps_adults$weight) | is.na(maps_adults$fat_content) | is.na(maps_adults$wing_chord))
maps_adults_qc <- maps_adults[-to.rm, ]


#only species with > 1000 data points
d_cnt <- plyr::count(maps_adults_qc, 'sci_name')
sp_p <- dplyr::filter(d_cnt, freq > 1000)[,1]

#subset of species
set.seed(1)
sp <- base::sample(sp_p, size = 6)
#sp <- sp_p

maps_f <- dplyr::filter(maps_adults_qc, sci_name %in% sp)

#factor for species_id
sp_f <- data.frame(sp = maps_f$sci_name, sp_factor = as.numeric(factor(maps_f$sci_name)))
usp <- unique(sp_f)[order(unique(sp_f$sp_factor)),]

#range(unique(maps_f$year))
#range(usp[,2])

maps_f$year_f <- as.numeric(factor(maps_f$year))
maps_f$sweight <- maps_f$weight/maps_f$wing_chord



# plot and lm -------------------------------------------------------------

# ggplot(maps_f, aes(year_f, sweight, col = sci_name)) +
#   geom_point(alpha = 0.2) +
#   theme(legend.position="none") +
#   geom_smooth(method='lm') +
#   ylim(c(0, 3))
# 
# summary(lm(maps_f$sweight ~ maps_f$year_f))




# model input -------------------------------------------------------------


DATA <- list(N = NROW(maps_f),
             J = NROW(usp),
             year = maps_f$year_f,
             sp = sp_f$sp_factor,
             mass = maps_f$sweight,
             fat = maps_f$fat_content)



# model -------------------------------------------------------------------

#varying slopes, varying intercept models
#need to accomodate missing data for years

wf_time <- '
data {
int<lower = 0> N;                                     // number of data points
int<lower = 0> J;                                     // number of species
int<lower = 1, upper = 26> year[N];                   // year id
int<lower = 1> sp[N];                                 // species id
real<lower = 0> mass[N];                              // mass standardized by wing chord
real<lower = 0> fat[N];                               // fat score
}

parameters {
real alpha_mass[J];                                   // intercept mass
real beta_mass[J];                                    // slope mass
// real alpha_fat[J];                                    // intercept fat
// real beta_fat[J];                                     // slope fat

real<lower = 0> sigma_mass;
// real<lower = 0> sigma_fat;

real mu_alpha_mass_raw;
// real mu_alpha_fat_raw;
real mu_beta_mass_raw;
// real mu_beta_fat_raw;

real<lower = 0> sigma_alpha_mass;
real<lower = 0> sigma_beta_mass;
// real<lower = 0> sigma_alpha_fat;
// real<lower = 0> sigma_beta_fat;
}

transformed parameters {
real mu_mass[N];
//real mu_fat[N];
real mu_alpha_mass;
// real mu_alpha_fat;
real mu_beta_mass;
// real mu_beta_fat;

mu_alpha_mass = mu_alpha_mass_raw * 10;
// mu_alpha_fat = mu_alpha_fat_raw * 10;

mu_beta_mass = mu_beta_mass_raw * 5;
// mu_beta_fat = mu_beta_fat_raw * 5;

for (i in 1:N)
{
  mu_mass[i] = alpha_mass[sp[i]] + beta_mass[sp[i]] * year[sp[i]];
  // mu_fat[i] = alpha_fat[sp[i]] + beta_fat[sp[i]] * year[sp[i]];
}

}

model {

// priors
mu_alpha_mass_raw ~ normal(0, 1);
// mu_alpha_fat_raw ~ normal(0, 1);
mu_beta_mass_raw ~ normal(0, 1);
// mu_beta_fat_raw ~ normal(0, 1);

sigma_alpha_mass ~ normal(0, 5);
sigma_beta_mass ~ normal(0, 5);
// sigma_alpha_fat ~ normal(0, 5);
// sigma_beta_fat ~ normal(0, 5);

sigma_mass ~ normal(0, 5);
// sigma_fat ~ normal(0, 5);


for (j in 1:J)
{
  alpha_mass[j] ~ normal(mu_alpha_mass, sigma_alpha_mass);
  beta_mass[j] ~ normal(mu_beta_mass, sigma_beta_mass);
  
  // alpha_fat[j] ~ normal(mu_alpha_fat, sigma_alpha_fat);
  // beta_fat[j] ~ normal(mu_beta_fat, sigma_beta_fat);
}


// model
mass ~ normal(mu_mass, sigma_mass);
// fat ~ normal(mu_fat, sigma_fat);

}

generated quantities {

// real mass_rep[N];
// real fat_rep[N];

// mass_rep = normal_rng(mu_mass, sigma_mass);
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
            pars = c('alpha_mass', #'alpha_fat',
                     'beta_mass', #'beta_fat',
                     'mu_alpha_mass', #'mu_alpha_fat',
                     'mu_beta_mass', #'mu_beta_fat',
                     'sigma_mass'))#, 'sigma_fat'),
                     #'mu_mass', 'mu_fat',
                     #'mass_rep', 'fat_rep'),
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

MCMCvis::MCMCsummary(fit, round = 2, n.eff = TRUE)
MCMCvis::MCMCplot(fit, params = 'beta_mass')
MCMCvis::MCMCplot(fit, params = 'beta_smass')
MCMCvis::MCMCplot(fit, params = 'beta_fat')

# library(shinystan)
# launch_shinystan(fit)







# process data for brms ---------------------------------------------------

brms_data <- maps_f

#factors for sp and fat
brms_data$sp_f <- factor(brms_data$sci_name)
brms_data$fat_f <- ordered(brms_data$fat_content)



# fat model -------------------------------------------------------------------

#varying slopes, varying intercept models

#random slopes, random intercepts
#fat ~ year + (year | sp)

DELTA <- 0.90
TREE_DEPTH <- 15
STEP_SIZE <- 0.05
CHAINS <- 4
ITER <- 2000

tt <- proc.time()
b_fit_fat <- brms::brm(formula = fat_f ~ year_f + (year_f | sp_f), 
                       data = brms_data, 
                       family = cumulative('logit'),
                       cores = CHAINS,
                       chains = CHAINS,
                       iter = ITER,
                       control = list(adapt_delta = DELTA,
                                      max_treedepth = TREE_DEPTH,
                                      stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
saveRDS(b_fit_fat, file = 'MAPS-fat-time-brms_output.rds')

#b_fit_fat$model

#https://groups.google.com/forum/#!msg/brms-users/2dH6EawVrtM/X7NPi0CACAAJ
#https://discourse.mc-stan.org/t/interpretation-brms-results/6374/17
#https://discourse.mc-stan.org/t/coef-versus-fixef-and-ranef/3914/8
#https://www.statisticssolutions.com/ordinal-regression-2/


#non bayes shows that there is a decline over time in fat content
polr_fit <- MASS::polr(formula = fat_f ~ year_f, 
                       data = brms_data, 
                       Hess = TRUE,
                       method = 'logistic')
summary(polr_fit)

#loglog, also known as negative loglog, used when lots of low cats
polr_fit <- MASS::polr(formula = fat_f ~ year_f, 
                       data = brms_data, 
                       Hess = TRUE,
                       method = 'loglog')
summary(polr_fit)

clm_fit <- ordinal::clmm(fat_f ~ year_f + (year_f | sp_f),  
                         data = brms_data,
                         link = 'loglog')
summary(clm_fit)





# weight model with brms --------------------------------------------------

DELTA <- 0.90
TREE_DEPTH <- 15
STEP_SIZE <- 0.05
CHAINS <- 4
ITER <- 2000

b_fit_sweight <- brms::brm(formula = sweight ~ year_f + (year_f | sp_f), 
                          data = brms_data, 
                          cores = CHAINS,
                          chains = CHAINS,
                          iter = ITER,
                          control = list(adapt_delta = DELTA,
                                         max_treedepth = TREE_DEPTH,
                                         stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
saveRDS(b_fit_sweight, file = 'MAPS-sweight-time-brms_output.rds')


#b_fit_weight$model
summary(b_fit_weight)
