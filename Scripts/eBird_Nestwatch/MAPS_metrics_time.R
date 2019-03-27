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
maps_f$sp_f <- factor(maps_f$sci_name)
usp <- unique(maps_f$sci_name)

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

brms_data <- dplyr::filter(maps_f, day <= 160, day >=130)

#factors for sp and fat
brms_data$fat_f <- ordered(brms_data$fat_content)
brms_data$st_f <- as.numeric(factor(brms_data$station))


# fat model -------------------------------------------------------------------

#varying slopes, varying intercept models

#random slopes, random intercepts
#fat ~ year + (year | sp)

DELTA <- 0.90
TREE_DEPTH <- 15
STEP_SIZE <- 0.05
CHAINS <- 1
ITER <- 10

tt <- proc.time()
b_fit_fat <- brms::brm(formula = fat_f ~ year_f, 
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

summary(b_fit_fat)

#b_fit_fat$model

#https://groups.google.com/forum/#!msg/brms-users/2dH6EawVrtM/X7NPi0CACAAJ
#https://discourse.mc-stan.org/t/interpretation-brms-results/6374/17
#https://discourse.mc-stan.org/t/coef-versus-fixef-and-ranef/3914/8
#https://www.statisticssolutions.com/ordinal-regression-2/


#non bayes shows that there is a decline over time in fat content
polr_fit <- MASS::polr(formula = fat_f ~ year_f + lat + poly(day, 2), 
                       data = brms_data, 
                       Hess = TRUE,
                       method = 'logistic')
summary(polr_fit)

lm_fit <- lm(sweight ~ year_f, data = brms_data)
plot(brms_data$year_f, brms_data$sweight, ylim = c(0, 0.6))
abline(lm_fit, col = 'red')
summary(lm_fit)
lm_fit <- lm(wing_chord ~ year_f, data = brms_data)
plot(brms_data$wing_chord, brms_data$fat_content)
lm_fit <- lm(fat_content ~ wing_chord, data = brms_data)
summary(lm_fit)
abline(lm_fit)


#There is an increase in wing_chord over time

#decreasing fat but increasing wing chord over time
species <- unique(brms_data$sp_f)
for (i in 1:length(species))
{
  #i <- 1
  temp <- dplyr::filter(brms_data, sci_name == species[i])
  
  # plot(temp$wing_chord, temp$weight, xlim = c(35, 70), ylim = c(5, 20))
  # plot(temp$lat, temp$wing_chord)
  # summary(lm(wing_chord ~ lat + weight, data = temp))
  # summary(lm(wing_chord ~ year_f + lat + year_f*lat, data = temp))
  # plot(temp$lat, temp$weight)
  # summary(lm(weight ~ lat, data = temp))
  #summary(lm(weight ~ wing_chord, data = temp))
  # boxplot(temp$fat_content ~ temp$year)
  # boxplot(temp$sweight ~ temp$year, ylim = c(0.1, 0.3))
  # boxplot(temp$weight ~ temp$year, ylim = c(7, 16))
  # boxplot(temp$wing_chord ~ temp$year, ylim = c(45, 60))
  # boxplot(temp$wing_chord ~ temp$day, ylim = c(45, 60))
  # plot(temp$wing_chord ~ temp$day)
  # plot(temp$wing_chord ~ temp$lat)
  # plot(temp$fat_content ~ temp$lat)
  # plot(temp$sweight ~ temp$lat)
  # summary(lm(temp$sweight ~ temp$lat))
  # summary(lm(temp$wing_chord ~ temp$lat))
  # plot(temp$fat_content ~ temp$lat)
  # plot(temp$fat_content ~ temp$wing_chord)
  # boxplot(temp$fat_content ~ temp$day)
  # boxplot(temp$day ~ temp$year)
  # yrs <- sort(unique(temp$year))
  # mf <- c()
  # plot(temp$day, temp$fat_content, col = 'white')
  # colors <- colorspace::sequential_hcl(length(yrs), alpha = 0.5)
  # plyr::count(temp, 'year')
  
  # for (j in length(yrs):1)
  # {
  #   #j <- 24
  #   t2 <- dplyr::filter(temp, year == yrs[j])
  #   mf <- c(mf, quantile(t2$fat_content, probs = c(0.6)))
  #   plot(t2$day, t2$fat_content, col = colors[j], lwd = 2, ylim = c(0, 7))
  #   #points(t2$day, t2$fat_content, col = colors[j], lwd = 2)
  #   #hist(sqrt(t2$fat_content), xlim = c(0, 7))
  #   #Sys.sleep(1)
  # }
  # 
  # plot(temp$year, temp$fat_content)
  # plot(temp$year, temp$day)
  
  polr_fit2 <- MASS::polr(formula = fat_f ~ year_f + poly(day, 2), 
                         data = temp, 
                         Hess = TRUE,
                         method = 'logistic')
  print(summary(polr_fit2)$coef[1,])
  
  lm_fit <- lm(formula = wing_chord ~ year_f, 
                         data = temp)
  print(summary(lm_fit)$coef[2,])
}


#loglog, also known as negative loglog, used when lots of low cats
polr_fit <- MASS::polr(formula = fat_f ~ year_f + poly(day, 2), 
                       data = brms_data, 
                       Hess = TRUE,
                       method = 'loglog')
summary(polr_fit)

clm_fit <- ordinal::clmm(fat_f ~ year_f + (year_f | sp_f),  
                         data = brms_data,
                         link = 'loglog')
summary(clm_fit)




# fat stan model ----------------------------------------------------------

#https://groups.google.com/forum/#!category-topic/stan-users/general/sgX2Edo8qiQ
#page 28 Stan users guide for accessing array/vector

#count # of obs for each year
plyr::count(brms_data, c('st_f', 'year_f'))

#fat over julian day - should increase 2nd order poly
plot(brms_data$day, brms_data$fat_content)
fd_fit <- lm(fat_content ~ day, data = brms_data)
summary(fd_fit)
abline(fd_fit, col = 'red', lwd = 3)
fd2_fit <- lm(fat_content ~ poly(day,2,raw=TRUE), data = brms_data)
summary(fd2_fit)
dd <- 120:160
lines(dd, predict(fd2_fit, data.frame(day = dd)), col="blue", lwd = 3)
AIC(fd_fit)
AIC(fd2_fit)

#lat over time

boxplot(brms_data$lat ~ brms_data$year)
boxplot(brms_data$fat_content ~ brms_data$year)
ly_fit <- lm(lat ~ year, data = brms_data)
abline(ly_fit, col = 'red', lwd = 3)

#fat over lat
plot(brms_data$lat, brms_data$fat_content)
fl_fit <- lm(fat_content ~ lat, data = brms_data)
abline(fl_fit, col = 'red', lwd = 3)

#julian day over time
plot(brms_data$year, brms_data$day)
summary(lm(day ~ year, data = brms_data))



DATA <- list(K = length(unique(brms_data$fat_content)),
             N = NROW(brms_data),
             Nsp = length(usp),
             Nst = length(brms_data$st_f),
             y = brms_data$fat_content + 1,
             sp = as.numeric(brms_data$sp_f), #species
             st = brms_data$st_f, #station
             year = brms_data$year_f,
             day = brms_data$day,
             lat = brms_data$lat)

stanmodel <- "
data {
  int<lower=0> K;                     // number of bins
  int<lower=0> N;                     // number of obs
  int<lower=0> Nsp;                   // number of species
  int<lower=0> Nst;                   // number of stations
  int<lower=1, upper=K> y[N];         // response
  int<lower=1, upper=Nsp> sp[N];      // groups
  int<lower=1, upper=Nst> st[N];      // stations
  vector<lower=0>[N] year;
  vector<lower=0>[N] day;
  vector<lower=0>[N] lat;
}

parameters {
  ordered[K-1] c;
  vector[Nsp] eta_raw;
  real gamma_raw;
  real nu_1_raw;
  real nu_2_raw;
  real beta_raw;
}

transformed parameters {
  vector[N] mu;
  vector[Nsp] eta;
  real gamma;
  real nu_1;
  real nu_2;
  real beta;

  eta = eta_raw * 3;              // eta ~ normal(0, 3)
  gamma = gamma_raw * 2;          // gamma ~ normal(0, 2)
  nu_1 = nu_1_raw * 2;            // nu_1 ~ normal(0, 2)
  nu_2 = nu_2_raw * 2;        // nu_2 ~ normal(0, 2)
  beta = beta_raw * 2;            // beta ~ normal(0, 2)

  for (i in 1:N)
  {
    mu[i] = eta[sp[i]] + gamma * lat[i] + nu_1 * day[i] + nu_2 * day[i]^2 + beta * year[i];
  }
}

model {
  eta_raw ~ normal(0, 1);
  gamma_raw ~ normal(0, 1);
  nu_1_raw ~ normal(0, 1);
  nu_2_raw ~ normal(0, 1);
  beta_raw ~ normal(0, 1);

  y ~ ordered_logistic(mu, c);
}
"

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.95
TREE_DEPTH <- 16
STEP_SIZE <- 0.01
CHAINS <- 4
ITER <- 2000

tt <- proc.time()
fit <- rstan::stan(model_code = stanmodel,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('eta',
                            'gamma',
                            'nu_1',
                            'nu_2',
                            'beta',
                            'c'), 
                   control = list(adapt_delta = DELTA,
                   max_treedepth = TREE_DEPTH,
                   stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60


setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
saveRDS(fit, file = 'MAPS-fat-time-stan_output.rds')


library(shinystan)
MCMCvis::MCMCsummary(fit, round = 2, n.eff = TRUE)


# fat latitude model ------------------------------------------------------

#varying slopes, varying intercept models

DELTA <- 0.95
TREE_DEPTH <- 15
STEP_SIZE <- 0.05
CHAINS <- 4
ITER <- 2000

tt <- proc.time()
b_fit_fat_lat <- brms::brm(formula = fat_f ~ lat, 
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
saveRDS(b_fit_fat_lat, file = 'MAPS-fat-lat-brms_output.rds')

summary(b_fit_fat_lat)

#PPC
pp_check(b_fit_fat_lat) #dens overlay
pp_check(b_fit_fat_lat, type = 'error_hist')
pp_check(b_fit_fat_lat, type = 'rootogram')

b_fit_fat_lat$model

# wing chord latitude model ------------------------------------------------------

#varying slopes, varying intercept models

DELTA <- 0.95
TREE_DEPTH <- 15
STEP_SIZE <- 0.05
CHAINS <- 4
ITER <- 2000

tt <- proc.time()
b_fit_wc_lat <- brms::brm(formula = wing_chord ~ lat + (lat | sp_f), 
                           data = brms_data, 
                           cores = CHAINS,
                           chains = CHAINS,
                           iter = ITER,
                           control = list(adapt_delta = DELTA,
                                          max_treedepth = TREE_DEPTH,
                                          stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
saveRDS(b_fit_wc_lat, file = 'MAPS-wc-lat-brms_output.rds')

summary(b_fit_wc_lat)

#PPC
pp_check(b_fit_wc_lat) #dens overlay
pp_check(b_fit_wc_lat, type = 'error_hist')
pp_check(b_fit_wc_lat, type = 'intervals')

# weight model with brms --------------------------------------------------

DELTA <- 0.92
TREE_DEPTH <- 16
STEP_SIZE <- 0.05
CHAINS <- 4
ITER <- 2000

tt <- proc.time()
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
summary(b_fit_sweight)










# wing chord --------------------------------------------------------------

maps_adults <- dplyr::filter(maps_data, age %in% c('1', '5', '6', '7', '8'))

#QC data - remove zeros
to.rm <- which(maps_adults$weight == 0 | maps_adults$wing_chord == 0 | is.na(maps_adults$weight) | is.na(maps_adults$weight) | is.na(maps_adults$fat_content) | is.na(maps_adults$wing_chord))
maps_adults_qc <- maps_adults[-to.rm, ]

#only species with > 1000 data points
d_cnt <- plyr::count(maps_adults_qc, 'sci_name')
sp_p <- dplyr::filter(d_cnt, freq > 1000)[,1]

#subset of species
set.seed(1)
#sp <- base::sample(sp_p, size = 30)
sp <- sp_p

maps_f <- dplyr::filter(maps_adults_qc, sci_name %in% sp)

#factor for species_id
maps_f$sp_f <- factor(maps_f$sci_name)
usp <- unique(maps_f$sci_name)

maps_f$year_f <- as.numeric(factor(maps_f$year))
maps_f$sweight <- maps_f$weight/maps_f$wing_chord

#brms_data <- dplyr::filter(maps_f, day <= 160, day >=130)

#only males age 2+
brms_data <- dplyr::filter(maps_f, sex == 'M', age >= 5)

#factors for sp and fat
brms_data$fat_f <- ordered(brms_data$fat_content)
brms_data$st_f <- as.numeric(factor(brms_data$station))





DATA <- list(K = length(unique(brms_data$fat_content)),
             N = NROW(brms_data),
             Nsp = length(usp),
             Nst = length(brms_data$st_f),
             y = brms_data$wing_chord,
             sp = as.numeric(brms_data$sp_f), #species
             st = brms_data$st_f, #station
             year = brms_data$year_f,
             day = brms_data$day,
             lat = brms_data$lat)

stanmodel2 <- "
data {
int<lower=0> K;                     // number of bins
int<lower=0> N;                     // number of obs
int<lower=0> Nsp;                   // number of species
int<lower=0> Nst;                   // number of stations
real<lower=0> y[N];                 // response
int<lower=1, upper=Nsp> sp[N];      // groups
int<lower=1, upper=Nst> st[N];      // stations
vector<lower=0>[N] year;
vector<lower=0>[N] day;
vector<lower=0>[N] lat;
}

parameters {
vector[Nsp] eta_raw;
real gamma_raw;
real beta_raw;
real rho_raw;
real<lower = 0> sigma_raw;
}

transformed parameters {
vector[N] mu;
vector[Nsp] eta;
real gamma;
real beta;
real rho;
real<lower = 0> sigma;

eta = eta_raw * 20;              // eta ~ normal(0, 3)
gamma = gamma_raw * 2;          // gamma ~ normal(0, 2)
beta = beta_raw * 2;            // beta ~ normal(0, 2)
rho = rho_raw * 2;              // rho ~ normal(0, 2)
sigma = sigma_raw * 5;

for (i in 1:N)
{
  mu[i] = eta[sp[i]] + gamma * lat[i] + beta * year[i] + rho * (year[i] * lat[i]);
}
}

model {
eta_raw ~ normal(0, 1);
gamma_raw ~ normal(0, 1);
beta_raw ~ normal(0, 1);
rho_raw ~ normal(0, 1);
sigma_raw ~ normal(0, 1);

y ~ normal(mu, sigma);
}

generated quantities {
real y_rep[N];

y_rep = normal_rng(mu, sigma);
}
"

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.93
TREE_DEPTH <- 16
STEP_SIZE <- 0.005
CHAINS <- 4
ITER <- 3000

tt <- proc.time()
fit2 <- rstan::stan(model_code = stanmodel2,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('eta',
                            'gamma',
                            'beta',
                            'rho',
                            'sigma',
                            'y_rep'), 
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
saveRDS(fit2, file = 'MAPS-wc-time-lat-stan_output.rds')

MCMCvis::MCMCsummary(fit2, n.eff = TRUE, round = 2, excl = 'y_rep')
MCMCvis::MCMCplot(fit2, excl = c('eta', 'lp__'))

library(shinystan)
launch_shinystan(fit2)
?ppc_dens_overlay
# y_rep <- MCMCvis::MCMCchains(fit2, params = 'y_rep')
# bayesplot::ppc_dens_overlay(DATA$y, y_rep)
# bayesplot::ppc_dens_overlay(y, yrep_poisson[1:50, ])





# wing chord --------------------------------------------------------------

maps_adults <- dplyr::filter(maps_data, age %in% c('1', '5', '6', '7', '8'))

#QC data - remove zeros
to.rm <- which(maps_adults$weight == 0 | maps_adults$wing_chord == 0 | is.na(maps_adults$weight) | is.na(maps_adults$weight) | is.na(maps_adults$fat_content) | is.na(maps_adults$wing_chord))
maps_adults_qc <- maps_adults[-to.rm, ]

#only species with > 1000 data points
d_cnt <- plyr::count(maps_adults_qc, 'sci_name')
sp_p <- dplyr::filter(d_cnt, freq > 1000)[,1]

#subset of species
set.seed(1)
sp <- base::sample(sp_p, size = 8)
#sp <- sp_p

maps_f <- dplyr::filter(maps_adults_qc, sci_name %in% sp)

#factor for species_id
maps_f$sp_f <- factor(maps_f$sci_name)
usp <- unique(maps_f$sci_name)

#numeric year
maps_f$year_f <- as.numeric(factor(maps_f$year))

#0/1 sex - M = 1
maps_f$sex_f <- 0
maps_f$sex_f[which(maps_f$sex == 'M')] <- 1

#only birds 2+ years of age (bc of changing wing chord from Age 1 to 2)
#age here is age code, with >= 5 denoting ASY ()
stan_data <- dplyr::filter(maps_f, age >= 5)

#filter by time period
#brms_data <- dplyr::filter(maps_f, day <= 160, day >=130)


DATA <- list(N = NROW(stan_data),
             Nsp = length(usp),
             Nst = length(stan_data$st_f),
             y = stan_data$wing_chord,
             sp = as.numeric(stan_data$sp_f), #species
             year = stan_data$year_f,
             day = stan_data$day,
             lat = stan_data$lat,
             sex = stan_data$sex_f)

#NEXT:
#hierarchical gamma (effect of lat)
#hierarchical nu (effect of sex)

stanmodel3 <- "
data {
int<lower=0> N;                     // number of obs
int<lower=0> Nsp;                   // number of species
real<lower=0> y[N];                 // response
int<lower=1, upper=Nsp> sp[N];      // species
vector<lower=0>[N] year;
vector<lower=0>[N] day;
vector<lower=0>[N] lat;
int<lower=0, upper=1> sex[N];
}

parameters {
real mu_alpha;
real mu_beta;
real<lower = 0> sigma_raw;
vector<lower = 0>[2] sigma_sp;                    // standard deviations
cholesky_factor_corr[2] L_Rho;                    // correlation matrix
matrix[2, Nsp] z;
real gamma_raw;
real eta_raw;
real nu_raw;
}

transformed parameters {
vector[N] mu;
matrix[Nsp, 2] alphabeta;                         // matrix for alpha and beta
matrix[2, 2] Rho;                                 // covariance matrix
real<lower = 0> sigma;
real gamma;
real eta;
real nu;
vector[Nsp] alpha;
vector[Nsp] beta;

sigma = sigma_raw * 5;
gamma = gamma_raw * 3;
eta = eta_raw * 3;
nu = nu_raw * 5;

alphabeta = (diag_pre_multiply(sigma_sp, L_Rho) * z)';  // cholesky factor of covariance matrix multiplied by z score
alpha = alphabeta[,1];
beta = alphabeta[,2];
Rho = L_Rho * L_Rho';

for (i in 1:N)
{
  mu[i] = (mu_alpha + alphabeta[sp[i], 1]) + 
          (mu_beta + alphabeta[sp[i], 2]) * year[i] + 
          gamma * lat[i] + eta * (year[i] * lat[i]) + 
          nu * sex[i];
}
}

model {
sigma_raw ~ normal(0, 1);
gamma_raw ~ normal(0, 1);
eta_raw ~ normal(0, 1);
nu_raw ~ normal(0, 1);

to_vector(z) ~ normal(0, 1);
mu_alpha ~ normal(0, 10);
mu_beta ~ normal(0, 10);
L_Rho ~ lkj_corr_cholesky(2);
sigma_sp ~ normal(0, 5);

y ~ normal(mu, sigma);
}

generated quantities {
real y_rep[N];

y_rep = normal_rng(mu, sigma);
}
"

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.90
TREE_DEPTH <- 16
STEP_SIZE <- 0.005
CHAINS <- 3
ITER <- 1000

tt <- proc.time()
fit3 <- rstan::stan(model_code = stanmodel3,
                    data = DATA,
                    chains = CHAINS,
                    iter = ITER,
                    cores = CHAINS,
                    pars = c('alpha',
                             'beta',
                             'mu_alpha',
                             'mu_beta',
                             'Rho',
                             'sigma_sp',
                             'sigma',
                             'gamma',
                             'eta',
                             'nu'), 
                    control = list(adapt_delta = DELTA,
                                   max_treedepth = TREE_DEPTH,
                                   stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
saveRDS(fit3, file = 'MAPS-wc-time-lat-chol-stan_output.rds')

MCMCvis::MCMCsummary(fit3, n.eff = TRUE, round = 2)
MCMCvis::MCMCplot(fit3, excl = c('eta', 'lp__'))

# library(shinystan)
# launch_shinystan(fit3)

# y_rep <- MCMCvis::MCMCchains(fit2, params = 'y_rep')
# bayesplot::ppc_dens_overlay(DATA$y, y_rep[1:25,])
# plot(DATA$y, y_rep[4,], pch = '.', xlim = c(0, 200), ylim = c(0, 200))
# abline(a = 0, b = 1, col = 'red', lty = 2)

