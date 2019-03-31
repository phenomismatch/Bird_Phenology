####################
# Changes in wing chord over time/lat
#
# Only age 2+ males
# 
# NEXT: fit female only model
####################


# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/'


# load packages -----------------------------------------------------------

library(dplyr)
library(rstan)


# load data ---------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))

maps_data <- readRDS('MAPS-age-filled.rds')



# wing chord --------------------------------------------------------------

#QC data - remove zeros
to.rm <- which(maps_data$weight == 0 | maps_data$wing_chord == 0 | 
                 is.na(maps_data$weight) | is.na(maps_data$weight) | 
                 is.na(maps_data$fat_content) | is.na(maps_data$wing_chord))
maps_data_qc <- maps_data[-to.rm, ]

#just age 2+ males caught on or before day 200
maps_ad <- dplyr::filter(maps_data_qc, sex == 'M', age >= 5, day <=200)

#remove WC values outside 3 sd of mean
usp_maps_ad <- unique(maps_ad$sci_name)
maps_ad_qc <- data.frame()
for (i in 1:length(usp_maps_ad))
{
  #i <- 7
  temp <- dplyr::filter(maps_ad, sci_name == usp_maps_ad[i])
  sd_wc <- sd(temp$wing_chord)
  mn_wc <- mean(temp$wing_chord)
  #outside 3 sds
  low <- mn_wc - 3*sd_wc
  high <-  mn_wc + 3*sd_wc
  to.rm <- which(temp$wing_chord < low | temp$wing_chord > high)
  if (length(to.rm) > 0)
  {
    tt <- temp[-to.rm,]
  } else {
    tt <- temp
  }

  maps_ad_qc <- rbind(maps_ad_qc, tt)
}


#only species with > 500 data points
d_cnt <- plyr::count(maps_ad_qc, 'sci_name')
sp_p <- dplyr::filter(d_cnt, freq > 500)[,1]

#subset of species
set.seed(1)
sp <- base::sample(sp_p, size = 50)
#sp <- sp_p

stan_data <- dplyr::filter(maps_ad_qc, sci_name %in% sp)

#factor for species_id and year
stan_data$sp_f <- factor(stan_data$sci_name)
usp <- sort(unique(stan_data$sci_name))
stan_data$year_f <- as.numeric(factor(stan_data$year))

DATA <- list(N = NROW(stan_data),
             Nsp = length(usp),
             y = stan_data$wing_chord,
             sp = as.numeric(stan_data$sp_f), #species
             year = stan_data$year_f,
             lat = stan_data$lat)


stanmodel3 <- "
data {
int<lower=0> N;                     // number of obs
int<lower=0> Nsp;                   // number of species
real<lower=0> y[N];                 // response
int<lower=1, upper=Nsp> sp[N];      // species ids
vector<lower=0>[N] year;
vector<lower=0>[N] lat;
}

parameters {
real mu_alpha_raw;
real mu_beta_raw;
real mu_gamma_raw;
real<lower = 0> sigma_raw;
vector<lower = 0>[3] sigma_sp_raw;                // standard deviations
cholesky_factor_corr[3] L_Rho;                    // correlation matrix
matrix[3, Nsp] z;
// real eta_raw;
real nu_raw;
}

transformed parameters {
vector[N] mu;
matrix[Nsp, 3] abg;                               // matrix for alpha, beta, and gamma
matrix[3, 3] Rho;                                 // covariance matrix
real<lower = 0> sigma;
// real eta;
vector[Nsp] alpha;
vector[Nsp] beta;
vector[Nsp] gamma;
vector<lower = 0>[3] sigma_sp;
real mu_alpha;
real mu_beta;
real mu_gamma;

sigma = sigma_raw * 5;
sigma_sp[1] = sigma_sp_raw[1] * 20;             // variance alpha
sigma_sp[2] = sigma_sp_raw[2] * 3;              // variance beta
sigma_sp[3] = sigma_sp_raw[3] * 3;              // variance gamma
// eta = eta_raw * 3;
mu_alpha = mu_alpha_raw * 10 + 50;
mu_beta = mu_beta_raw * 3;
mu_gamma = mu_gamma_raw * 3;


// cholesky factor of covariance matrix multiplied by z score
abg = (diag_pre_multiply(sigma_sp, L_Rho) * z)';
alpha = abg[,1];
beta = abg[,2];
gamma = abg[,3];
Rho = L_Rho * L_Rho';

for (i in 1:N)
{
  mu[i] = (mu_alpha + abg[sp[i], 1]) + 
  (mu_beta + abg[sp[i], 2]) * year[i] + 
  (mu_gamma + abg[sp[i], 3]) * lat[i];
}
}

model {
sigma_raw ~ normal(0, 1);         // sigma ~ halfnormal(0, 5)
// eta_raw ~ normal(0, 1);           // eta ~ normal(0, 3)

to_vector(z) ~ normal(0, 1);
mu_alpha_raw ~ normal(0, 1);      // mu_alpha ~ normal(50, 10)
mu_beta_raw ~ normal(0, 1);       // mu_beta ~ normal(0, 3)
mu_gamma_raw ~ normal(0, 1);      // mu_gamma ~ normal(0, 3)
L_Rho ~ lkj_corr_cholesky(1);
sigma_sp_raw ~ normal(0, 1);      // sigma_sp[1] ~ halfnormal(0, 20); sigma_sp[2] ~ halfnormal(0, 3); sigma_sp[3] ~ halfnormal(0, 3); 

y ~ normal(mu, sigma);
}

generated quantities {
real y_rep[N];

y_rep = normal_rng(mu, sigma);
}
"

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.85
TREE_DEPTH <- 16
STEP_SIZE <- 0.005
CHAINS <- 4
ITER <- 1000

tt <- proc.time()
fit3 <- rstan::stan(model_code = stanmodel3,
                    data = DATA,
                    chains = CHAINS,
                    iter = ITER,
                    cores = CHAINS,
                    pars = c('alpha',
                             'beta',
                             'gamma',
                             'mu_alpha',
                             'mu_beta',
                             'mu_gamma',
                             'Rho',
                             'L_Rho',
                             'sigma_sp',
                             'sigma',
                             #'eta',
                             'z',
                             'y_rep'), 
                    control = list(adapt_delta = DELTA,
                                   max_treedepth = TREE_DEPTH,
                                   stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
saveRDS(fit3, file = 'MAPS-wc-time-lat-chol-stan_output-vary-gamma-50.rds')
#fit3 <- readRDS('MAPS-wc-time-lat-chol-stan_output-vary-gamma.rds')

MCMCvis::MCMCsummary(fit3, n.eff = TRUE, round = 2, excl = 'y_rep')
MCMCvis::MCMCsummary(fit3, n.eff = TRUE, round = 2, params = 'mu', ISB = FALSE)
#MCMCvis::MCMCplot(fit3, excl = c('eta', 'lp__', 'y_rep'))
MCMCvis::MCMCplot(fit3, params = c('gamma'))
MCMCvis::MCMCplot(fit3, params = c('beta'))
MCMCvis::MCMCplot(fit3, params = c('z'))
MCMCvis::MCMCplot(fit3, params = c('L_Rho'))

# library(shinystan)
# launch_shinystan(fit3)

# y_rep <- MCMCvis::MCMCchains(fit3, params = 'y_rep')
# bayesplot::ppc_dens_overlay(DATA$y, y_rep[1:25,])
# plot(DATA$y, y_rep[1,], pch = '.', xlim = c(0, 200), ylim = c(0, 200))
# abline(a = 0, b = 1, col = 'red', lty = 2)



# PPO ---------------------------------------------------------------------

# 
# 
# sigma_raw ~ normal(0, 1);         // sigma ~ halfnormal(0, 5)
# eta_raw ~ normal(0, 1);           // eta ~ normal(0, 3)
# 
# to_vector(z) ~ normal(0, 1);
# mu_alpha_raw ~ normal(0, 1);      // mu_alpha ~ normal(30, 10)
# mu_beta_raw ~ normal(0, 1);       // mu_beta ~ normal(0, 3)
# mu_gamma_raw ~ normal(0, 1);      // mu_gamma ~ normal(0, 3)
# L_Rho ~ lkj_corr_cholesky(3);
# sigma_sp_raw ~ normal(0, 1);      // sigma_sp[1] ~ halfnormal(0, 10); sigma_sp[2] ~ halfnormal(0, 3); sigma_sp[3] ~ halfnormal(0, 3); 



#mu_alpha ~ N(50, 10)
PR <- rnorm(10000, 50, 10)
MCMCvis::MCMCtrace(fit3,
                   params = 'mu_alpha',
                   priors = PR,
                   pdf = FALSE)

#mu_beta ~ N(0, 3)
PR <- rnorm(10000, 0, 1)
MCMCvis::MCMCtrace(fit3,
                   params = 'mu_beta',
                   priors = PR,
                   pdf = FALSE)

#mu_gamma ~ N(0, 3)
PR <- rnorm(10000, 0, 3)
MCMCvis::MCMCtrace(fit3,
                   params = 'mu_gamma',
                   priors = PR,
                   pdf = FALSE)

#sigma_sp[1] ~ HN(0, 10)
PR_p <- rnorm(10000, 0, 10)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit3,
                   params = 'sigma_sp\\[1',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = FALSE)

#sigma_sp[2] ~ HN(0, 3)
PR_p <- rnorm(10000, 0, 3)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit3,
                   params = 'sigma_sp\\[2',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = FALSE)

#sigma_sp[3] ~ HN(0, 3)
PR_p <- rnorm(10000, 0, 3)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit3,
                   params = 'sigma_sp\\[3',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = FALSE)

#eta ~ N(0, 5)
PR <- rnorm(10000, 0, 3)
MCMCvis::MCMCtrace(fit3,
                   params = 'eta',
                   priors = PR,
                   pdf = FALSE)

#sigma ~ HN(0, 5)
PR_p <- rnorm(10000, 0, 5)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit3,
                   params = 'sigma',
                   priors = PR,
                   pdf = FALSE, 
                   post_zm = FALSE)


head(stan_data)
head(usp)

plyr::count(stan_data, 'sci_name')
