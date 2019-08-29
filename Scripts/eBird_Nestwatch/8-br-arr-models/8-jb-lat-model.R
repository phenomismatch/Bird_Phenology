######################
# 8 - juv - brood patch ~ cell lat
#
# not much there, it seems
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/UCHC/LABS/Tingley/phenomismatch/'



# model dir ------------------------------------------------------------

#juveniles MAPS - date input data processed
bp_date <- '2019-08-26'

#IAR data
arr_master_dir <- 'arrival_master_2019-05-26'



# Load packages -----------------------------------------------------------

library(rstan)
library(ggplot2)
library(dplyr)
library(MCMCvis)
library(dggridR)



# Filter data ------------------------------------------------------------------

#juveniles hitting nets - MAPS

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/br_arr_', bp_date))

#read in 
bp_master <- readRDS(paste0('bp-output-', bp_date, '.rds'))

juv_master <- readRDS(paste0('juv-output-', bp_date, '.rds'))

#master arrival data (from IAR output) - for cell lat
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', arr_master_dir))

arr_master <- readRDS(paste0(arr_master_dir, '.rds'))

#merge data sets
mrg <- dplyr::inner_join(bp_master, juv_master, by = c('species', 'cell', 'year'))

#create hex grid
hexgrid6 <- dggridR::dgconstruct(res = 6)

#get lat/lon
ll <- dggridR::dgSEQNUM_to_GEO(hexgrid6, in_seqnum = mrg$cell)

#add to df
mrg$cell_lat <- ll$lat_deg

#only species/cells/years with data for juvs
mrg_f <- dplyr::filter(mrg, !is.na(bp_mean), !is.na(juv_mean))

#only species that have at least 5 data points
cnt_arr <- plyr::count(mrg_f, 'species')
sp_f <- dplyr::filter(cnt_arr, freq >= 5)$species

#only species that cover more than 2 degrees lat
mrg_f2 <- dplyr::filter(mrg_f, species %in% sp_f)
lat_rng <- aggregate(cell_lat ~ species, data = mrg_f2, FUN = function(x) diff(range(x)))
sp_f2 <- lat_rng$species[which(lat_rng$cell_lat > 2)]

#final filter
mrg_f3 <- dplyr::filter(mrg_f2, species %in% sp_f2)

sp_idx <- as.numeric(factor(mrg_f3$species))

#difference between arrival and juveniles hitting nets
#for difference of normals, subtract means, add variances - http://mathworld.wolfram.com/NormalDifferenceDistribution.html
mrg_f3$jb <- mrg_f3$juv_mean - mrg_f3$bp_mean
mrg_f3$sd_jb <- sqrt(mrg_f3$juv_sd^2 + mrg_f3$bp_sd^2)



# Stan model --------------------------------------------------------------

DATA <- list(y = mrg_f3$jb,
             sd_y = mrg_f3$sd_jb,
             N = NROW(mrg_f3),
             sp = sp_idx,
             cell_lat = mrg_f3$cell_lat,
             Nsp = length(unique(sp_idx)))

stanmodel1 <- "
data {
int<lower=0> N;                      // number of data points
vector[N] y;
vector<lower=0>[N] sd_y;
int<lower=0> sp[N];              
int<lower=0> Nsp;
vector<lower=0>[N] cell_lat;
}

parameters {
real<lower = 0> sigma_raw;
real mu_y_raw[N];
real mu_alpha_raw;
real mu_beta_raw;
vector<lower = 0>[2] sigma_sp_raw;
cholesky_factor_corr[2] L_Rho;             // correlation matrix
matrix[2, Nsp] z;                          // z-scores
}

transformed parameters {
real<lower = 0> sigma;
vector[N] mu_y;
real mu[N];
matrix[Nsp, 2] ab;                              // matrix for alpha, beta, gamma, and theta
matrix[2, 2] Rho;                                 // covariance matrix
vector[Nsp] alpha;
vector[Nsp] beta;
vector[Nsp] alpha_c;
vector[Nsp] beta_c;
vector<lower = 0>[2] sigma_sp;
real mu_alpha;
real mu_beta;

mu_alpha = mu_alpha_raw * 300;
mu_beta = mu_beta_raw * 2;
sigma = sigma_raw * 20;
sigma_sp[1] = sigma_sp_raw[1] * 100;
sigma_sp[2] = sigma_sp_raw[2] * 10;


// cholesky factor of covariance matrix multiplied by z score
ab = (diag_pre_multiply(sigma_sp, L_Rho) * z)';
alpha = ab[,1];
beta = ab[,2];
Rho = L_Rho * L_Rho';

for (i in 1:N)
{
  mu[i] = (mu_alpha + ab[sp[i], 1]) +
  (mu_beta + ab[sp[i], 2]) * cell_lat[i];
  
  // implies mu_y ~ normal(mu, sigma)
  mu_y[i] = mu_y_raw[i] * sigma + mu[i];
}

alpha_c = mu_alpha + alpha;
beta_c = mu_beta + beta;

}

model {
sigma_raw ~ std_normal();
sigma_sp_raw ~ std_normal();
mu_alpha_raw ~ std_normal();
mu_beta_raw ~ std_normal();
mu_y_raw ~ std_normal();

to_vector(z) ~ std_normal();
L_Rho ~ lkj_corr_cholesky(1);

// observation model
y ~ normal(mu_y, sd_y);
}

generated quantities {
real y_rep[N];

y_rep = normal_rng(mu_y, sd_y);
}
"

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.99
TREE_DEPTH <- 1
STEP_SIZE <- 0.0001
CHAINS <- 4
ITER <- 4000

tt <- proc.time()
fit <- rstan::stan(model_code = stanmodel1,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha_c',
                            'beta_c',
                            'alpha',
                            'beta',
                            'mu_alpha',
                            'mu_beta',
                            'sigma_sp',
                            'Rho',
                            'sigma',
                            'mu_y',
                            'y_rep'), 
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60


setwd(paste0(dir, 'Bird_Phenology/Data/Processed/br_arr_', bp_date))
saveRDS(fit, file = paste0('jb-lat-stan-output-', bp_date, '.rds'))
#fit <- readRDS(paste0('jb-lat-stan-output-', bp_date, '.rds'))


num_diverge <- rstan::get_num_divergent(fit)
num_tree <- rstan::get_num_max_treedepth(fit)
num_BFMI <- length(rstan::get_low_bfmi_chains(fit))



# Calc diagnostics ---------------------------------------------------

# library(shinystan)  
# launch_shinystan(fit)

sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
mn_stepsize <- sapply(sampler_params, 
                      function(x) mean(x[, 'stepsize__']))
mn_treedepth <- sapply(sampler_params, 
                       function(x) mean(x[, 'treedepth__']))
accept_stat <- sapply(sampler_params, 
                      function(x) mean(x[, 'accept_stat__']))



# Summaries ---------------------------------------------------------------

#get summary of model output
model_summary <- MCMCvis::MCMCsummary(fit, Rhat = TRUE, n.eff = TRUE, round = 2, excl = c('y_rep', 'mu_y'))

#extract Rhat and neff values
rhat_output <- as.vector(model_summary[, grep('Rhat', colnames(model_summary))])
neff_output <- as.vector(model_summary[, grep('n.eff', colnames(model_summary))])



# write model results to file ---------------------------------------------

options(max.print = 50000)
sink(paste0('jb-lat-stan-results-', DATE, '.txt'))
cat(paste0('Total minutes: ', round(run_time, digits = 2), ' \n'))
cat(paste0('Adapt delta: ', DELTA, ' \n'))
cat(paste0('Max tree depth: ', TREE_DEPTH, ' \n'))
cat(paste0('Step size: ', STEP_SIZE, ' \n'))
cat(paste0('Number of divergences: ', num_diverge, ' \n'))
cat(paste0('Number of tree exceeds: ', num_tree, ' \n'))
cat(paste0('Number chains low BFMI: ', num_BFMI, ' \n'))
cat(paste0('Mean stepsize: ', round(mean(mn_stepsize), 5), ' \n'))
cat(paste0('Mean treedepth: ', round(mean(mn_treedepth), 1), ' \n'))
cat(paste0('Mean accept stat: ', round(mean(accept_stat), 2), ' \n'))
cat(paste0('Max Rhat: ', max(rhat_output, na.rm = TRUE), ' \n'))
cat(paste0('Min n.eff: ', min(neff_output, na.rm = TRUE), ' \n'))
print(model_summary)
sink()




# check output ------------------------------------------------------------


# MCMCvis::MCMCsummary(fit, n.eff = TRUE, round = 2, params = 'alpha')
# MCMCvis::MCMCplot(fit, params = c('alpha'), labels = usp, main = 'alpha', rank = TRUE)
# 
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, round = 2, params = 'beta')
# MCMCvis::MCMCplot(fit, params = c('beta_c'), labels = usp,
#                   main = 'beta - change over time')
#  
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, round = 2, params = 'gamma')
# MCMCvis::MCMCplot(fit, params = c('gamma_c'),
#                   labels = usp, main = 'gamma - change over latitude')
# 
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, round = 2, params = 'theta')
# MCMCvis::MCMCplot(fit, params = c('theta_c'),
#                   labels = usp, main = 'theta - change over elevation')
# 
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, round = 2, params = 'sigma')
# 
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, round = 5, params = c('mu_beta', 'mu_gamma', 'mu_theta'), ISB = FALSE)
# MCMCvis::MCMCplot(fit, params = 'mu', ISB = FALSE)

# sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
# sapply(sampler_params, function(x) mean(x[, 'stepsize__']))
# sapply(sampler_params, function(x) mean(x[, 'divergent__']))
# sapply(sampler_params, function(x) mean(x[, 'treedepth__']))
# sapply(sampler_params, function(x) mean(x[, 'accept_stat__']))
# sapply(sampler_params, function(x) mean(x[, 'n_leapfrog__']))
# sapply(sampler_params, function(x) mean(x[, 'energy__']))



# # PPC ---------------------------------------------------------------------


y_val <- DATA$y
y_rep <- MCMCvis::MCMCchains(fit, params = 'y_rep')
bayesplot::ppc_stat(DATA$y, y_rep, stat = 'mean')
bayesplot::ppc_dens_overlay(DATA$y, y_rep[1:100,])
PPC_fun <- function(FUN, YR = y_rep, D = DATA$y)
{
  out <- sum(apply(YR, 1, FUN) > FUN(D)) / NROW(YR)
  print(out)
}
PPC_fun(mean)
PPC_fun(min)
PPC_fun(max)

# 
# 
# # PPO ---------------------------------------------------------------------

# mu_alpha = mu_alpha_raw * 300;
# mu_beta = mu_besta_raw * 2;
# sigma = sigma_raw * 20;
# sigma_sp[1] = sigma_sp_raw[1] * 100;
# sigma_sp[2] = sigma_sp_raw[2] * 10;


#mu_alpha ~ N(0, 300)
PR <- rnorm(10000, 0, 300)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_alpha',
                   priors = PR,
                   ind = TRUE,
                   pdf = TRUE,
                   filename = paste0('jb-lat-', bp_date, '-trace_mu_alpha.pdf'))

#mu_beta ~ N(0, 2)
PR <- rnorm(10000, 0, 2)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_beta',
                   priors = PR,
                   ind = TRUE,
                   pdf = TRUE,
                   filename = paste0('jb-lat-', bp_date, '-trace_mu_beta.pdf'))

#sigma_sp[1] ~ HN(0, 100)
PR_p <- rnorm(10000, 0, 100)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_sp\\[1',
                   ISB = 'FALSE',
                   priors = PR,
                   ind = TRUE,
                   pdf = TRUE,
                   filename = paste0('jb-lat-', bp_date, '-trace_sigma_sp[1].pdf'))

#sigma_sp[2] ~ HN(0, 10)
PR_p <- rnorm(10000, 0, 10)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_sp\\[2',
                   ISB = 'FALSE',
                   priors = PR,
                   ind = TRUE,
                   pdf = TRUE,
                   filename = paste0('jb-lat-', bp_date, '-trace_sigma_sp[2].pdf'))


#sigma ~ HN(0, 20)
PR_p <- rnorm(10000, 0, 20)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma',
                   priors = PR,
                   ind = TRUE,
                   pdf = TRUE,
                   filename = paste0('jb-lat-', bp_date, '-trace_sigma.pdf'))



# Plot fit ----------------------------------------------------------------


data_vis_fun <- function(SPECIES = 'all')
{
  #SPECIES <- 'Geothlypis_trichas'
  
  #extract posterior estimates for true states for y and x
  y_true_mean <- MCMCvis::MCMCpstr(fit, params = 'mu_y', type = 'summary', 
                                   func = mean)[[1]]
  y_true_LCI <- MCMCvis::MCMCpstr(fit, params = 'mu_y', type = 'summary', 
                                  func = function(x) quantile(x, probs = c(0.025)))[[1]]
  y_true_UCI <- MCMCvis::MCMCpstr(fit, params = 'mu_y', type = 'summary', 
                                  func = function(x) quantile(x, probs = c(0.975)))[[1]]
  
  
  DATA_PLOT <- data.frame(mean_y = y_true_mean,
                          mean_y_l = y_true_LCI,
                          mean_y_u = y_true_UCI,
                          y_obs = DATA$y,
                          y_obs_l = DATA$y - (1.96 * DATA$sd_y),
                          y_obs_u = DATA$y + (1.96 * DATA$sd_y),
                          x_obs = DATA$cell_lat,
                          sp_id = factor(mrg_f3$species))
  
  if (SPECIES == 'all')
  {
    #model fit for mu_beta and mu_alpha
    a_ch <- MCMCvis::MCMCchains(fit, params = 'mu_alpha')[,1]
    b_ch <- MCMCvis::MCMCchains(fit, params = 'mu_beta')[,1]
    
    DATA_PLOT2 <- DATA_PLOT
  } else {
    
    idx <- which(unique(mrg_f3$species) == SPECIES)
    if (length(idx) > 0)
    {
      a_ch <- MCMCvis::MCMCchains(fit, ISB = FALSE, params = paste0('alpha_c\\[', idx, '\\]'))[,1]
      b_ch <- MCMCvis::MCMCchains(fit, ISB = FALSE, params = paste0('beta_c\\[', idx, '\\]'))[,1]
      
      DATA_PLOT2 <- dplyr::filter(DATA_PLOT, sp_id == SPECIES)
    } else {
      stop(paste0('Species: ', SPECIES, ' not found!'))
    }
  }
  
  sim_x <- seq(min(DATA_PLOT2$x_obs) - 1, max(DATA_PLOT2$x_obs) + 1, length = 100)
  
  mf <- matrix(nrow = length(a_ch), ncol = 100)
  for (i in 1:length(sim_x))
  {
    mf[,i] <- a_ch + b_ch * sim_x[i]
  }
  
  med_mf <- apply(mf, 2, median)
  LCI_mf <- apply(mf, 2, function(x) quantile(x, probs = 0.025))
  UCI_mf <- apply(mf, 2, function(x) quantile(x, probs = 0.975))
  
  FIT_PLOT <- data.frame(MN = med_mf,
                         MN_X = sim_x,
                         LCI = LCI_mf,
                         UCI = UCI_mf)
  
  p <- ggplot(data = DATA_PLOT2, aes(x_obs, mean_y)) +
    #model fit
    geom_ribbon(data = FIT_PLOT,
                aes(x = MN_X, ymin = LCI, ymax = UCI),
                fill = 'grey', alpha = 0.7,
                inherit.aes = FALSE) +
    geom_line(data = FIT_PLOT, aes(MN_X, MN), color = 'red',
              alpha = 0.9,
              inherit.aes = FALSE,
              size = 1.4) +
    #latent state
    geom_errorbar(data = DATA_PLOT2, 
                  aes(ymin = mean_y_l, ymax = mean_y_u), width = 0.3,
                  color = 'black', alpha = 0.2) +
    geom_point(data = DATA_PLOT2, aes(x_obs, mean_y), color = 'black',
               inherit.aes = FALSE, size = 3, alpha = 0.3) +
    # geom_point(data = DATA_PLOT2, aes(mean_x, mean_y, color = sp_id),
    #            inherit.aes = FALSE, size = 3, alpha = 0.3) +
    #observed data
    # geom_errorbar(data = DATA_PLOT2,
    #               aes(x = x_obs, ymin = y_obs_l, ymax = y_obs_u), width = 0.3,
    #               color = 'red', alpha = 0.2) +
    # geom_point(data = DATA_PLOT2, aes(x_obs, y_obs), color = 'red',
    #            inherit.aes = FALSE, size = 3, alpha = 0.3) +
    theme_bw() +
    #scale_x_discrete(limits = c(seq(18,30, by = 2))) +
    ylab('BP - Arr') +
    xlab('Latitude') +
    ggtitle(paste0('Species: ', SPECIES)) +
    theme(
      plot.title = element_text(size = 22),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
      axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
      axis.ticks.length= unit(0.2, 'cm')) #length of axis tick
  
  print(p)
}jb

#all species together
pdf('jb-lat-ALL.pdf')
data_vis_fun(SPECIES = 'all')
dev.off()

#each species individually
sps <- unique(mrg_f3$species)
for (i in 1:length(sps))
{
  #i <- 3
  pdf(paste0('jb-lat-', sps[i], '.pdf'))
  data_vis_fun(SPECIES = sps[i])
  dev.off()
}

