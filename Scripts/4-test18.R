######################
# 4 - IAR model - NO LUMPED SPATIAL/NON-SPATIAL NO HIERARCHICAL SIGMA_PHI
#
# Fit IAR model
#
# VVV MODEL VVV
# i = year
# j = cell
# y_{obs[i,j]} \sim N(y_{true[i,j]}, \sigma_{y[i,j]})
# y_{true[i,j]} \sim N(\mu_{[i,j]}, \sigma_{y_true})
# \mu_{[i,j]} = \beta_{0[i]} + \gamma_{[j]} + \phi_{[i,j]} * sigma_phi
# \beta_{0[i]} \sim N(0, \sigma_{\beta_{0}})
# \gamma_{[j]} \sim N(\mu_{\gamma[j]}, \sigma_{\gamma})
# \mu_{\gamma[j]} = \alpha_{\gamma} + \beta_{\gamma} * lat_{[j]}
# \sigma_{\phi} \sim HN(0, 5)
# \phi_{[i,j]} \sim N(0, [D - W]^{-1})
# \sum_{j}{} \phi_{[i,j]} = 0
######################

#Stan resources:
#http://mc-stan.org/users/documentation/case-studies/icar_stan.html
#http://mc-stan.org/users/documentation/case-studies/divergences_and_bias.html
#https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html
#https://chi-feng.github.io/mcmc-demo/app.html#HamiltonianMC,standard
#https://groups.google.com/forum/#!msg/stan-users/zOjAeJC4x_E/OyCOfJo8AwAJ (regarding non-centered parameterization on sd)


# Top-level dir -----------------------------------------------------------

#desktop/laptop
#dir <- '~/Google_Drive/R/'

#Xanadu
dir <- '/labs/Tingley/phenomismatch/'



# db/hm query dir ------------------------------------------------------------

IAR_in_dir <- 'IAR_input_2019-05-03'
IAR_out_dir <- '4-test'



# Load packages -----------------------------------------------------------

library(rstan)
library(geosphere)
library(ggplot2)
library(maps)
library(dplyr)
library(dggridR)
library(MCMCvis)
#Also need to be installed, but not loaded: rgeos, maptools, mapproj



# Set wd ------------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_in_dir))

IAR_in_date <- substr(IAR_in_dir, start = 11, stop = 20)
IAR_out_date <- '2020-03-04'



# species arg -----------------------------------------------------

#args <- commandArgs(trailingOnly = TRUE)
#args <- as.character('Icterus_spurius')
#args <- as.character('Catharus_minimus')
#args <- as.character('Empidonax_virescens')
#args <- as.character('Vireo_olivaceus')
args <- as.character('Agelaius_phoeniceus')


# Filter data by species/years ------------------------------------------------------

#read in master df
df_master <- readRDS(paste0('IAR_input-', IAR_in_date, '.rds'))

#filter by species and year to be modeled
f_out <- dplyr::filter(df_master, species == args & MODEL == TRUE)

#define cells and years to be modeled
cells <- unique(f_out$cell)
years <- unique(f_out$year)
nyr <- length(years)
ncell <- length(cells)



# create adjacency matrix -------------------------------------------------

#make hexgrid
hexgrid6 <- dggridR::dgconstruct(res = 6)

#get hexgrid cell centers
cellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid6, cells)

#add lat col to df
f_out$lat <- cellcenters$lat_deg

#create adjacency matrix - 1 if adjacent to cell, 0 if not
adjacency_matrix <- matrix(data = NA, nrow = ncell, ncol = ncell)

for (i in 1:ncell)
{
  #i <- 1
  for (j in i:ncell)
  {
    #j <- 69
    dists <- geosphere::distm(c(cellcenters$lon_deg[i], cellcenters$lat_deg[i]),
                              c(cellcenters$lon_deg[j], cellcenters$lat_deg[j]))
    adjacency_matrix[i,j] <- as.numeric((dists/1000) > 0 & (dists/1000) < 311)
  }
}

#indices for 1s
ninds <- which(adjacency_matrix == 1, arr.ind = TRUE)


#if a cell doesn't border any other cells, drop it and redefine objects
DROP <- FALSE

s_cols <- apply(adjacency_matrix, 2, function(x) sum(x, na.rm = TRUE))
s_rows <- apply(adjacency_matrix, 1, function(x) sum(x, na.rm = TRUE))
to.rm.ind <- which((s_cols + s_rows) == 0)

if (length(to.rm.ind) > 0)
{
  DROP <- cells[to.rm.ind]
  
  cells <- cells[-to.rm.ind]
  ncell <- length(cells)
  f_out <- dplyr::filter(f_out, cell %in% cells)
  cellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid6, cells)
  
  #create adjacency matrix - 1 if adjacent to cell, 0 if not
  adjacency_matrix <- matrix(data = NA, nrow = ncell, ncol = ncell)
  
  for (i in 1:ncell)
  {
    #i <- 1
    for (j in i:ncell)
    {
      #j <- 4
      dists <- geosphere::distm(c(cellcenters$lon_deg[i], cellcenters$lat_deg[i]),
                                c(cellcenters$lon_deg[j], cellcenters$lat_deg[j]))
      adjacency_matrix[i,j] <- as.numeric((dists/1000) > 0 & (dists/1000) < 311)
    }
  }
  
  #indices for 1s
  ninds <- which(adjacency_matrix == 1, arr.ind = TRUE)
}



# create Stan data object -------------------------------------------------

#create and fill sds, obs, data for PPC, and temp
sigma_y_in <- matrix(nrow = nyr, ncol = ncell)
y_obs_in <- matrix(nrow = nyr, ncol = ncell)
y_PPC <- rep(NA, nyr * ncell)

#number of observation and NAs for each year
len_y_obs_in <- rep(NA, nyr)
len_y_mis_in <- rep(NA, nyr)

#indices for observed and missing
ii_obs_in <- matrix(NA, nrow = nyr, ncol = ncell)
ii_mis_in <- matrix(NA, nrow = nyr, ncol = ncell)

#counter to fill y_PPC
counter <- 1
for (i in 1:nyr)
{
  #j <- 16
  temp_yr <- dplyr::filter(f_out, year == years[i])
  
  #don't need to manipulate position of sigmas
  sigma_y_in[i,] <- temp_yr$HM_sd
  
  for (j in 1:ncell)
  {
    #n <- 1
    #matrix with observed values with NAs
    y_PPC[counter] <- temp_yr$HM_mean[j]
    counter <- counter + 1
  }
  
  #which are not NA
  no_na <- temp_yr$HM_mean[which(!is.na(temp_yr$HM_mean))]
  
  #pad end with NAs
  if (length(no_na) < ncell)
  {
    num_na <- ncell - length(no_na)
    
    #add NAs to end
    t_y_obs_in <- c(no_na, rep(NA, num_na))
    t_obs_in <- c(which(!is.na(temp_yr$HM_mean)), rep(NA, num_na)) 
    t_mis_in <- c(which(is.na(temp_yr$HM_mean)), rep(NA, length(no_na)))
    
    #fill objects
    ii_obs_in[i,] <- t_obs_in
    ii_mis_in[i,] <- t_mis_in
    y_obs_in[i,] <- t_y_obs_in
  } else {
    #no NAs to end (no mimssing values)
    y_obs_in[i,] <- no_na
    ii_mis_in[i,] <- which(!is.na(temp_yr$HM_mean))
    y_obs_in[i,] <- which(is.na(temp_yr$HM_mean))
  }
  
  #length of data/miss for each year
  len_y_obs_in[i] <- length(no_na)
  len_y_mis_in[i] <- ncell - length(no_na)
}


#fill 0 where NA in y_obs - Stan does not like NA and zeros are not being used to estimate any param (y_obs is used to fill y)
y_obs_in[which(is.na(y_obs_in), arr.ind = TRUE)] <- 0
#see script here showing this value fo sigma_y_in has no impact on y_true: Archive/test_sigma_insert.R
sigma_y_in[which(is.na(sigma_y_in), arr.ind = TRUE)] <- 0.1
ii_obs_in[which(is.na(ii_obs_in), arr.ind = TRUE)] <- 0
ii_mis_in[which(is.na(ii_mis_in), arr.ind = TRUE)] <- 0


#create data list for Stan
DATA <- list(J = ncell,
             cells = cells,
             N = nyr, 
             NJ = nyr * ncell,
             N_obs = len_y_obs_in,
             N_mis = len_y_mis_in,
             N_edges = nrow(ninds), 
             node1 = ninds[,1],
             node2 = ninds[,2],
             y_obs = y_obs_in,
             sigma_y = sigma_y_in,
             ii_obs = ii_obs_in,
             ii_mis = ii_mis_in,
             lat = cellcenters$lat_deg,
             y_PPC = y_PPC)



# Stan model --------------------------------------------------------------

#years in rows
#cells in columns

IAR <- '
data {
int<lower = 0> N;                                     // number of years
int<lower = 0> J;                                     // number of cells
int<lower = 0> NJ;                                    // number of cell/years
int<lower = 0> N_obs[N];                              // number of non-missing for each year
int<lower = 0> N_mis[N];                              // number missing for each year
int<lower = 0> N_edges;                               // number of edges in adjacency matrix
int<lower = 1, upper = J> node1[N_edges];             // node1[i] adjacent to node2[i]
int<lower = 1, upper = J> node2[N_edges];             // and node1[i] < node2[i]
vector[J] y_obs[N];                                   // observed response data (add NAs to end)
vector<lower = 0>[J] sigma_y[N];                      // observed sd of data (observation error)
int<lower = 0> ii_obs[N, J];                          // indices of observed data
int<lower = 0> ii_mis[N, J];                          // indices of missing data
vector<lower = 24, upper = 90>[J] lat;
}

parameters {
vector[J] y_mis[N];
real alpha_gamma;
real beta_gamma;                                  // effect of latitude
real<lower = 0> sigma_gamma;
vector[J] gamma_raw;
vector[J] phi[N];                                     // sptial error componenet
real<lower = 0> sigma_phi;
vector[N] beta0_raw;
real<lower = 0> sigma_beta0;
vector[J] y_true_raw[N];                           // J vectors (years in rows) of length N (cells in cols)
real<lower = 0> sigma_y_true;
}

transformed parameters {
vector[J] y[N];
vector[J] gamma;
vector[J] mu_gamma;
vector[J] y_true[N];
vector[N] beta0;

mu_gamma = alpha_gamma + beta_gamma * lat;
// implies gamma ~ N(mu_gamma, sigma_gamma)
gamma = gamma_raw * sigma_gamma + mu_gamma;
beta0 = beta0_raw * sigma_beta0;

// for each year, vectorize over cells
for (i in 1:N)
{
  // implies y_true ~ N(beta0 + gamma + phi * sigma_phi, sigma_y_true)
  y_true[i] = y_true_raw[i] * sigma_y_true + beta0[i] + gamma + phi[i] * sigma_phi;
  
  // indexing to avoid NAs  
  y[i, ii_obs[i, 1:N_obs[i]]] = y_obs[i, 1:N_obs[i]];
  y[i, ii_mis[i, 1:N_mis[i]]] = y_mis[i, 1:N_mis[i]];
}
}

model {
// priors for non-centered parameters
gamma_raw ~ std_normal();
beta0_raw ~ std_normal();

// priors for centered parameters
alpha_gamma ~ normal(0, 100);
// beta_gamma = 3 represents 37 km / day (which matches speeds reported in lit and derived in La Sorte et al. 2013 Ecology)
beta_gamma ~ normal(3, 3);
sigma_gamma ~ normal(0, 10);
sigma_beta0 ~ normal(0, 10);
sigma_y_true ~ normal(0, 10);
sigma_phi ~ normal(0, 10);

for (i in 1:N)
{
  y_true_raw[i] ~ std_normal();
  // index array first (each year), then vector (for cells)
  target += -0.5 * dot_self(phi[i, node1] - phi[i, node2]);
  // soft sum to 0 constraint - J is number of phis per year
  sum(phi[i]) ~ normal(0, 0.001 * J);
  
  // observation model for y
  y[i] ~ normal(y_true[i], sigma_y[i]);
}
}

generated quantities {

vector[NJ] y_rep;
int<lower = 0> counter;

counter = 1;
for (n in 1:N)
{
  for (j in 1:J)
  {
  y_rep[counter] = normal_rng(y_true[n,j], sigma_y[n,j]);
  counter = counter + 1;
  }
}
}'



# Run model ---------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.90
TREE_DEPTH <- 18
STEP_SIZE <- 0.0005
CHAINS <- 6
ITER <- 5000

tt <- proc.time()
fit <- rstan::stan(model_code = IAR,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('sigma_beta0', 
                            'beta0',
                            'alpha_gamma', 
                            'beta_gamma', 
                            'sigma_gamma', 
                            'gamma',
                            'phi',
                            'sigma_phi',
                            'sigma_y_true',
                            'y_true', 
                            'y_rep'),
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60


# Summaries ---------------------------------------------------------------

#get summary of model output
model_summary <- MCMCvis::MCMCsummary(fit, Rhat = TRUE, n.eff = TRUE, 
                                      round = 2, excl = 'y_rep')

#extract Rhat and neff values
rhat_output <- as.vector(model_summary[, grep('Rhat', colnames(model_summary))])
neff_output <- as.vector(model_summary[, grep('n.eff', colnames(model_summary))])



# rerun if necessary ------------------------------------------------------

#double iterations and run again
while ((max(rhat_output) > 1.05 | min(neff_output) < 400) & ITER < 10001)
{
  
  ITER <- ITER * 2
  
  fit <- rstan::stan(model_code = IAR,
                     data = DATA,
                     chains = CHAINS,
                     iter = ITER,
                     cores = CHAINS,
                     pars = c('sigma_beta0', 
                              'beta0',
                              'alpha_gamma', 
                              'beta_gamma', 
                              'sigma_gamma', 
                              'gamma',
                              'phi',
                              'sigma_phi',
                              'sigma_y_true',
                              'y_true', 
                              'y_rep'),
                     control = list(adapt_delta = DELTA,
                                    max_treedepth = TREE_DEPTH,
                                    stepsize = STEP_SIZE))
  run_time <- (proc.time()[3] - tt[3]) / 60
  
  #get updated summary
  model_summary <- MCMCvis::MCMCsummary(fit, Rhat = TRUE, n.eff = TRUE, 
                                        round = 2, excl = 'y_rep')
  
  #extract Rhat and neff values
  rhat_output <- as.vector(model_summary[, grep('Rhat', colnames(model_summary))])
  neff_output <- as.vector(model_summary[, grep('n.eff', colnames(model_summary))])
}


#save to RDS
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_out_dir))
saveRDS(fit, file = paste0(args, '-iar-stan_output-18-', IAR_out_date, '.rds'))

#save data to RDS (has which cells are modeled)
saveRDS(DATA, file = paste0(args, '-iar-stan_input-18-',  IAR_out_date, '.rds'))


# Calc diagnostics ---------------------------------------------------

# library(shinystan)
# launch_shinystan(fit)

num_diverge <- rstan::get_num_divergent(fit)
num_tree <- rstan::get_num_max_treedepth(fit)
num_BFMI <- length(rstan::get_low_bfmi_chains(fit))

sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
mn_stepsize <- sapply(sampler_params,
                      function(x) mean(x[, 'stepsize__']))
mn_treedepth <- sapply(sampler_params,
                       function(x) mean(x[, 'treedepth__']))
accept_stat <- sapply(sampler_params,
                      function(x) mean(x[, 'accept_stat__']))



# write model results to file ---------------------------------------------

options(max.print = 5e6)
sink(paste0(args, '-iar-stan_results-18-', IAR_out_date, '.txt'))
cat(paste0('IAR results ', args, ' \n'))
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
cat(paste0('Cell drop: ', DROP, ' \n'))
cat(paste0('Max Rhat: ', max(rhat_output), ' \n'))
cat(paste0('Min n.eff: ', min(neff_output), ' \n'))
print(model_summary)
sink()


