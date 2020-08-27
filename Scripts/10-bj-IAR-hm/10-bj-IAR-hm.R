######################
# 10 - bj IAR model - halfmax
# 
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

breed_date <- '2020-06-04'
juv_date <- '2020-06-04'
run_date <- '2020-08-27'
bj_IAR_out_dir <- paste0('bj_IAR_', run_date)


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

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/breeding_master_', breed_date))
br_master <- readRDS(paste0('breeding_master_', breed_date, '.rds'))

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/juv_master_', juv_date))
juv_master <- readRDS(paste0('juv_master_', juv_date, '.rds'))


# species arg -----------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
#args <- c('Dumetella_carolinensis', 5000)
#args <- c('Geothlypis_trichas', 5000)


# Filter data by species/years ------------------------------------------------------

br_f <- dplyr::select(br_master, species, year, cell, br_GAM_mean, 
                      br_GAM_sd, VALID)

juv_f <- dplyr::select(juv_master, species, year, cell, juv_GAM_mean, 
                      juv_GAM_sd, breed_cell, other_cell, VALID, 
                      per_ovr, cell_lat, cell_lng)

#join
mrg1 <- dplyr::full_join(br_f, juv_f, by = c('species', 'year', 'cell'))

br_na <- which(mrg1$VALID.x == FALSE)
mrg1$br_GAM_mean[br_na] <- NA
mrg1$br_GAM_sd[br_na] <- NA
juv_na <- which(mrg1$VALID.y == FALSE)
mrg1$juv_GAM_mean[juv_na] <- NA
mrg1$juv_GAM_sd[juv_na] <- NA


#filter by year and species
mrg2 <- dplyr::filter(mrg1, year >= 2002, year <= 2017, per_ovr >= 0.05, 
                      species == args[1], breed_cell == TRUE, other_cell == FALSE)
if (NROW(mrg2) < 3)
{
  stop('Not enough data')
}

agg_br <- aggregate(br_GAM_mean ~ year, data = mrg2, function(x) sum(!is.na(x)))
agg_juv <- aggregate(juv_GAM_mean ~ year, data = mrg2, function(x) sum(!is.na(x)))

#filter for valid years
agg_mrg <- dplyr::full_join(agg_br, agg_juv, by = 'year')
agg_mrg$j <- apply(agg_mrg[,2:3], 1, function(x) sum(x, na.rm = TRUE))
vyrs <- dplyr::filter(agg_mrg, j >=3)$year

mrg3 <- dplyr::filter(mrg2, year %in% vyrs)

#stop if species has fewer than 3 valid years
v_idx <- which(!is.na(mrg3$br_GAM_mean) | !is.na(mrg3$juv_GAM_mean))
df <- mrg3[v_idx,]
if (length(unique(df$year)) < 3 & NROW(df) < 3)
{
  stop('Not enough data')
}

# #number of cell/years with overlapping data
novr <- NROW(dplyr::filter(mrg3, !is.na(br_GAM_mean), !is.na(juv_GAM_mean)))

#define cells and years to be modeled
cells <- unique(mrg3$cell)
years <- unique(mrg3$year)
nyr <- length(years)
ncell <- length(cells)


# create adjacency matrix -------------------------------------------------

#make hexgrid
hexgrid6 <- dggridR::dgconstruct(res = 6)

#get hexgrid cell centers
cellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid6, cells)

#add lat col to df
mrg3$lat <- cellcenters$lat_deg

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
  mrg3 <- dplyr::filter(mrg3, cell %in% cells)
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

#create and fill sds, obs
sigma_br <- matrix(nrow = nyr, ncol = ncell)
br_obs <- matrix(nrow = nyr, ncol = ncell)
sigma_juv <- matrix(nrow = nyr, ncol = ncell)
juv_obs <- matrix(nrow = nyr, ncol = ncell)
br_PPC <- rep(NA, nyr * ncell)
juv_PPC <- rep(NA, nyr * ncell)

#number of observation and NAs for each year
len_br_obs <- rep(NA, nyr)
len_br_mis <- rep(NA, nyr)
len_juv_obs <- rep(NA, nyr)
len_juv_mis <- rep(NA, nyr)

#indices for observed and missing
ii_br_obs <- matrix(NA, nrow = nyr, ncol = ncell)
ii_br_mis <- matrix(NA, nrow = nyr, ncol = ncell)
ii_juv_obs <- matrix(NA, nrow = nyr, ncol = ncell)
ii_juv_mis <- matrix(NA, nrow = nyr, ncol = ncell)

#counter to fill
counter <- 1
for (i in 1:nyr)
{
  #i <- 1
  temp_yr <- dplyr::filter(mrg3, year == years[i])
  
  #don't need to manipulate position of sigmas
  sigma_br[i,] <- temp_yr$br_GAM_sd
  sigma_juv[i,] <- temp_yr$juv_GAM_sd
  
  for (j in 1:ncell)
  {
    #n <- 1
    #matrix with observed values with NAs
    br_PPC[counter] <- temp_yr$br_GAM_mean[j]
    juv_PPC[counter] <- temp_yr$juv_GAM_mean[j]
    counter <- counter + 1
  }
  
  #which are not NA
  no_na_br <- temp_yr$br_GAM_mean[which(!is.na(temp_yr$br_GAM_mean))]
  no_na_juv <- temp_yr$juv_GAM_mean[which(!is.na(temp_yr$juv_GAM_mean))]
  
  #br - pad end with NAs
  if (length(no_na_br) < ncell)
  {
    num_na_br <- ncell - length(no_na_br)
    
    #add NAs to end
    t_br_obs <- c(no_na_br, rep(NA, num_na_br))
    t_br_obs_ind <- c(which(!is.na(temp_yr$br_GAM_mean)), rep(NA, num_na_br))
    t_br_mis_ind <- c(which(is.na(temp_yr$br_GAM_mean)), rep(NA, length(no_na_br)))
    
    #fill objects
    ii_br_obs[i,] <- t_br_obs_ind
    ii_br_mis[i,] <- t_br_mis_ind
    br_obs[i,] <- t_br_obs
  } else {
    #no NAs to end (no mimssing values)
    br_obs[i,] <- no_na_br
    ii_br_mis[i,] <- which(!is.na(temp_yr$br_GAM_mean))
    br_obs[i,] <- which(is.na(temp_yr$br_GAM_mean))
  }
  
  #br - length of data/miss for each year
  len_br_obs[i] <- length(no_na_br)
  len_br_mis[i] <- ncell - length(no_na_br)
  
  #juv - pad end with NAs
  if (length(no_na_juv) < ncell)
  {
    num_na_juv <- ncell - length(no_na_juv)
    
    #add NAs to end
    t_juv_obs <- c(no_na_juv, rep(NA, num_na_juv))
    t_juv_obs_ind <- c(which(!is.na(temp_yr$juv_GAM_mean)), rep(NA, num_na_juv))
    t_juv_mis_ind <- c(which(is.na(temp_yr$juv_GAM_mean)), rep(NA, length(no_na_juv)))
    
    #fill objects
    ii_juv_obs[i,] <- t_juv_obs_ind
    ii_juv_mis[i,] <- t_juv_mis_ind
    juv_obs[i,] <- t_juv_obs
  } else {
    #no NAs to end (no mimssing values)
    juv_obs[i,] <- no_na_juv
    ii_juv_mis[i,] <- which(!is.na(temp_yr$juv_GAM_mean))
    juv_obs[i,] <- which(is.na(temp_yr$juv_GAM_mean))
  }
  
  #juv - length of data/miss for each year
  len_juv_obs[i] <- length(no_na_juv)
  len_juv_mis[i] <- ncell - length(no_na_juv)
}


#fill 0 where NA in y_obs - Stan does not like NA and zeros are not being used to estimate any param (y_obs is used to fill y)
br_obs[which(is.na(br_obs), arr.ind = TRUE)] <- 0
sigma_br[which(is.na(sigma_br), arr.ind = TRUE)] <- 0.1
ii_br_obs[which(is.na(ii_br_obs), arr.ind = TRUE)] <- 0
ii_br_mis[which(is.na(ii_br_mis), arr.ind = TRUE)] <- 0

juv_obs[which(is.na(juv_obs), arr.ind = TRUE)] <- 0
sigma_juv[which(is.na(sigma_juv), arr.ind = TRUE)] <- 0.1
ii_juv_obs[which(is.na(ii_juv_obs), arr.ind = TRUE)] <- 0
ii_juv_mis[which(is.na(ii_juv_mis), arr.ind = TRUE)] <- 0


#create data list for Stan
DATA <- list(J = ncell,
             cells = cells,
             N = nyr, 
             NJ = nyr * ncell,
             N_br_obs = len_br_obs,
             N_br_mis = len_br_mis,
             N_juv_obs = len_juv_obs,
             N_juv_mis = len_juv_mis,
             N_edges = nrow(ninds), 
             node1 = ninds[,1],
             node2 = ninds[,2],
             br_obs = br_obs,
             sigma_br = sigma_br,
             juv_obs = juv_obs,
             sigma_juv = sigma_juv,
             ii_br_obs = ii_br_obs,
             ii_br_mis = ii_br_mis,
             ii_juv_obs = ii_juv_obs,
             ii_juv_mis = ii_juv_mis,
             lat = cellcenters$lat_deg,
             br_PPC = br_PPC,
             juv_PPC = juv_PPC)


# Stan model --------------------------------------------------------------

#years in rows
#cells in columns

IAR <- '
data {
int<lower = 0> N;                                     // number of years
int<lower = 0> J;                                     // number of cells
int<lower = 0> NJ;                                    // number of cell/years
int<lower = 0> N_br_obs[N];                           // number of non-missing for each year
int<lower = 0> N_br_mis[N];                           // number missing for each year
int<lower = 0> N_juv_obs[N];                          
int<lower = 0> N_juv_mis[N];                          
int<lower = 0> N_edges;                               // number of edges in adjacency matrix
int<lower = 1, upper = J> node1[N_edges];             // node1[i] adjacent to node2[i]
int<lower = 1, upper = J> node2[N_edges];             // and node1[i] < node2[i]
vector[J] br_obs[N];                                   // observed response data (add NAs to end)
vector[J] juv_obs[N];
vector<lower = 0>[J] sigma_br[N];                      // observed sd of data (observation error)
vector<lower = 0>[J] sigma_juv[N];
int<lower = 0> ii_br_obs[N, J];                          // indices of observed data
int<lower = 0> ii_br_mis[N, J];                          // indices of missing data
int<lower = 0> ii_juv_obs[N, J];
int<lower = 0> ii_juv_mis[N, J];
vector<lower = 24, upper = 90>[J] lat;
}

parameters {
vector[J] br_mis[N];
vector[J] juv_mis[N];
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
real alpha;                                       // offset for juv to estimate back to lay date
}

transformed parameters {
vector[J] br[N];
vector[J] juv[N];
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
  br[i, ii_br_obs[i, 1:N_br_obs[i]]] = br_obs[i, 1:N_br_obs[i]];
  br[i, ii_br_mis[i, 1:N_br_mis[i]]] = br_mis[i, 1:N_br_mis[i]];
  juv[i, ii_juv_obs[i, 1:N_juv_obs[i]]] = juv_obs[i, 1:N_juv_obs[i]];
  juv[i, ii_juv_mis[i, 1:N_juv_mis[i]]] = juv_mis[i, 1:N_juv_mis[i]];
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
// offset = 29 represents 29 days between hatch and fledge
alpha ~ normal(29, 3);

for (i in 1:N)
{
  y_true_raw[i] ~ std_normal();
  // index array first (each year), then vector (for cells)
  target += -0.5 * dot_self(phi[i, node1] - phi[i, node2]);
  // soft sum to 0 constraint - J is number of phis per year
  sum(phi[i]) ~ normal(0, 0.001 * J);
  
  // observation model for y
  br[i] ~ normal(y_true[i], sigma_br[i]);
  juv[i] ~ normal(y_true[i] + alpha, sigma_juv[i]);
}
}

generated quantities {

vector[NJ] br_rep;
vector[NJ] juv_rep;
int<lower = 0> counter;

counter = 1;
for (n in 1:N)
{
  for (j in 1:J)
  {
    br_rep[counter] = normal_rng(y_true[n,j], sigma_br[n,j]);
    juv_rep[counter] = normal_rng(y_true[n,j] + alpha, sigma_juv[n,j]);
    counter = counter + 1;
  }
}
}'


# Run model ---------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.90
TREE_DEPTH <- 16
STEP_SIZE <- 0.0005
CHAINS <- 4
ITER <- as.numeric(args[2])

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
                     'alpha',
                     'phi',
                     'sigma_phi',
                     'sigma_y_true',
                     'y_true', 
                     'br_rep',
                     'juv_rep'),
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
while ((max(rhat_output) >= 1.02 | min(neff_output) < (CHAINS * 100)) & ITER < 20001)
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
                              'alpha',
                              'phi',
                              'sigma_phi',
                              'sigma_y_true',
                              'y_true', 
                              'br_rep',
                              'juv_rep'),
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
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', bj_IAR_out_dir))
saveRDS(fit, file = paste0(args[1], '-bj-iar-hm-stan_output-', run_date, '.rds'))

#save data to RDS (has which cells are modeled)
saveRDS(DATA, file = paste0(args[1], '-bj-iar-hm-stan_input-',  run_date, '.rds'))


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


# Checks -------------------------------------------------------------

# for PPC extract y_rep and transpose (so iter are rows as required by shiny stan)
br_rep_ch <- MCMCvis::MCMCpstr(fit, params = 'br_rep', type = 'chains')[[1]]
t_br_rep <- t(br_rep_ch)
juv_rep_ch <- MCMCvis::MCMCpstr(fit, params = 'juv_rep', type = 'chains')[[1]]
t_juv_rep <- t(juv_rep_ch)

#remove NA vals
na.br.rm <- which(is.na(DATA$br_PPC))
n_br_PPC <- DATA$br_PPC[-na.br.rm]
n_br_rep <- t_br_rep[, -na.br.rm]
na.juv.rm <- which(is.na(DATA$juv_PPC))
n_juv_PPC <- DATA$juv_PPC[-na.juv.rm]
n_juv_rep <- t_juv_rep[, -na.juv.rm]

#density overlay plot - first 100 iter
#modified bayesplot::ppc_dens_overlay function
tdata_br <- bayesplot::ppc_data(n_br_PPC, n_br_rep[1:100,])
tdata_juv <- bayesplot::ppc_data(n_juv_PPC, n_juv_rep[1:100,])


annotations <- data.frame(xpos = c(-Inf, -Inf),
                          ypos = c(Inf, Inf),
                          annotateText = c(paste0('Max Rhat: ', max(rhat_output)),
                                           paste0('Min n.eff: ', min(neff_output))),
                          hjustvar = c(0, 0),
                          vjustvar = c(4, 6))

p_br <- ggplot(tdata_br) +
  aes_(x = ~value) +
  stat_density(aes_(group = ~rep_id, color = "yrep"),
               data = function(x) dplyr::filter(x, !tdata_br$is_y),
               geom = "line", position = "identity", size = 0.25,
               alpha = 0.3, trim = FALSE, bw = 'nrd0', adjust = 1,
               kernel = 'gaussian', n = 1024) +
  stat_density(aes_(color = "y"), data = function(x) dplyr::filter(x, tdata_br$is_y),
               geom = "line", position = "identity", lineend = "round", size = 1, trim = FALSE,
               bw = 'nrd0', adjust = 1, kernel = 'gaussian', n = 1024) +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(size = 24)) +
  labs(colour = '') +
  scale_color_manual(values = c('red', 'black')) +
  ggtitle(paste0(args[1])) +
  geom_text(data = annotations, aes(x = xpos, y = ypos,
                                    hjust = hjustvar, vjust = vjustvar,
                                    label = annotateText),
            size = 3, col = 'black')

p_juv <- ggplot(tdata_juv) +
  aes_(x = ~value) +
  stat_density(aes_(group = ~rep_id, color = "yrep"),
               data = function(x) dplyr::filter(x, !tdata_juv$is_y),
               geom = "line", position = "identity", size = 0.25,
               alpha = 0.3, trim = FALSE, bw = 'nrd0', adjust = 1,
               kernel = 'gaussian', n = 1024) +
  stat_density(aes_(color = "y"), data = function(x) dplyr::filter(x, tdata_juv$is_y),
               geom = "line", position = "identity", lineend = "round", size = 1, trim = FALSE,
               bw = 'nrd0', adjust = 1, kernel = 'gaussian', n = 1024) +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(size = 24)) +
  labs(colour = '') +
  scale_color_manual(values = c('red', 'black')) +
  ggtitle(paste0(args[1])) +
  geom_text(data = annotations, aes(x = xpos, y = ypos,
                                    hjust = hjustvar, vjust = vjustvar,
                                    label = annotateText),
            size = 3, col = 'black')

ggsave(paste0(args[1], '_br_dens_overlay.pdf'), p_br)
ggsave(paste0(args[1], '_juv_dens_overlay.pdf'), p_juv)


# write model results to file ---------------------------------------------

options(max.print = 5e6)
sink(paste0(args[1], '-bj-iar-hm-stan_results-', run_date, '.txt'))
cat(paste0('bj IAR hm results ', args[1], ' \n'))
cat(paste0('Total minutes: ', round(run_time, digits = 2), ' \n'))
cat(paste0('Iterations: ', ITER, ' \n'))
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
cat(paste0('Num data points: ', length(v_idx), ' \n'))
cat(paste0('Num cell/year overlap: ', novr, ' \n'))
cat(paste0('Max Rhat: ', max(rhat_output), ' \n'))
cat(paste0('Min n.eff: ', min(neff_output), ' \n'))
print(model_summary)
sink()


# Plot pre-IAR/post_IAR halfmax estimates ------------------------------------------

#estimated half-max in grey, sd in white (derived from GAM)

#extract median and sd estimates for mu params
med_fit <- MCMCvis::MCMCpstr(fit, params = 'y_true', func = median)[[1]]
sd_fit <- MCMCvis::MCMCpstr(fit, params = 'y_true', func = sd)[[1]]

#transform cells to grid
cell_grid <- dggridR::dgcellstogrid(hexgrid6, cells)
cell_grid$cell <- as.numeric(cell_grid$cell)
cell_centers <- dggridR::dgSEQNUM_to_GEO(hexgrid6, cells)
ll_df <- data.frame(cell = cells,
                    lon_deg = cell_centers$lon_deg,
                    lat_deg = cell_centers$lat_deg)

#load maps
worldmap <- data.frame(maps::map("world", plot = FALSE)[c("x", "y")])

#min/max for plotting using input data
#f_rng <- c(range(f_out$arr_GAM_hm_mean, na.rm = TRUE), range(med_fit, na.rm = TRUE))
#MIN <- round(min(f_rng))
#MAX <- round(max(f_rng))

#min/max for plotting using output data
MIN <- (round(min(med_fit)) - 1)
MAX <- (round(max(med_fit)) + 1)


#create output image dir if it doesn't exist
ifelse(!dir.exists(paste0(dir, 'Bird_Phenology/Figures/pre_post_IAR_maps/', IAR_out_dir)),
       dir.create(paste0(dir, 'Bird_Phenology/Figures/pre_post_IAR_maps/', IAR_out_dir)),
       FALSE)


setwd(paste0(dir, 'Bird_Phenology/Figures/pre_post_IAR_maps/', IAR_out_dir))

#loop plots for each year
for (i in 1:length(years))
{
  #i <- 1
  
  #filter data for year[i]
  f_out_filt <- dplyr::filter(mrg3, year == years[i])
  
  #merge hex spatial data with GAM data
  to_plt <- dplyr::inner_join(f_out_filt, cell_grid, by = 'cell')
  to_plt2 <- dplyr::inner_join(to_plt, ll_df, by = 'cell')
  
  #pre-IAR
  p <- ggplot() +
    geom_path(data = worldmap,
              aes(x = x, y = y), color = 'black') +
    coord_map("ortho", orientation = c(35, -80, 0),
              xlim = c(-100, -55), ylim = c(23, 66)) +
    geom_polygon(data = to_plt2, aes(x = long, y = lat.y, group = group, fill = arr_GAM_hm_mean),
                 alpha = 0.4) +
    geom_path(data = to_plt2, aes(x = long, y = lat.y, group = group),
              alpha = 0.4, color = 'black') +
    scale_fill_gradientn(colors = c('red', 'blue'),
                         limits = c(MIN, MAX)) +
    labs(fill = 'Estimated Arrival') +
    annotate('text', x = to_plt2$lon_deg, y = to_plt2$lat_deg + 0.5,
             label = round(to_plt2$arr_GAM_hm_mean, digits = 0), col = 'black', alpha = 0.2,
             size = 4) +
    annotate('text', x = to_plt2$lon_deg, y = to_plt2$lat_deg - 0.5,
             label = round(to_plt2$arr_GAM_hm_sd, digits = 0), col = 'white', alpha = 0.3,
             size = 3) +
    ggtitle(paste0(f_out_filt$species[1], ' - ', f_out_filt$year[1], ' - Pre-IAR')) +
    theme_bw() +
    xlab('Longitude') +
    ylab('Latitude')
  
  ggsave(plot = p,
         filename = paste0(f_out_filt$species[1], '_', f_out_filt$year[1], '-pre_IAR.pdf'))
  
  
  #post-IAR
  t_med_fit <- med_fit[i,]
  t_sd_fit <- sd_fit[i,]
  
  #median of mu and sd of mu
  m_fit <- data.frame(med_mu = t_med_fit, sd_mu = t_sd_fit, cell = cells)
  
  #merge hex spatial data with GAM data
  to_plt_post <- dplyr::inner_join(m_fit, cell_grid, by = 'cell')
  to_plt2_post <- dplyr::inner_join(to_plt_post, ll_df, by = 'cell')
  
  
  #plot
  p_post <- ggplot() +
    geom_path(data = worldmap,
              aes(x = x, y = y), color = 'black') +
    coord_map("ortho", orientation = c(35, -80, 0),
              xlim = c(-100, -55), ylim = c(23, 66)) +
    geom_polygon(data = to_plt2_post, aes(x = long, y = lat, group = group, fill = med_mu),
                 alpha = 0.4) +
    geom_path(data = to_plt2_post, aes(x = long, y = lat, group = group),
              alpha = 0.4, color = 'black') +
    scale_fill_gradientn(colors = c('red', 'blue'),
                         limits = c(MIN, MAX)) +
    labs(fill = 'Estimated Arrival') +
    annotate('text', x = to_plt2_post$lon_deg, y = to_plt2_post$lat_deg + 0.5,
             label = round(to_plt2_post$med_mu, digits = 0), col = 'black', alpha = 0.2,
             size = 4) +
    annotate('text', x = to_plt2_post$lon_deg, y = to_plt2_post$lat_deg - 0.5,
             label = round(to_plt2_post$sd_mu, digits = 0), col = 'white', alpha = 0.3,
             size = 3) +
    ggtitle(paste0(f_out_filt$species[1], ' - ', f_out_filt$year[1], ' - Post-IAR')) +
    theme_bw() +
    xlab('Longitude') +
    ylab('Latitude')
  
  ggsave(plot = p_post,
         filename = paste0(f_out_filt$species[1], '_', f_out_filt$year[1], '-post_IAR.pdf'))
}


# Trace plots with PPO ----------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', bj_IAR_out_dir))

#alpha_gamma ~ normal(0, 100)
PR <- rnorm(10000, 0, 100)
MCMCvis::MCMCtrace(fit,
                   params = 'alpha_gamma',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args[1], '-trace_alpha_gamma-', run_date, '.pdf'))

#beta_gamma ~ normal(3, 3)
PR <- rnorm(10000, 3, 3)
MCMCvis::MCMCtrace(fit,
                   params = 'beta_gamma',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args[1], '-trace_beta_gamma-', run_date, '.pdf'))

#sigma_gamma ~ halfnormal(0, 10)
PR_p <- rnorm(10000, 0, 10)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_gamma',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args[1], '-trace_sigma_gamma-', run_date, '.pdf'))

#sigma_beta0 ~ HN(0, 10)
PR_p <- rnorm(10000, 0, 10)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_beta0',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args[1], '-trace_sigma_beta0-', run_date, '.pdf'))

#sigma_y_true ~ HN(0, 10)
PR_p <- rnorm(10000, 0, 10)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_y_true',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args[1], '-trace_sigma_y_true-', run_date, '.pdf'))

#sigma_phi ~ HN(0, 10)
PR_p <- rnorm(10000, 0, 10)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_phi',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args[1], '-trace_sigma_phi-', run_date, '.pdf'))

#alpha ~ N(29, 3)
PR <- rnorm(10000, 29, 3)
MCMCvis::MCMCtrace(fit,
                   params = 'alpha',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args[1], '-trace_alpha-', run_date, '.pdf'))

if ('Rplots.pdf' %in% list.files())
{
  file.remove('Rplots.pdf')
}

print('I completed!')

