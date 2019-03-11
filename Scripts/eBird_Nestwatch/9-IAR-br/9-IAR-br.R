######################
# 9 - IAR model breeding
#
# Fit IAR model
#
# VVV MODEL VVV
# i = cell
# j = year
# y_{obs[i,j]} \sim N(y_{true[i,j]}, \sigma_{y[i,j]})
# y_{true[i,j]} = \beta_{0[j]} + \gamma_{[i]} + \nu_{[i,j]} * \sigma_{\nu[j]}
# \beta_{0[j]} \sim N(0, \sigma_{\beta_{0}})
# \gamma_{[i]} \sim N(\mu_{\gamma[i]}, \sigma_{\gamma})
# \mu_{\gamma[i]} = \alpha_{\gamma} + \beta_{\gamma} * lat_{[i]}
# \nu_{[i,j]} = \sqrt{1 - \rho} * \theta[i,j] + \sqrt{\frac{\rho}{sf}} * \phi_{[i,j]}
# \sigma_{\nu[j]} \sim LN(\mu_{\sigma_{\nu}}, 0.7)
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
dir <- '/UCHC/LABS/Tingley/phenomismatch/'



# db/hm query dir ------------------------------------------------------------

IAR_in_dir <- 'IAR_br_input_2019-02-02'
IAR_out_dir <- 'IAR_br_output_2019-02-24'



# Load packages -----------------------------------------------------------

library(rstan)
library(INLA)
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
IAR_out_date <- substr(IAR_out_dir, start = 12, stop = 21)



# species arg -----------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
#args <- as.character('Catharus_minimus')
#args <- as.character('Empidonax_virescens')
#args <- as.character('Vireo_olivaceus')


#species for which one cell had to be dropped, bc it did not border any others
#args <- as.character('Cistothorus_palustris')
#args <- as.character('Setophaga_dominica')
#args <- as.character('Melospiza_lincolnii')
#args <- as.character('Zonotrichia_leucophrys')
#args <- as.character('Bombycilla_cedrorum')


# Filter data by species/years ------------------------------------------------------

#read in master df
df_master <- readRDS(paste0('IAR_br_input-', IAR_in_date, '.rds'))

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
    #j <- 4
    dists <- geosphere::distm(c(cellcenters$lon_deg[i], cellcenters$lat_deg[i]),
                              c(cellcenters$lon_deg[j], cellcenters$lat_deg[j]))
    adjacency_matrix[i,j] <- as.numeric((dists/1000) > 0 & (dists/1000) < 311)
  }
}

#indices for 1s
ninds <- which(adjacency_matrix == 1, arr.ind = TRUE)


#if a cell doesn't border any other cells, drop it and redefine objects
DROP <- FALSE
if (max(ninds) < ncell)
{
  s_cols <- apply(adjacency_matrix, 2, function(x) sum(x, na.rm = TRUE))
  s_rows <- apply(adjacency_matrix, 1, function(x) sum(x, na.rm = TRUE))
  to.rm.ind <- which((s_cols + s_rows) == 0)
  
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



# Estimate scaling factor for BYM2 model with INLA ------------------------

#Build the adjacency matrix using INLA library functions
adj.matrix <- Matrix::sparseMatrix(i = ninds[,1], j = ninds[,2], x = 1, symmetric = TRUE)

#The IAR precision matrix (note! This is singular)
Q <- Matrix::Diagonal(ncell, Matrix::rowSums(adj.matrix)) - adj.matrix
#Add a small jitter to the diagonal for numerical stability (optional but recommended)
Q_pert <- Q + Matrix::Diagonal(ncell) * max(diag(Q)) * sqrt(.Machine$double.eps)

# Compute the diagonal elements of the covariance matrix subject to the 
# constraint that the entries of the ICAR sum to zero.
# See the inla.qinv function help for further details.
Q_inv <- INLA::inla.qinv(Q_pert, 
                         constr = list(A = matrix(1, 1, ncell), e = 0))

#Compute the geometric mean of the variances, which are on the diagonal of Q.inv
scaling_factor <- exp(mean(log(diag(Q_inv))))



# create Stan data object -------------------------------------------------

#create and fill sds, obs, data for PPC, and temp
sigma_y_in <- matrix(nrow = ncell, ncol = nyr)
y_obs_in <- matrix(nrow = ncell, ncol = nyr)
y_PPC <- rep(NA, ncell * nyr)

#number of observation and NAs for each year
len_y_obs_in <- rep(NA, nyr)
len_y_mis_in <- rep(NA, nyr)

#indices for observed and missing
ii_obs_in <- matrix(NA, nrow = ncell, ncol = nyr)
ii_mis_in <- matrix(NA, nrow = ncell, ncol = nyr)

#counter to fill y_PPC
counter <- 1
for (j in 1:nyr)
{
  #j <- 16
  temp_yr <- dplyr::filter(f_out, year == years[j])
  
  #don't need to manipulate position of sigmas
  sigma_y_in[,j] <- temp_yr$HM_sd
  
  for (n in 1:ncell)
  {
    #n <- 1
    #matrix with observed values with NAs
    y_PPC[counter] <- temp_yr$HM_mean[n]
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
    ii_obs_in[,j] <- t_obs_in
    ii_mis_in[,j] <- t_mis_in
    y_obs_in[,j] <- t_y_obs_in
  } else {
    #no NAs to end (no mimssing values)
    y_obs_in[,j] <- no_na
    ii_mis_in[,j] <- which(!is.na(temp_yr$HM_mean))
    y_obs_in[,j] <- which(is.na(temp_yr$HM_mean))
  }
  
  #length of data/miss for each year
  len_y_obs_in[j] <- length(no_na)
  len_y_mis_in[j] <- ncell - length(no_na)
}


#fill 0 where NA in y_obs - Stan does not like NA and zeros are not being used to estimate any param (y_obs is used to fill y)
y_obs_in[which(is.na(y_obs_in), arr.ind = TRUE)] <- 0
#see script here showing this value fo sigma_y_in has no impact on y_true: Archive/test_sigma_insert.R
sigma_y_in[which(is.na(sigma_y_in), arr.ind = TRUE)] <- 0.1
ii_obs_in[which(is.na(ii_obs_in), arr.ind = TRUE)] <- 0
ii_mis_in[which(is.na(ii_mis_in), arr.ind = TRUE)] <- 0


#create data list for Stan
DATA <- list(J = nyr,
             N = ncell, 
             NJ = nyr * ncell,
             N_obs = len_y_obs_in,
             N_mis = len_y_mis_in,
             N_edges = nrow(ninds), 
             node1 = ninds[,1],
             node2 = ninds[,2],
             y_obs = y_obs_in,
             sigma_y = sigma_y_in,
             scaling_factor = scaling_factor,
             ii_obs = ii_obs_in,
             ii_mis = ii_mis_in,
             lat = cellcenters$lat_deg)



# Stan model --------------------------------------------------------------


# Model in latex
# ==============
# y_{obs[i,j]} \sim N(y_{true[i,j]}, \sigma_{y[i,j]})
# y_{true[i,j]} = \beta_{0[j]} + \gamma_{[i]} + \nu_{[i,j]} * \sigma_{\nu[j]}
# \beta_{0[j]} \sim N(0, \sigma_{\beta_{0}})
# \gamma_{[i]} \sim N(\mu_{\gamma[i]}, \sigma_{\gamma})
# \mu_{\gamma[i]} = \alpha_{\gamma} + \beta_{\gamma} * lat_{[i]}
# \nu_{[i,j]} = \sqrt{1 - \rho} * \theta[i,j] + \sqrt{\frac{\rho}{sf}} * \phi_{[i,j]}
# \sigma_{\nu[j]} \sim LN(\mu_{\sigma_{\nu}}, 0.7)
# \rho \sim Beta(0.5, 0.5)
# \theta \sim N(0, 1)
# \sigma_{\beta_{0}} \ sim HN(0, 5)
# \alpha_{\gamma} \sim N(0, 30)
# \beta_{\gamma} \sim N(2, 3)
# \sigma_{\gamma} \sim HN(0, 5)
# \mu_{\sigma_{\nu}} \sim N(0, 1.5)
# \phi_{[i,j]} \sim N(0, [D - W]^{-1})
# \forall j \in \left \{1, ..., J  \right \}; \sum_{i}{} \phi_{[i,j]} = 0


#cells in rows
#years in columns

IAR_2 <- '
data {
int<lower = 0> J;                                     // number of years      
int<lower = 0> N;                                     // number of cells
int<lower = 0> NJ;                                    // number of cell/years
int<lower = 0> N_obs[J];                              // number of non-missing for each year
int<lower = 0> N_mis[J];                              // number missing for each year
int<lower = 0> N_edges;                               // number of edges in adjacency matrix
int<lower = 1, upper = N> node1[N_edges];             // node1[i] adjacent to node2[i]
int<lower = 1, upper = N> node2[N_edges];             // and node1[i] < node2[i]
real<lower = 0, upper = 200> y_obs[N, J];             // observed response data (add NAs to end)
real<lower = 0> sigma_y[N, J];                        // observed sd of data (observation error)
int<lower = 0> ii_obs[N, J];                          // indices of observed data
int<lower = 0> ii_mis[N, J];                          // indices of missing data
real<lower = 0> scaling_factor;                       // scales variances of spatial effects (estimated from INLA)
real<lower = 26, upper = 90> lat[N];
}

parameters {
real<lower = 1, upper = 200> y_mis[N, J];             // missing response data
real alpha_gamma_raw;
real beta_gamma_raw;                                       // effect of latitude
real<lower = 0> sigma_gamma_raw;
vector[N] gamma_raw;
matrix[N, J] phi;                                     // spatial error component (scaled to N(0,1))
matrix[N, J] theta;                                   // non-spatial error component (scaled to N(0,1))
real<lower = 0, upper = 1> rho;                       // proportion unstructured vs spatially structured variance
vector<lower = 0>[J] sigma_nu_raw;
real beta0_raw[J];
real mu_sn_raw;
real<lower = 0> sigma_beta0_raw;
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

for (i in 1:N)
{
  mu_gamma[i] = alpha_gamma + beta_gamma * lat[i];
  gamma[i] = gamma_raw[i] * sigma_gamma + mu_gamma[i];
}

for (j in 1:J)
{
  nu[,j] = sqrt(1 - rho) * theta[,j] + sqrt(rho / scaling_factor) * phi[,j]; // combined spatial/non-spatial
  sigma_nu[j] = exp(sigma_nu_raw[j] * 0.7 + mu_sn);               //implies sigma_nu[j] ~ lognormal(mu_sn, 0.7) 
  beta0[j] = beta0_raw[j] * sigma_beta0;
  y_true[,j] = beta0[j] + gamma + nu[,j] * sigma_nu[j];
  
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
alpha_gamma_raw ~ normal(0, 1);
beta_gamma_raw ~ normal(0, 1);
sigma_gamma_raw ~ normal(0, 1);
sigma_nu_raw ~ normal(0, 1);
rho ~ beta(0.5, 0.5);
gamma_raw ~ normal(0, 1);
beta0_raw ~ normal(0, 1);
mu_sn_raw ~ normal(0, 1);
sigma_beta0_raw ~ normal(0, 1);

for (j in 1:J)
{
  theta[,j] ~ normal(0, 1);
  target += -0.5 * dot_self(phi[node1, j] - phi[node2, j]);
  sum(phi[,j]) ~ normal(0, 0.001 * N);
  
  y[,j] ~ normal(y_true[,j], sigma_y[,j]);
}

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



#save to RDS
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_out_dir))
#setwd("~/Google_Drive/R/Bird_Phenology/Data/Processed/Empidonax_virescens_test_no_slope/Ev_ns_ye")
saveRDS(fit, file = paste0(args, '-', IAR_out_date, '-iar-stan_output.rds'))





# Calc diagnostics and rerun if needed ------------------------------------

num_diverge <- rstan::get_num_divergent(fit)
num_tree <- rstan::get_num_max_treedepth(fit)
num_BFMI <- length(rstan::get_low_bfmi_chains(fit))


# #rerun model if things didn't go well
# while (sum(c(num_diverge, num_tree, num_BFMI)) > 0 & DELTA <= 0.99)
# {
#   DELTA <- DELTA + 0.01
#   TREE_DEPTH <- TREE_DEPTH + 1
#   STEP_SIZE <- STEP_SIZE * 0.75
# 
#   tt <- proc.time()
#   fit <- rstan::stan(model_code = IAR_2,
#                      data = DATA,
#                      chains = CHAINS,
#                      iter = ITER,
#                      cores = CHAINS,
#                      pars = c('sigma_beta0', 'beta0',
#                               'alpha_gamma', 'beta_gamma', 'sigma_gamma', 'gamma',
#                               'sigma_nu', 'mu_sn', 'rho', 'nu', 'theta', 'phi', 
#                               'y_true', 'y_rep'),
#                      control = list(adapt_delta = DELTA,
#                                     max_treedepth = TREE_DEPTH,
#                                     stepsize = STEP_SIZE))
#   run_time <- (proc.time()[3] - tt[3]) / 60
# 
#   num_diverge <- rstan::get_num_divergent(fit)
#   num_tree <- rstan::get_num_max_treedepth(fit)
#   num_BFMI <- length(rstan::get_low_bfmi_chains(fit))
# }






# Checks -------------------------------------------------------------

# setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_out_dir))
# fit <- readRDS('Catharus_minimus-2019-02-21-iar-stan_output.rds')
# fit <- readRDS('Empidonax_virescens-2019-02-21-iar-stan_output.rds')
# fit <- readRDS('Vireo_olivaceus-2019-02-13-IAR_stan-test-3.rds')

# for PPC extract y_rep and transpose (so iter are rows as required by shiny stan)
y_rep_ch <- MCMCvis::MCMCpstr(fit, params = 'y_rep', type = 'chains')[[1]]
t_y_rep <- t(y_rep_ch)

#for each iter see if y_rep value is greater or less than true value
tsum <- 0
for (i in 1:NROW(t_y_rep))
{
  #i <- 1
  temp <- sum(t_y_rep[i,] > y_PPC, na.rm = TRUE)
  tsum <- tsum + temp
}

#number of true y values that are not NA
l_PPC <- sum(!is.na(y_PPC))
PPC_p <- tsum / (l_PPC * NROW(t_y_rep))


# #shiny stan
#for shiny stan PPC
# na.y.rm <- which(is.na(y_PPC))
# n_y_PPC <- y_PPC[-na.y.rm]
# n_t_y_rep <- t_y_rep[,-na.y.rm]

# library(shinystan)
# launch_shinystan(fit)

#PPC
# bayesplot::ppc_stat(n_y_PPC, n_t_y_rep, stat = 'mean')
# bayesplot::ppc_stat(n_y_PPC, n_t_y_rep, stat = 'max')
# bayesplot::ppc_stat(n_y_PPC, n_t_y_rep, stat = 'min')



# write model results to file ---------------------------------------------


options(max.print = 50000)
sink(paste0(args, '-', IAR_out_date, '-iar-stan_results.txt'))
cat(paste0('IAR results ', args, ' \n'))
cat(paste0('Total minutes: ', round(run_time, digits = 2), ' \n'))
cat(paste0('Adapt delta: ', DELTA, ' \n'))
cat(paste0('Max tree depth: ', TREE_DEPTH, ' \n'))
cat(paste0('Step size: ', STEP_SIZE, ' \n'))
cat(paste0('Number of divergences: ', num_diverge, ' \n'))
cat(paste0('Number of tree exceeds: ', num_tree, ' \n'))
cat(paste0('Number chains low BFMI: ', num_BFMI, ' \n'))
cat(paste0('PPC p-val: ', round(PPC_p, 3), ' \n'))
cat(paste0('Cell drop: ', DROP, ' \n'))
print(MCMCvis::MCMCsummary(fit, Rhat = TRUE, n.eff = TRUE, round = 2))
sink()




# Plot pre-IAR/post_IAR halfmax estimates ------------------------------------------

#estimated half-max in grey, sd in white (derived from logit cubic)

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
usamap <- data.frame(maps::map("world", "USA", plot = FALSE)[c("x", "y")])
canadamap <- data.frame(maps::map("world", "Canada", plot = FALSE)[c("x", "y")])
mexicomap <- data.frame(maps::map("world", "Mexico", plot = FALSE)[c("x", "y")])

#determine min/max for plotting
f_rng <- c(range(f_out$HM_mean, na.rm = TRUE), range(med_fit, na.rm = TRUE))
MIN <- round(min(f_rng))
MAX <- round(max(f_rng))


#read in breeding/migration range shp file
setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/shapefiles/'))
sp_rng <- rgdal::readOGR(f_out$shp_fname[1], verbose = FALSE)

#filter by breeding (2) and migration (4) range - need to convert spdf to sp
nrng <- sp_rng[which(sp_rng$SEASONAL == 2 | sp_rng$SEASONAL == 4),]
nrng_sp <- sp::SpatialPolygons(nrng@polygons)

#filter by resident (1) and over winter (3) range - need to convert spdf to sp
nrng_rm <- sp_rng[which(sp_rng$SEASONAL == 1 | sp_rng$SEASONAL == 3),]
nrng_rm_sp <- sp::SpatialPolygons(nrng_rm@polygons)


#plotting species range
nrng@data$id <- rownames(nrng@data)
nrng.points <- ggplot2::fortify(nrng, region = "id")
nrng.df <- plyr::join(nrng.points, nrng@data, by = "id")

nrng_rm@data$id <- rownames(nrng_rm@data)
nrng_rm.points <- ggplot2::fortify(nrng_rm, region = "id")
nrng_rm.df <- plyr::join(nrng_rm.points, nrng_rm@data, by = "id")


#create output image dir if it doesn't exist

ifelse(!dir.exists(paste0(dir, 'Bird_Phenology/Figures/pre_post_IAR_maps/', IAR_out_date)),
       dir.create(paste0(dir, 'Bird_Phenology/Figures/pre_post_IAR_maps/', IAR_out_date)),
       FALSE)


setwd(paste0(dir, 'Bird_Phenology/Figures/pre_post_IAR_maps/', IAR_out_date))

#loop plots for each year
for (i in 1:length(years))
{
  #i <- 1
  
  #filter data for year[i]
  f_out_filt <- filter(f_out, year == years[i])
  
  #merge hex spatial data with HM data
  to_plt <- dplyr::inner_join(f_out_filt, cell_grid, by = 'cell')
  to_plt2 <- dplyr::inner_join(to_plt, ll_df, by = 'cell')
  
  #pre-IAR
  p <- ggplot() +
    geom_path(data = usamap,
              aes(x = x, y = y), color = 'black') +
    geom_path(data = canadamap,
              aes(x = x, y = y), color = 'black') +
    geom_path(data = mexicomap,
              aes(x = x, y = y), color = 'black') +
    # geom_polygon(data = nrng.df,
    #           aes(x = long, y = lat, group=group), fill = 'green', alpha = 0.4) +
    # geom_polygon(data = nrng_rm.df,
    #              aes(x = long, y = lat, group=group), fill = 'orange', alpha = 0.4) +
    coord_map("ortho", orientation = c(35, -80, 0),
              xlim = c(-100, -55), ylim = c(25, 66)) +
    geom_polygon(data = to_plt2, aes(x = long, y = lat.y, group = group, fill = HM_mean),
                 alpha = 0.4) +
    geom_path(data = to_plt2, aes(x = long, y = lat.y, group = group),
              alpha = 0.4, color = 'black') +
    scale_fill_gradientn(colors = c('red', 'blue'),
                         limits = c(MIN, MAX)) +
    labs(fill = 'Estimated Arrival') +
    annotate('text', x = to_plt2$lon_deg, y = to_plt2$lat_deg + 0.5,
             label = round(to_plt2$HM_mean, digits = 0), col = 'black', alpha = 0.2,
             size = 4) +
    annotate('text', x = to_plt2$lon_deg, y = to_plt2$lat_deg - 0.5,
             label = round(to_plt2$HM_sd, digits = 0), col = 'white', alpha = 0.3,
             size = 3) +
    ggtitle(paste0(f_out_filt$species[1], ' - ', f_out_filt$year[1], ' - Pre-IAR')) +
    theme_bw() +
    xlab('Longitude') +
    ylab('Latitude')
  
  ggsave(plot = p,
         filename = paste0(f_out_filt$species[1], '_', f_out_filt$year[1], '-pre_IAR.pdf'))
  
  
  #post-IAR
  t_med_fit <- med_fit[,i]
  t_sd_fit <- sd_fit[,i]
  
  #median of mu and sd of mu
  m_fit <- data.frame(med_mu = t_med_fit, sd_mu = t_sd_fit, cell = cells)
  
  #merge hex spatial data with HM data
  to_plt_post <- dplyr::inner_join(m_fit, cell_grid, by = 'cell')
  to_plt2_post <- dplyr::inner_join(to_plt_post, ll_df, by = 'cell')
  
  
  #plot
  p_post <- ggplot() +
    geom_path(data = usamap,
              aes(x = x, y = y), color = 'black') +
    geom_path(data = canadamap,
              aes(x = x, y = y), color = 'black') +
    geom_path(data = mexicomap,
              aes(x = x, y = y), color = 'black') +
    # geom_polygon(data = nrng.df,
    #           aes(x = long, y = lat, group=group), fill = 'green', alpha = 0.4) +
    # geom_polygon(data = nrng_rm.df,
    #              aes(x = long, y = lat, group=group), fill = 'orange', alpha = 0.4) +
    coord_map("ortho", orientation = c(35, -80, 0),
              xlim = c(-100, -55), ylim = c(25, 66)) +
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

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_out_dir))


#alpha_gamma ~ normal(0, 30)
PR <- rnorm(10000, 0, 30)
MCMCvis::MCMCtrace(fit,
                   params = 'alpha_gamma',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-', IAR_out_date, '-trace_alpha_gamma.pdf'))

#beta_gamma ~ normal(2, 3)
PR <- rnorm(10000, 2, 3)
MCMCvis::MCMCtrace(fit,
                   params = 'beta_gamma',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-', IAR_out_date, '-trace_beta_gamma.pdf'))

#sigma_gamma ~ halfnormal(0, 5)
PR_p <- rnorm(10000, 0, 5)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_gamma',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-', IAR_out_date, '-trace_sigma_gamma.pdf'))

#mu_sn ~ halfnormal(0, 1.5)
PR <- rnorm(10000, 0, 1.5)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_sn',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-', IAR_out_date, '-trace_mu_sn.pdf'))

#rho ~ beta(0.5, 0.5)
PR <- rbeta(10000, 0.5, 0.5)
MCMCvis::MCMCtrace(fit,
                   params = 'rho',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-', IAR_out_date, '-trace_rho.pdf'))

#sigma_beta0 ~ HN(0, 5)
PR_p <- rnorm(10000, 0, 5)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_beta0',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-', IAR_out_date, '-trace_sigma_beta0.pdf'))


if ('Rplots.pdf' %in% list.files())
{
  file.remove('Rplots.pdf')
}


print('I completed!')

