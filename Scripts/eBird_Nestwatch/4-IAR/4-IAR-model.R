######################
# 4 - IAR model
#
# Fit IAR model - true arrival dates for each species-cell-year are modeled as latent states, with the observed 
# state derived from 2-logit-cubic.R. 
#
# Model fits different spatial (phi)/non-spatial (theta) params for each year but ONE rho (the degree of spatial vs. non-spatial error) and ONE sigma (scale factor for how large those errors are)
######################

#Stan resources:
#http://mc-stan.org/users/documentation/case-studies/icar_stan.html
#http://mc-stan.org/users/documentation/case-studies/divergences_and_bias.html
#https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html
#https://chi-feng.github.io/mcmc-demo/app.html#HamiltonianMC,standard
#https://groups.google.com/forum/#!msg/stan-users/zOjAeJC4x_E/OyCOfJo8AwAJ (non-centered parameterization on sd should be fine)

# Top-level dir -----------------------------------------------------------

#desktop/laptop
#dir <- '~/Google_Drive/R/'

#Xanadu
dir <- '/UCHC/LABS/Tingley/phenomismatch/'



# db/hm query dir ------------------------------------------------------------

db_dir <- 'eBird_query_2018-10-15'
hm_dir <- 'halfmax_species_2018-10-16'
IAR_in_dir <- 'IAR_input_2018-11-12'
IAR_out_dir <- 'IAR_output_2019-01-16-supp'

#create output dir if it doesn't exist - need to create before running script bc STD out and STD error are written there
# ifelse(!dir.exists(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_out_dir)), 
#        dir.create(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_out_dir)), 
#        FALSE)

# runtime -----------------------------------------------------------------

tt <- proc.time()



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
#args <- as.character('Vireo_olivaceus')
#args <- as.character('Catharus_minimus')
#args <- as.character('Empidonax_virescens')

#species for which one cell had to be dropped, bc it did not border any others
#args <- as.character('Cistothorus_palustris')
#args <- as.character('Setophaga_dominica')
#args <- as.character('Melospiza_lincolnii')
#args <- as.character('Zonotrichia_leucophrys')
#args <- as.character('Bombycilla_cedrorum')


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


# filter cells ------------------------------------------------------------

# #cells ID'd as a 'test set'
# cells_to_keep <- c(622, 595, 567, 594, 621, 648, 649, 676, 675, 647, 619, 
#                    620, 621, 593, 592, 619, 618, 591, 564, 565, 566, 646, 
#                    674, 702, 703)



# create adjacency matrix -------------------------------------------------

#make hexgrid
hexgrid6 <- dggridR::dgconstruct(res = 6)

#get hexgrid cell centers
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

#create and fill sds and obs
sigma_y_in <- matrix(nrow = ncell, ncol = nyr)
y_obs_in <- matrix(nrow = ncell, ncol = nyr)

#number of observation and NAs for each year
len_y_obs_in <- rep(NA, nyr)
len_y_mis_in <- rep(NA, nyr)

#indices for observed and missing
ii_obs_in <- matrix(NA, nrow = ncell, ncol = nyr)
ii_mis_in <- matrix(NA, nrow = ncell, ncol = nyr)

for (j in 1:nyr)
{
  #j <- 16
  temp_yr_p <- dplyr::filter(f_out, year == years[j])
  temp_yr <- temp_yr_p[order(temp_yr_p$cell),]
  
  sigma_y_in[,j] <- temp_yr$HM_sd
  
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
sigma_y_in[which(is.na(sigma_y_in), arr.ind = TRUE)] <- 0.1
ii_obs_in[which(is.na(ii_obs_in), arr.ind = TRUE)] <- 0
ii_mis_in[which(is.na(ii_mis_in), arr.ind = TRUE)] <- 0



#create data list for Stan
DATA <- list(J = nyr,
             N = ncell, 
             N_obs = len_y_obs_in,
             N_mis = len_y_mis_in,
             N_edges = nrow(ninds), 
             node1 = ninds[,1],
             node2 = ninds[,2],
             y_obs = y_obs_in,
             sigma_y = sigma_y_in,
             scaling_factor = scaling_factor,
             ii_obs = ii_obs_in,
             ii_mis = ii_mis_in)



# Stan model --------------------------------------------------------------

#spatial and non-spatial component (reparameterized) - also model obs error
#uses scaling factor from INLA to set priors more easily
#matrix[N,J] -> cells in rows, years in columns


IAR_bym2 <- '
data {
int<lower = 0> J;                                     // number of years      
int<lower = 0> N;                                     // number of cells
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
}

parameters {
real<lower = 1, upper = 200> y_mis[N, J];             // missing response data
real beta0[J];                                        // intercept
matrix[N, J] theta;                                   // non-spatial error component (centered on 0)
matrix[N, J] phi;                                     // spatial error component (centered on 0)
real<lower = 0> sigma_raw[J];                                // scaling factor for spatial and non-spatial components
real<lower = 0, upper = 1> rho;                       // proportion unstructured vs spatially structured variance
real<lower = 0> mu_sigma_raw;
}

transformed parameters {
real<lower = 0, upper = 200> y[N, J];                 // response data to be modeled

real<lower = 0> mu_sigma;
real<lower = 0> sigma[J];

matrix[N, J] convolved_re;                            // spatial and non-spatial component
matrix[N, J] mu;                                      // latent true halfmax values

mu_sigma = 3 * mu_sigma_raw;                          // non-centered parameterization

for (j in 1:J)
{
  sigma[j] = mu_sigma + sigma_raw[j] * 3;   // non-centered parameterization
}

for (j in 1:J)
{
  convolved_re[,j] = sqrt(1 - rho) * theta[,j] + sqrt(rho / scaling_factor) * phi[,j];
  mu[,j] = beta0[j] + convolved_re[,j] * sigma[j];
}

// indexing to avoid NAs
for (j in 1:J)
{
  y[ii_obs[1:N_obs[j], j], j] = y_obs[1:N_obs[j], j];
  y[ii_mis[1:N_mis[j], j], j] = y_mis[1:N_mis[j], j];
}
}

model {
// Separate sets of phis/thetas (no pool) - same rho, sigma (diff betas, phis, thetas)
for (j in 1:J)
{
  y[,j] ~ normal(mu[,j], sigma_y[,j]);
  target += -0.5 * dot_self(phi[node1, j] - phi[node2, j]);
  sum(phi[,j]) ~ normal(0, 0.001 * N);
  theta[,j] ~ normal(0, 1);
  beta0[j] ~ normal(120, 10);
  sigma_raw[j] ~ normal(0, 1); // implies sigma[j] ~ normal(mu_sigma, 3)
}

rho ~ beta(0.5, 0.5);
mu_sigma_raw ~ normal(0, 1); // implies mu_sigma ~ halfnormal(0, 3)
}'



# Run model ---------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.95
TREE_DEPTH <- 18
STEP_SIZE <- 0.001

tt <- proc.time()
fit <- stan(model_code = IAR_bym2,
            data = DATA,
            chains = 4,
            iter = 8000,
            cores = 4,
            pars = c('sigma', 'mu_sigma', 
                     'rho', 'beta0', 'theta', 'phi', 'mu'),
            control = list(max_treedepth = TREE_DEPTH, adapt_delta = DELTA, stepsize = STEP_SIZE)) # modified control parameters based on warnings
run_time <- (proc.time()[3] - tt[3]) / 60



# Calc diagnostics and rerun if needed ------------------------------------

num_diverge <- get_num_divergent(fit)
num_tree <- sum(get_max_treedepth_iterations(fit))
num_BFMI <- length(get_low_bfmi_chains(fit))


#rerun model if things didn't go well
while (sum(c(num_diverge, num_tree, num_BFMI)) > 0 & DELTA <= 0.99)
{
  DELTA <- DELTA + 0.1
  TREE_DEPTH <- TREE_DEPTH + 1
  STEP_SIZE <- STEP_SIZE * 0.75

  tt <- proc.time()
  fit <- stan(model_code = IAR_bym2,
              data = DATA,
              chains = 4,
              iter = 8000,
              cores = 4,
              pars = c('sigma', 'mu_sigma', 
                       'rho', 'beta0', 'theta', 'phi', 'mu'),
              control = list(max_treedepth = TREE_DEPTH, adapt_delta = DELTA, stepsize = STEP_SIZE)) # modified control parameters based on warnings
  run_time <- (proc.time()[3] - tt[3]) / 60
  
  num_diverge <- get_num_divergent(fit)
  num_tree <- sum(get_max_treedepth_iterations(fit))
  num_BFMI <- length(get_low_bfmi_chains(fit))
}






# Checks -------------------------------------------------------------

# MCMCtrace(fit)
# MCMCsummary(fit, params = c('sigma', 'rho', 'beta0', 'mu_sigma'), n.eff = TRUE)
# MCMCsummary(fit, params = c('theta', 'phi'), n.eff = TRUE)

# print(fit, pars = c('sigma', 'rho'))

# #shiny stan
# library(shinystan)
# launch_shinystan(fit)
# fit <- readRDS('IAR_stan_Bombycilla_cedrorum-2019-01-16.rds')


# write model results to file ---------------------------------------------

  
#save to RDS
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_out_dir))
saveRDS(fit, file = paste0('IAR_stan_', args, '-', IAR_out_date, '.rds'))


options(max.print = 50000)
sink(paste0('IAR_results_', args, '.txt'))
cat(paste0('IAR results ', args, ' \n'))
cat(paste0('Total minutes: ', round(run_time, digits = 2), ' \n'))
cat(paste0('Adapt delta: ', DELTA, ' \n'))
cat(paste0('Max tree depth: ', TREE_DEPTH, ' \n'))
cat(paste0('Step size: ', STEP_SIZE, ' \n'))
cat(paste0('Number of divergences: ', num_diverge, ' \n'))
cat(paste0('Number of tree exceeds: ', num_tree, ' \n'))
cat(paste0('Number chains low BFMI: ', num_BFMI, ' \n'))
cat(paste0('Cell drop: ', DROP, ' \n'))
print(fit)
sink()




# Plot pre-IAR/post_IAR halfmax estimates ------------------------------------------

#estimated half-max in grey, sd in white (derived from logit cubic)

#extract median and sd estimates for mu params
med_fit <- MCMCpstr(fit, params = 'mu', func = median)[[1]]
sd_fit <- MCMCpstr(fit, params = 'mu', func = sd)[[1]]


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
  #i <- 4
  
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
    geom_polygon(data = to_plt2, aes(x = long, y = lat, group = group, fill = HM_mean), 
                 alpha = 0.4) +
    geom_path(data = to_plt2, aes(x = long, y = lat, group = group), 
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
         filename = paste0(f_out_filt$species[1], '_', f_out_filt$year[1], '_pre_IAR.pdf'))
  
  
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
         filename = paste0(f_out_filt$species[1], '_', f_out_filt$year[1], '_post_IAR.pdf'))
}



# Trace plots with PPO ----------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_out_dir))


#beta0[j] ~ normal(120, 10)
#rho ~ beta(0.5, 0.5);
#mu_sigma ~ normal(0, 3);
#sigma_sigma ~ uniform(0, 5);

#beta0
PR <- rnorm(10000, 120, 10)
MCMCtrace(fit, 
          params = 'beta0',
          priors = PR,
          open_pdf = FALSE,
          filename = paste0('trace_beta0_', args, '-', IAR_out_date, '.pdf'))

#rho
PR <- rbeta(10000, 0.5, 0.5)
MCMCtrace(fit, 
          params = 'rho',
          priors = PR,
          open_pdf = FALSE,
          filename = paste0('trace_rho_', args, '-', IAR_out_date, '.pdf'))

#mu_sigma
PR <- rnorm(10000, 0, 3)
MCMCtrace(fit, 
          params = 'mu_sigma',
          priors = PR,
          open_pdf = FALSE,
          filename = paste0('trace_mu_sigma_', args, '-', IAR_out_date, '.pdf'))


if ('Rplots.pdf' %in% list.files())
{
  file.remove('Rplots.pdf')
}




print('I completed!')

