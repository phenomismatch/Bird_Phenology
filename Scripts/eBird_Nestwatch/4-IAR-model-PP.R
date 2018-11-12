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


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/home/CAM/cyoungflesh/phenomismatch/'



# db/hm query dir ------------------------------------------------------------

db_dir <- 'db_query_2018-10-15'
hm_dir <- 'halfmax_species_2018-10-16'
IAR_dir <- 'IAR_2018-10-26'


# runtime -----------------------------------------------------------------

tt <- proc.time()



# Load packages -----------------------------------------------------------

library(rstan)
library(dplyr)
library(dggridR)
library(geosphere)
library(ggplot2)
library(MCMCvis)
library(maps)
library(INLA)


# Set wd ------------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/'))



# import eBird species list -----------------------------------------------------

species_list_i <- read.table('eBird_species_list.txt', stringsAsFactors = FALSE)
species_list <- species_list_i[,1]
nsp <- length(species_list)


# read in data files ------------------------------------------------------

#DATA ONLY VALID THROUGH 2017 (2018 data only goes to ~ jday 60 as of 2018-10-15 query)

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_dir))

IAR_date <- substr(IAR_dir, start = 5, stop = 15)

diagnostics_frame <- readRDS(paste0('diagnostics_frame-', IAR_date,'.rds'))
cells_frame <- readRDS(paste0('cells_frame-', IAR_date, '.rds'))
yrs_frame <- readRDS(paste0('yrs_frame-', IAR_date, '.rds'))


#which species/year has the most cells - to model 'data rich species'
# DR_sp <- as.character(yrs_frame[which.max(yrs_frame[,1:3]$n_cells),1])
# DR_filt <- dplyr::filter(yrs_frame, species == DR_sp)[,1:3]
# DR_yr <- 2002:2017

aggregate(n_cells ~ species, data = yrs_frame, FUN = max)
aggregate(n_cells ~ species, data = yrs_frame, FUN = mean)

#species with very little data
DR_sp <- 'Ammodramus_nelsoni'
DR_yr <- 2015:2017
nyr <- length(DR_yr)

# Filter data by species/years ------------------------------------------------------

#MAKE SURE THERE ARE DATA FOR THAT YEAR (see script 3-....R) - can't accomodate no values at all for a given year


#filter by species/year here
df_filt <- filter(diagnostics_frame, species == DR_sp, year %in% DR_yr)


#explore filter for good data
sum(df_filt$min_n.eff < 500, na.rm = TRUE)
sum(df_filt$max_Rhat > 1.1, na.rm = TRUE)
sum(df_filt$nphen_bad > 100, na.rm = TRUE)



# filter cells ------------------------------------------------------------

# #cells ID'd as a 'test set'
# cells_to_keep <- c(622, 595, 567, 594, 621, 648, 649, 676, 675, 647, 619, 
#                    620, 621, 593, 592, 619, 618, 591, 564, 565, 566, 646, 
#                    674, 702, 703)
# 
# diagnostics_frame <- diagnostics_frame_p[which(diagnostics_frame_p$cell %in% cells_to_keep),]


#cells chosen based on species range - if cell overlaps species range -> include it
#create hex cell grid
hexgrid6 <- dggridR::dgconstruct(res = 6)

#get boundaries of all cells over earth - ADD TO PREVIOUS SCRIPT
setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/shapefiles/'))
dggridR::dgearthgrid(hexgrid6, savegrid = 'global_hex.shp')

#load species range map
#sp_rng <- rgdal::readOGR('Vireo_olivaceus_22705243.shp', verbose = FALSE)
sp_rng <- rgdal::readOGR('Ammodramus_nelsoni_22728393.shp', verbose = FALSE)

#filter by breeding (2) and migration (4) range - need to convert spdf to sp
nrng <- sp_rng[which(sp_rng$SEASONAL == 2 | sp_rng$SEASONAL == 4),]
nrng_sp <- sp::SpatialPolygons(nrng@polygons)
sp::proj4string(nrng_sp) <- sp::CRS(sp::proj4string(nrng))

hge <- rgdal::readOGR('global_hex.shp', verbose = FALSE)
ptsreg <- sp::spsample(nrng, 50000, type = "regular")
ovrlp <- as.numeric(which(!is.na(sp::over(hge, ptsreg))))
overlap_cells <- ovrlp

#get cell centers
cell_centers <- dggridR::dgSEQNUM_to_GEO(hexgrid6, overlap_cells)
cc_df <- data.frame(cell = overlap_cells, lon = cell_centers$lon_deg, 
                    lat = cell_centers$lat_deg)

#cells only within the range that ebird surveys were filtered to
n_cc_df <- cc_df[which(cc_df$lon > -100 & cc_df$lon < -50 & cc_df$lat > 26),]
cells <- n_cc_df$cell

#retain rows that match selected cells
df_filt2 <- df_filt[which(df_filt$cell %in% cells),]

'%ni%' <- Negate('%in%')

#create rows for cells that were missing in ebird data
missing_cells <- cells[which(cells %ni% df_filt2$cell)]

temp_dff <- df_filt2[1,]
temp_dff[,2:20] <- NA

nmc <- length(missing_cells)
nyrs <- length(DR_yr)
nreps <- nmc * nyrs

temp_dff2 <- temp_dff[rep(row.names(temp_dff), nreps),]
rownames(temp_dff2) <- NULL

temp_dff2$year <- rep(DR_yr, nmc)
temp_dff2$cell <- rep(missing_cells, each = nyrs)


#combine filtered data with missing cells
f_out <- rbind(df_filt2, temp_dff2)




# create adjacency matrix -------------------------------------------------

#unique cells
ncel <- length(cells)

#get hexgrid cell centers
cellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid6, cells)

#create adjacency matrix - 1 if adjacent to cell, 0 if not
adjacency_matrix <- matrix(data = NA, nrow = length(cells), ncol = length(cells))

for (i in 1:length(cells))
{
  #i <- 1
  for (j in i:length(cells))
  {
    #j <- 4
    dists <- geosphere::distm(c(cellcenters$lon_deg[i], cellcenters$lat_deg[i]),
                              c(cellcenters$lon_deg[j], cellcenters$lat_deg[j]))
    adjacency_matrix[i,j] <- as.numeric((dists/1000) > 0 & (dists/1000) < 311)
  }
}

#indices for 1s
ninds <- which(adjacency_matrix == 1, arr.ind = TRUE)




# Estimate scaling factor for BYM2 model with INLA ------------------------

#Build the adjacency matrix using INLA library functions
adj.matrix <- sparseMatrix(i = ninds[,1], j = ninds[,2], 
                           x = 1, symmetric = TRUE)

#The ICAR precision matrix (note! This is singular)
Q <- Diagonal(ncel, rowSums(adj.matrix)) - adj.matrix
#Add a small jitter to the diagonal for numerical stability (optional but recommended)
Q_pert <- Q + Diagonal(ncel) * max(diag(Q)) * sqrt(.Machine$double.eps)

# Compute the diagonal elements of the covariance matrix subject to the 
# constraint that the entries of the ICAR sum to zero.
# See the inla.qinv function help for further details.
Q_inv <- INLA::inla.qinv(Q_pert, 
                         constr = list(A = matrix(1, 1, ncel), e = 0))

#Compute the geometric mean of the variances, which are on the diagonal of Q.inv
scaling_factor <- exp(mean(log(diag(Q_inv))))




# create Stan data object -------------------------------------------------

#create and fill sds and obs
sigma_y_in <- matrix(nrow = ncel, ncol = nyr)
y_obs_in <- matrix(nrow = ncel, ncol = nyr)

#number of observation and NAs for each year
len_y_obs_in <- rep(NA, nyr)
len_y_mis_in <- rep(NA, nyr)

#indices for observed and missing
ii_obs_in <- matrix(NA, nrow = ncel, ncol = nyr)
ii_mis_in <- matrix(NA, nrow = ncel, ncol = nyr)

for (j in 1:nyr)
{
  #j <- 16
  temp_yr_p <- filter(f_out, year == DR_yr[j])
  temp_yr <- temp_yr_p[order(temp_yr_p$cell),]
  
  sigma_y_in[,j] <- temp_yr$HM_sd
  
  #which are not NA
  no_na <- temp_yr$HM_mean[which(!is.na(temp_yr$HM_mean))]
  
  #pad end with NAs
  if (length(no_na) < ncel)
  {
    num_na <- ncel - length(no_na)
    
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
  len_y_mis_in[j] <- ncel - length(no_na)
}


#fill 0 where NA in y_obs - Stan does not like NA and zeros are not being used to estimate any param (y_obs is used to fill y)
y_obs_in[which(is.na(y_obs_in), arr.ind = TRUE)] <- 0
sigma_y_in[which(is.na(sigma_y_in), arr.ind = TRUE)] <- 0.1
ii_obs_in[which(is.na(ii_obs_in), arr.ind = TRUE)] <- 0
ii_mis_in[which(is.na(ii_mis_in), arr.ind = TRUE)] <- 0



#create data list for Stan
DATA <- list(J = nyr,
             N = ncel, 
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
real<lower = 0> sigma;                                // scaling factor for spatial and non-spatial components
real<lower = 0, upper = 1> rho;                       // proportion unstructured vs spatially structured variance
}

transformed parameters {
real<lower = 0, upper = 200> y[N, J];                 // response data to be modeled
matrix[N, J] convolved_re;                            // spatial and non-spatial component
matrix[N, J] mu;                                      // latent true halfmax values
for (j in 1:J)
{
  convolved_re[,j] = sqrt(1 - rho) * theta[,j] + sqrt(rho / scaling_factor) * phi[,j];
  mu[,j] = beta0[j] + convolved_re[,j] * sigma;
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
}

rho ~ beta(0.5, 0.5);
sigma ~ normal(0, 5);
}'



# Run model ---------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

tt <- proc.time()
fit <- stan(model_code = IAR_bym2,
            data = DATA,
            chains = 3,
            iter = 1000,
            cores = 3,
            pars = c('sigma', 'rho', 'beta0', 'theta', 'phi', 'mu'),
            control = list(max_treedepth = 20, adapt_delta = 0.90, stepsize = 0.01)) # modified control parameters based on warnings
proc.time() - tt


#save to RDS
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_dir))
saveRDS(fit, 'stan_bym2_allyr_allcells_sep_phis_1000_A_nelsoni.rds')
# fit <- readRDS('stan_bym2_allyr_allcells_sep_phis_500.rds')

#diagnostics
# pairs(fit, pars = c('sigma', 'rho'))

sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
max_treedepth_by_chain <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
get_elapsed_time(fit)

MCMCtrace(fit)
MCMCsummary(fit, params = c('sigma', 'rho', 'beta0'), n.eff = TRUE)
MCMCsummary(fit, params = c('theta', 'phi'), n.eff = TRUE)

print(fit, pars = c('sigma', 'rho'))

#shiny stan
library(shinystan)
launch_shinystan(fit)



# put copy of script in dir -----------------------------------------------

system(paste0('cp ', dir, 'Bird_Phenology/Scripts/ebird_Nestwatch/4-IAR-model.R ', 
              IAR_dir, '/4-ICAR-model-', Sys.Date(), '.R'))



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

#plotting species range
nrng@data$id <- rownames(nrng@data)
nrng.points <-fortify(nrng, region="id")
nrng.df <- plyr::join(nrng.points, nrng@data, by="id")

setwd(paste0(dir, 'Bird_Phenology/Figures/pre_post_IAR_maps'))


#loop plots for each year
for (i in 1:length(DR_yr))
{
  #i <- 1
  
  #filter data for year[i]
  f_out_filt <- filter(f_out, year == DR_yr[i])
  
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
    geom_polygon(data = nrng.df, 
              aes(x = long, y = lat, group=group), fill = 'green', alpha = 0.5) + 
    coord_map("ortho", orientation = c(35, -80, 0), 
              xlim = c(-100, -55), ylim = c(20, 60)) + 
    geom_polygon(data = to_plt2, aes(x = long, y = lat, group = group, fill = HM_mean), 
                 alpha = 0.4) +
    geom_path(data = to_plt2, aes(x = long, y = lat, group = group), 
              alpha = 0.4, color = 'black') + 
    scale_fill_gradientn(colors = c('red', 'blue'),
                         limits = c(MIN, MAX)) +
    labs(fill = 'Estimated Arrival') +
    annotate('text', x = to_plt2$lon_deg, y = to_plt2$lat_deg + 0.5, 
             label = round(to_plt2$HM_mean, digits = 0), col = 'black', alpha = 0.3,
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
    geom_polygon(data = nrng.df, 
                 aes(x = long, y = lat, group=group), fill = 'green', alpha = 0.5) + 
    coord_map("ortho", orientation = c(35, -80, 0), 
              xlim = c(-100, -55), ylim = c(20, 60)) + 
    geom_polygon(data = to_plt2_post, aes(x = long, y = lat, group = group, fill = med_mu), 
                 alpha = 0.4) +
    geom_path(data = to_plt2_post, aes(x = long, y = lat, group = group), 
              alpha = 0.4, color = 'black') + 
    scale_fill_gradientn(colors = c('red', 'blue'),
                         limits = c(MIN, MAX)) +
    labs(fill = 'Estimated Arrival') +
    annotate('text', x = to_plt2_post$lon_deg, y = to_plt2_post$lat_deg + 0.5, 
             label = round(to_plt2_post$med_mu, digits = 0), col = 'black', alpha = 0.1,
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
