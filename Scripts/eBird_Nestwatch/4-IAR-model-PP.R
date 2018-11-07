######################
# 4 - ICA model PARTIAL POOLING
######################

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




#DATA ONLY VALID THROUGH 2017 (2018 data only goes to ~ jday 60 as of 2018-10-15 query)
years <- 2002:2017
nyr <- length(years)




# read in data files ------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_dir))

IAR_date <- substr(IAR_dir, start = 5, stop = 15)

diagnostics_frame <- readRDS(paste0('diagnostics_frame-', IAR_date,'.rds'))
cells_frame <- readRDS(paste0('cells_frame-', IAR_date, '.rds'))
yrs_frame <- readRDS(paste0('yrs_frame-', IAR_date, '.rds'))


#which species/year has the most cells - to model 'data rich species'
DR_sp <- as.character(yrs_frame[which.max(yrs_frame[,1:3]$n_cells),1])
DR_filt <- dplyr::filter(yrs_frame, species == DR_sp)[,1:3]
DR_yr <- 2002:2017




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

#get boundaries of all cells over earth
hge <- dggridR::dgearthgrid(hexgrid6)

#load species range map
setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/shapefiles/'))
sp_rng <- rgdal::readOGR('Vireo_olivaceus_22705243.shp', verbose = FALSE)

#filter by breeding range - need to convert spdf to sp
nrng <- subset(sp_rng, SEASONAL == 2)
nrng_sp <- sp::SpatialPolygons(nrng@polygons)
sp::proj4string(nrng_sp) <- sp::CRS(sp::proj4string(nrng))

#convert hex cell vetices to spatial points
ecells_sp <- sp::SpatialPoints(cbind(hge$long, hge$lat), 
                               proj4string = sp::CRS(sp::proj4string(nrng)))

#which vertices overlap species breeding range
ovrlp <- which(!is.na(sp::over(ecells_sp, nrng_sp)))
overlap_cells <- unique(as.numeric(hge[ovrlp,]$cell))

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

#taken from bym2 example here: http://mc-stan.org/users/documentation/case-studies/icar_stan.html

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
  temp_yr <- filter(f_out, year == DR_yr[j])
  
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
  
  #assign(paste0('ii_obs', j, '_in'), which(!is.na(temp_yr$HM_mean)))
  #assign(paste0('ii_mis', j, '_in'), which(is.na(temp_yr$HM_mean)))
}


#fill 0 where NA in y_obs - Stan does not like NA and zeros are not being used to estimate any param (y_obs is used to fill y)
y_obs_in[which(is.na(y_obs_in), arr.ind = TRUE)] <- 0
sigma_y_in[which(is.na(sigma_y_in), arr.ind = TRUE)] <- 0
ii_obs_in[which(is.na(ii_obs_in), arr.ind = TRUE)] <- 0
ii_mis_in[which(is.na(ii_mis_in), arr.ind = TRUE)] <- 0



#create data list for Stan
DATA <- list(J = length(unique(f_out$year)),
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

# #add observation indices to DATA list
# i <- 1
# while (i <= length(DR_yr))
# {
#   DATA[[paste0('ii_obs', i)]] <- get(paste0('ii_obs', i, '_in'))
#   DATA[[paste0('ii_mis', i)]] <- get(paste0('ii_mis', i, '_in'))
#   
#   i <- i + 1
# }




# Stan model --------------------------------------------------------------

#spatial and non-spatial component (reparameterized) - also model obs error
#uses scaling factor fomr INLA to set priors more easily
#matrix[N,J] -> cells in rows, years in columns


IAR_bym2_PP <- '
data {
int<lower = 0> J;                                     // number of years      
int<lower = 0> N;                                     // number of cells

int<lower = 0> N_obs[J];                              // number of non-missing for each year
int<lower = 0> N_mis[J];                              // number missing for each year

int<lower = 0> N_edges;                               // number of edges in adjacency matrix
int<lower = 1, upper = N> node1[N_edges];             // node1[i] adjacent to node2[i]
int<lower = 1, upper = N> node2[N_edges];             // and node1[i] < node2[i]

real<lower = 0, upper = 200> y_obs[N, J];            // observed response data (add NAs to end)
real<lower = 0> sigma_y[N, J];                           // observed sd of data (observation error)
int<lower = 0> ii_obs[N, J];
int<lower = 0> ii_mis[N, J];

real<lower = 0> scaling_factor;                       // scales variances of spatial effects

int<lower = 1, upper = N> ii_obs1[N_obs[1]];
int<lower = 1, upper = N> ii_mis1[N_mis[1]];
int<lower = 1, upper = N> ii_obs2[N_obs[2]];
int<lower = 1, upper = N> ii_mis2[N_mis[2]];
int<lower = 1, upper = N> ii_obs3[N_obs[3]];
int<lower = 1, upper = N> ii_mis3[N_mis[3]];
int<lower = 1, upper = N> ii_obs4[N_obs[4]];
int<lower = 1, upper = N> ii_mis4[N_mis[4]];
int<lower = 1, upper = N> ii_obs5[N_obs[5]];
int<lower = 1, upper = N> ii_mis5[N_mis[5]];
int<lower = 1, upper = N> ii_obs6[N_obs[6]];
int<lower = 1, upper = N> ii_mis6[N_mis[6]];
int<lower = 1, upper = N> ii_obs7[N_obs[7]];
int<lower = 1, upper = N> ii_mis7[N_mis[7]];
int<lower = 1, upper = N> ii_obs8[N_obs[8]];
int<lower = 1, upper = N> ii_mis8[N_mis[8]];
int<lower = 1, upper = N> ii_obs9[N_obs[9]];
int<lower = 1, upper = N> ii_mis9[N_mis[9]];
int<lower = 1, upper = N> ii_obs10[N_obs[10]];
int<lower = 1, upper = N> ii_mis10[N_mis[10]];
int<lower = 1, upper = N> ii_obs11[N_obs[11]];
int<lower = 1, upper = N> ii_mis11[N_mis[12]];
int<lower = 1, upper = N> ii_obs12[N_obs[12]];
int<lower = 1, upper = N> ii_mis12[N_mis[13]];
int<lower = 1, upper = N> ii_obs13[N_obs[13]];
int<lower = 1, upper = N> ii_mis13[N_mis[14]];
int<lower = 1, upper = N> ii_obs14[N_obs[14]];
int<lower = 1, upper = N> ii_mis14[N_mis[15]];
int<lower = 1, upper = N> ii_obs15[N_obs[15]];
int<lower = 1, upper = N> ii_mis15[N_mis[16]];
int<lower = 1, upper = N> ii_obs16[N_obs[16]];
int<lower = 1, upper = N> ii_mis16[N_mis[16]];
}



parameters {
real<lower = 1, upper = 200> y_mis[N, J];             // missing response data
real beta0[J];                                        // intercept
matrix[N, J] theta;                                   // non-spatial error component (centered on 0)
matrix[N, J] phi;
real<lower = 0> sigma;                                // spatial and non-spatial sd
real<lower = 0, upper = 1> rho;                       // proportion unstructure vs spatially structured variance

// matrix[N, J] phi_raw;                              // hierarchical
// vector[N] mu_phi;                                  //hierarchical
// vector[N] sigma_phi;                               //hierarchical
}

transformed parameters {
real<lower = 0, upper = 200> y[N, J];
matrix[N, J] convolved_re;
matrix[N, J] mu;                                      // latent true halfmax values

for (j in 1:J)
{
  // phi[,j] = mu_phi + phi_raw[,j] .* sigma_phi;      //hierarchical: .* for elementwise multiplication of two vectors
  convolved_re[,j] = sqrt(1 - rho) * theta[,j] + sqrt(rho / scaling_factor) * phi[,j];
  mu[,j] = beta0[j] + convolved_re[,j] * sigma;          // scaling by sigma_phi rather than phi ~ N(0, sigma_phi)
}


y[ii_obs1, 1] = y_obs[1:N_obs[1], 1];
y[ii_mis1, 1] = y_mis[1:N_mis[1], 1];
y[ii_obs2, 2] = y_obs[1:N_obs[2], 2];
y[ii_mis2, 2] = y_mis[1:N_mis[2], 2];
y[ii_obs3, 3] = y_obs[1:N_obs[3], 3];
y[ii_mis3, 3] = y_mis[1:N_mis[3], 3];
y[ii_obs4, 4] = y_obs[1:N_obs[4], 4];
y[ii_mis4, 4] = y_mis[1:N_mis[4], 4];
y[ii_obs5, 5] = y_obs[1:N_obs[5], 5];
y[ii_mis5, 5] = y_mis[1:N_mis[5], 5];
y[ii_obs6, 6] = y_obs[1:N_obs[6], 6];
y[ii_mis6, 6] = y_mis[1:N_mis[6], 6];
y[ii_obs7, 7] = y_obs[1:N_obs[7], 7];
y[ii_mis7, 7] = y_mis[1:N_mis[7], 7];

y[ii_obs8, 8] = y_obs[1:N_obs[8], 8];
y[ii_mis8, 8] = y_mis[1, 8];                        
y[ii_obs9, 9] = y_obs[1:N_obs[9], 9];
y[ii_mis9, 9] = y_mis[1, 9];                        
y[ii_obs10, 10] = y_obs[1:N_obs[10], 10];
y[ii_mis10, 10] = y_mis[1, 10];                        // had to change due to dimsion mismatch (integer != length 1 vector)

y[ii_obs11, 11] = y_obs[1:N_obs[11], 11];
y[ii_obs12, 12] = y_obs[1:N_obs[12], 12];
y[ii_obs13, 13] = y_obs[1:N_obs[13], 13];
y[ii_obs14, 14] = y_obs[1:N_obs[14], 14];
y[ii_obs15, 15] = y_obs[1:N_obs[15], 15];
y[ii_obs16, 16] = y_obs[1:N_obs[16], 16];

}


model {

// 1) one set of phis/thetas (complete pool)
// 2) separate sets of phis (no pool)
// 3) partial pool

// #1
// One set of phis/thetas (complete pool) - same rho, phis, thetas, sigma (diff betas)
/* 
for (j in 1:J)
{
  y[,j] ~ normal(mu[,j], sigma_y[,j]);
  beta0[j] ~ normal(120, 10);
}

target += -0.5 * dot_self(phi[node1] - phi[node2]);
sum(phi) ~ normal(0, 0.001 * N);
theta ~ normal(0, 1);
rho ~ beta(0.5, 0.5);
sigma ~ normal(0, 5);
*/


// #2
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


// #3
// Hierarchical phis (partial pool) - same rho, sigma (partial pool phis, diff betas, thetas)
/*
for (j in 1:J)
{
  y[,j] ~ normal(mu[,j], sigma_y[,j]);
  target += -0.5 * dot_self(phi[node1, j] - phi[node2, j]);
  phi_raw[,j] ~ normal(0, 1);                   // reparameterize non-centered to optimize
  beta0[j] ~ normal(120, 10);
  theta[j] ~ normal(0,1);
}

mu_phi ~ normal(0, 1);
sigma_phi ~ normal(0,1);
sum(mu_phi) ~ normal(0, 0.001 * N);
rho ~ beta(0.5, 0.5);
sigma ~ normal(0, 5);
*/

}'



# Run model ---------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

tt <- proc.time()
fit_PP <- stan(model_code = IAR_bym2_PP,              # Model
            data = DATA,                           # Data
            chains = 3,                            # Number chains
            iter = 4000,                           # Iterations per chain
            cores = 3,                             # Number cores to use
            control = list(max_treedepth = 20, adapt_delta = 0.80)) # modified control parameters based on warnings;
# see http://mc-stan.org/misc/warnings.html
proc.time() - tt


#save to RDS
# saveRDS(fit_PP, 'stan_bym2_PP_fit_allyr_sep.rds')
# fit_PP <- readRDS('stan_bym2_PP_fit_3yr_sep.rds')

#diagnostics
# pairs(fit_PP, pars = c('sigma', 'rho', 'beta0'))


#https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html
sampler_params <- get_sampler_params(fit_PP, inc_warmup = FALSE)
mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
max_treedepth_by_chain <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
get_elapsed_time(fit_PP)


MCMCsummary(fit_PP, params = c('sigma', 'rho', 'beta0'), n.eff = TRUE)
MCMCsummary(fit_PP, params = c('theta', 'phi'), n.eff = TRUE)

print(fit_PP, pars = c('sigma', 'rho'))

#shiny stan
library(shinystan)
launch_shinystan(fit_PP)



# put copy of script in dir -----------------------------------------------

system(paste0('cp ', dir, 'Bird_Phenology/Scripts/ebird_Nestwatch/4-IAR-model.R ', 
              IAR_dir, '/4-ICAR-model-', Sys.Date(), '.R'))



# Plot pre-IAR/post_IAR halfmax estimates ------------------------------------------

#estimated half-max in grey, sd in white (derived from logit cubic)

#extract median and sd estimates for mu params
med_fit <- MCMCpstr(fit_PP, params = 'mu', func = median)[[1]]
sd_fit <- MCMCpstr(fit_PP, params = 'mu', func = sd)[[1]]


#transform cells to grid
cell_grid <- dggridR::dgcellstogrid(hexgrid6, cells)
cell_grid$cell <- as.numeric(cell_grid$cell)
cell_centers <- dggridR::dgSEQNUM_to_GEO(hexgrid6, cells)
ll_df <- data.frame(cell = cells, lon_deg = cell_centers$lon_deg, lat_deg = cell_centers$lat_deg)

#load maps
usamap <- data.frame(maps::map("world", "USA", plot = FALSE)[c("x", "y")])
canadamap <- data.frame(maps::map("world", "Canada", plot = FALSE)[c("x", "y")])
mexicomap <- data.frame(maps::map("world", "Mexico", plot = FALSE)[c("x", "y")])

#determine min/max for plotting
f_rng <- c(range(f_out$HM_mean, na.rm = TRUE), range(med_fit, na.rm = TRUE))
MIN <- round(min(f_rng))
MAX <- round(max(f_rng))



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
             label = round(to_plt2$HM_mean, digits = 0), col = 'black', alpha = 0.1,
             size = 4) +
    annotate('text', x = to_plt2$lon_deg, y = to_plt2$lat_deg - 0.5, 
             label = round(to_plt2$HM_sd, digits = 0), col = 'white', alpha = 0.3,
             size = 3) +
    ggtitle(paste0(f_out_filt$species[1], ' - ', f_out_filt$year[1], ' - Pre-IAR')) +
    theme_bw() +
    xlab('Longitude') +
    ylab('Latitude')

  ggsave(plot = p, filename = paste0(f_out_filt$species[1], '_', f_out_filt$year[1], '_pre_IAR.pdf'))


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
  
  ggsave(plot = p_post, filename = paste0(f_out_filt$species[1], '_', f_out_filt$year[1], '_post_IAR.pdf'))
}
