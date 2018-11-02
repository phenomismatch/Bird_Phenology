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

diagnostics_frame_p <- readRDS(paste0('diagnostics_frame-', IAR_date,'.rds'))
cells_frame <- readRDS(paste0('cells_frame-', IAR_date, '.rds'))
yrs_frame <- readRDS(paste0('yrs_frame-', IAR_date, '.rds'))


#which species/year has the most cells - to model 'data rich species'
DR_sp <- as.character(yrs_frame[which.max(yrs_frame[,1:3]$n_cells),1])
DR_filt <- dplyr::filter(yrs_frame, species == DR_sp)[,1:3]
DR_yr <- c(2002, 2010, 2017)


# filter cells ------------------------------------------------------------

#cells identified manually from plot (plot_code.R) to filter out
cells_to_rm <- c(208, 845, 289, 316, 343, 345, 346, 319, 320, 373, 
                 401, 402, 375, 376, 402, 403, 404, 370, 770, 796, 
                 766, 765, 763, 735, 734, 731, 758, 813, 3669, 397,
                 398)

'%ni%' <- Negate('%in%') 
diagnostics_frame <- diagnostics_frame_p[which(diagnostics_frame_p$cell %ni% cells_to_rm),]




# Species list to model ----------------------------------------------------


NC <- 3

#species list to model (greater than or = to 'NC' cells for 2015:2017)
m_species_list <- c()
#for (i in 1:nsp)
#{
i <- 80
t_yrsf <- dplyr::filter(yrs_frame, species == species_list[i])

if (sum(t_yrsf$n_cells[which(t_yrsf$year %in% 2015:2017)] >= NC) == 3)
{
  m_species_list <- c(m_species_list, species_list[i])
}
#}



# create adjacency matrix -------------------------------------------------

#unique cells
cells <- unique(diagnostics_frame$cell)
ncel <- length(cells)

#get hexgrid cell centers
hexgrid6 <- dggridR::dgconstruct(res = 6)
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



# Filter data by species ------------------------------------------------------

#filter by species/year here
f_out <- filter(diagnostics_frame, species == DR_sp, year %in% DR_yr)


#explore filter for good data
sum(f_out$min_n.eff < 500, na.rm = TRUE)
sum(f_out$max_Rhat > 1.1, na.rm = TRUE)
sum(f_out$nphen_bad > 100, na.rm = TRUE)



# if (m_crit == TRUE)
# {
#   #proceed
# } else {
#   #stop
# }

#check to make sure cells order == order in f_out
all.equal(unique(f_out$cell), cells)



# Estimate scaling factor for BYM2 model with INLA ------------------------

#taken from bym2 example here: http://mc-stan.org/users/documentation/case-studies/icar_stan.html
library(INLA)

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


#create data list for Stan
DATA <- list(N = ncel, 
             N_obs = sum(!is.na(f_out$HM_mean)),
             N_mis = sum(is.na(f_out$HM_mean)),
             N_edges = nrow(ninds), 
             node1 = ninds[,1], 
             node2 = ninds[,2],
             y_obs = f_out$HM_mean[which(!is.na(f_out$HM_mean))], 
             sigma_y = f_out$HM_sd,
             ii_obs = which(!is.na(f_out$HM_mean)),
             ii_mis = which(is.na(f_out$HM_mean)),
             scaling_factor = scaling_factor)


DATA$sigma_y[which(is.na(DATA$sigma_y))] <- 0.01




# Stan model --------------------------------------------------------------

#spatial and non-spatial component (reparameterized) - also model obs error
#uses scaling factor fomr INLA to set priors more easily
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

for (j in 1:J)
{
  real<lower = 0, upper = 200> y_obs[N_obs[j], J];            // observed response data (excluding NAs)
}

int<lower = 1, upper = N_obs + N_mis> ii_obs[N_obs];  // indices of observed values
int<lower = 1, upper = N_obs + N_mis> ii_mis[N_mis];  // indices of unobserved values

real<lower = 0> sigma_y[N, J];                           // observed sd of data (observation error)

real<lower = 0> scaling_factor;                       // scales variances of spatial effects
}

parameters {
real<lower = 1, upper = 200> y_mis[N_mis];            // missing response data

real beta0[J];                                           // intercept

real<lower = 0> sigma[J];                             // spatial and non-spatial sd
real<lower = 0, upper = 1> rho;                       // proportion unstructure vs spatially structured variance

vector[N] phi;                                        // spatial error component (centered on 0)
vector[N] theta;                                      // non-spatial error component (centered on 0)
}

transformed parameters {
real<lower = 0, upper = 200> y[N];

vector[N] convolved_re;
matrix[N, J] mu;                                         // latent true halfmax values

// variance of each component should be approx equal to 1
convolved_re = sqrt(1 - rho) * theta + sqrt(rho / scaling_factor) * phi;

for (j in 1:J)
{
  mu[,j] = beta0[j] + convolved_re * sigma[j];                    // scaling by sigma_phi rather than phi ~ N(0, sigma_phi)

  y[ii_obs, J] = y_obs;                                    // transform to accomodate missing data in Stan
  y[ii_mis, J] = y_mis;                                    // transform to accomodate missing data in Stan
}


}


model {

for (j in 1:J)
{
  y[,j] ~ normal(mu[,j], sigma_y[,j]);
  
  // the following computes the prior on phi on the unit scale with sd = 1
  target += -0.5 * dot_self(phi[node1] - phi[node2]);
  
  // soft sum-to-zero constraint on phi) - equivalent to mean(phi) ~ normal(0,0.001)
  sum(phi) ~ normal(0, 0.001 * N);
  
  beta0[j] ~ normal(120, 10);
  theta ~ normal(0, 1);
  sigma[j] ~ normal(0, sigma_sigma[j]);
  rho ~ beta(0.5, 0.5);
}

// previously used sigma ~ N(0, 5)
sigma_sigma ~ uniform(0, 5)

}

}'



# Run model ---------------------------------------------------------------

# Reference Ver Hoef et al. 2018 Eco Monographs for nice overview of spatial autoregressive models
# Model structure to accomodate missing data derived from section 11.3 of Stan reference manual 2.17.0 (p. 182)
# No need to worry about Jacobian warning - the 'transformation' is linear (just using an indexing trick to accomodate missing data in Y)
# No need to worry about Parser warning about 'Unknown variable: sum' - known error in stan (as of Oct 26, 2018)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

tt <- proc.time()
fit <- stan(model_code = IAR_bym2,                 # Model
            data = DATA,                           # Data
            chains = 5,                            # Number chains
            iter = 5000,                           # Iterations per chain
            cores = 5,                             # Number cores to use
            control = list(max_treedepth = 18, adapt_delta = 0.90)) # modified control parameters based on warnings;
# see http://mc-stan.org/misc/warnings.html
proc.time() - tt


#save to RDS
# saveRDS(fit, 'stan_bym2_fit.rds')


#diagnostics
# pairs(fit, pars = c('sigma', 'rho', 'beta0'))


#https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html
sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
max_treedepth_by_chain <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
get_elapsed_time(fit)



#shiny stan
library(shinystan)
launch_shinystan(fit)



# put copy of script in dir -----------------------------------------------

system(paste0('cp ', dir, 'Bird_Phenology/Scripts/ebird_Nestwatch/4-IAR-model.R ', IAR_dir, '/4-ICAR-model-', Sys.Date(), '.R'))



# Plot pre-IAR halfmax estimates ------------------------------------------


#transform cells to grid
cell_grid <- dggridR::dgcellstogrid(hexgrid6, cells)
cell_grid$cell <- as.numeric(cell_grid$cell)
cell_centers <- dggridR::dgSEQNUM_to_GEO(hexgrid6, cells)
ll_df <- data.frame(cell = cells, lon_deg = cell_centers$lon_deg, lat_deg = cell_centers$lat_deg)

#merge hex spatial data with HM data
to_plt <- dplyr::inner_join(f_out, cell_grid, by = 'cell')
to_plt2 <- dplyr::inner_join(to_plt, ll_df, by = 'cell')

#load maps
usamap <- data.frame(maps::map("world", "USA", plot = FALSE)[c("x", "y")])
canadamap <- data.frame(maps::map("world", "Canada", plot = FALSE)[c("x", "y")])
mexicomap <- data.frame(maps::map("world", "Mexico", plot = FALSE)[c("x", "y")])

#plot
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
                       limits = c(75, 180)) +
  labs(fill = 'Estimated Arrival') +
  annotate('text', x = to_plt2$lon_deg, y = to_plt2$lat_deg + 0.5, 
           label = round(to_plt2$HM_mean, digits = 0), col = 'black', alpha = 0.1,
           size = 4) +
  annotate('text', x = to_plt2$lon_deg, y = to_plt2$lat_deg - 0.5, 
           label = round(to_plt2$HM_sd, digits = 0), col = 'white', alpha = 0.3,
           size = 3) +
  ggtitle(paste0(f_out$species[1], ' - ', f_out$year[1], ' - Pre-IAR')) +
  theme_bw() +
  xlab('Longitude') +
  ylab('Latitude')


#estimated half-max in grey, sd in white (derived from logit cubic)
p



# Plot post-IAR halfmax estimates -----------------------------------------

med_fit <- MCMCpstr(fit, params = 'mu', func = median)[[1]]
sd_fit <- MCMCpstr(fit, params = 'mu', func = sd)[[1]]

m_fit <- data.frame(mean = med_fit, sd = sd_fit, cell = cells)

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
  geom_polygon(data = to_plt2_post, aes(x = long, y = lat, group = group, fill = mean), 
               alpha = 0.4) +
  geom_path(data = to_plt2_post, aes(x = long, y = lat, group = group), 
            alpha = 0.4, color = 'black') + 
  scale_fill_gradientn(colors = c('red', 'blue'),
                       limits = c(75, 180)) +
  labs(fill = 'Estimated Arrival') +
  annotate('text', x = to_plt2_post$lon_deg, y = to_plt2_post$lat_deg + 0.5, 
           label = round(to_plt2_post$mean, digits = 0), col = 'black', alpha = 0.1,
           size = 4) +
  annotate('text', x = to_plt2_post$lon_deg, y = to_plt2_post$lat_deg - 0.5, 
           label = round(to_plt2_post$sd, digits = 0), col = 'white', alpha = 0.3,
           size = 3) +
  ggtitle(paste0(f_out$species[1], ' - ', f_out$year[1], ' - Post-IAR')) +
  theme_bw() +
  xlab('Longitude') +
  ylab('Latitude')


#estimated half-max in grey, sd in white (derived from logit cubic)
p_post
