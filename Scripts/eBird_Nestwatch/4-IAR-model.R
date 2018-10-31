######################
# 4 - ICAR model (only spatial component)
#
# Fit ICAR model - true arrival dates for each species-cell-year are modeled as latent states, with the observed 
# state derived from 2-logit-cubic.R. 
#
# Formerly ICAR_parallel.R
######################

# This script takes the output of NA_birdPhen1.R, which gives posterior distributions for the half-max
# parameter by species-cell-years, and uses these as the basis for an ICAR model of the half-max
# parameter.  The model assumes that the posterior distributions for the half-max parameter from 
# NA_birdPhen1.R are normally distributed. **I have spot-checked this assumption and it seems reasonable, 
# but have not exhaustively evaluated it**

# Then, the model assumes that the posterior means are drawn from normal distributions
# parameterized by a latent true half-max and the variance in the posterior estimate.

# The latent true half-maxima (LTHMs) are given a prior with a spatial component corresponding to 
# an intrinsic conditional autoregressive (ICAR) model. The model is fit across cells within a 
# species-year; there is no data sharing across years or species. The ICAR model includes only a spatial
# variance term for the LTHMs. Thus, the LTHM from a cell is drawn from a normal distribution centered 
# around the average of the values of the neighboring cells.

# Note than in an ICAR model with a nonspatial component, the LTHM from a cell would be drawn from
# a normal distribution centered around around some latent spatial mean, which itself is drawn from 
# a normal distribution centered on the average of the spatial means of the cell's neighbors. So 
# the difference between the spatial variance and the nonspatial variance is that the spatial 
# residual propagates to spatially to influence the spatial mean of the neighbors, whereas the 
# nonspatial residual does not do this. 

# In NA_birdPhen2.R (now `3b-ICAR-model-ns.R`), I fit models with a nonspatial variance term, but I found that the estimated
# nonspatial variance is small, difficult to estimate, and leads to diagnostic warnings from Stan.
# I experimented with strongly informative priors for that value, but here I simply use a model
# with no nonspatial component. 



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
DR_yr <- as.numeric(yrs_frame[which.max(yrs_frame[,1:3]$n_cells),2])


# filter cells ------------------------------------------------------------

#cells identified manually from plot (plot_code.R) to filter out
cells_to_rm <- c(208, 845, 289, 316, 343, 345, 346, 319, 320, 373, 
                 401, 402, 375, 376, 402, 403, 404, 370, 770, 796, 
                 766, 765, 763, 735, 734, 731, 758, 813, 3669, 397,
                 398)

'%ni%' <- Negate('%in%') 
diagnostics_frame <- diagnostics_frame_p[which(diagnostics_frame_p$cell %ni% cells_to_rm),]




#explore filter for good data
sum(diagnostics_frame$min_n.eff < 500, na.rm = TRUE)
sum(diagnostics_frame$max_Rhat > 1.1, na.rm = TRUE)
sum(diagnostics_frame$nphen_bad > 100, na.rm = TRUE)

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
#cells <- unique(diagnostics_frame$cell)
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
f_out <- filter(diagnostics_frame, species == DR_sp, year == DR_yr)


# if (m_crit == TRUE)
# {
#   #proceed
# } else {
#   #stop
# }

#check to make sure cells order == order in f_out
all.equal(f_out$cell, cells)



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


IAR_no_obs_error <- '
data {
int<lower = 0> N;                                     // number of cells (= number of observations, including NAs)
int<lower = 0> N_obs;                                 // number of non-missing
int<lower = 0> N_mis;                                 // number missing
int<lower = 0> N_edges;                               // number of edges in adjacency matrix
int<lower = 1, upper = N> node1[N_edges];             // node1[i] adjacent to node2[i]
int<lower = 1, upper = N> node2[N_edges];             // and node1[i] < node2[i]

real<lower = 0, upper = 200> y_obs[N_obs];            // observed response data (excluding NAs)

int<lower = 1, upper = N_obs + N_mis> ii_obs[N_obs];  // indices of observed values
int<lower = 1, upper = N_obs + N_mis> ii_mis[N_mis];  // indices of unobserved values

// real<lower = 0> sds_obs[N_obs];
// real<lower = 0> sds[N];
}

parameters {
real<lower = 1, upper = 200> y_mis[N_mis];            // missing response data
real<lower = 0> sigma_y;                              // non-spatial error component
real beta0;                                           // intercept
vector[N] phi;                                        // spatial error component (centered on 0)
real<lower = 0> sigma_phi;                            // to scale spatial error component

// real<lower = 0> sds_mis[N_mis];                    
}

transformed parameters {
real<lower = 0, upper = 200> y[N];

vector[N] mu;                                         // latent true halfmax values
mu = beta0 + phi * sigma_phi;                         // scaling by sigma_phi rather than phi ~ N(0, sigma_phi)

y[ii_obs] = y_obs;                                    // transform to accomodate missing data in Stan
y[ii_mis] = y_mis;                                    // transform to accomodate missing data in Stan

// real<lower = 0> sd_y[N];
// sds[ii_obs] = sds_obs;
// sds[ii_mis] = sds_mis;
}

  
model {
// the way this is coded, it assumes y is observed without error
// check notation - formulation from Herarchical modeling and analysis for spatial data Eq. 3.16, Ver Hoef et al. 2018, and Morris online Stan IAR
// http://mc-stan.org/users/documentation/case-studies/icar_stan.html

// y[i] sim N(mu[i], sigma_y)
// mu[i] = beta_0 + phi[i]
// phi sim N(0, [tau(D-W)]^{-1})
// which rewrites to the *pairwise difference* formulation:
// phi[i] sim N(0, sigma_{phi})
// sigma_{phi} = exp(-frac{1}{2} sum_{i sim j}{(phi[i] - phi[j])^{2}})

// y[i] is response data
// mu[i] is latent half max
// sigma_y is non-spatial error component
// phi[i] is spatial error component - must be centered on 0 to fit properly
// sigma_{phi} is the pairwise difference formulation 
// D is n x n diagonal matrix, where d[i,j] = # of neighbors for cell i
// W is weights matrix, w[i,i] = 0, w[i,j] = 1 if i is neighbor of j and w[i,j] = 0 otherwise


y ~ normal(mu, sigma_y);

// the following computes the prior on phi on the unit scale with sd = 1
target += -0.5 * dot_self(phi[node1] - phi[node2]);

// soft sum-to-zero constraint on phi) - equivalent to mean(phi) ~ normal(0,0.001)
sum(phi) ~ normal(0, 0.001 * N);
}'






#spatial and non-spatial component - also model obs error
#setting fair prior in this model is difficult

# 'fair prior' following Carlin and Perez 2000 p. 16
# K = 0.7
# n_bar = avg number of neighbors
# tau_h = K^2 * n_bar * tau_c
# 
# tau_h = 0.62^2 * 4.71 * tau_c
# 
# 
# tau_c ~ gamma(a_c, b_c) and tau_h ~ gamma(a_h, b_h)
# 
# Set tau_c = 1:
#   E(tau_c) = 1.0 and Var(tau_c) = sigma^2_c
# a_c = b_c = 1/sigma^2_c
# E(tau_h) = 1.81 [from above] and Var(tau_h) = sigma^2_h
# a_h = (1.81)^2/sigma^2_h and b_h = 1.81/sigma^2_h


IAR_bym <- '
data {
int<lower = 0> N;                                     // number of cells (= number of observations, including NAs)
int<lower = 0> N_obs;                                 // number of non-missing
int<lower = 0> N_mis;                                 // number missing
int<lower = 0> N_edges;                               // number of edges in adjacency matrix
int<lower = 1, upper = N> node1[N_edges];             // node1[i] adjacent to node2[i]
int<lower = 1, upper = N> node2[N_edges];             // and node1[i] < node2[i]

real<lower = 0, upper = 200> y_obs[N_obs];            // observed response data (excluding NAs)

int<lower = 1, upper = N_obs + N_mis> ii_obs[N_obs];  // indices of observed values
int<lower = 1, upper = N_obs + N_mis> ii_mis[N_mis];  // indices of unobserved values

real<lower = 0> sigma_y[N];
}

parameters {
real<lower = 1, upper = 200> y_mis[N_mis];            // missing response data

real beta0;                                           // intercept
vector[N] phi;                                        // spatial error component (centered on 0)
real<lower = 0> tau_phi;                            // to scale spatial error component

vector[N] theta;                                        // non-spatial error component (centered on 0)
real<lower = 0> tau_theta;                            // to scale non-spatial error component
}

transformed parameters {
real<lower = 0, upper = 200> y[N];

real<lower = 0> sigma_theta = inv(sqrt(tau_theta));
real<lower = 0> sigma_phi = inv(sqrt(tau_phi));

vector[N] mu;                                         // latent true halfmax values
mu = beta0 + phi * sigma_phi + theta * sigma_theta;                         // scaling by sigma_phi rather than phi ~ N(0, sigma_phi)

y[ii_obs] = y_obs;                                    // transform to accomodate missing data in Stan
y[ii_mis] = y_mis;                                    // transform to accomodate missing data in Stan

}


model {

y ~ normal(mu, sigma_y);

// the following computes the prior on phi on the unit scale with sd = 1
target += -0.5 * dot_self(phi[node1] - phi[node2]);

// soft sum-to-zero constraint on phi) - equivalent to mean(phi) ~ normal(0,0.001)
sum(phi) ~ normal(0, 0.001 * N);

beta0 ~ normal(120, 10);
tau_phi ~ gamma(1, 1);
theta ~ normal(0, 1);
tau_theta ~ gamma(3.27 ,1.81 );

}'





#spatial and non-spatial component (reparameterized) - also model obs error
#uses scaling factor fomr INLA to set priors more easily

IAR_bym2 <- '
data {
int<lower = 0> N;                                     // number of cells (= number of observations, including NAs)
int<lower = 0> N_obs;                                 // number of non-missing
int<lower = 0> N_mis;                                 // number missing
int<lower = 0> N_edges;                               // number of edges in adjacency matrix
int<lower = 1, upper = N> node1[N_edges];             // node1[i] adjacent to node2[i]
int<lower = 1, upper = N> node2[N_edges];             // and node1[i] < node2[i]

real<lower = 0, upper = 200> y_obs[N_obs];            // observed response data (excluding NAs)

int<lower = 1, upper = N_obs + N_mis> ii_obs[N_obs];  // indices of observed values
int<lower = 1, upper = N_obs + N_mis> ii_mis[N_mis];  // indices of unobserved values

real<lower = 0> sigma_y[N];                           // observed sd of data (observation error)

real<lower = 0> scaling_factor;                       // scales variances of spatial effects
}

parameters {
real<lower = 1, upper = 200> y_mis[N_mis];            // missing response data

real beta0;                                           // intercept

real<lower = 0> sigma;                                // spatial and non-spatial sd
real<lower = 0, upper = 1> rho;                       // proportion unstructure vs spatially structured variance

vector[N] phi;                                        // spatial error component (centered on 0)
vector[N] theta;                                      // non-spatial error component (centered on 0)
}

transformed parameters {
real<lower = 0, upper = 200> y[N];

vector[N] convolved_re;
vector[N] mu;                                         // latent true halfmax values

// variance of each component should be approx equal to 1
convolved_re = sqrt(1 - rho) * theta + sqrt(rho / scaling_factor) * phi;

mu = beta0 + convolved_re * sigma;                    // scaling by sigma_phi rather than phi ~ N(0, sigma_phi)

y[ii_obs] = y_obs;                                    // transform to accomodate missing data in Stan
y[ii_mis] = y_mis;                                    // transform to accomodate missing data in Stan
}


model {

y ~ normal(mu, sigma_y);

// the following computes the prior on phi on the unit scale with sd = 1
target += -0.5 * dot_self(phi[node1] - phi[node2]);

// soft sum-to-zero constraint on phi) - equivalent to mean(phi) ~ normal(0,0.001)
sum(phi) ~ normal(0, 0.001 * N);

beta0 ~ normal(120, 10);
theta ~ normal(0, 1);
sigma ~ normal(0, 5);
rho ~ beta(0.5, 0.5);

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
            chains = 3,                            # Number chains
            iter = 2000,                           # Iterations per chain
            cores = 3,                             # Number cores to use
            control = list(max_treedepth = 25, adapt_delta = 0.85)) # modified control parameters based on warnings;
# see http://mc-stan.org/misc/warnings.html
proc.time() - tt


#save to RDS
# saveRDS(fit, 'stan_bym2_fit.rds')



#diagnostics
pairs(fit, pars = c('sigma', 'rho', 'beta0'))


#https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html
sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
max_treedepth_by_chain <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
get_elapsed_time(fit)



#plotting diverged params
c_dark <- c("#8F272780")
green <- c("#00FF0080")

partition <- partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

par(mar = c(4, 4, 0.5, 0.5))
plot(nondiv_params$'sigma_y[1]', nondiv_params$'sigma_phi',
     col=c_dark, pch=16, cex=0.8, xlab="sigma_theta", ylab="sigma_phi")
points(div_params$'sigma_y[1]', div_params$'sigma_phi',
       col=green, pch=16, cex=0.8)


#shiny stan
library(shinystan)
launch_shinystan(fit)



# put copy of script in dir -----------------------------------------------

system(paste0('cp ', dir, 'Bird_Phenology/Scripts/ebird_Nestwatch/4-IAR-model.R ', IAR_dir, '/4-ICAR-model-', Sys.Date(), '.R'))


proc.time() - tt

MCMCsummary(fit, params = 'theta')
ttt <- MCMCsummary(fit, n.eff = TRUE, params = c('beta0', 'sigma', 'rho', 'phi', 'theta'))
min(ttt[,7])
fit

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
  scale_fill_gradientn(colors = c('red', 'blue')) +
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
  scale_fill_gradientn(colors = c('red', 'blue')) +
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

