######################
# 4 - IAR model - MORE FLEX SIGMA
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
# \sigma_{\nu[j]} \sim LN(\mu_{\sigma_{\nu}}, \sigma_{\sigma_{\nu}})
# \alpha_{\gamma} \sim N(0, 50)
# \beta_{\gamma} \sim N(2, 3)
# \sigma_{\gamma} \sim HN(0, 5)
# \mu_{\sigma_{\nu}} \sim N(0, 1.5)
# \sigma_{\sigma_{\nu}} \sim N(0, 0.5)
# \phi_{[i,j]} \sim N(0, [D - W]^{-1})
# \forall j \in \left \{1, ..., J  \right \}; \sum_{i}{} \phi_{[i,j]} = 0


#Joint distribution
#[\beta_{0_{y}}, \gamma_{c}, \rho, \theta_{cy}, \phi_{cy}, \sigma_{\nu_{y}}, \sigma_{\beta_{0}}, \alpha_{\gamma}, \beta_{\gamma}, \mu_{\sigma_{\nu}}, \sigma_{\sigma_{\nu}} | y_{obs_{cy}}, \sigma_{y_{cy}}] \propto \prod [y_{obs_{cy}}, \sigma_{y_{cy}} | \beta_{0_{y}}, \gamma_{c}, \rho, \theta_{cy}, \phi_{cy}, \sigma_{\nu_{y}}]  [\beta_{0_{y}} | \sigma_{\beta_{0}}] [\gamma_{c} | \alpha_{\gamma}, \beta_{\gamma}] [\sigma_{\nu_{y}} | \mu_{\sigma_{\nu}}, \sigma_{\sigma_{\nu}}][\sigma_{\beta_{0}}] [\alpha_{\gamma}] [\beta_{\gamma}][\rho][\theta_{cy}][\phi_{cy}][\mu_{\sigma_{\nu}}][\sigma_{\sigma_{\nu}}]
######################

#Stan resources:
#http://mc-stan.org/users/documentation/case-studies/icar_stan.html
#http://mc-stan.org/users/documentation/case-studies/divergences_and_bias.html
#https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html
#https://chi-feng.github.io/mcmc-demo/app.html#HamiltonianMC,standard
#https://groups.google.com/forum/#!msg/stan-users/zOjAeJC4x_E/OyCOfJo8AwAJ (regarding non-centered parameterization on sd)


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/labs/Tingley/phenomismatch/'



# db/hm query dir ------------------------------------------------------------

IAR_in_dir <- 'IAR_input_2019-05-03'
IAR_out_dir <- 'BYM_output_2019-11-12'



# Load packages -----------------------------------------------------------

library(rstan)
library(INLA)
library(geosphere)
library(ggplot2)
library(maps)
library(dplyr)
library(dggridR)
library(MCMCvis)
library(Matrix)
#Also need to be installed, but not loaded: rgeos, maptools, mapproj



# Set wd ------------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_in_dir))

IAR_in_date <- substr(IAR_in_dir, start = 11, stop = 20)
IAR_out_date <- substr(IAR_out_dir, start = 12, stop = 21)



# species arg -----------------------------------------------------

#args <- commandArgs(trailingOnly = TRUE)
#args <- as.character('Icterus_spurius')
#args <- as.character('Catharus_minimus')
#args <- as.character('Empidonax_virescens')
args <- as.character('Vireo_olivaceus')


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
             scaling_factor = scaling_factor,
             ii_obs = ii_obs_in,
             ii_mis = ii_mis_in,
             lat = cellcenters$lat_deg,
             y_PPC = y_PPC)



# Stan model --------------------------------------------------------------

#years in rows
#cells in columns

IAR_2 <- '
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
real<lower = 0> scaling_factor;                       // scales variances of spatial effects (estimated from INLA)
vector<lower = 24, upper = 90>[N] lat;
}

parameters {
// real<lower = 1, upper = 200> y_mis[N, J];             // missing response data
vector[J] y_mis[N];
real alpha_gamma_raw;
real beta_gamma_raw;                                       // effect of latitude
real<lower = 0> sigma_gamma_raw;
vector[j] gamma_raw;
vector[J] phi[N];                                     // spatial error component (scaled to N(0,1))
vector[J] theta[N];                                   // non-spatial error component (scaled to N(0,1))
real<lower = 0, upper = 1> rho;                       // proportion unstructured vs spatially structured variance
vector[N] sigma_nu_raw;
vector[N] beta0_raw;
real mu_sn_raw;
real<lower = 0> sigma_beta0_raw;
real<lower = 0> sigma_sn_raw;
}

transformed parameters {
// real<lower = 0, upper = 200> y[N, J];                 // response data to be modeled
vector[J] y[N];
vector[J] gamma;
real alpha_gamma;
real beta_gamma;
real<lower = 0> sigma_gamma;
vector[J] mu_gamma;
vector<lower = 0>[N] sigma_nu;
vector[J] y_true[N];
vector[J] nu[N];                                      // spatial and non-spatial component
real mu_sn;
real<lower = 0> sigma_beta0;
vector[N] beta0;
real<lower = 0> sigma_sn;

alpha_gamma = alpha_gamma_raw * 60;
beta_gamma = beta_gamma_raw * 3 + 2;
sigma_gamma = sigma_gamma_raw * 5;
sigma_beta0 = sigma_beta0_raw * 5;
mu_sn = mu_sn_raw * 1.5;
sigma_sn = sigma_sn_raw * 0.5;

mu_gamma = alpha_gamma + beta_gamma * lat;
gamma = gamma_raw * sigma_gamma + mu_gamma;
beta0 = beta0_raw * sigma_beta0;
sigma_nu = exp(sigma_nu_raw * sigma_sn + mu_sn);    //implies sigma_nu[i] ~ lognormal(mu_sn, sigma_sn) 

for (i in 1:N)
{
  nu[i] = sqrt(1 - rho) * theta[i] + sqrt(rho / scaling_factor) * phi[i]; // combined spatial/non-spatial
  
  y_true[j] = beta0[i] + gamma + nu[i] * sigma_nu[i];

  // indexing to avoid NAs  
  y[i, ii_obs[i, 1:N_obs[i]]] = y_obs[i, 1:N_obs[i]];
  y[i, ii_mis[i, 1:N_mis[i]]] = y_mis[i, 1:N_mis[i]];
}
}

model {

// priors
alpha_gamma_raw ~ std_normal();       // faster than normal(0, 1)
beta_gamma_raw ~ std_normal();
sigma_gamma_raw ~ std_normal();
gamma_raw ~ std_normal();
beta0_raw ~ std_normal();
sigma_beta0_raw ~ std_normal();
rho ~ beta(0.5, 0.5);
sigma_nu_raw ~ std_normal();
mu_sn_raw ~ std_normal();
sigma_sn_raw ~ std_normal();


for (i in 1:N)
{
  theta[j] ~ std_normal();
  target += -0.5 * dot_self(phi[i, node1] - phi[i, node2]);
  sum(phi[i]) ~ normal(0, 0.001 * N);
  
  y[i] ~ normal(y_true[i], sigma_y[i]);
}

}

generated quantities {

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
TREE_DEPTH <- 18
STEP_SIZE <- 0.0003
CHAINS <- 6
ITER <- 10000

tt <- proc.time()
fit <- rstan::stan(model_code = IAR_2,
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
                            'sigma_nu', 
                            'mu_sn', 
                            'sigma_sn', 
                            'rho', 
                            'nu', 
                            'theta', 
                            'phi',
                            'y_true', 
                            'y_rep'),
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60



#save to RDS
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_out_dir))
saveRDS(fit, file = paste0(args, '-', IAR_out_date, '-iar-stan_output.rds'))

#save data to RDS (has which cells are modeled)
saveRDS(DATA, file = paste0(args, '-', IAR_out_date, '-iar-stan_input.rds'))


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


# Summaries ---------------------------------------------------------------

#get summary of model output
model_summary <- MCMCvis::MCMCsummary(fit, Rhat = TRUE, n.eff = TRUE, round = 2, excl = 'y_rep')

#extract Rhat and neff values
rhat_output <- as.vector(model_summary[, grep('Rhat', colnames(model_summary))])
neff_output <- as.vector(model_summary[, grep('n.eff', colnames(model_summary))])



# Checks -------------------------------------------------------------

# for PPC extract y_rep and transpose (so iter are rows as required by shiny stan)
y_rep_ch <- MCMCvis::MCMCpstr(fit, params = 'y_rep', type = 'chains')[[1]]
t_y_rep <- t(y_rep_ch)

#remove NA vals
na.y.rm <- which(is.na(y_PPC))
n_y_PPC <- y_PPC[-na.y.rm]
n_y_rep <- t_y_rep[, -na.y.rm]

#mean resid (pred - actual) for each datapoint
ind_resid <- rep(NA, length(n_y_PPC))
for (i in 1:length(n_y_PPC))
{
  #i <- 1
  ind_resid[i] <- mean(n_y_rep[,i] - n_y_PPC[i])
}

#density overlay plot - first 100 iter
#modified bayesplot::ppc_dens_overlay function
tdata <- bayesplot::ppc_data(n_y_PPC, n_y_rep[1:100,])

annotations <- data.frame(xpos = c(-Inf, -Inf),
                          ypos = c(Inf, Inf),
                          annotateText = c(paste0('Max Rhat: ', max(rhat_output)),
                                           paste0('Min n.eff: ', min(neff_output))),
                          hjustvar = c(0, 0),
                          vjustvar = c(4, 6))

p <- ggplot(tdata) +
  aes_(x = ~value) +
  stat_density(aes_(group = ~rep_id, color = "yrep"),
               data = function(x) dplyr::filter(x, !tdata$is_y),
               geom = "line", position = "identity", size = 0.25,
               alpha = 0.3, trim = FALSE, bw = 'nrd0', adjust = 1,
               kernel = 'gaussian', n = 1024) +
  stat_density(aes_(color = "y"), data = function(x) dplyr::filter(x, tdata$is_y),
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
  ggtitle(paste0(args)) +
  geom_text(data = annotations, aes(x = xpos, y = ypos,
                                    hjust = hjustvar, vjust = vjustvar,
                                    label = annotateText),
            size = 3, col = 'black')

ggsave(paste0(args, '_dens_overlay.pdf'), p)

#bayesplot::ppc_intervals(n_y_PPC, y_rep[1:200,])
#bayesplot::ppc_ribbon(n_y_PPC, y_rep[1:200,])

#average yrep for each pnt
yrm <- apply(n_y_rep, 2, mean)
pdf(paste0(args, '_pred_true.pdf'))
plot(n_y_PPC, yrm, pch = 19, col = rgb(0,0,0,0.4),
     xlim = range(n_y_PPC, yrm), ylim = range(n_y_PPC, yrm),
     xlab = 'y', ylab = 'y_rep', main = paste0(args))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = 'red')
dev.off()


#histrogram of residuals
pdf(paste0(args, '_hist_resid.pdf'))
hist(ind_resid, main = paste0(args),
     xlab = 'Residuals (predicted - true)')
abline(v = 0, lty = 2, lwd = 3, col = 'red')
dev.off()

#PPC
# bayesplot::ppc_stat(n_y_PPC, n_t_y_rep, stat = 'mean')
# bayesplot::ppc_stat(n_y_PPC, n_t_y_rep, stat = 'max')
# bayesplot::ppc_stat(n_y_PPC, n_t_y_rep, stat = 'min')
# bayesplot::ppc_dens_overlay(n_y_PPC, n_t_y_rep[1:500,])


# write model results to file ---------------------------------------------

options(max.print = 5e6)
sink(paste0(args, '-iar-stan_results-', IAR_out_date, '.txt'))
cat(paste0('BYM results ', args, ' \n'))
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
usamap <- data.frame(maps::map("world", "USA", plot = FALSE)[c("x", "y")])
canadamap <- data.frame(maps::map("world", "Canada", plot = FALSE)[c("x", "y")])
mexicomap <- data.frame(maps::map("world", "Mexico", plot = FALSE)[c("x", "y")])

#min/max for plotting using input data
#f_rng <- c(range(f_out$HM_mean, na.rm = TRUE), range(med_fit, na.rm = TRUE))
#MIN <- round(min(f_rng))
#MAX <- round(max(f_rng))

#min/max for plotting using output data
MIN <- (round(min(med_fit)) - 1)
MAX <- (round(max(med_fit)) + 1)

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
ifelse(!dir.exists(paste0(dir, 'Bird_Phenology/Figures/pre_post_IAR_maps/', IAR_out_dir)),
       dir.create(paste0(dir, 'Bird_Phenology/Figures/pre_post_IAR_maps/', IAR_out_dir)),
       FALSE)


setwd(paste0(dir, 'Bird_Phenology/Figures/pre_post_IAR_maps/', IAR_out_dir))

#loop plots for each year
for (i in 1:length(years))
{
  #i <- 1
  
  #filter data for year[i]
  f_out_filt <- dplyr::filter(f_out, year == years[i])
  
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
              xlim = c(-100, -55), ylim = c(23, 66)) +
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
  t_med_fit <- med_fit[i,]
  t_sd_fit <- sd_fit[i,]
  
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

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_out_dir))


#alpha_gamma ~ normal(0, 60)
PR <- rnorm(10000, 0, 60)
MCMCvis::MCMCtrace(fit,
                   params = 'alpha_gamma',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-trace-alpha_gamma-', IAR_out_date, '.pdf'))

#beta_gamma ~ normal(2, 3)
PR <- rnorm(10000, 2, 3)
MCMCvis::MCMCtrace(fit,
                   params = 'beta_gamma',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-trace-beta_gamma-', IAR_out_date, '.pdf'))

#sigma_gamma ~ halfnormal(0, 5)
PR_p <- rnorm(10000, 0, 5)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_gamma',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-trace-sigma_gamma-', IAR_out_date, '.pdf'))

#mu_sn ~ halfnormal(0, 1.5)
PR <- rnorm(10000, 0, 1.5)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_sn',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-trace-mu_sn-', IAR_out_date, '.pdf'))

#sigma_sn ~ halfnormal(0, 1.5)
PR <- rnorm(10000, 0, 05)
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_sn',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-trace-sigma_sn-', IAR_out_date, '.pdf'))

#rho ~ beta(0.5, 0.5)
PR <- rbeta(10000, 0.5, 0.5)
MCMCvis::MCMCtrace(fit,
                   params = 'rho',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-trace-rho-', IAR_out_date, '.pdf'))

#sigma_beta0 ~ HN(0, 5)
PR_p <- rnorm(10000, 0, 5)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_beta0',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-trace-sigma_beta0-', IAR_out_date, '.pdf'))


if ('Rplots.pdf' %in% list.files())
{
  file.remove('Rplots.pdf')
}


print('I completed!')