######################
# 4 - arr IAR model
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

IAR_in_dir <- 'IAR_input_2020-05-07'
IAR_out_dir <- 'IAR_output_2020-06-08'


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
IAR_out_date <- substr(IAR_out_dir, start = 12, stop = 21)


# species arg -----------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)


# Filter data by species/years ------------------------------------------------------

#read in master df
df_master <- readRDS(paste0('IAR_input-', IAR_in_date, '.rds'))

#filter by species and year to be modeled
f_out <- dplyr::filter(df_master, species == args[1] & MODEL == TRUE)

#fill invalid rows with NA
f_idx <- which(f_out$VALID == FALSE)
f_out$arr_GAM_mean[f_idx] <- NA
f_out$arr_GAM_sd[f_idx] <- NA

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
  sigma_y_in[i,] <- temp_yr$arr_GAM_sd
  
  for (j in 1:ncell)
  {
    #n <- 1
    #matrix with observed values with NAs
    y_PPC[counter] <- temp_yr$arr_GAM_mean[j]
    counter <- counter + 1
  }
  
  #which are not NA
  no_na <- temp_yr$arr_GAM_mean[which(!is.na(temp_yr$arr_GAM_mean))]
  
  #pad end with NAs
  if (length(no_na) < ncell)
  {
    num_na <- ncell - length(no_na)
    
    #add NAs to end
    t_y_obs_in <- c(no_na, rep(NA, num_na))
    t_obs_in <- c(which(!is.na(temp_yr$arr_GAM_mean)), rep(NA, num_na)) 
    t_mis_in <- c(which(is.na(temp_yr$arr_GAM_mean)), rep(NA, length(no_na)))
    
    #fill objects
    ii_obs_in[i,] <- t_obs_in
    ii_mis_in[i,] <- t_mis_in
    y_obs_in[i,] <- t_y_obs_in
  } else {
    #no NAs to end (no mimssing values)
    y_obs_in[i,] <- no_na
    ii_mis_in[i,] <- which(!is.na(temp_yr$arr_GAM_mean))
    y_obs_in[i,] <- which(is.na(temp_yr$arr_GAM_mean))
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
while ((max(rhat_output) > 1.02 | min(neff_output) < (CHAINS * 100)) & ITER < 10001)
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
saveRDS(fit, file = paste0(args[1], '-iar-stan_output-', IAR_out_date, '.rds'))

#save data to RDS (has which cells are modeled)
saveRDS(DATA, file = paste0(args[1], '-iar-stan_input-',  IAR_out_date, '.rds'))


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
y_rep_ch <- MCMCvis::MCMCpstr(fit, params = 'y_rep', type = 'chains')[[1]]
t_y_rep <- t(y_rep_ch)

#remove NA vals
na.y.rm <- which(is.na(DATA$y_PPC))
n_y_PPC <- DATA$y_PPC[-na.y.rm]
n_y_rep <- t_y_rep[, -na.y.rm]

#PPC
PPC_fun <- function(FUN, YR = n_y_rep, D = n_y_PPC)
{
  out <- sum(apply(YR, 1, FUN) > FUN(D)) / NROW(YR)
  print(out)
}
PPC_mn <- PPC_fun(mean)

pdf(paste0(args[1], '_PPC_mn.pdf'))
bayesplot::ppc_stat(n_y_PPC, n_y_rep, stat = 'mean')
dev.off()

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
  ggtitle(paste0(args[1])) +
  geom_text(data = annotations, aes(x = xpos, y = ypos,
                                    hjust = hjustvar, vjust = vjustvar,
                                    label = annotateText),
            size = 3, col = 'black')

ggsave(paste0(args[1], '_dens_overlay.pdf'), p)


#average yrep for each pnt
yrm <- apply(n_y_rep, 2, mean)
pdf(paste0(args[1], '_pred_true.pdf'))
plot(n_y_PPC, yrm, pch = 19, col = rgb(0,0,0,0.4),
     xlim = range(n_y_PPC, yrm), ylim = range(n_y_PPC, yrm),
     xlab = 'y', ylab = 'y_rep', main = paste0(args[1]))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = 'red')
dev.off()


#histrogram of residuals
pdf(paste0(args[1], '_hist_resid.pdf'))
hist(ind_resid, main = paste0(args[1]),
     xlab = 'Residuals (predicted - true)')
abline(v = 0, lty = 2, lwd = 3, col = 'red')
dev.off()


# write model results to file ---------------------------------------------

options(max.print = 5e6)
sink(paste0(args[1], '-iar-stan_results-', IAR_out_date, '.txt'))
cat(paste0('IAR results ', args[1], ' \n'))
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
cat(paste0('Max Rhat: ', max(rhat_output), ' \n'))
cat(paste0('Min n.eff: ', min(neff_output), ' \n'))
cat(paste0('PPC (mean): ', round(PPC_mn, 3), ' \n'))
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
#f_rng <- c(range(f_out$arr_GAM_mean, na.rm = TRUE), range(med_fit, na.rm = TRUE))
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
  
  #merge hex spatial data with GAM data
  to_plt <- dplyr::inner_join(f_out_filt, cell_grid, by = 'cell')
  to_plt2 <- dplyr::inner_join(to_plt, ll_df, by = 'cell')
  
  #pre-IAR
  p <- ggplot() +
    geom_path(data = worldmap,
              aes(x = x, y = y), color = 'black') +
    # geom_polygon(data = nrng.df,
    #           aes(x = long, y = lat, group=group), fill = 'green', alpha = 0.4) +
    # geom_polygon(data = nrng_rm.df,
    #              aes(x = long, y = lat, group=group), fill = 'orange', alpha = 0.4) +
    coord_map("ortho", orientation = c(35, -80, 0),
              xlim = c(-100, -55), ylim = c(23, 66)) +
    geom_polygon(data = to_plt2, aes(x = long, y = lat.y, group = group, fill = arr_GAM_mean),
                 alpha = 0.4) +
    geom_path(data = to_plt2, aes(x = long, y = lat.y, group = group),
              alpha = 0.4, color = 'black') +
    scale_fill_gradientn(colors = c('red', 'blue'),
                         limits = c(MIN, MAX)) +
    labs(fill = 'Estimated Arrival') +
    annotate('text', x = to_plt2$lon_deg, y = to_plt2$lat_deg + 0.5,
             label = round(to_plt2$arr_GAM_mean, digits = 0), col = 'black', alpha = 0.2,
             size = 4) +
    annotate('text', x = to_plt2$lon_deg, y = to_plt2$lat_deg - 0.5,
             label = round(to_plt2$arr_GAM_sd, digits = 0), col = 'white', alpha = 0.3,
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

#alpha_gamma ~ normal(0, 100)
PR <- rnorm(10000, 0, 100)
MCMCvis::MCMCtrace(fit,
                   params = 'alpha_gamma',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args[1], '-trace_alpha_gamma-', IAR_out_date, '.pdf'))

#beta_gamma ~ normal(3, 3)
PR <- rnorm(10000, 3, 3)
MCMCvis::MCMCtrace(fit,
                   params = 'beta_gamma',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args[1], '-trace_beta_gamma-', IAR_out_date, '.pdf'))

#sigma_gamma ~ halfnormal(0, 10)
PR_p <- rnorm(10000, 0, 10)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_gamma',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args[1], '-trace_sigma_gamma-', IAR_out_date, '.pdf'))

#sigma_beta0 ~ HN(0, 10)
PR_p <- rnorm(10000, 0, 10)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_beta0',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args[1], '-trace_sigma_beta0-', IAR_out_date, '.pdf'))

#sigma_y_true ~ HN(0, 10)
PR_p <- rnorm(10000, 0, 10)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_y_true',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args[1], '-trace_sigma_y_true-', IAR_out_date, '.pdf'))

#sigma_phi ~ HN(0, 10)
PR_p <- rnorm(10000, 0, 10)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_phi',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args[1], '-trace_sigma_phi-', IAR_out_date, '.pdf'))


if ('Rplots.pdf' %in% list.files())
{
  file.remove('Rplots.pdf')
}

print('I completed!')

