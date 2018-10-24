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
#dir <- '~/Google_Drive/R/'

#Xanadu
dir <- '/home/CAM/cyoungflesh/phenomismatch/'



# db/hm query dir ------------------------------------------------------------

db_dir <- 'db_query_2018-10-15'
hm_dir <- 'halfmax_species_2018-10-23'
ICAR_dir <- 'ICAR_2018-10-24'


# runtime -----------------------------------------------------------------

tt <- proc.time()



# Load packages -----------------------------------------------------------

library(rstan)
library(dplyr)
library(dggridR)
library(geosphere)


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

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', ICAR_dir))

diagnostics_frame <- readRDS('diagnostics_frame-2018-10-24.rds')
cells_frame <- readRDS('cells_frame-2018-10-24.rds')
yrs_frame <- readRDS('yrs_frame-2018-10-24.rds')




# Species list to model ----------------------------------------------------


NC <- 3

#species list to model (greater than or = to 'NC' cells for 2015:2017)
m_species_list <- c()
for (i in 1:nsp)
{
  #i <- 80
  t_yrsf <- dplyr::filter(yrs_frame, species == species_list[i])
  
  if (sum(t_yrsf$n_cells[which(t_yrsf$year %in% 2015:2017)] >= NC) == 3)
  {
    m_species_list <- c(m_species_list, species_list[i])
  }
}



# create adjacency matrix -------------------------------------------------

#unique cells (exclude NA as debuggin done on subset of species)
cells <- unique(diagnostics_frame$cell[!is.na(diagnostics_frame$cell)])
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
    #j <- 1
    dists <- geosphere::distm(c(cellcenters$lon_deg[i], cellcenters$lat_deg[i]),
                              c(cellcenters$lon_deg[j], cellcenters$lat_deg[j]))
    adjacency_matrix[i,j] <- as.numeric((dists/1000) > 0 & (dists/1000) < 311)
  }
}

#indices for 1s
ninds <- which(adjacency_matrix == 1, arr.ind = TRUE)




# Create Stan data object ------------------------------------------------------

#filter by species/year here
i <- 1
j <- 15
f_out <- filter(diagnostics_frame, species == m_species_list[i], year == years[j])

if (m_crit == TRUE)
{
  #proceed
} else {
  #stop
}

#check to make sure cells order == order in f_out
all.equal(f_out$cell, cells)


#create data list for Stan
DATA <- list(N = length(cells), 
             N_edges = nrow(ninds), 
             node1 = ninds[,1], 
             node2 = ninds[,2],
             obs = f_out$HM_mean[which(!is.na(f_out$HM_mean))], 
             ii_obs = which(!is.na(f_out$HM_mean)),
             ii_mis = which(is.na(f_out$HM_mean)),
             N_obs = sum(!is.na(f_out$HM_mean)),
             N_mis = sum(is.na(f_out$HM_mean)),
             sds = f_out$HM_sd)

#insert 0.01 as sd for missing vals
DATA$sds[which(is.na(DATA$sds))] <- 0.01


# Stan model --------------------------------------------------------------


stan_ICAR_no_nonspatial <- '
data {
int<lower = 0> N; // number of cells (= number of observations, including NAs)
int<lower = 0> N_obs; // number of non-missing
int<lower = 0> N_mis; // number missing
int<lower = 0> N_edges; // number of edges in adjacency matrix
int<lower = 1, upper = N> node1[N_edges];  // node1[i] adjacent to node2[i]
int<lower = 1, upper = N> node2[N_edges];  // and node1[i] < node2[i]
real<lower = 0, upper = 200> obs[N_obs];    // observed data (excluding NAs)
int<lower = 1, upper = N_obs + N_mis> ii_obs[N_obs];
int<lower = 1, upper = N_obs + N_mis> ii_mis[N_mis];
real<lower = 0> sds[N];                 // sds for ALL data (observed and unobserved)
}

parameters {
real<lower = 1, upper = 200> y_mis[N_mis];         // missing data
real beta0;                // intercept
real<lower = 0> sigma_phi;     // sd of spatial effects
vector[N] theta;       // heterogeneous effects - DOESNT APPEAR TO BE ELSEWHERE IN THE MODEL
vector[N] phi;         // spatial effects
}

transformed parameters {
real<lower = 0, upper = 200> y[N];
vector[N] latent; // latent true halfmax values
y[ii_obs] = obs;
y[ii_mis] = y_mis;
latent = beta0 + phi * sigma_phi; 
}


model {
y ~ normal(latent, sds);

// the following computes the prior on phi on the unit scale with sd = 1
target += -0.5 * dot_self(phi[node1] - phi[node2]);

// soft sum-to-zero constraint on phi)
sum(phi) ~ normal(0, 0.001 * N);  // equivalent to mean(phi) ~ normal(0,0.001)

beta0 ~ normal(0, 100);
theta ~ normal(0, 30);
}'



# Run model ---------------------------------------------------------------

#don't worry about compiler error messages (https://discourse.mc-stan.org/t/boost-and-rcppeigen-warnings-for-r-package-using-stan/3478) - need to include control line to avoid jaconian warning messages

out <- stan(model_code = stan_ICAR_no_nonspatial,  # Stan program
            data = DATA,                           # named list of data
            chains = 3,                            # number of Markov chains
            iter = 6000,                           # total number of iterations per chain
            cores = 3,                             # number of cores
            control = list(max_treedepth = 20, adapt_delta = .9)) # modified control parameters based on warnings;
            # see http://mc-stan.org/misc/warnings.html




# put copy of script in dir -----------------------------------------------

system(paste0('cp ', dir, 'Bird_Phenology/Scripts/ebird_Nestwatch/3-ICAR-model.R ', ICAR_dir_path, '/3-ICAR-model-', Sys.Date(), '.R'))


proc.time() - tt

