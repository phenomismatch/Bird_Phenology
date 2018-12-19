####################
# 6 - nesting date (from Nestwatch) ~ nesting date (from IAR model)
#
# *How well do arrival dates (as determined from the IAR model) predict nesting dates (as determined frmo Nestwatch)
####################



# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/'


# Load packages -----------------------------------------------------------

library(dplyr)
library(dggridR)
library(ggplot2)
library(rstan)
library(MCMCvis)


# import IAR data ---------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))

IAR_out_dir <- 'IAR_output_2018-11-12'
IAR_out_date <- substr(IAR_out_dir, start = 12, stop = 21)

IAR_data <- readRDS(paste0('master_arrival_', IAR_out_date, '.rds'))



# import IAR species list -----------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/'))

species_list_i <- read.table('IAR_species_list.txt', stringsAsFactors = FALSE)

#remove underscore and coerce to vector
species_list_i2 <- as.vector(apply(species_list_i, 2, function(x) gsub("_", " ", x)))


# create hex grid ---------------------------------------------------------

hexgrid6 <- dggridR::dgconstruct(res = 6) 



# Process Nestwatch data ---------------------------------------------------------------

#ebird specices codes
setwd(paste0(dir, 'Bird_Phenology/Data/'))

ebird_tax <- read.csv('eBird_Taxonomy_v2018_14Aug2018.csv')
ebird_tax2 <- dplyr::select(ebird_tax, SPECIES_CODE, PRIMARY_COM_NAME, SCI_NAME)
#replace space with underscore
ebird_tax2$SCI_NAME <- gsub(" ", "_", ebird_tax2$SCI_NAME)


#nestwatch data
setwd(paste0(dir, 'Bird_Phenology/Data/Nestwatch'))

nw_data <- read.csv('Nestwatch_2018_1026.csv')

#add grid cells
nw_data$CELL <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                         in_lon_deg = nw_data$LONGITUDE, 
                                         in_lat_deg = nw_data$LATITUDE)[[1]]

#add year
nw_data$YEAR <- substr(nw_data$PROJ_PERIOD_ID, start = 4, stop = 7)

#merge with ebird taxonomy
nw_data2 <- dplyr::left_join(nw_data, ebird_tax2, by = 'SPECIES_CODE')

#only keep entries that have dates - convert to julian day
nw_data2$FIRST_LAY_DT <- format(as.Date(substr(nw_data2$FIRST_LAY_DT, 1, 9), format = '%d%B%Y'), '%j')
nw_data2$HATCH_DT <- format(as.Date(substr(nw_data2$HATCH_DT, 1, 9), format = '%d%B%Y'), '%j')
nw_data2$FLEDGE_DT <- format(as.Date(substr(nw_data2$FLEDGE_DT, 1, 9), format = '%d%B%Y'), '%j')
to.rm <- which(is.na(nw_data2$FIRST_LAY_DT))
nw_data3 <- nw_data2[-to.rm, ]

#filter by species used in IAR
nw_data4 <- dplyr::filter(nw_data3, SCI_NAME %in% species_list_i[,1])

#number of Nestwatch observations for each species
plyr::count(nw_data4, "SCI_NAME")


species <- sort(unique(nw_data4$SCI_NAME))

#determine if Nestwatch data is sufficicent for that species/year
#combine nest watch data with IAR data
nw_IAR <- data.frame()
for (i in 1:length(species))
{
  #i <- 34

  #filter by species
  sp <- species[i]
  t_IAR <- dplyr::filter(IAR_data, species == sp)
  t_nw <- dplyr::filter(nw_data4, SCI_NAME == sp)
  
  #for Nestwatch, only keep cell/years that are in the IAR data
  #merge by multiple metrics
  t_mrg <- left_join(t_nw, t_IAR, by = c('CELL' = 'cell', 'YEAR' = 'year'))

  #determine how many obs there are for each cell/year
  t_nr <- plyr::count(t_nw, c('CELL', 'YEAR'))

  #make sure there are at least three data points per cell/year??
  if (t_nr > 3)
  {
    #USE THOSE CELL/YEARS FOR THAT SPECIES
  }
}





#hierarchical model to relate arrival date (from IAR model) to nesting date (from Nest Watch)

#one model per species






#############################
#Process data for model input
#############################






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

#Might be alright - data needs to be formatted properly though
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

nest_arrival <- '
data {
int<lower = 0> N;                                     // number of cell/years
int<lower = 0> J;                                     // max number of obs in a cell/year for Nestwatch
real<lower = 0, upper = 200> obs_IAR[N];              // mean output IAR model
real sigma_IAR[N];                                    // sd output IAR model
real p_obs_NW[J,N];                                   // nestwatch data
int<lower = 0> ii_obs[J, N];                          // indices of observed data
int<lower = 0> ii_mis[J, N];                          // indices of missing data
int<lower = 0> N_obs[N];                              // number of non-missing for each cell/year
int<lower = 0> N_mis[N];                              // number missing for each cell/year
}


parameters {
real<lower = 0, upper = 200> p_obs_NW_mis[J, N];      // missing nestwatch data
real<lower = 0> true_IAR[N];                          //true arrival date
real<lower = 0> true_NW[N];                           //true nesting date
real<lower = 0> sigma_NW;                             //sd for cell/year averaging
real mu[N];
real<lower = 0> sigma;
real alpha;
real beta;
}

transformed parameters {
real<lower = 0, upper = 200> obs_NW[J, N];            // nestwatch data to be modeled

// indexing to avoid NAs
for (n in 1:N)
{
  obs_NW[ii_obs[1:N_obs[n], n], n] = p_obs_NW[1:N_obs[n], n];
  obs_NW[ii_mis[1:N_mis[n], n], n] = p_obs_NW_mis[1:N_mis[n], n];
}
}

model {

// observation model - modeling true state as a function of some observed state

obs_IAR ~ normal(true_IAR, sigma_IAR); //sigma_IAR is given (one for each observation)

for (n in 1:N) //n is cell/year number
{
  obs_NW[,n] ~ normal(true_NW[n], sigma_NW); //sigma_NW is estimated (one per species)
}

  true_NW ~ normal(mu, sigma); //one sigma, one mu for each cell/year
  mu = alpha + beta * true_IAR; //one alpha, one beta
}'




# Run model ---------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

tt <- proc.time()
fit <- stan(model_code = nest_arrival,
            data = DATA,
            chains = 4,
            iter = 2,
            cores = 4,
            pars = c('true_IAR', 'true_NW', 'sigma_NW', 'alpha', 'beta', 'mu', 'sigma'),
            control = list(max_treedepth = 25, adapt_delta = 0.95, stepsize = 0.005)) # modified control parameters based on warnings
run_time <- (proc.time() - tt[3]) / 60







#save to RDS
# setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_out_dir))
# saveRDS(fit, file = paste0('IAR_stan_', args, '-', IAR_out_date, '.rds'))
# fit <- readRDS('IAR_stan_Catharus_minimus-2018-11-12.rds')




# diagnostics -------------------------------------------------------------


# pairs(fit, pars = c('sigma', 'rho'))

# sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
# mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
# max_treedepth_by_chain <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
# get_elapsed_time(fit)

# MCMCtrace(fit)
# MCMCsummary(fit, params = c('sigma', 'rho', 'beta0'), n.eff = TRUE)
# MCMCsummary(fit, params = c('theta', 'phi'), n.eff = TRUE)

# print(fit, pars = c('sigma', 'rho'))

# #shiny stan
# library(shinystan)
# launch_shinystan(fit)




# write model results to file ---------------------------------------------

# options(max.print = 50000)
# sink(paste0('IAR_results_', args, '.txt'))
# cat(paste0('IAR results ', args, ' \n'))
# cat(paste0('Total minutes: ', round(run_time, digits = 2), ' \n'))
# print(fit)
# sink()








