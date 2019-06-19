#######
# calculate (interannual) variance in arrival date for each cell for a single species
# does not appear to be much differences among cells
# may need to model all species together?
#######


# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/UCHC/LABS/Tingley/phenomismatch/'


# species -----------------------------------------------------------------

#args <- commandArgs(trailingOnly = TRUE)
#args <- as.character('Catharus_minimus')
#args <- as.character('Empidonax_virescens')
args <- as.character('Vireo_olivaceus')



# other dir ---------------------------------------------------------------

IAR_in_dir <- paste0(dir, 'Bird_Phenology/Data/Processed/IAR_input_2019-05-03')
#IAR_out_dir <- paste0(dir, 'Bird_Phenology/Data/Processed/IAR_output_2019-05-26')
IAR_out_dir <- '~/Desktop/Bird_Phenology_Offline/Data/Processed/IAR_output_2019-05-26'
trends_out_dir <- '~/Desktop/Bird_Phenology_Offline/Data/Processed/trends_output_2019-05-26'


# Load packages -----------------------------------------------------------

library(dplyr)
library(ggplot2)
library(rstan)
library(MCMCvis)
library(dggridR)



# setwd -------------------------------------------------------------------

setwd(IAR_in_dir)

#CHANGE
nc_in_dir <- nchar(IAR_in_dir)
nc_out_dir <- nchar(IAR_out_dir)
IAR_in_date <- substr(IAR_in_dir, start = (nc_in_dir - 9), stop = nc_in_dir)
IAR_out_date <- substr(IAR_out_dir, start = (nc_out_dir - 9), stop = nc_out_dir)



# Filter data by species/years ------------------------------------------------------

#read in master df
df_master <- readRDS(paste0('IAR_input-', IAR_in_date, '.rds'))

#switch to out dir
setwd(IAR_out_dir)


#if that species RDS object exists in dir
pro_data <- data.frame()
if (length(grep(paste0(args, '-', IAR_out_date, '-iar-stan_output.rds'), list.files())) > 0)
{
  f_in_p <- dplyr::filter(df_master, species == args & MODEL == TRUE)
  
  #read in IAR model output and input
  t_fit <- readRDS(paste0(args, '-', IAR_out_date, '-iar-stan_output.rds'))
  t_data <- readRDS(paste0(args, '-', IAR_out_date, '-iar-stan_input.rds'))
  
  #only cells and years that were modeled (to account for any lone cells that were dropped in 4-IAR-arr.R)
  f_in <- dplyr::filter(f_in_p, cell %in% t_data$cells)
  
  t_cells <- unique(f_in$cell)
  t_years <- unique(f_in$year)
  
  #make hexgrid
  hexgrid6 <- dggridR::dgconstruct(res = 6)
  
  #get hexgrid cell centers
  cellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid6, t_cells)
  
  #extract median and sd for IAR arrival dates
  mean_fit <- MCMCvis::MCMCpstr(t_fit, params = 'y_true', func = mean)[[1]]
  med_fit <- MCMCvis::MCMCpstr(t_fit, params = 'y_true', func = median)[[1]]
  sd_fit <- MCMCvis::MCMCpstr(t_fit, params = 'y_true', func = sd)[[1]]
  
  #loop through years
  for (j in 1:length(t_years))
  {
    #j <- 1
    print(paste0('species: ', args, ', ', 
                 'year: ', t_years[j]))
    
    t_f_in <- dplyr::filter(f_in, year == t_years[j])
    
    t_full <- data.frame(t_f_in[,c('species','cell')], 
                         cell_lat = round(cellcenters$lat_deg, digits = 2), 
                         cell_lon = round(cellcenters$lon_deg, digits = 2),
                         t_f_in[,c('year', 'HM_mean', 'HM_sd')],
                         mean_post_IAR = mean_fit[,j], sd_post_IAR = sd_fit[,j])
    
    colnames(t_full)[6:7] <- c('mean_pre_IAR', 'sd_pre_IAR')
    
    pro_data <- rbind(pro_data, t_full)
  } #end year loop
} else {
  stop(paste0('.rds file for ', sp, ' not found in directory'))
}



# Process data ------------------------------------------------------------

#cell years with input data
#data_f <- pro_data[which(!is.na(pro_data$mean_pre_IAR)),]
data_f <- pro_data

#cells with more than three years of data
cnts <- plyr::count(data_f, 'cell')
u_cells <- cnts[which(cnts[,2] > 3),1]

data_f2 <- dplyr::filter(data_f, cell %in% u_cells)

t_cl <- unique(data_f2[,c('cell', 'cell_lat')])
ot_cl <- t_cl[order(t_cl[,1]),]

#years to generate low, mid, high lat phenology
sim_year <- min(data_f2$year - 2001):max(data_f2$year - 2001)

#create data list for Stan
DATA <- list(N = NROW(data_f2),
             y_obs = data_f2$mean_post_IAR,
             y_sd = data_f2$sd_post_IAR,
             cn_id = as.numeric(factor(data_f2$cell)),
             NC = length(u_cells),
             year = (data_f2$year - 2001),
             lat = ot_cl$cell_lat)


#y_true[i] ~ normal(mu[cell[i]], sigma[cell[i]])
#y_obs[i] ~ normal(y_true[i], y_sd[i])


# Stan model --------------------------------------------------------------

model <- "
data {
int<lower = 0> N;                                     // number of obs
vector<lower = 0, upper = 200>[N] y_obs;
vector<lower = 0>[N] y_sd;
int<lower = 1> cn_id[N];                              // cell ids
int<lower = 0> NC;                                    // number of cells
vector[NC] lat;
}

parameters {
vector[N] y_true_raw;
vector[NC] mu_raw;
vector[NC] sigma_raw;
real<lower = 0> sigma_mu_raw;
real mu_sigma_raw;
real<lower = 0> sigma_sigma_raw;
real alpha;
real beta;
}

transformed parameters {

vector[N] y_true;
vector<lower = 0>[NC] mu;
vector<lower = 0>[NC] sigma;
vector<lower = 0>[NC] mu_mu;
real<lower = 0> sigma_mu;
real mu_sigma;
real<lower = 0> sigma_sigma;

sigma_mu = sigma_mu_raw * 5;
mu_mu = alpha + beta * lat;
mu = mu_raw * sigma_mu + mu_mu;

mu_sigma = mu_sigma_raw * 1;
sigma_sigma = sigma_sigma_raw * 0.5;
sigma = exp(sigma_raw * sigma_sigma + mu_sigma);

for (i in 1:N)
{
  y_true[i] = y_true_raw[i] * sigma[cn_id[i]] + mu[cn_id[i]];
}
}

model {

y_true_raw ~ std_normal();
sigma_mu_raw ~ std_normal();
mu_raw ~ std_normal();
mu_sigma_raw ~ std_normal();
sigma_sigma_raw ~ std_normal();
sigma_raw ~ std_normal();
alpha ~ normal(125, 30);
beta ~ normal(0, 3);

y_obs ~ normal(y_true, y_sd);
}

generated quantities {

real y_rep[N];

// PPC
y_rep = normal_rng(y_true, y_sd);
}
"


# Run model ---------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.95
TREE_DEPTH <- 15
STEP_SIZE <- 0.001
CHAINS <- 4
ITER <- 6000

tt <- proc.time()
fit <- rstan::stan(model_code = model,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('mu', 'alpha', 'beta', 'sigma_mu', 
                            'sigma', 'mu_sigma', 'sigma_sigma',
                            'y_true', 'y_rep'),
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60

#save to RDS
#setwd(trends_out_dir)
setwd('~/Desktop')
saveRDS(fit, file = paste0(args, '-', IAR_out_date, '-variance.rds'))

MCMCvis::MCMCsummary(fit, excl = c('y_true', 'y_rep'))

MCMCvis::MCMCplot(fit, params = 'sigma')