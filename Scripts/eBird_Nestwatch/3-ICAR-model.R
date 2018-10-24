######################
# 3 - ICAR model (only spatial component)
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

# In NA_birdPhen2.R, I fit models with a nonspatial variance term, but I found that the estimated
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


# runtime -----------------------------------------------------------------

tt <- proc.time()



# Load packages -----------------------------------------------------------

library(rstan)
library(dplyr)



# Set wd ------------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/'))



# import eBird species list -----------------------------------------------------

species_list_i <- read.table('eBird_species_list.txt', stringsAsFactors = FALSE)
species_list <- species_list_i[,1]
nsp <- length(species_list)


#DATA ONLY VALID THROUGH 2017 (2018 data only goes to ~ jday 60 as of 2018-10-15 query)
years <- 2002:2017
nyr <- length(years)



# combine logit cubic results and diagnostic info -----------------------------------------------------------------

counter <- 0
for (i in 1:nsp)
{
  #i <- 80
  
  #import presence absence ebird data for each specices
  setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', db_dir))
  spdata <- readRDS(paste0('ebird_NA_phen_proc_', species_list[i], '.rds'))
  
  #import halfmax estimates and diagnostics from logit cubic model
  setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', hm_dir))
  temp_halfmax <- readRDS(paste0('halfmax_matrix_list_', species_list[i], '.rds'))
  temp_diag <- readRDS(paste0('halfmax_fit_diag_', species_list[i], '.rds'))
  
  if (i == 1)
  {
    #get number of unique cells
    cells <- unique(spdata$cell6)
    ncel <- length(cells)
    
    #create data.frame to fill
    diagnostics_frame <- as.data.frame(matrix(data = NA, nrow = nsp*ncel*nyr, ncol = 19))
    names(diagnostics_frame) <- c("species", "cell", "year", "n1", "n1W", "n0", "n0i", "njd1", "njd0", "njd0i",
                                  "nphen_bad", "min_n.eff", "max_Rhat", "HM_n.eff", "HM_Rhat",
                                  "HM_mean", "HM_sd", "HM_LCI", "HM_UCI")
  }
  
  #loop through years
  for (j in 1:nyr)
  {
    #j <- 10
    print(paste(i,j))
    ysdata <- dplyr::filter(spdata, year == years[j])
  
    for (k in 1:ncel)
    {
      #k <- 1
      counter <- counter + 1
      diagnostics_frame$species[counter] <- species_list[i]
      diagnostics_frame$year[counter] <- years[j]
      diagnostics_frame$cell[counter] <- cells[k]
      
      cysdata <- dplyr::filter(ysdata, cell6 == cells[k])
      
      #number of surveys where species was detected
      diagnostics_frame$n1[counter] <- sum(cysdata$detect)
      #number of surveys where species was not detected
      diagnostics_frame$n0[counter] <- sum(cysdata$detect == 0)
      #number of detections that came before jday 60
      diagnostics_frame$n1W[counter] <- sum(cysdata$detect*as.numeric(cysdata$day < 60))
      
      if (diagnostics_frame$n1[counter] > 0)
      {
        #number of non-detections before first detection
        diagnostics_frame$n0i[counter] <- length(which(cysdata$detect == 0 & 
                                                         cysdata$day < min(cysdata$day[which(cysdata$detect == 1)])))
        #number of unique days with detections
        diagnostics_frame$njd1[counter] <- length(unique(cysdata$day[which(cysdata$detect == 1)]))
        #number of unique days of non-detections before first detection
        diagnostics_frame$njd0i[counter] <- length(unique(cysdata$day[which(cysdata$detect == 0 & 
                                                                                cysdata$day < min(cysdata$day[which(cysdata$detect == 1)]))]))
      }
      
      #number of unique days with non-detection
      diagnostics_frame$njd0[counter] <- length(unique(cysdata$day[which(cysdata$detect == 0)]))
      
      if (diagnostics_frame$n1[counter] > 29 & 
          diagnostics_frame$n1W[counter] < (diagnostics_frame$n1[counter] / 50) &
          diagnostics_frame$n0[counter] > 29)
      {
        diagnostics_frame$min_n.eff[counter] <- temp_diag[[j]][[k]]$mineffsSize
        diagnostics_frame$max_Rhat[counter] <- temp_diag[[j]][[k]]$maxRhat
        
        halfmax_posterior <- as.vector(temp_halfmax[[j]][k,])
        
        #convert to mcmc.list and calc n_eff and Rhat using coda (DIFFERENT THAN STAN ESTIMATES)
        halfmax_mcmcList <- coda::mcmc.list(coda::as.mcmc(halfmax_posterior[1:500]), 
                                            coda::as.mcmc(halfmax_posterior[501:1000]),
                                            coda::as.mcmc(halfmax_posterior[1001:1500]), 
                                            coda::as.mcmc(halfmax_posterior[1501:2000]))
        
        diagnostics_frame$HM_n.eff[counter] <- round(coda::effectiveSize(halfmax_mcmcList), digits = 0)
        diagnostics_frame$HM_Rhat[counter] <- round(coda::gelman.diag(halfmax_mcmcList)$psrf[1], digits = 2)
        
        #determine how many estimates are 1 and not 1 (estimates of 1 are bogus)
        diagnostics_frame$nphen_bad[counter] <- sum(halfmax_posterior == 1)
        halfmax_posterior2 <- halfmax_posterior[which(halfmax_posterior != 1)]
        
        #calculate posterior mean and sd
        diagnostics_frame$HM_mean[counter] <- mean(halfmax_posterior)
        diagnostics_frame$HM_sd[counter] <- sd(halfmax_posterior)
        
        diagnostics_frame$HM_LCI[counter] <- quantile(halfmax_posterior, probs = 0.025)
        diagnostics_frame$HM_UCI[counter] <- quantile(halfmax_posterior, probs = 0.975)
        
        # #fit normal and cauchy distributions to data
        # normfit <- fitdistrplus::fitdist(halfmax_posterior2, "norm")
        # diagnostics_frame$Nloglik[counter] <- normfit$loglik
        # 
        # cfit <- NA
        # #cfit <- tryCatch(fitdistrplus::fitdist(halfmax_posterior2,"cauchy"), error=function(e){return(NA)})
        # if(!is.na(cfit))
        # {
        #   diagnostics_frame$HM_c_loc[counter] <- cfit$estimate[1]
        #   diagnostics_frame$HM_c_scale[counter] <- cfit$estimate[2]
        #   diagnostics_frame$Cloglik[counter] <- cfit$loglik
        # }
      }
    } # k -cell
  } # j - year
} # i - species


# #how many have crazy estimates for halfmax (1 for some iter)
# sum(diagnostics_frame$nphen_bad > 0, na.rm = TRUE)

#add species_cell column
diagnostics_frame$spCel <- paste(diagnostics_frame$species, diagnostics_frame$cell, sep="_")



# write to RDS ------------------------------------------------------------

ICAR_dir_path <- paste0(dir, 'Bird_phenology/Data/Processed/ICAR_', Sys.Date())

dir.create(ICAR_dir_path)
setwd(ICAR_dir_path)

saveRDS(diagnostics_frame, paste0('diagnostics_frame-', Sys.Date(), '.rds'))




# initial data filtering based on diagnostics -------------------------------------------------------------

#filter based on priors metrics 
diag_frame_p <- filter(diagnostics_frame, n1 > 29, n0 > 29, 
                     n1W < (n1/50), njd0i > 29, njd1 > 19)

#number of winter obs is greater than (number of detections/50) - winter obs make up more than 2% of obs
winter_spcels <- unique(diagnostics_frame$spCel[which(diagnostics_frame$n1W > (diagnostics_frame$n1/50))])

#alias for  !%in%
'%ni%' <- Negate('%in%')

#exclude 'winter pixels'
diag_frame <- diag_frame_p[which(diag_frame_p$spCel %ni% winter_spcels), ]


#how many cells have data for each year
for (i in min(years):max(years))
{
  #i <- 2002
  n_cell_yr <- length(unique(diag_frame$spCel[which(diag_frame$year==i)]))
  print(paste0('Cells for ', i, ': ', n_cell_yr))
}




# #determine which species/cells/years there are good data for -----------------

#create vector of year names
yrs_vec <- c()
for (i in min(years):max(years))
{
  yrs_vec <- c(yrs_vec, paste0('yr_', i))
}

#create df with species/cells/n_yrs per sp_cell
cells_frame <- data.frame(species = rep(species_list, each = length(cells)), 
                          cell = rep(cells, length(species_list)), 
                          n_yrs = NA)

#add years
cells_frame[yrs_vec] <- NA

#fill cells_frame
for (i in 1:nsp)
{
  #i <- 80
  t_sp <- filter(cells_frame, species == species_list[i])
  
  for (k in 1:ncel)
  {
    #k <- 30
    t_cell <- filter(t_sp, cell == cells[k])
    
    #filter diag frame by species and cell
    diag_temp <- filter(diag_frame, species == species_list[i], cell == cells[k])
    
    #get index for species/cell row
    idx <- which(cells_frame$species == species_list[i] & 
                               cells_frame$cell == cells[k])
    
    #insert number of yrs with data
    cells_frame[idx,'n_yrs'] <- NROW(diag_temp)
    
    #find out which years have data and insert TRUE into df - leave NA otherwise
    diag_temp_yrs <- paste0('yr_', diag_temp$year)
    cells_frame[idx, which(colnames(cells_frame) %in% diag_temp_yrs)] <- TRUE
    # f_yrs <- which(colnames(cells_frame) %ni% diag_temp_yrs)[-c(1:3)]
    # cells_frame[idx, f_yrs] <- 0
  }
}

#look at just one species
filter(cells_frame, species == species_list[80])



# Which species meet data threshold ---------------------------------------

# Which species-years are worth modeling? At a bare minimum:
#     Species with at least 'NC' cells in all three years from 2015-2017
#     Species-years with at least 'NC' cells for those species

#threshold
NC <- 3

out <- data.frame()

#which species to model
spYs <- data.frame(species = species_list)
for (i in 1:length(species_list))
{
  #i <- 80
  #filter by species
  t_sp <- dplyr::filter(cells_frame, species == species_list[i])
  
  #get number of cells with data for each year in 2015-2017
  d_yrs <- apply(t_sp[c('yr_2015', 'yr_2016', 'yr_2017')], 2, function(x) sum(x, na.rm = TRUE))
  
  #if all three years have more than 'NC' cells of data, figure out which years have at least 'NC' cells
  if (sum(d_yrs >= NC) == 3)
  {
    #calc nunber of cells for each year (apply function only to year columns)
    n_cells_yr <- apply(t_sp[-which(c(names(t_sp) %in% c('species', 'cell', 'n_yrs')))], 2, 
                        function(x) sum(x, na.rm = TRUE))
    
    #get which years meet NC threshold
    cells_yr <- as.numeric(substr(names(which(n_cells_yr >= NC)), start = 4, stop = 7))
    
    #filter diag_frame for species
    diag_frame_t_sp <- filter(diag_frame, species == species_list[i])
    
    #only years with > NC cells of data
    t_out <- diag_frame_t_sp[which(diag_frame$year %in% cells_yr),]
    
    #add rows that meet thresholds to output df
    out <- rbind(out, t_out)
  }
}




# ICAR model data prep ----------------------------------------------------


# ICAR models

# General data preparation

ncel <- length(cells)
cellcenters <- dgSEQNUM_to_GEO(hexgrid6, as.numeric(cells))
adjacency_matrix <- matrix(data=NA, nrow = length(cells), ncol=length(cells))

for(i in 1:length(cells))
{
  for(j in i:length(cells))
  {
    dists <- geosphere::distm(c(cellcenters$lon_deg[i], cellcenters$lat_deg[i]),
                              c(cellcenters$lon_deg[j], cellcenters$lat_deg[j]))
    adjacency_matrix[i,j] <- as.numeric((dists/1000) > 0 & (dists/1000) < 311)
  }
}

ninds <- which(adjacency_matrix==1, arr.ind = T)

stan_ICAR_infoPrior <- '
data {
int<lower=0> N; // number of cells (= number of observations, including NAs)
int<lower=0> N_obs; // number of non-missing
int<lower=0> N_mis; // number missing
int<lower=0> N_edges; // number of edges in adjacency matrix
int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]
real<lower=0, upper=200> obs[N_obs];    // observed data (excluding NAs)
int<lower = 1, upper = N_obs + N_mis> ii_obs[N_obs];
int<lower = 1, upper = N_obs + N_mis> ii_mis[N_mis];
real<lower=0> sds[N];                 // sds for ALL data (observed and unobserved)
}

parameters {
real<lower=1, upper=200> y_mis[N_mis];         // missing data
real beta0;                // intercept

real<lower=0> sigma_theta;   // sd of heterogeneous effects
real<lower=0> sigma_phi;     // sd of spatial effects

vector[N] theta;       // heterogeneous effects
vector[N] phi;         // spatial effects
vector[N] latent;      // latent true halfmax values
}

transformed parameters {
real<lower=0, upper=200> y[N];
y[ii_obs] = obs;
y[ii_mis] = y_mis;
}


model {
y ~ normal(latent, sds);
latent ~ normal(beta0 + phi * sigma_phi, sigma_theta);

// the following computes the prior on phi on the unit scale with sd = 1
target += -0.5 * dot_self(phi[node1] - phi[node2]);

// soft sum-to-zero constraint on phi)
sum(phi) ~ normal(0, 0.001 * N);  // equivalent to mean(phi) ~ normal(0,0.001)

beta0 ~ normal(0, 100);
theta ~ normal(0, 30);
sigma_theta ~ normal(.7, .05);
}'

stan_ICAR_no_nonspatial <- '
data {
int<lower=0> N; // number of cells (= number of observations, including NAs)
int<lower=0> N_obs; // number of non-missing
int<lower=0> N_mis; // number missing
int<lower=0> N_edges; // number of edges in adjacency matrix
int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]
real<lower=0, upper=200> obs[N_obs];    // observed data (excluding NAs)
int<lower = 1, upper = N_obs + N_mis> ii_obs[N_obs];
int<lower = 1, upper = N_obs + N_mis> ii_mis[N_mis];
real<lower=0> sds[N];                 // sds for ALL data (observed and unobserved)
}

parameters {
real<lower=1, upper=200> y_mis[N_mis];         // missing data
real beta0;                // intercept
real<lower=0> sigma_phi;     // sd of spatial effects
vector[N] theta;       // heterogeneous effects
vector[N] phi;         // spatial effects
}

transformed parameters {
real<lower=0, upper=200> y[N];
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

diag_frame2$spyr <- paste(diag_frame2$species, diag_frame2$year, sep="_")
spyrs <- unique(diag_frame2$spyr)
nfits <- length(spyrs)

# CHANGE LINE BELOW TO REFLECT AVAILABLE NUMBER OF CORES
ncores <- 3

Ldatalist <- list()
counter <- 0
for(i in 1:ncores)
{
  Ldatalist[[i]] <- list()
  
  for(j in 1:ceiling(nfits/ncores))
  {
    counter <- counter + 1
    
    if(counter <= nfits)
    {
      Ldatalist[[i]][[j]] <- 
        diag_frame2[which(diag_frame2$spyr == spyrs[counter]), ]
    }
  }
}

fit_funct <- function(datalist, model)
{
  output_list <- list()
  
  for(j in 1:length(datalist))
  {
    dd <- datalist[[j]]
    SY <- data.frame(cell = cells, halfmax = NA, sd = NA)
    
    for(i in 1:length(cells))
    {
      d1 <- which(dd$cell == as.integer(cells[i]))
      
      if(length(d1) > 0)
      {
        SY$halfmax[i] <- dd$phen_mean[d1]
        SY$sd[i] <- dd$phen_sd[d1]
      }else{
        SY$sd[i] <- .01
      }
    }
    
    SYL <- list(N = length(cells), 
                N_edges = nrow(ninds), 
                node1 = ninds[,1], 
                node2 = ninds[,2],
                obs = SY$halfmax[which(!is.na(SY$halfmax))], 
                ii_obs = which(!is.na(SY$halfmax)),
                ii_mis = which(is.na(SY$halfmax)),
                N_obs = sum(!is.na(SY$halfmax)),
                N_mis = sum(is.na(SY$halfmax)),
                sds = SY$sd)
    
    
    
    output_list[[j]] <- stan(
      model_code = model,  # Stan program
      data = SYL,    # named list of data
      chains = 3,             # number of Markov chains
      iter = 6000,            # total number of iterations per chain
      cores = 1,              # number of cores
      control = list(max_treedepth = 20, adapt_delta = .9) # modified control parameters based on warnings;
      # see http://mc-stan.org/misc/warnings.html
    )
  }
  return(output_list)
}

cl <- makeCluster(ncores)
registerDoParallel(cl)

all_output <- 
  foreach(i=1:ncores, .packages = 'rstan') %dopar% 
  fit_funct(Ldatalist[[i]], stan_ICAR_no_nonspatial)









# put copy of script in dir -----------------------------------------------

system(paste0('cp ', dir, 'Bird_Phenology/Scripts/ebird_Nestwatch/3-ICAR-model.R ', ICAR_dir_path, '/3-ICAR-model-', Sys.Date(), '.R'))


proc.time() - tt

