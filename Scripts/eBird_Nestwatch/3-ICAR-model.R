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



cy_dir <- '~/Google_Drive/R/'


# Load packages -----------------------------------------------------------

library(dggridR)
library(coda)
library(ggplot2)
library(viridis)
library(rstan)
library(doParallel)
library(foreach)



# Set wd ------------------------------------------------------------------

setwd(paste0(cy_dir, 'Bird_Phenology/Data/'))



# import eBird species list -----------------------------------------------------

species_list_i <- read.table('eBird_species_list.txt')
species_list <- species_list_i[,1]
nsp <- length(species_list)



# import data ------------------------------------------------------------

'%ni%' <- Negate('%in%')


hexgrid6 <- dgconstruct(res=6) # Construct geospatial hexagonal grid
years <- 2002:2016
nyr <- length(years)







load('/Users/Tingleylab/Dropbox/Work/Phenomismatch/data_NA_birdPhen.Rdata')
load('/Users/Jacob/Dropbox/Work/Phenomismatch/data_NA_birdPhen.Rdata')

data <- data_NA_birdPhen[which(data_NA_birdPhen$PRIMARY_CHECKLIST_FLAG==1),]
data$EFFORT_HRS <- as.numeric(as.character(data$EFFORT_HRS))
data <- data[which(data$EFFORT_HRS > .1 & data$EFFORT_HRS<24), ]
data <- data[which(data$YEAR > 2001), ]
data <- data[-which(data$EFFORT_DISTANCE_KM == "?"),]
data <- data[-which(data$TIME+data$EFFORT_HRS <6), ]
data <- data[-which(data$TIME > 16), ]
data <- data[which(data$LONGITUDE > -100 & data$LONGITUDE < -50 & data$LATITUDE > 26), ]
data$EFFORT_DISTANCE_KM <- as.numeric(as.character(data$EFFORT_DISTANCE_KM))
data <- data[which(data$EFFORT_DISTANCE_KM >= 0 & data$EFFORT_DISTANCE_KM < 100),]



data$cell6 <- dgGEO_to_SEQNUM(hexgrid6, data$LONGITUDE, data$LATITUDE)[[1]]
cells <- unique(data$cell6)
ncel <- length(cells)
save(cells, file = "/Users/Tingleylab/Dropbox/Work/Phenomismatch/data_NA_birdPhen_cells.Rdata")

for(i in 17:130){
  print(i)
  ttt <- data[, i]
  ttt[ttt=="X"] <- 1
  ttt <- as.numeric(as.character(ttt))
  ttt[ttt > 0] <- 1
  if(min(ttt < 0)){stop()}
  data[, i] <- ttt
}

data$sjday <- as.vector(scale(as.numeric(as.character(data$DAY))))
data$sjday2 <- data$sjday^2
data$sjday3 <- data$sjday^3
data$shr <- as.vector(scale(data$EFFORT_HRS))







Mu.day <- mean(as.numeric(as.character(data$DAY)))
sd.day <- sd(as.numeric(as.character(data$DAY)))

predictDays = (c(1:200)-Mu.day)/sd.day
newdata <- data.frame(sjday = predictDays, sjday2 = predictDays^2, sjday3 = predictDays^3, shr = 1)

load("/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/halfmax_matrix_list.Rdata")
load("/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/fit_diag.Rdata")

diagnostics_frame <- as.data.frame(matrix(data=NA, nrow=nsp*ncel*nyr, ncol=21))
names(diagnostics_frame) <- c("species", "cell", "year", "n1", "n1W", "n0", "n0i", "njd1", "njd0", "njd0i",
                              "nphen_bad", "min_effSize", "max_Rhat", "phen_effSize", "phen_Rhat",
                              "phen_mean", "phen_sd", "phen_c_loc", "phen_c_scale", "Nloglik", "Cloglik")

counter <- 0
for(i in 1:nsp){
  sdata <- data[,c("YEAR","DAY","sjday","sjday2","sjday3","shr","cell6",species_list[i])]
  names(sdata)[8] <- "detect"
  for(j in 1:nyr){
    print(paste(i,j))
    ysdata <- sdata[which(sdata$YEAR == years[j]), ]
    for(k in 1:ncel){
      counter <- counter + 1
      diagnostics_frame$species[counter] <- species_list[i]
      diagnostics_frame$year[counter] <- years[j]
      diagnostics_frame$cell[counter] <- cells[k]
      cysdata <- ysdata[which(ysdata$cell6 == cells[k]), ]
      diagnostics_frame$n1[counter] <- sum(cysdata$detect)
      diagnostics_frame$n1W[counter] <- sum(cysdata$detect*as.numeric(cysdata$DAY < 60))
      diagnostics_frame$n0[counter] <- sum(cysdata$detect == 0)
      if(diagnostics_frame$n1[counter] > 0){
        diagnostics_frame$n0i[counter] <- length(which(cysdata$detect == 0 & cysdata$sjday < min(cysdata$sjday[which(cysdata$detect==1)])))
        diagnostics_frame$njd1[counter] <- length(unique(cysdata$sjday[which(cysdata$detect == 1)]))
        diagnostics_frame$njd0i[counter] <- length(unique(cysdata$sjday[which(cysdata$detect == 0 & cysdata$sjday < min(cysdata$sjday[which(cysdata$detect==1)]))]))
      }
      
      diagnostics_frame$njd0[counter] <- length(unique(cysdata$sjday[which(cysdata$detect == 0)]))
      
      if(diagnostics_frame$n1[counter] > 29 &
         diagnostics_frame$n1W[counter] < (diagnostics_frame$n1[counter]/50) &
         diagnostics_frame$n0[counter] > 29){
        
        diagnostics_frame$min_effSize[counter] <- fit_diag[[i]][[j]][[k]]$mineffsSize
        diagnostics_frame$max_Rhat[counter] <- fit_diag[[i]][[j]][[k]]$maxRhat
        
        halfmax_posterior <- as.vector(halfmax_matrix_list[[i]][[j]][k,])
        diagnostics_frame$nphen_bad[counter] <- sum(halfmax_posterior == 1)
        
        halfmax_mcmcList <- mcmc.list(as.mcmc(halfmax_posterior[1:500]), as.mcmc(halfmax_posterior[501:1000]),
                                      as.mcmc(halfmax_posterior[1001:1500]), as.mcmc(halfmax_posterior[1501:2000]))
        diagnostics_frame$phen_effSize[counter] <- effectiveSize(halfmax_mcmcList)
        diagnostics_frame$phen_Rhat[counter] <- gelman.diag(halfmax_mcmcList)$psrf[1]
        
        halfmax_posterior2 <- halfmax_posterior[which(halfmax_posterior != 1)]
        
        diagnostics_frame$phen_mean[counter] <- mean(halfmax_posterior2)
        diagnostics_frame$phen_sd[counter] <- sd(halfmax_posterior2)
        
        normfit <- fitdistrplus::fitdist(halfmax_posterior2,"norm")
        diagnostics_frame$Nloglik[counter] <- normfit$loglik
        
        cfit <- NA
        #cfit <- tryCatch(fitdistrplus::fitdist(halfmax_posterior2,"cauchy"), error=function(e){return(NA)})
        if(!is.na(cfit)){
          diagnostics_frame$phen_c_loc[counter] <- cfit$estimate[1]
          diagnostics_frame$phen_c_scale[counter] <- cfit$estimate[2]
          diagnostics_frame$Cloglik[counter] <- cfit$loglik
        }
      }
    }
  }
}

diagnostics_frame$spCel <- paste(diagnostics_frame$species, diagnostics_frame$cell, sep="_")
birdPhenRaw <- diagnostics_frame

save(birdPhenRaw, file="/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/birdPhenRaw.Rdata")







# old code ----------------------------------------------------------------


# CHANGE LINE BELOW TO CORRECT FILE PATH
load("/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/birdPhenRaw.Rdata")
load("/Users/Jacob/Dropbox/Work/Phenomismatch/NA_birdPhen/birdPhenRaw.Rdata")

diag_frame <- birdPhenRaw[which(birdPhenRaw$n1 > 29 & birdPhenRaw$n0 > 29 &
                                  birdPhenRaw$n1W < birdPhenRaw$n1/50 &
                                  birdPhenRaw$njd0i > 29 & birdPhenRaw$njd1 > 19), ]

winter_spcels <- unique(birdPhenRaw$spCel[which(birdPhenRaw$n1W > birdPhenRaw$n1/50)])



diag_frame <- diag_frame[which(diag_frame$spCel %ni% winter_spcels), ]
length(unique(diag_frame$spCel))
for(i in 2002:2016)
{
  print(length(unique(diag_frame$spCel[which(diag_frame$year==i)])))
}
hist(diag_frame$phen_mean)

# CHANGE LINE BELOW TO CORRECT FILE PATH
load("/Users/Jacob/Dropbox/Work/Phenomismatch/data_NA_birdPhen_cells.Rdata")
ncel <- length(cells)

cells_frame <- data.frame(cell=cells, n=NA)
for(i in 1:ncel)
{
  cells_frame$n[i] <- sum(diag_frame$cell == cells[i])
}

# Which species-years are worth modeling? At a bare minimum:
#     Species with at least 3 cells in all three years from 2014-2016
#     Species-years with at least 3 cells for those species

spYs <- list()
for(i in 2014:2016)
{
  framei <- diag_frame[which(diag_frame$year == i), ]
  spis <- table(framei$species)
  spYs[[i-2013]] <- names(spis[which(spis>2)])
}

SPs <- spYs[[1]]
for(i in 2:3)
{
  SPs <- SPs[which(SPs %in% spYs[[i]])]
}

species_years_use <- matrix(data = NA, nrow=15, ncol=length(SPs))
for(i in 1:15)
{
  for(j in 1:length(SPs))
  {
    species_years_use[i,j] <- length(which(diag_frame$year == years[i] & diag_frame$species == SPs[j]))
  }
}
length(which(species_years_use > 2))

diag_frame2 <- diag_frame
for(i in nrow(diag_frame):1)
{
  if(length(which(diag_frame2$year == diag_frame2$year[i] & diag_frame2$species == diag_frame2$species[j])) < 3)
  {
    diag_frame2 <- diag_frame2[-i, ]
  }
}


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


