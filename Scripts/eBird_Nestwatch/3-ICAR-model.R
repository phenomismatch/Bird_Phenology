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

library(dggridR)
library(coda)
library(ggplot2)
library(viridis)
library(rstan)
library(doParallel)
library(foreach)
'%ni%' <- Negate('%in%')


hexgrid6 <- dgconstruct(res=6) # Construct geospatial hexagonal grid
years <- 2002:2016
nyr <- length(years)
species_list <- c("Empidonax_virescens", "Myiarchus_crinitus", "Contopus_virens", "Vireo_olivaceus",
                  "Vireo_solitarius", "Vireo_gilvus", "Vireo_flavifrons", "Catharus_fuscescens",
                  "Dumetella_carolinensis", "Setophaga_dominica", "Limnothlypis_swainsonii",
                  "Setophaga_citrina", "Geothlypis_formosa", "Parkesia_motacilla", "Parkesia_noveboracensis",
                  "Mniotilta_varia", "Setophaga_americana", "Setophaga_ruticilla",
                  "Setophaga_virens", "Setophaga_caerulescens", "Protonotaria_citrea", "Setophaga_cerulea",
                  "Seiurus_aurocapilla", "Cardellina_canadensis", "Piranga_olivacea", "Piranga_rubra",
                  "Pheucticus_ludovicianus", "Icterus_galbula", "Empidonax_traillii", "Empidonax_alnorum", #30
                  "Empidonax_minimus", "Tyrannus_tyrannus", "Vireo_bellii", "Vireo_griseus", "Tachycineta_bicolor",
                  "Stelgidopteryx_serripennis","Hirundo_rustica","Riparia_riparia","Petrochelidon_pyrrhonota",
                  "Progne_subis", "Vermivora_cyanoptera","Vermivora_chrysoptera","Oreothlypis_ruficapilla", #43
                  "Setophaga_pensylvanica", "Setophaga_petechia", "Setophaga_discolor", "Geothlypis_philadelphia",
                  "Pooecetes_gramineus", "Ammodramus_nelsoni", "Passerina_cyanea", #50
                  "Passerina_caerulea","Spiza_americana","Icterus_spurius","Dolichonyx_oryzivorus","Contopus_cooperi", #55
                  "Empidonax_flaviventris","Regulus_satrapa","Regulus_calendula","Vireo_philadelphicus", #59
                  "Troglodytes_hiemalis","Catharus_guttatus","Catharus_ustulatus","Catharus_bicknelli", #63
                  "Catharus_minimus","Setophaga_fusca","Setophaga_striata","Setophaga_tigrina","Oreothlypis_peregrina", #68
                  "Setophaga_castanea","Setophaga_palmarum","Oreothlypis_celata","Cardellina_pusilla", #72
                  "Oporornis_agilis","Setophaga_magnolia","Setophaga_coronata","Passerella_iliaca",
                  "Melospiza_lincolnii","Spizelloides_arborea","Junco_hyemalis","Zonotrichia_leucophrys", #80
                  "Zonotrichia_albicollis","Euphagus_carolinus","Ictinia_mississippiensis",
                  "Elanoides_forficatus","Archilochus_colubris","Coccyzus_erythropthalmus","Coccyzus_americanus",
                  "Antrostomus_vociferus","Antrostomus_carolinensis","Sphyrapicus_varius","Sayornis_phoebe",
                  "Bombycilla_cedrorum","Cyanocitta_cristata","Turdus_migratorius","Polioptila_caerulea",
                  "Setophaga_pinus","Peucaea_aestivalis","Spinus_tristis",
                  "Geothlypis_trichas","Icteria_virens","Mimus_polyglottos","Lanius_ludovicianus",
                  "Melospiza_georgiana","Quiscalus_quiscula","Passerculus_sandwichensis","Ammodramus_caudacutus",
                  "Spizella_passerina","Pipilo_erythrophthalmus","Troglodytes_aedon","Sturnella_magna",
                  "Sialia_sialis","Cistothorus_palustris","Corvus_ossifragus","Agelaius_phoeniceus")
nsp <- length(species_list)

# CHANGE LINE BELOW TO CORRECT FILE PATH
load("/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/birdPhenRaw.Rdata")
load("/Users/Jacob/Dropbox/Work/Phenomismatch/NA_birdPhen/birdPhenRaw.Rdata")

diag_frame <- birdPhenRaw[which(birdPhenRaw$n1 > 29 & birdPhenRaw$n0 > 29 &
                                  birdPhenRaw$n1W < birdPhenRaw$n1/50 &
                                  birdPhenRaw$njd0i > 29 & birdPhenRaw$njd1 > 19), ]
winter_spcels <- unique(birdPhenRaw$spCel[which(birdPhenRaw$n1W > birdPhenRaw$n1/50)])

diag_frame <- diag_frame[which(diag_frame$spCel %ni% winter_spcels), ]
length(unique(diag_frame$spCel))
for(i in 2002:2016){print(length(unique(diag_frame$spCel[which(diag_frame$year==i)])))}
hist(diag_frame$phen_mean)

# CHANGE LINE BELOW TO CORRECT FILE PATH
load("/Users/Jacob/Dropbox/Work/Phenomismatch/data_NA_birdPhen_cells.Rdata")
ncel <- length(cells)

cells_frame <- data.frame(cell=cells, n=NA)
for(i in 1:ncel){
  cells_frame$n[i] <- sum(diag_frame$cell == cells[i])
}

# Which species-years are worth modeling? At a bare minimum:
#     Species with at least 3 cells in all three years from 2014-2016
#     Species-years with at least 3 cells for those species
spYs <- list()
for(i in 2014:2016){
  framei <- diag_frame[which(diag_frame$year == i), ]
  spis <- table(framei$species)
  spYs[[i-2013]] <- names(spis[which(spis>2)])
}
SPs <- spYs[[1]]
for(i in 2:3){
  SPs <- SPs[which(SPs %in% spYs[[i]])]
}
species_years_use <- matrix(data = NA, nrow=15, ncol=length(SPs))
for(i in 1:15){
  for(j in 1:length(SPs)){
    species_years_use[i,j] <- length(which(diag_frame$year == years[i] & diag_frame$species == SPs[j]))
  }
}
length(which(species_years_use > 2))

diag_frame2 <- diag_frame
for(i in nrow(diag_frame):1){
  if(length(which(diag_frame2$year == diag_frame2$year[i] & diag_frame2$species == diag_frame2$species[j])) < 3){
    diag_frame2 <- diag_frame2[-i, ]
  }
}


# ICAR models

# General data preparation

ncel <- length(cells)
cellcenters <- dgSEQNUM_to_GEO(hexgrid6, as.numeric(cells))
adjacency_matrix <- matrix(data=NA, nrow = length(cells), ncol=length(cells))
for(i in 1:length(cells)){
  for(j in i:length(cells)){
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
for(i in 1:ncores){
  Ldatalist[[i]] <- list()
  for(j in 1:ceiling(nfits/ncores)){
    counter <- counter + 1
    if(counter <= nfits){
      Ldatalist[[i]][[j]] <- 
        diag_frame2[which(diag_frame2$spyr == spyrs[counter]), ]
    }
    
  }
}

fit_funct <- function(datalist, model){
  output_list <- list()
  for(j in 1:length(datalist)){
    dd <- datalist[[j]]
    SY <- data.frame(cell = cells, halfmax = NA, sd = NA)
    for(i in 1:length(cells)){
      d1 <- which(dd$cell == as.integer(cells[i]))
      if(length(d1) > 0){
        SY$halfmax[i] <- dd$phen_mean[d1]
        SY$sd[i] <- dd$phen_sd[d1]
      }else{
        SY$sd[i] <- .01
      }
    }
    
    SYL <- list(N = length(cells), N_edges = nrow(ninds), node1 = ninds[,1], node2 = ninds[,2],
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


