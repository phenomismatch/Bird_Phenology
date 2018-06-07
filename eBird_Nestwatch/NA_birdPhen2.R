# This script takes the output of NA_birdPhen1.R, which gives posterior distributions for the half-max
# parameter by species-cell-years, and uses these as the basis for spatial models of the half-max
# parameter.  In a nutshell, the spatial models operate by assuming that the posterior distributions
# for the half-max parameter from NA_birdPhen1.R are normally distributed. **I have spot-checked this
# assumption and it seems reasonable, but have not exhaustively evaluated it**

# Then, the models assume that the posterior means are in fact drawn from normal distributions
# parameterized by a latent true half-max and the variance in the posterior estimate.

# The latent true half-maxima (LTHMs) are given a prior with a spatial component, corresponding either to 
# an intrinsic conditional autoregressive (ICAR) model, or to a Gaussian process (GP) model. For now,
# both models are fit across cells within a species-year; there is no data sharing across years or
# species. For the ICAR model, I assume both spatial and nonspatial variance in the LTHMs. Thus, the
# LTHM from a cell is drawn from a normal distribution centered around some latent spatial mean, which
# itself is drawn from a normal distribution centered on the average of the spatial means of the cell's
# neighbors. So the difference between the spatial variance and the nonspatial variance is that the 
# spatial residual propagates to spatially to influence the spatial mean of the neighbors, whereas the 
# nonspatial residual does not do this. In practice, for the few models that I have fit so far, I find
# that the nonspatial variance is small, difficult to estimate, and leads to diagnostic warnings from 
# Stan. For the time being, I have given that parameter a strongly informative prior centered slightly 
# higher than the typical (suspect, due to diagnostic warnings) values that I saw first.


# NOTE TO SELF: DEAL WITH CHOICE OF WHICH CELLS TO INCLUDE BASED ON HALFMAX UNCERTAINTY
library(dggridR)
library(coda)
library(ggplot2)
library(viridis)
library(rstan)
'%ni%' <- Negate('%in%')

# Function to help evaluate stan fit diagnostics
pairs_stan <- function(chain, stan_model, pars) {
  energy <- as.matrix(sapply(get_sampler_params(stan_model, inc_warmup = F), 
                             function(x) x[,"energy__"]))
  pars <- extract(stan_model, pars = pars, permuted = F)
  df <- data.frame(energy[,chain], pars[,chain,])
  names(df)[1] <- "energy"
  GGally::ggpairs(df, title = paste0("Chain", chain), 
                  lower = list(continuous = GGally::wrap("points", alpha = 0.2)))                    
}

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

# load('/Users/Tingleylab/Dropbox/Work/Phenomismatch/data_NA_birdPhen.Rdata')
# 
# data <- data_NA_birdPhen[which(data_NA_birdPhen$PRIMARY_CHECKLIST_FLAG==1),]
# data$EFFORT_HRS <- as.numeric(as.character(data$EFFORT_HRS))
# data <- data[which(data$EFFORT_HRS > .1 & data$EFFORT_HRS<24), ]
# data <- data[which(data$YEAR > 2001), ]
# data <- data[-which(data$EFFORT_DISTANCE_KM == "?"),]
# data <- data[-which(data$TIME+data$EFFORT_HRS <6), ]
# data <- data[-which(data$TIME > 16), ]
# data <- data[which(data$LONGITUDE > -100 & data$LONGITUDE < -50 & data$LATITUDE > 26), ]
# data$EFFORT_DISTANCE_KM <- as.numeric(as.character(data$EFFORT_DISTANCE_KM))
# data <- data[which(data$EFFORT_DISTANCE_KM >= 0 & data$EFFORT_DISTANCE_KM < 100),]
# data$cell6 <- dgGEO_to_SEQNUM(hexgrid6, data$LONGITUDE, data$LATITUDE)[[1]]
# cells <- unique(data$cell6)
# ncel <- length(cells)
# save(cells, file = "/Users/Tingleylab/Dropbox/Work/Phenomismatch/data_NA_birdPhen_cells.Rdata")
# 
# for(i in 17:130){
#   print(i)
#   ttt <- data[, i]
#   ttt[ttt=="X"] <- 1
#   ttt <- as.numeric(as.character(ttt))
#   ttt[ttt > 0] <- 1
#   if(min(ttt < 0)){stop()}
#   data[, i] <- ttt
# }
# 
# data$sjday <- as.vector(scale(as.numeric(as.character(data$DAY))))
# data$sjday2 <- data$sjday^2
# data$sjday3 <- data$sjday^3
# data$shr <- as.vector(scale(data$EFFORT_HRS))
# 
# Mu.day <- mean(as.numeric(as.character(data$DAY)))
# sd.day <- sd(as.numeric(as.character(data$DAY)))
# 
# predictDays = (c(1:200)-Mu.day)/sd.day
# newdata <- data.frame(sjday = predictDays, sjday2 = predictDays^2, sjday3 = predictDays^3, shr = 1)
# 
# load("/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/halfmax_matrix_list.Rdata")
# load("/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/fit_diag.Rdata")
# 
# diagnostics_frame <- as.data.frame(matrix(data=NA, nrow=nsp*ncel*nyr, ncol=21))
# names(diagnostics_frame) <- c("species", "cell", "year", "n1", "n1W", "n0", "n0i", "njd1", "njd0", "njd0i",
#                               "nphen_bad", "min_effSize", "max_Rhat", "phen_effSize", "phen_Rhat",
#                               "phen_mean", "phen_sd", "phen_c_loc", "phen_c_scale", "Nloglik", "Cloglik")
# 
# counter <- 0
# for(i in 1:nsp){
#   sdata <- data[,c("YEAR","DAY","sjday","sjday2","sjday3","shr","cell6",species_list[i])]
#   names(sdata)[8] <- "detect"
#   for(j in 1:nyr){
#     print(paste(i,j))
#     ysdata <- sdata[which(sdata$YEAR == years[j]), ]
#     for(k in 1:ncel){
#       counter <- counter + 1
#       diagnostics_frame$species[counter] <- species_list[i]
#       diagnostics_frame$year[counter] <- years[j]
#       diagnostics_frame$cell[counter] <- cells[k]
#       cysdata <- ysdata[which(ysdata$cell6 == cells[k]), ]
#       diagnostics_frame$n1[counter] <- sum(cysdata$detect)
#       diagnostics_frame$n1W[counter] <- sum(cysdata$detect*as.numeric(cysdata$DAY < 60))
#       diagnostics_frame$n0[counter] <- sum(cysdata$detect == 0)
#       if(diagnostics_frame$n1[counter] > 0){
#         diagnostics_frame$n0i[counter] <- length(which(cysdata$detect == 0 & cysdata$sjday < min(cysdata$sjday[which(cysdata$detect==1)])))
#         diagnostics_frame$njd1[counter] <- length(unique(cysdata$sjday[which(cysdata$detect == 1)]))
#         diagnostics_frame$njd0i[counter] <- length(unique(cysdata$sjday[which(cysdata$detect == 0 & cysdata$sjday < min(cysdata$sjday[which(cysdata$detect==1)]))]))
#       }
#       
#       diagnostics_frame$njd0[counter] <- length(unique(cysdata$sjday[which(cysdata$detect == 0)]))
#       
#       if(diagnostics_frame$n1[counter] > 29 & 
#          diagnostics_frame$n1W[counter] < (diagnostics_frame$n1[counter]/50) & 
#          diagnostics_frame$n0[counter] > 29){
#         
#         diagnostics_frame$min_effSize[counter] <- fit_diag[[i]][[j]][[k]]$mineffsSize
#         diagnostics_frame$max_Rhat[counter] <- fit_diag[[i]][[j]][[k]]$maxRhat
#         
#         halfmax_posterior <- as.vector(halfmax_matrix_list[[i]][[j]][k,])
#         diagnostics_frame$nphen_bad[counter] <- sum(halfmax_posterior == 1)
#         
#         halfmax_mcmcList <- mcmc.list(as.mcmc(halfmax_posterior[1:500]), as.mcmc(halfmax_posterior[501:1000]),
#                                         as.mcmc(halfmax_posterior[1001:1500]), as.mcmc(halfmax_posterior[1501:2000]))
#         diagnostics_frame$phen_effSize[counter] <- effectiveSize(halfmax_mcmcList)
#         diagnostics_frame$phen_Rhat[counter] <- gelman.diag(halfmax_mcmcList)$psrf[1]
#         
#         halfmax_posterior2 <- halfmax_posterior[which(halfmax_posterior != 1)]
#         
#         diagnostics_frame$phen_mean[counter] <- mean(halfmax_posterior2)
#         diagnostics_frame$phen_sd[counter] <- sd(halfmax_posterior2)
#         
#         normfit <- fitdistrplus::fitdist(halfmax_posterior2,"norm")
#         diagnostics_frame$Nloglik[counter] <- normfit$loglik
#         
#         cfit <- NA
#         #cfit <- tryCatch(fitdistrplus::fitdist(halfmax_posterior2,"cauchy"), error=function(e){return(NA)})
#         if(!is.na(cfit)){
#           diagnostics_frame$phen_c_loc[counter] <- cfit$estimate[1]
#           diagnostics_frame$phen_c_scale[counter] <- cfit$estimate[2]
#           diagnostics_frame$Cloglik[counter] <- cfit$loglik
#         }
#       }
#     }
#   }
# }
# 
# diagnostics_frame$spCel <- paste(diagnostics_frame$species, diagnostics_frame$cell, sep="_")
# birdPhenRaw <- diagnostics_frame
# 
# save(birdPhenRaw, file="/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/birdPhenRaw.Rdata")

load("/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/birdPhenRaw.Rdata")

diag_frame <- birdPhenRaw[which(birdPhenRaw$n1 > 29 & birdPhenRaw$n0 > 29 &
                                        birdPhenRaw$n1W < birdPhenRaw$n1/50 &
                                        birdPhenRaw$njd0i > 29 & birdPhenRaw$njd1 > 19), ]
winter_spcels <- unique(birdPhenRaw$spCel[which(birdPhenRaw$n1W > birdPhenRaw$n1/50)])

diag_frame <- diag_frame[which(diag_frame$spCel %ni% winter_spcels), ]
length(unique(diag_frame$spCel))
for(i in 2002:2016){print(length(unique(diag_frame$spCel[which(diag_frame$year==i)])))}
hist(diag_frame$phen_mean)

cells_frame <- data.frame(cell=cells, n=NA)
for(i in 1:ncel){
  cells_frame$n[i] <- sum(diag_frame$cell == cells[i])
}

# visualize how many species-years are available per cell
grid <- dgcellstogrid(hexgrid6, cells_frame$cell, frame=TRUE,wrapcells=TRUE)
grid <- merge(grid, cells_frame, by = "cell")
grid$logn <- log(grid$n)

usa <- map_data("state") 

p<- ggplot() + 
  geom_polygon(data=usa, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid,      aes(x=long, y=lat, group=group, fill=logn), alpha=0.9)    +
  geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  #  geom_point  (aes(x=cellcenters$lon_deg, y=cellcenters$lat_deg)) +
  scale_fill_gradient()
p+coord_map("ortho", orientation = c(39, -80, 0))+
  xlab('')+ylab('')+
  theme(axis.ticks.x=element_blank())+
  theme(axis.ticks.y=element_blank())+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y=element_blank())

grid <- grid[which(grid$n > 8), ]
p<- ggplot() + 
  geom_polygon(data=usa, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid,      aes(x=long, y=lat, group=group, fill=logn), alpha=0.9)    +
  geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  #  geom_point  (aes(x=cellcenters$lon_deg, y=cellcenters$lat_deg)) +
  scale_fill_gradient()
p+coord_map("ortho", orientation = c(39, -80, 0))+
  xlab('')+ylab('')+
  theme(axis.ticks.x=element_blank())+
  theme(axis.ticks.y=element_blank())+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y=element_blank())

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

# ICAR models

# General data preparation
cells <- unique(grid$cell)
save(cells, file = "/Users/TingleyLab/Dropbox/Work/Phenomismatch/data_NA_birdPhen_cells.Rdata")
load("/Users/TingleyLab/Dropbox/Work/Phenomismatch/data_NA_birdPhen_cells.Rdata")
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

# Stan code
# 1st attempt, tested on REVI for 2005 and 2016 yields warning about low estimated fraction of Bayesian
# missing information. A pairs plot diagnoses the problem as 
stan_ICAR <- '
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
}'

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

# Data preparation: Red-eyed Vireo in 2016
frame_16 <- diag_frame[which(diag_frame$year == 2016), ]
REVI_16 <- frame_16[which(frame_16$species == "Vireo_olivaceus"), ]
REVI_data2016 <- data.frame(cell = cells, halfmax = NA, sd = NA)

for(i in 1:length(cells)){
  d1 <- which(diag_frame$species == "Vireo_olivaceus" & diag_frame$year == 2016 & diag_frame$cell == cells[i])
  if(length(d1) > 0){
    if(length(d1) != 1){
      print(i)
      stop("length d1 != 1")
    }
    REVI_data2016$halfmax[i] <- diag_frame$phen_mean[d1]
    REVI_data2016$sd[i] <- diag_frame$phen_sd[d1]
  }else{
    REVI_data2016$sd[i] <- .01
  }
}

REVI2016 <- list(N = length(cells), N_edges = nrow(ninds), node1 = ninds[,1], node2 = ninds[,2],
                 obs = REVI_data2016$halfmax[which(!is.na(REVI_data2016$halfmax))], 
                 ii_obs = which(!is.na(REVI_data2016$halfmax)),
                 ii_mis = which(is.na(REVI_data2016$halfmax)),
                 N_obs = sum(!is.na(REVI_data2016$halfmax)),
                 N_mis = sum(is.na(REVI_data2016$halfmax)),
                 sds = REVI_data2016$sd)



REVI2016_fit1 <- stan(
  model_code = stan_ICAR,  # Stan program
  data = REVI2016,    # named list of data
  chains = 3,             # number of Markov chains
  iter = 6000,            # total number of iterations per chain
  cores = 3,              # number of cores
  control = list(max_treedepth = 20, adapt_delta = .9) # modified control parameters based on warnings;
  # see http://mc-stan.org/misc/warnings.html
)
save(REVI2016_fit1, file="/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/REVI2016_fit1.Rdata")


REVI2016_fit2 <- stan(
  model_code = stan_ICAR_infoPrior,  # Stan program
  data = REVI2016,    # named list of data
  chains = 3,             # number of Markov chains
  iter = 6000,            # total number of iterations per chain
  cores = 3,              # number of cores
  control = list(max_treedepth = 20, adapt_delta = .9)
)
save(REVI2016_fit2, file="/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/REVI2016_fit2.Rdata")

load("/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/REVI2016_fit2.Rdata")
pairs_stan(chain = 1, stan_model = REVI2016_fit2, pars=c("sigma_theta", "sigma_phi", "beta0", "lp__"))
check_energy(REVI2016_fit2)
summary(REVI2016_fit2)



REVI2016_fit3 <- stan(
  model_code = stan_ICAR_no_nonspatial,  # Stan program
  data = REVI2016,    # named list of data
  chains = 3,             # number of Markov chains
  iter = 6000,            # total number of iterations per chain
  cores = 3,              # number of cores
  control = list(max_treedepth = 20, adapt_delta = .9)
)
save(REVI2016_fit3, file="/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/REVI2016_fit3.Rdata")

load("/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/REVI2016_fit3.Rdata")
pairs_stan(chain = 1, stan_model = REVI2016_fit3, pars=c("sigma_theta", "sigma_phi", "beta0", "lp__"))
check_energy(REVI2016_fit3)
summary(REVI2016_fit3)




### Models for year 2005
REVI_data <- data.frame(cell = cells, halfmax = NA, sd = NA)
for(i in 1:length(cells)){
  d1 <- which(diag_frame$species == "Vireo_olivaceus" & diag_frame$year == 2005 & diag_frame$cell == cells[i])
  if(length(d1) > 0){
    if(length(d1) != 1){
      print(i)
      stop("length d1 != 1")
    }
    REVI_data$halfmax[i] <- diag_frame$phen_mean[d1]
    REVI_data$sd[i] <- diag_frame$phen_sd[d1]
  }else{
    REVI_data$sd[i] <- .01
  }
}
REVI2005 <- list(N = length(cells), N_edges = nrow(ninds), node1 = ninds[,1], node2 = ninds[,2],
                 obs = REVI_data$halfmax[which(!is.na(REVI_data$halfmax))], 
                 ii_obs = which(!is.na(REVI_data$halfmax)),
                 ii_mis = which(is.na(REVI_data$halfmax)),
                 N_obs = sum(!is.na(REVI_data$halfmax)),
                 N_mis = sum(is.na(REVI_data$halfmax)),
                 sds = REVI_data$sd)


REVI2005_fit2 <- stan(
  model_code = stan_ICAR,  # Stan program
  data = REVI2005,    # named list of data
  chains = 3,             # number of Markov chains
  iter = 6000,            # total number of iterations per chain
  cores = 3,              # number of cores
  control = list(max_treedepth = 20, adapt_delta = .9)
)

pairs_stan(chain = 1, stan_model = REVI2005_fit2, pars=c("sigma_theta", "sigma_phi", "beta0", "lp__"))
check_energy(REVI2005_fit2)

save(REVI2005_fit2, file="/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/REVI2005_fit2.Rdata")
load("/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/REVI2005_fit2.Rdata")
summary(REVI2005_fit2)






######### ATTEMPTING A GAUSSIAN PROCESS MODEL (AS WRITTEN THIS DOESN'T REALLY WORK)

distance_matrix <- matrix(data=NA, nrow = length(cells), ncol=length(cells))
for(i in 1:length(cells)){
  for(j in 1:length(cells)){
    distance_matrix[i,j] <- geosphere::distm(c(cellcenters$lon_deg[i], cellcenters$lat_deg[i]),
                              c(cellcenters$lon_deg[j], cellcenters$lat_deg[j]))
  }
}
distance_matrix <- distance_matrix/10^6

# Stan code
stanCode3 <- '
data {
int<lower=0> N; // number of cells (= number of observations, including NAs)
int<lower=0> N_obs; // number of non-missing
int<lower=0> N_mis; // number missing
real<lower=-20, upper=20> obs[N_obs];    // observed data (excluding NAs)
int<lower = 1, upper = N_obs + N_mis> ii_obs[N_obs];
int<lower = 1, upper = N_obs + N_mis> ii_mis[N_mis];
real<lower=0> sds[N];                 // sds for ALL data (observed and unobserved)
matrix[N,N] dist; // distance matrix
}

parameters {
real<lower=-20, upper=20> y_mis[N_mis];         // missing data
real beta;                // intercept

real<lower=0> sigma_sq;   // sd of heterogeneous effects
real<lower=0> phi;     // sd of spatial effects
vector[N] latent;      // latent true halfmax values
}

transformed parameters {
real<lower=-20, upper=20> y[N];
matrix[N,N] Sigma;
matrix[N,N] L;
vector[N] mu;

y[ii_obs] = obs;
y[ii_mis] = y_mis;

for(i in 1:(N-1)){
for(j in (i+1):N){
Sigma[i,j] = exp((-1)*phi*dist[i,j]);
Sigma[j,i] = Sigma[i,j];
}
}
for(i in 1:N){
 Sigma[i,i] = sigma_sq;
 mu[i] = beta;
}

L = cholesky_decompose(Sigma);

}


model {
y ~ normal(latent, sds);
sigma_sq ~ normal(0,50);
phi ~ normal(0,50);
latent ~ multi_normal_cholesky(mu, L);
beta ~ normal(1000, 100);
}'

REVI2005 <- list(N = length(cells), 
                 obs = (REVI_data$halfmax[which(!is.na(REVI_data$halfmax))]-mean(REVI_data$halfmax, na.rm = T))/sd(REVI_data$halfmax, na.rm = T), 
                 ii_obs = which(!is.na(REVI_data$halfmax)),
                 ii_mis = which(is.na(REVI_data$halfmax)),
                 N_obs = sum(!is.na(REVI_data$halfmax)),
                 N_mis = sum(is.na(REVI_data$halfmax)),
                 sds = REVI_data$sd,
                 dist=distance_matrix)

GP_REVI2005 <- stan(
  model_code = stanCode,  # Stan program
  data = REVI2005,    # named list of data
  chains = 3,             # number of Markov chains
  iter = 6000,            # total number of iterations per chain
  cores = 3,              # number of cores
  control = list(max_treedepth = 20, adapt_delta = .9)
)

pairs_stan(chain = 1, stan_model = GP_REVI2005, pars=c("sigma_sq", "phi", "beta", "latent[1]", "lp__"))
check_energy(GP_REVI2005)

save(GP_REVI2005, file="/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/GP_REVI2005.Rdata")

summary(GP_REVI2005)




