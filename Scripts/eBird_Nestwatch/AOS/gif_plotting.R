######################
# Plots for post IAR gif
######################

# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'



# db/hm query dir ------------------------------------------------------------

IAR_in_dir <- 'IAR_input_2019-05-03'
IAR_out_dir <- 'IAR_output_2019-06-10'



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

args <- commandArgs(trailingOnly = TRUE)
#args <- as.character('Icterus_spurius')
#args <- as.character('Catharus_minimus')
#args <- as.character('Empidonax_virescens')
#args <- as.character('Vireo_olivaceus')


#species for which one cell had to be dropped, bc it did not border any others
#args <- as.character('Cistothorus_palustris')
#args <- as.character('Setophaga_dominica')
#args <- as.character('Melospiza_lincolnii')
#args <- as.character('Zonotrichia_leucophrys')
#args <- as.character('Bombycilla_cedrorum')
#args <- as.character('Coccyzus_americanus')
#args <- as.character('Vireo_bellii')

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
sigma_y_in <- matrix(nrow = ncell, ncol = nyr)
y_obs_in <- matrix(nrow = ncell, ncol = nyr)
y_PPC <- rep(NA, ncell * nyr)

#number of observation and NAs for each year
len_y_obs_in <- rep(NA, nyr)
len_y_mis_in <- rep(NA, nyr)

#indices for observed and missing
ii_obs_in <- matrix(NA, nrow = ncell, ncol = nyr)
ii_mis_in <- matrix(NA, nrow = ncell, ncol = nyr)

#counter to fill y_PPC
counter <- 1
for (j in 1:nyr)
{
  #j <- 16
  temp_yr <- dplyr::filter(f_out, year == years[j])
  
  #don't need to manipulate position of sigmas
  sigma_y_in[,j] <- temp_yr$HM_sd
  
  for (n in 1:ncell)
  {
    #n <- 1
    #matrix with observed values with NAs
    y_PPC[counter] <- temp_yr$HM_mean[n]
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
    ii_obs_in[,j] <- t_obs_in
    ii_mis_in[,j] <- t_mis_in
    y_obs_in[,j] <- t_y_obs_in
  } else {
    #no NAs to end (no mimssing values)
    y_obs_in[,j] <- no_na
    ii_mis_in[,j] <- which(!is.na(temp_yr$HM_mean))
    y_obs_in[,j] <- which(is.na(temp_yr$HM_mean))
  }
  
  #length of data/miss for each year
  len_y_obs_in[j] <- length(no_na)
  len_y_mis_in[j] <- ncell - length(no_na)
}


#fill 0 where NA in y_obs - Stan does not like NA and zeros are not being used to estimate any param (y_obs is used to fill y)
y_obs_in[which(is.na(y_obs_in), arr.ind = TRUE)] <- 0
#see script here showing this value fo sigma_y_in has no impact on y_true: Archive/test_sigma_insert.R
sigma_y_in[which(is.na(sigma_y_in), arr.ind = TRUE)] <- 0.1
ii_obs_in[which(is.na(ii_obs_in), arr.ind = TRUE)] <- 0
ii_mis_in[which(is.na(ii_mis_in), arr.ind = TRUE)] <- 0


#create data list for Stan
DATA <- list(J = nyr,
             cells = cells,
             N = ncell, 
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



# read in IAR fit ---------------------------------------------------------

setwd("~/Desktop/Bird_Phenology_Offline/Data/Processed/IAR_output_2019-05-26/")
fit <- readRDS(paste0(args, '-2019-05-26-iar-stan_output.rds'))

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



setwd('~/Desktop/gif_maps')
#loop plots for each year
#for (i in 1:length(years))
#{
i <- 1

#filter data for year[i]
f_out_filt <- dplyr::filter(f_out, year == years[i])

#merge hex spatial data with HM data
to_plt <- dplyr::inner_join(f_out_filt, cell_grid, by = 'cell')
to_plt2 <- dplyr::inner_join(to_plt, ll_df, by = 'cell')

#post-IAR
t_med_fit <- med_fit[,i]
t_sd_fit <- sd_fit[,i]

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
  # annotate('text', x = to_plt2_post$lon_deg, y = to_plt2_post$lat_deg + 0.5,
  #          label = round(to_plt2_post$med_mu, digits = 0), col = 'black', alpha = 0.2,
  #          size = 4) +
  # annotate('text', x = to_plt2_post$lon_deg, y = to_plt2_post$lat_deg - 0.5,
  #          label = round(to_plt2_post$sd_mu, digits = 0), col = 'white', alpha = 0.3,
  #          size = 3) +
  # ggtitle(paste0(f_out_filt$species[1], ' - ', f_out_filt$year[1], ' - Post-IAR')) +
  theme_bw() +
  xlab('') +
  ylab('') +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none')


ggsave(plot = p_post,
       filename = paste0(f_out_filt$species[1], '_', f_out_filt$year[1], '-post_IAR.pdf'))