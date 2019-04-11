######################
# 2 - halfmax model arrival
#
# Fit model (GAM logistic regression) to eBird data to get half-max parameter (bird arrival) for each species-cell-year
#
# y ~ bernouli(p)
# logit(p) = s(day)
#
# Halfmax is derived from model fit
#
# Species name should be given as an arg to this script. The model will then be fit to that species only.
# Runtime: Up to 7 days on Xanadu
######################  


# top-level dir -----------------------------------------------------------

#dir <- '~/Google_Drive/R/'

#Xanadu
dir <- '/UCHC/LABS/Tingley/phenomismatch/'


# db query dir ------------------------------------------------------------

db_dir <- 'eBird_query_2019-03-29'
RUN_DATE <- '2019-03-29'



# model settings ----------------------------------------------------------

#number of iterations each model should be run
ITER <- 1500
#ITER <- 10
CHAINS <- 4



# runtime -----------------------------------------------------------------

tt <- proc.time()



# Load packages -----------------------------------------------------------

library(rstanarm)
library(rstan)
library(dplyr)
library(dggridR)
library(sp)
library(raster)



# Set wd ------------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/'))



# Get args fed to script --------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
#args <- 'Empidonax_virescens'
#args <- 'Mimus_polyglottos'
#args <- 'Agelaius_phoeniceus'
#args <- 'Vireo_olivaceus'
#args <- 'Catharus_fuscescens'
#args <- 'Ammospiza_nelsoni'


# import processed data ---------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', db_dir))

#import data for species
spdata <- readRDS(paste0('ebird_NA_phen_proc_', args, '.rds'))



# create grid -------------------------------------------------------------

hexgrid6 <- dggridR::dgconstruct(res = 6)

#get boundaries of all cells over earth
setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/shapefiles/'))

#dggridR::dgearthgrid(hexgrid6, savegrid = 'global_hex.shp')
#read in grid
hge <- rgdal::readOGR('global_hex.shp', verbose = FALSE)



# filter cells by range  ---------------------------------------------------

'%ni%' <- Negate('%in%')

#reference key for species synonyms
setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/metadata/'))
sp_key <- read.csv('species_filenames_key.csv')

#change dir to shp files
setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/shapefiles/'))

#filter by breeding/migration cells
#match species name to shp file name
g_ind <- grep(args, sp_key$file_names_2016)

#check for synonyms if there are no matches
if (length(g_ind) == 0)
{
  g_ind2 <- grep(args, sp_key$BL_Checklist_name)
} else {
  g_ind2 <- g_ind
}

#get filename and read in
fname <- as.character(sp_key[g_ind2,]$filenames[grep('.shp', sp_key[g_ind2, 'filenames'])])
sp_rng <- rgdal::readOGR(fname, verbose = FALSE)
#crop to area of interest
sp_rng2 <- raster::crop(sp_rng, raster::extent(-100, -50, 26, 90))

#filter by breeding (2) and migration (4) range - need to convert spdf to sp
nrng <- sp_rng2[which(sp_rng2$SEASONAL == 2 | sp_rng2$SEASONAL == 4),]

#filter by resident (1) and non-breeding (3) to exclude hex cells that contain 2/4 and 1/3
nrng_rm <- sp_rng2[which(sp_rng2$SEASONAL == 1 | sp_rng2$SEASONAL == 3),]

#remove unneeded objects
rm(sp_rng)
rm(sp_rng2)
rm(fname)


#if there is a legitimate range
if (NROW(nrng@data) > 0)
{
  #good cells
  nrng_sp <- sp::SpatialPolygons(nrng@polygons)
  sp::proj4string(nrng_sp) <- sp::CRS(sp::proj4string(nrng))
  ptsreg <- sp::spsample(nrng, 50000, type = "regular")
  br_mig_cells <- as.numeric(which(!is.na(sp::over(hge, ptsreg))))
  
  #bad cells
  if (length(nrng_rm) > 0)
  {
    nrng_rm_sp <- sp::SpatialPolygons(nrng_rm@polygons)
    sp::proj4string(nrng_rm_sp) <- sp::CRS(sp::proj4string(nrng_rm))
    ptsreg_rm <- sp::spsample(nrng_rm_sp, 50000, type = "regular")
    res_ovr_cells <- as.numeric(which(!is.na(sp::over(hge, ptsreg_rm))))

    #remove cells that appear in resident and overwinter range that also appear in breeding range
    cell_mrg <- c(br_mig_cells, res_ovr_cells)
    to_rm <- cell_mrg[duplicated(cell_mrg)]
      
    rm(nrng_rm)
    rm(nrng_rm_sp)  
    rm(ptsreg_rm)
    rm(res_ovr_cells)
    
  } else {
    cell_mrg <- br_mig_cells
    to_rm <- NULL
  }
  
  #remove unneeded objects
  rm(hge)
  rm(nrng)
  rm(nrng_sp)
  rm(ptsreg)
  
  if (length(to_rm) > 0)
  {
    overlap_cells <- br_mig_cells[-which(br_mig_cells %in% to_rm)]
  } else {
    overlap_cells <- br_mig_cells
  }
  
  #remove unneeded objects
  rm(br_mig_cells)
  rm(cell_mrg)
  rm(to_rm)
  
  #get cell centers
  cell_centers <- dggridR::dgSEQNUM_to_GEO(hexgrid6, overlap_cells)
  cc_df <- data.frame(cell = overlap_cells, lon = cell_centers$lon_deg, 
                      lat = cell_centers$lat_deg)
  
  #remove unneeded objects
  rm(hexgrid6)
  rm(cell_centers)
  rm(overlap_cells)
  
  #cells only within the range that ebird surveys were filtered to
  n_cc_df <- cc_df[which(cc_df$lon > -100 & cc_df$lon < -50 & cc_df$lat > 26),]
  cells <- n_cc_df$cell
  
  #retain rows that match selected cells
  spdata2 <- spdata[which(spdata$cell %in% cells),]
  
  #remove unneeded objects
  rm(cc_df)
  rm(n_cc_df)
  rm(spdata)
  
} else {
  #write blank .rds file
  t_mat <- matrix(data = NA, nrow = 1, ncol = ((ITER/2)*CHAINS))
  colnames(t_mat) <- paste0('iter_', 1:((ITER/2)*CHAINS))
  halfmax_df <- data.frame(species = args, 
                           year = NA, 
                           cell = NA, 
                           max_Rhat = NA,
                           min_neff = NA,
                           num_diverge = NA,
                           num_tree = NA,
                           num_BFMI = NA,
                           delta = NA,
                           tree_depth = NA,
                           n1 = NA,
                           n1W = NA,
                           n0 = NA,
                           n0i = NA,
                           njd1 = NA,
                           njd0 = NA,
                           njd0i = NA,
                           t_mat)

  #save to rds object
  setwd(paste0(dir, '/Bird_Phenology/Data/Processed/halfmax_species_', RUN_DATE))
  saveRDS(halfmax_df, file = paste0('halfmax_df_arrival_', args, '.rds'))
  
  stop('Range not suitable for modeling!')
}



# process data ------------------------------------------------------------

ncell <- length(cells)

years <- min(spdata2$year):max(spdata2$year)
nyr <- length(years)



# fit model ---------------------------------------------------------

t_mat <- matrix(data = NA, nrow = ncell*nyr, ncol = ((ITER/2)*CHAINS))
colnames(t_mat) <- paste0('iter_', 1:((ITER/2)*CHAINS))
halfmax_df <- data.frame(species = args, 
                         year = rep(years, each = ncell), 
                         cell = rep(cells, nyr), 
                         max_Rhat = NA,
                         min_neff = NA,
                         num_diverge = NA,
                         num_tree = NA,
                         num_BFMI = NA,
                         delta = NA,
                         tree_depth = NA,
                         n1 = NA,
                         n1W = NA,
                         n0 = NA,
                         n0i = NA,
                         njd1 = NA,
                         njd0 = NA,
                         njd0i = NA,
                         t_mat)


#create dir for figs if doesn't exist
ifelse(!dir.exists(paste0(dir, 'Bird_Phenology/Figures/halfmax/arrival_', RUN_DATE)), 
       dir.create(paste0(dir, 'Bird_Phenology/Figures/halfmax/arrival_', RUN_DATE)), 
       FALSE)

setwd(paste0(dir, 'Bird_Phenology/Figures/halfmax/arrival_', RUN_DATE))


#loop through each species, year, cell and extract half-max parameter
counter <- 1
for (j in 1:nyr)
{
  #j <- 1
  yspdata <- spdata2[which(spdata2$year == years[j]), ]
  
  for (k in 1:ncell)
  {
    #k <- 1
    print(paste0('species: ', args, ', year: ', j, ', cell: ', k))
    
    cyspdata <- yspdata[which(yspdata$cell == cells[k]), ]
    
    #number of surveys where species was detected
    n1 <- sum(cyspdata$detect)
    #number of surveys where species was not detected
    n0 <- sum(cyspdata$detect == 0)
    #number of detections that came before jday 60
    n1W <- sum(cyspdata$detect * as.numeric(cyspdata$jday < 60))
    #number of unique days with detections
    njd1 <- length(unique(cyspdata$jday[which(cyspdata$detect == 1)]))
    #number of unique days with non-detection
    njd0 <- length(unique(cyspdata$jday[which(cyspdata$detect == 0)]))
    
    
    if (n1 > 0)
    {
      #number of unique days of non-detections before first detection
      njd0i <- length(unique(cyspdata$jday[which(cyspdata$detect == 0 & cyspdata$jday < 
                                                 min(cyspdata$jday[which(cyspdata$detect == 1)]))]))
      #number of non-detections before first detection
      n0i <- length(which(cyspdata$detect == 0 & 
                            cyspdata$jday < min(cyspdata$jday[which(cyspdata$detect == 1)])))
    } else {
      njd0i <- 0
      n0i <- 0
    }
    
    halfmax_df$n1[counter] <- n1
    halfmax_df$n1W[counter] <- n1W
    halfmax_df$n0[counter] <- n0
    halfmax_df$n0i[counter] <- n0i
    halfmax_df$njd1[counter] <- njd1
    halfmax_df$njd0[counter] <- njd0
    halfmax_df$njd0i[counter] <- njd0i
    
    #defaults for rstanarm are 0.95 and 15
    DELTA <- 0.95
    TREE_DEPTH <- 15
    
    if (n1 > 29 & n1W < (n1/50) & n0 > 29 & njd0i > 29 & njd1 > 19)
    {
      fit2 <- rstanarm::stan_gamm4(detect ~ s(jday) + shr, 
                                 data = cyspdata,
                                 family = binomial(link = "logit"),
                                 algorithm = 'sampling',
                                 iter = ITER,
                                 chains = CHAINS,
                                 cores = CHAINS,
                                 adapt_delta = DELTA,
                                 control = list(max_treedepth = TREE_DEPTH))
      
      #calculate diagnostics
      num_diverge <- rstan::get_num_divergent(fit2$stanfit)
      num_tree <- rstan::get_num_max_treedepth(fit2$stanfit)
      num_BFMI <- length(rstan::get_low_bfmi_chains(fit2$stanfit))
      
      #rerun model if things didn't go well
      while (sum(c(num_diverge, num_BFMI)) > 0 & DELTA <= 0.98)
      {
        DELTA <- DELTA + 0.01
        TREE_DEPTH <- TREE_DEPTH + 1

        fit2 <- rstanarm::stan_gamm4(detect ~ s(jday) + shr,
                                   data = cyspdata,
                                   family = binomial(link = "logit"),
                                   algorithm = 'sampling',
                                   iter = ITER,
                                   chains = CHAINS,
                                   cores = CHAINS,
                                   adapt_delta = DELTA,
                                   control = list(max_treedepth = TREE_DEPTH))

        num_diverge <- rstan::get_num_divergent(fit2$stanfit)
        num_tree <- rstan::get_num_max_treedepth(fit2$stanfit)
        num_BFMI <- length(rstan::get_low_bfmi_chains(fit2$stanfit))
      }
      
      halfmax_df$num_diverge[counter] <- num_diverge
      halfmax_df$num_tree[counter] <- num_tree
      halfmax_df$num_BFMI[counter] <- num_BFMI
      halfmax_df$delta[counter] <- DELTA
      halfmax_df$tree_depth[counter] <- TREE_DEPTH
      
      #generate predict data
      predictDays <- range(cyspdata$jday)[1]:range(cyspdata$jday)[2]
      newdata <- data.frame(jday = predictDays, shr = 0)
      
      #predict response
      dfit <- rstanarm::posterior_linpred(fit2, newdata = newdata, transform = T)
      halfmax_fit <- rep(NA, ((ITER/2)*CHAINS))
      
      for (L in 1:((ITER/2)*CHAINS))
      {
        #L <- 1
        rowL <- as.vector(dfit[L,])
        halfmax_fit[L] <- predictDays[min(which(rowL > (max(rowL)/2)))]
      }
      
      ########################
      #PLOT MODEL FIT AND DATA

      #summary(fit2)
      mn_dfit <- apply(dfit, 2, mean)
      LCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.025))
      UCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.975))
      mn_hm <- mean(halfmax_fit)
      LCI_hm <- quantile(halfmax_fit, probs = 0.025)
      UCI_hm <- quantile(halfmax_fit, probs = 0.975)

      pdf(paste0(args, '_', years[j], '_', cells[k], '_arrival.pdf'))
      plot(predictDays, UCI_dfit, type = 'l', col = 'red', lty = 2, lwd = 2,
           ylim = c(0, max(UCI_dfit)),
           main = paste0(args, ' - ', years[j], ' - ', cells[k]),
           xlab = 'Julian Day', ylab = 'Detection Probability')
      lines(predictDays, LCI_dfit, col = 'red', lty = 2, lwd = 2)
      lines(predictDays, mn_dfit, lwd = 2)
      cyspdata$detect[which(cyspdata$detect == 1)] <- max(UCI_dfit)
      points(cyspdata$jday, cyspdata$detect, col = rgb(0,0,0,0.25))
      abline(v = mn_hm, col = rgb(0,0,1,0.5), lwd = 2)
      abline(v = LCI_hm, col = rgb(0,0,1,0.5), lwd = 2, lty = 2)
      abline(v = UCI_hm, col = rgb(0,0,1,0.5), lwd = 2, lty = 2)
      legend('topleft',
             legend = c('Cubic fit', 'CI fit', 'Half max', 'CI HM'),
             col = c('black', 'red', rgb(0,0,1,0.5), rgb(0,0,1,0.5)),
             lty = c(1,2,1,2), lwd = c(2,2,2,2), cex = 1.3)
      dev.off()
      ########################
      
      iter_ind <- grep('iter', colnames(halfmax_df))
      halfmax_df[counter,iter_ind] <- halfmax_fit
      halfmax_df$max_Rhat[counter] <- round(max(summary(fit2)[, "Rhat"]), 2)
      halfmax_df$min_neff[counter] <- min(summary(fit2)[, "n_eff"])
    }
    counter <- counter + 1
  } #k
} #j


#save to rds object
setwd(paste0(dir, '/Bird_Phenology/Data/Processed/halfmax_species_', RUN_DATE))
saveRDS(halfmax_df, file = paste0('halfmax_df_arrival_', args, '.rds'))



# runtime -----------------------------------------------------------------

time <- proc.time() - tt
rtime <- round(time[3]/60, 2)
paste0('Runtime (minutes): ', rtime)


print('I completed!')
