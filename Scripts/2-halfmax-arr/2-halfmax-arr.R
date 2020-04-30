######################
# 2 - halfmax model arrival
#
# Fit model (GAM logistic regression) to eBird data to get half-max parameter (bird arrival) for each species-cell-year
#
# y ~ bernouli(p)
# logit(p) = s(day)
#
#
# GAM notes
# ---------
# https://stats.stackexchange.com/questions/352995/does-there-exist-theory-behind-how-many-knots-one-should-use-in-a-stan-gamm4-m?noredirect=1&lq=1
# https://github.com/milkha/Splines_in_Stan/blob/master/splines_in_stan.pdf
# 
# Species name should be given as an arg to this script. The model will then be fit to that species only.
# Runtime: Up to 21 days on Xanadu (very long tail here, depends on data volume, etc.)
######################  


# top-level dir -----------------------------------------------------------

#dir <- '~/Google_Drive/R/'

#Xanadu
dir <- '/labs/Tingley/phenomismatch/'


# db query dir ------------------------------------------------------------

db_dir <- 'eBird_arrival_query_2020-02-25'

RUN_DATE <- '2020-02-26'


# model settings ----------------------------------------------------------

#number of iterations each model should be run
ITER <- 1500

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
library(rgeos)
library(rgdal)


# Set wd ------------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/'))


# Get args fed to script --------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)


# import processed data ---------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', db_dir))

#import data for species
spdata <- readRDS(paste0('ebird_arrival_query_', args, '.rds'))


# create grid -------------------------------------------------------------

#make hexgrid
hexgrid6 <- dggridR::dgconstruct(res = 6)

setwd(paste0(dir, 'Bird_Phenology/Data/hex_grid_crop/'))

#read in grid
hge <- rgdal::readOGR('hex_grid_crop.shp', verbose = FALSE)
hge_cells <- as.numeric(as.character(hge@data[,1]))


# filter cells by range  ---------------------------------------------------

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
sp_rng <- rgdal::readOGR(fname[1], verbose = FALSE)

#filter by breeding (2) and migration (4) range - need to convert spdf to sp
nrng <- sp_rng[which(sp_rng$SEASONAL == 2 | sp_rng$SEASONAL == 4),]

#filter by resident (1) and non-breeding (3) to exclude hex cells that contain 2/4 and 1/3
nrng_rm <- sp_rng[which(sp_rng$SEASONAL == 1 | sp_rng$SEASONAL == 3),]

#remove unneeded objects
rm(sp_rng)
rm(fname)


#if there is a legitimate range
if (NROW(nrng@data) > 0 & raster::extent(nrng)@xmax > -95)
{
  #good cells
  nrng_sp <- sp::SpatialPolygons(nrng@polygons)
  sp::proj4string(nrng_sp) <- sp::CRS(sp::proj4string(nrng))
  #find intersections with code from here: https://gis.stackexchange.com/questions/140504/extracting-intersection-areas-in-r
  poly_int <- rgeos::gIntersects(hge, nrng_sp, byid = TRUE)
  tpoly <- which(poly_int == TRUE, arr.ind = TRUE)[,2]
  br_mig_cells <- hge_cells[as.numeric(tpoly[!duplicated(tpoly)])]
  
  #bad cells
  if (length(nrng_rm) > 0)
  {
    nrng_rm_sp <- sp::SpatialPolygons(nrng_rm@polygons)
    sp::proj4string(nrng_rm_sp) <- sp::CRS(sp::proj4string(nrng_rm))
    poly_int_rm <- rgeos::gIntersects(hge, nrng_rm_sp, byid = TRUE)
    tpoly_rm <- which(poly_int_rm == TRUE, arr.ind = TRUE)[,2]
    res_ovr_cells <- hge_cells[as.numeric(tpoly_rm[!duplicated(tpoly_rm)])]
    
    #remove cells that appear in resident and overwinter range that also appear in breeding range
    cell_mrg <- c(br_mig_cells, res_ovr_cells)
    to_rm <- c(cell_mrg[duplicated(cell_mrg)], 812, 813, 841)
    
    rm(nrng_rm)
    rm(nrng_rm_sp)
    rm(res_ovr_cells)
    
  } else {
    cell_mrg <- br_mig_cells
    to_rm <- c(812, 813, 841)
  }
  
  #remove unneeded objects
  rm(hge)
  rm(nrng)
  rm(nrng_sp)
  
  c_rm <- which(br_mig_cells %in% to_rm)
  if (length(c_rm) > 0)
  {
    overlap_cells <- br_mig_cells[-c_rm]  
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
  
  cells <- cc_df$cell
  
  #retain rows that match selected cells
  spdata2 <- spdata[which(spdata$cell %in% cells),]
  
  #remove unneeded objects
  rm(cc_df)
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
                           mlmax = NA,
                           plmax = NA,
                           num_diverge = NA,
                           num_tree = NA,
                           num_BFMI = NA,
                           delta = NA,
                           tree_depth = NA,
                           t_iter = NA,
                           n1 = NA,
                           n1W = NA,
                           n0 = NA,
                           n0i = NA,
                           njd = NA,
                           njd1 = NA,
                           njd0 = NA,
                           njd0i = NA,
                           t_mat)

  #save to rds object
  setwd(paste0(dir, 'Bird_Phenology/Data/Processed/halfmax_arrival_', RUN_DATE))
  saveRDS(halfmax_df, file = paste0('halfmax_arrival_', args, '.rds'))
  
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
                         mlmax = NA,
                         plmax = NA,
                         num_diverge = NA,
                         num_tree = NA,
                         num_BFMI = NA,
                         delta = NA,
                         tree_depth = NA,
                         t_iter = NA,
                         n1 = NA,
                         n1W = NA,
                         n0 = NA,
                         n0i = NA,
                         njd = NA,
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
  yspdata <- dplyr::filter(spdata2, year == years[j])
  
  for (k in 1:ncell)
  {
    #k <- 16
    print(paste0('species: ', args, ', year: ', j, ', cell: ', k))
    
    cyspdata <- dplyr::filter(yspdata, cell == cells[k])
    
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
    #number of total unique days
    njd <- length(unique(cyspdata$jday))
    
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
    halfmax_df$njd[counter] <- njd
    halfmax_df$njd1[counter] <- njd1
    halfmax_df$njd0[counter] <- njd0
    halfmax_df$njd0i[counter] <- njd0i
    
    #center effort (to make sure supplying value of 0 in post-model prediction is meaningful)
    cyspdata$shr <- scale(cyspdata$duration_minutes, scale = FALSE)[,1]
    
    #defaults for rstanarm are 0.95 and 15
    DELTA <- 0.95
    TREE_DEPTH <- 15
    
    if (n1 > 19 & n1W < (n1/50) & n0 > 19 & njd0i > 19 & njd1 > 9)
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
      
      #rerun model if divergences occurred
      while (sum(c(num_diverge, num_BFMI)) > 0 & DELTA <= 0.98)
      {
        DELTA <- DELTA + 0.01

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
      
      max_Rhat <- round(max(summary(fit2)[, 'Rhat']), 2)
      min_neff <- min(summary(fit2)[, 'n_eff'])
      
      halfmax_df$num_diverge[counter] <- num_diverge
      halfmax_df$num_tree[counter] <- num_tree
      halfmax_df$num_BFMI[counter] <- num_BFMI
      halfmax_df$delta[counter] <- DELTA
      halfmax_df$tree_depth[counter] <- TREE_DEPTH
      halfmax_df$t_iter[counter] <- ITER
      halfmax_df$max_Rhat[counter] <- max_Rhat
      halfmax_df$min_neff[counter] <- min_neff
      
      #generate predict data
      predictDays <- range(cyspdata$jday)[1]:range(cyspdata$jday)[2]
      newdata <- data.frame(jday = predictDays, shr = 0)
      
      #predict response
      dfit <- rstanarm::posterior_linpred(fit2, newdata = newdata, transform = T)
      halfmax_fit <- rep(NA, ((ITER/2)*CHAINS))
      tlmax <- rep(NA, ((ITER/2)*CHAINS))
      #day at which probability of occurence is half local maximum value
      for (L in 1:((ITER/2)*CHAINS))
      {
        rowL <- as.vector(dfit[L,])
        #first detection
        fd <- min(cyspdata$jday[which(cyspdata$detect == 1)])
        #local maximum(s)
        #from: stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
        lmax_idx <- which(diff(sign(diff(rowL))) == -2) + 1
        lmax <- predictDays[lmax_idx]
        #first local max to come after first detection
        flm <- which(lmax > fd)
        if (length(flm) > 0)
        {
          #first local max to come after first detection
          lmax2_idx <- lmax_idx[min(flm)]
          lmax2 <- lmax[min(flm)]
          tlmax[L] <- TRUE
        } else {
          #no local max
          lmax2_idx <- which.max(rowL)
          lmax2 <- predictDays[which.max(rowL)]
          tlmax[L] <- FALSE
        }
        #position of min value before max - typically, where 0 is
        lmin_idx <- which.min(rowL[1:lmax2_idx])
        #value at local max - value at min (typically 0)
        dmm <- rowL[lmax2_idx] - rowL[lmin_idx]
        #all positions less than or equal to half diff between max and min + value min
        tlm <- which(rowL <= ((dmm/2) + rowL[lmin_idx]))
        #which of these come before max and after min
        vgm <- which(tlm < lmax2_idx & tlm > lmin_idx)
        #insert halfmax (first day for situations where max is a jday = 1)
        if (length(vgm) > 0)
        {
          halfmax_fit[L] <- predictDays[max(vgm)]
        } else {
          halfmax_fit[L] <- predictDays[1]
        }
      }
      
      #number of iterations that had local max
      halfmax_df$plmax[counter] <- round(sum(tlmax)/((ITER/2)*CHAINS), 3)
      
      #model fit
      mn_dfit <- apply(dfit, 2, mean)
      LCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.025))
      UCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.975))
      
      #check whether local max exists for mean model fit
      mlmax <- sum(diff(sign(diff(mn_dfit))) == -2)
      if (mlmax > 0)
      {
        halfmax_df$mlmax[counter] <- TRUE
      } else {
        halfmax_df$mlmax[counter] <- FALSE
      }
      
      #estimated halfmax
      mn_hm <- mean(halfmax_fit)
      LCI_hm <- quantile(halfmax_fit, probs = 0.025)
      UCI_hm <- quantile(halfmax_fit, probs = 0.975)
      
      #fill df with halfmax iter
      cndf <- colnames(halfmax_df)
      iter_ind <- grep('iter', cndf)
      to.rm <- which(cndf == 't_iter')
      #remove t_iter col
      halfmax_df[counter, iter_ind[-which(iter_ind == to.rm)]] <- halfmax_fit
      
      ########################
      #PLOT MODEL FIT AND DATA
      
      pdf(paste0(args, '_', years[j], '_', cells[k], '_arrival.pdf'))
      plot(predictDays, UCI_dfit, type = 'l', col = 'red', lty = 2, lwd = 2,
           ylim = c(0, max(UCI_dfit)),
           main = paste0(args, ' - ', years[j], ' - ', cells[k]),
           xlab = 'Julian Day', ylab = 'Probability of occurrence')
      lines(predictDays, LCI_dfit, col = 'red', lty = 2, lwd = 2)
      lines(predictDays, mn_dfit, lwd = 2)
      dd <- cyspdata$detect
      dd[which(dd == 1)] <- max(UCI_dfit)
      points(cyspdata$jday, dd, col = rgb(0,0,0,0.25))
      abline(v = mn_hm, col = rgb(0,0,1,0.5), lwd = 2)
      abline(v = LCI_hm, col = rgb(0,0,1,0.5), lwd = 2, lty = 2)
      abline(v = UCI_hm, col = rgb(0,0,1,0.5), lwd = 2, lty = 2)
      legend('topleft',
             legend = c('Model fit', 'CI fit', 'Half max', 'CI HM'),
             col = c('black', 'red', rgb(0,0,1,0.5), rgb(0,0,1,0.5)),
             lty = c(1,2,1,2), lwd = c(2,2,2,2), cex = 1.3)
      dev.off()
      
      #PLOT CODE FOR PRESENTATION
      # mn_dfit <- apply(dfit, 2, mean)
      # LCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.025))
      # UCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.975))
      # mn_hm <- mean(halfmax_fit)
      # LCI_hm <- quantile(halfmax_fit, probs = 0.025)
      # UCI_hm <- quantile(halfmax_fit, probs = 0.975)
      # 
      # pdf(paste0(args, '_', years[j], '_', cells[k], '_arrival.pdf'))
      # plot(predictDays, UCI_dfit, type = 'l', col = 'white', lty = 2, lwd = 5,
      #      ylim = c(0, max(UCI_dfit)),
      #      #main = paste0(args, ' - ', years[j], ' - ', cells[k]),
      #      #xlab = 'Julian Day', ylab = 'Probability of occurrence',
      #      tck = 0, ann = FALSE, xaxt = 'n', yaxt = 'n')
      # # lines(predictDays, LCI_dfit, col = 'red', lty = 2, lwd = 5)
      # # lines(predictDays, mn_dfit, lwd = 5)
      # cyspdata$detect[which(cyspdata$detect == 1)] <- max(UCI_dfit)
      # points(cyspdata$jday, cyspdata$detect, col = rgb(0,0,0,0.25), cex = 2)
      # # abline(v = mn_hm, col = rgb(0,0,1,0.5), lwd = 5)
      # # abline(v = LCI_hm, col = rgb(0,0,1,0.5), lwd = 5, lty = 2)
      # # abline(v = UCI_hm, col = rgb(0,0,1,0.5), lwd = 5, lty = 2)
      # # legend('topleft',
      # #        legend = c('Model fit', 'CI fit', 'Half max', 'CI HM'),
      # #        col = c('black', 'red', rgb(0,0,1,0.5), rgb(0,0,1,0.5)),
      # #        lty = c(1,2,1,2), lwd = c(2,2,2,2), cex = 1.3)
      # dev.off()
      
      # #alternative visualization
      # pdf(paste0(args, '_', years[j], '_', cells[k], '_arrival_realizations.pdf'))
      # plot(NA, xlim = c(range(cyspdata$jday)[1]:range(cyspdata$jday)[2]), ylim = c(0, quantile(dfit, 0.999)),
      #      xlab = 'Julian Day', ylab = 'Probabiity of occurrence')
      # for (L in 1:((ITER/2)*CHAINS))
      # {
      #   lines(range(cyspdata$jday)[1]:range(cyspdata$jday)[2], as.vector(dfit[L,]), type = 'l', col = rgb(0,0,0,0.025))
      # }
      # for (L in 1:((ITER/2)*CHAINS))
      # {
      #   abline(v = halfmax_fit[L], col = rgb(1,0,0,0.025))
      # }
      # dev.off()
      ########################
    }
    counter <- counter + 1
  } #k
} #j


#order by year then cell
OUT <- halfmax_df[order(halfmax_df[,'year'], halfmax_df[,'cell']),]

#save to rds object
setwd(paste0(dir, '/Bird_Phenology/Data/Processed/halfmax_arrival_', RUN_DATE))
saveRDS(OUT, file = paste0('halfmax_arrival_', args, '.rds'))


# runtime -----------------------------------------------------------------

time <- proc.time() - tt
rtime <- round(time[3]/60, 2)
paste0('Runtime (minutes): ', rtime)


print('I completed!')
