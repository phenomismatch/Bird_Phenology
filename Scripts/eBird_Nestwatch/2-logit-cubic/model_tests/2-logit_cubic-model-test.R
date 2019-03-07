#########################
#Determine complexity of polynomial models
#
#calculate WAIC and SSE
#produce diagnostic plots with PPC on binned data
#########################

#Species, years, cells investigated (inital)
#Agelaius phoeniceus - 2014 - 536
#Agelaius phoeniceus - 2014 - 538
#Agelaius phoeniceus - 2014 - 564
#Agelaius phoeniceus - 2014 - 566


# top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/UCHC/LABS/Tingley/phenomismatch/'


# db query dir ------------------------------------------------------------

db_dir <- 'eBird_query_2019-02-02'
RUN_DATE <- '2019-02-02'



# model settings ----------------------------------------------------------

#number of iterations each model should be run
ITER <- 750
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
library(loo)

# Set wd ------------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/'))



# Get args fed to script --------------------------------------------------

#args <- commandArgs(trailingOnly = TRUE)
#args <- 'Empidonax_virescens'
#args <- 'Mimus_polyglottos'
#args <- 'Vireo_olivaceus'
args <- 'Agelaius_phoeniceus'


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
sp_rng2 <- raster::crop(sp_rng, extent(-100, -50, 26, 90))

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

#cells <- c(536, 538, 564, 566)
ncell <- length(cells)

years <- min(spdata2$year):max(spdata2$year)
#years <- 2014
nyr <- length(years)


formulas <- list(f2 = 'detect ~ sjday + sjday2 + shr',
                 f3 = 'detect ~ sjday + sjday2 + sjday3 + shr',
                 f4 = 'detect ~ sjday + sjday2 + sjday3 + sjday4 + shr',
                 f5 = 'detect ~ sjday + sjday2 + sjday3 + sjday4 + sjday5 + shr',
                 f6 = 'detect ~ sjday + sjday2 + sjday3 + sjday4 + sjday5 + sjday6 + shr')

#number of models
nm <- length(formulas)

# fit logit cubic ---------------------------------------------------------

# t_mat <- matrix(data = NA, nrow = ncell*nyr, ncol = ((ITER/2)*CHAINS))
# colnames(t_mat) <- paste0('iter_', 1:((ITER/2)*CHAINS))
# halfmax_df <- data.frame(species = args, 
#                          year = rep(years, each = ncell), 
#                          cell = rep(cells, nyr), 

#t_mat <- matrix(data = NA, nrow = ncell * nyr * nm, ncol = ((ITER/2)*CHAINS))
#colnames(t_mat) <- paste0('iter_', 1:((ITER/2)*CHAINS))

halfmax_df <- data.frame(species = rep(args, times = (ncell * nyr * nm)), 
                         year = NA, 
                         cell = NA,
                         poly = NA,
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
                         hf_mean = NA,
                         hf_sd = NA,
                         waic = NA,
                         waic_se = NA,
                         ppc_b1 = NA,
                         ppc_b2 = NA,
                         ppc_b3 = NA,
                         ppc_b4 = NA,
                         ppc_b5 = NA,
                         ppc_b6 = NA,
                         ppc_b7 = NA,
                         ppc_b8 = NA,
                         ppc_b9 = NA,
                         ppc_b10 = NA,
                         SSE_b1 = NA,
                         SSE_b2 = NA,
                         SSE_b3 = NA,
                         SSE_b4 = NA,
                         SSE_b5 = NA,
                         SSE_b6 = NA,
                         SSE_b7 = NA,
                         SSE_b8 = NA,
                         SSE_b9 = NA,
                         SSE_b10 = NA)


# #add other polynomials
# spdata2$ssjday <- scale(spdata2$sjday)
# spdata2$ssjday2 <- scale(spdata2$sjday^2)
# spdata2$ssjday3 <- scale(spdata2$sjday^3)
# spdata2$ssjday4 <- scale(spdata2$sjday^4)
# spdata2$ssjday5 <- scale(spdata2$sjday^5)
# spdata2$ssjday6 <- scale(spdata2$sjday^6)
# spdata2$ssjday7 <- scale(spdata2$sjday^7)

#add other polynomials
spdata2$sjday3 <- spdata2$sjday^3
spdata2$sjday4 <- spdata2$sjday^4
spdata2$sjday5 <- spdata2$sjday^5
spdata2$sjday6 <- spdata2$sjday^6
spdata2$sjday7 <- spdata2$sjday^7

# formulas <- list(f3 = 'detect ~ ssjday + ssjday2 + ssjday3 + shr',
#                  f4 = 'detect ~ ssjday + ssjday2 + ssjday3 + ssjday4 + shr',
#                  f5 = 'detect ~ ssjday + ssjday2 + ssjday3 + ssjday4 + ssjday5 + shr',
#                  f6 = 'detect ~ ssjday + ssjday2 + ssjday3 + ssjday4 + ssjday5 + ssjday6 + shr',
#                  f7 = 'detect ~ ssjday + ssjday2 + ssjday3 + ssjday4 + ssjday5 + ssjday6 + ssjday7 + shr')


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
    n1W <- sum(cyspdata$detect * as.numeric(cyspdata$day < 60))
    #number of unique days with detections
    njd1 <- length(unique(cyspdata$day[which(cyspdata$detect == 1)]))
    #number of unique days with non-detection
    njd0 <- length(unique(cyspdata$day[which(cyspdata$detect == 0)]))
    
    
    if (n1 > 0)
    {
      #number of unique days of non-detections before first detection
      njd0i <- length(unique(cyspdata$day[which(cyspdata$detect == 0 & cyspdata$day < 
                                                  min(cyspdata$day[which(cyspdata$detect == 1)]))]))
      #number of non-detections before first detection
      n0i <- length(which(cyspdata$detect == 0 & 
                            cyspdata$day < min(cyspdata$day[which(cyspdata$detect == 1)])))
    } else {
      njd0i <- 0
      n0i <- 0
    }
    
    
    #defaults for rstanarm are 0.95 and 15
    DELTA <- 0.95
    TREE_DEPTH <- 15
    
    if (n1 > 29 & n1W < (n1/50) & n0 > 29 & njd0i > 29 & njd1 > 19)
    {
      #generate predict data
      predictDays <- range(cyspdata$sjday)[1]:range(cyspdata$sjday)[2]
      predictDays2 <- predictDays^2
      predictDays3 <- predictDays^3
      predictDays4 <- predictDays^4
      predictDays5 <- predictDays^5
      predictDays6 <- predictDays^6
      predictDays7 <- predictDays^7
      
      newdata <- data.frame(sjday = predictDays, sjday2 = predictDays2,
                            sjday3 = predictDays3, sjday4 = predictDays4,
                            sjday5 = predictDays5, sjday6 = predictDays6,
                            sjday7 = predictDays7, shr = 0)
      
      # #generate predict data (scaled)
      # predictDays <- scale(range(cyspdata$sjday)[1]:range(cyspdata$sjday)[2])
      # predictDays2 <- scale(predictDays^2)
      # predictDays3 <- scale(predictDays^3)
      # predictDays4 <- scale(predictDays^4)
      # predictDays5 <- scale(predictDays^5)
      # predictDays6 <- scale(predictDays^6)
      # predictDays7 <- scale(predictDays^7)
      # 
      # newdata <- data.frame(ssjday = predictDays, ssjday2 = predictDays2,
      #                       ssjday3 = predictDays3, ssjday4 = predictDays4,
      #                       ssjday5 = predictDays5, ssjday6 = predictDays6,
      #                       ssjday7 = predictDays7, shr = 0)
      
      for (m in 1:length(formulas))
      {
        #m <- 1
        
        halfmax_df$year[counter] <- years[j]
        halfmax_df$cell[counter] <- cells[k]
        halfmax_df$poly[counter] <- m + 1
        
        halfmax_df$n1[counter] <- n1
        halfmax_df$n1W[counter] <- n1W
        halfmax_df$n0[counter] <- n0
        halfmax_df$n0i[counter] <- n0i
        halfmax_df$njd1[counter] <- njd1
        halfmax_df$njd0[counter] <- njd0
        halfmax_df$njd0i[counter] <- njd0i
        
        tfit <- rstanarm::stan_glm(formulas[[m]], 
                                   data = cyspdata,
                                   family = binomial(link = "logit"),
                                   algorithm = 'sampling',
                                   iter = ITER,
                                   chains = CHAINS,
                                   cores = CHAINS,
                                   adapt_delta = DELTA,
                                   control = list(max_treedepth = TREE_DEPTH))
      
        name <- paste0(args, '_', years[j], '_', cells[k], '-poly_', m + 1)
        assign(name, tfit)
        
        #calculate diagnostics
        halfmax_df$num_diverge[counter] <- rstan::get_num_divergent(tfit$stanfit)
        halfmax_df$num_tree[counter] <- rstan::get_num_max_treedepth(tfit$stanfit)
        halfmax_df$num_BFMI[counter] <- length(rstan::get_low_bfmi_chains(tfit$stanfit))
        halfmax_df$delta[counter] <- DELTA
        halfmax_df$tree_depth[counter] <- TREE_DEPTH
        halfmax_df$max_Rhat[counter] <- round(max(summary(tfit)[, "Rhat"]), 2)
        halfmax_df$min_neff[counter] <- min(summary(tfit)[, "n_eff"])
        
        
        #predict response
        dfit <- rstanarm::posterior_linpred(tfit, newdata = newdata, transform = T)
        halfmax_fit <- rep(NA, ((ITER/2)*CHAINS))
      
        for (L in 1:((ITER/2)*CHAINS))
        {
          #L <- 1
          rowL <- as.vector(dfit[L,])
          halfmax_fit[L] <- predictDays[min(which(rowL > (max(rowL)/2)))]
        }
      
        #draw from posterior predictive distribution
        pp <- rstanarm::posterior_predict(tfit)
        
        #assign bin numbers to each point
        bin_number <- 1 + cyspdata$day %/% 20

        ########################
        #PLOT MODEL FIT AND DATA
        
        setwd('~/Desktop/Agelaius_phoeniceus_tests')
        cyspdata2 <- cyspdata
        
        #summary(fit2)
        mn_dfit <- apply(dfit, 2, mean)
        LCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.025))
        UCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.975))
        mn_hm <- mean(halfmax_fit)
        LCI_hm <- quantile(halfmax_fit, probs = 0.025)
        UCI_hm <- quantile(halfmax_fit, probs = 0.975)
        
        pdf(paste0(args, '_', years[j], '_', cells[k], '-poly_', m+1, '_arrival.pdf'))
        plot(predictDays, UCI_dfit, type = 'l', col = 'red', lty = 2, lwd = 3,
             ylim = c(-(max(UCI_dfit)/5), max(UCI_dfit)),
             main = paste0(args, ' - ', years[j], ' - ', cells[k], '- poly ', m + 1),
             xlab = 'Julian Day', ylab = 'Detection Probability')
        lines(predictDays, LCI_dfit, col = 'red', lty = 2, lwd = 3)
        lines(predictDays, mn_dfit, lwd = 3)
        cyspdata2$detect[which(cyspdata2$detect == 1)] <- max(UCI_dfit)
        points(cyspdata2$day, cyspdata2$detect, col = rgb(0,0,0,0.25))
        
        ds <- sort(unique(cyspdata$day))
        prop_d <- rep(NA, length(ds))
        for (p in 1:length(ds))
        {
          #p <- 2
          td <- filter(cyspdata, day == ds[p])
          prop_d[p] <- mean(td$detect)
        }
        
        ma <- function(x, n = 5)
        {
          stats::filter(x, rep(1/n, n), sides = 2)
        }
        ma_prop_d <- ma(prop_d) * max(UCI_dfit)
        
        lines(ds, ma_prop_d, col = rgb(0.5,0.5,0.5,0.9), lwd = 2)
        
        segments(x0 = mn_hm, x1 = mn_hm, y0 = 0, y1 = 1, col = rgb(0,0,1,0.5), lwd = 3)
        segments(x0 = LCI_hm, x1 = LCI_hm, y0 = 0, y1 = 1, col = rgb(0,0,1,0.5), lwd = 3, lty = 2)
        segments(x0 = UCI_hm, x1 = UCI_hm, y0 = 0, y1 = 1, col = rgb(0,0,1,0.5), lwd = 3, lty = 2)

        rect(xleft = 0, ybottom = 0, xright = 20, ytop = max(UCI_dfit),
             col = rgb(0,1,1,0.2), border = FALSE)
        rect(xleft = 40, ybottom = 0, xright = 60, ytop = max(UCI_dfit),
             col = rgb(0,1,1,0.2), border = FALSE)
        rect(xleft = 80, ybottom = 0, xright = 100, ytop = max(UCI_dfit),
             col = rgb(0,1,1,0.2), border = FALSE)
        rect(xleft = 120, ybottom = 0, xright = 140, ytop = max(UCI_dfit),
             col = rgb(0,1,1,0.2), border = FALSE)
        rect(xleft = 160, ybottom = 0, xright = 180, ytop = max(UCI_dfit),
             col = rgb(0,1,1,0.2), border = FALSE)
        legend('topleft',
               legend = c('Cubic fit', 'CI fit', 'Half max', 'CI HM', '5d MA'),
               col = c('black', 'red', rgb(0,0,1,0.5), rgb(0,0,1,0.5), rgb(0.5,0.5,0.5,0.7)),
               lty = c(1,2,1,2), lwd = c(2,2,2,2), cex = 1.3)
        for (bin in 1:10) 
        {
          #bin <- 3
          #for each draw, how many detections are in a particular bin
          if (sum(bin_number == bin) > 1)
          {
            binpredict <- rowSums(pp[, bin_number == bin])
            truenumber <- sum(cyspdata$detect[bin_number == bin])
          
            #what proportion of time are more detections predicted
            halfmax_df[counter, paste0('ppc_b', bin)] <- sum(binpredict > truenumber) / ((ITER/2)*CHAINS)
          
            #sum of squared errors for that bin
            halfmax_df[counter, paste0('SSE_b', bin)] <- sum((truenumber - binpredict)^2)
          
            # sum(binpredict == truenumber) / 1000
            # sum(binpredict > truenumber) / 1000
            # sum(binpredict < truenumber) / 1000
          
            h_plt <- function()
            {
              hist(binpredict, xaxt = 'n', yaxt = 'n', main = NULL,
                   xlab = NULL, ylab = NULL)
              abline(v = truenumber, col = 'red', lwd = 3)
            }
            Hmisc::subplot(h_plt, x = ((bin * 20) - 8), y = -(max(UCI_dfit) / 8), size = c(0.5, 0.5))
          }
        }
        
        dev.off()
        ########################
        
        #halfax mean and sd
        halfmax_df[counter, 'hf_mean'] <- mean(halfmax_fit)
        halfmax_df[counter, 'hf_sd'] <- sd(halfmax_fit)
        
        #waic and se
        halfmax_df[counter, 'waic'] <- waic(tfit)$estimates[3, 1]
        halfmax_df[counter, 'waic_se'] <- waic(tfit)$estimates[3, 2]
        
        counter <- counter + 1
      }
    }
  } #k
} #j


saveRDS(halfmax_df, 'Agelaius_phoeniceus_poly_tests.rds')
  
