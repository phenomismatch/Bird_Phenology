######################
# 7 - GAM model juv MAPS
#
# run for all cells (in and out of study area)
######################


# runtime -----------------------------------------------------------------

tt <- proc.time()


# Top-level dir -----------------------------------------------------------

#desktop/laptop
#dir <- '~/Google_Drive/R/'

#Xanadu
dir <- '/labs/Tingley/phenomismatch/'


# other dir ------------------------------------------------------------

#run date
RUN_DATE <- '2020-06-04'

#MAPS query date
MAPS_date <- '2020-03-03'


# load packages -----------------------------------------------------------

library(rstanarm)
library(rstan)
library(dplyr)
library(dggridR)
library(sp)
library(raster)
library(rgeos)
library(rgdal)


# Get args fed to script --------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
#args <- 'Vireo_olivaceus'

# model settings ----------------------------------------------------------

ITER <- 3000
CHAINS <- 4


# read in data ------------------------------------------------------------

#read in - RDS create with 1-query-db.R in wing_chord_changes project
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/MAPS/'))

data_p <- readRDS(paste0('MAPS-age-filled-', MAPS_date, '.rds'))
colnames(data_p)[grep('sci_name', colnames(data_p))] <- 'species'
#add underscore to species names
data_p$species <- gsub(' ', '_', data_p$species)

data <- dplyr::filter(data_p, species == args)


# create grid -------------------------------------------------------------

#make hexgrid
hexgrid6 <- dggridR::dgconstruct(res = 6)
data$cell <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                      in_lon_deg = data$lng, 
                                      in_lat_deg = data$lat)[[1]]

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

#filter by breeding (2) - need to convert spdf to sp
nrng <- sp_rng[which(sp_rng$SEASONAL == 2),]

#do not exclude resident, non-breeding, and migration cells as breeding at station is a requirement for modeling

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
  poly_int <- rgeos::gIntersects(hge, nrng_sp, byid=TRUE)
  tpoly <- which(poly_int == TRUE, arr.ind = TRUE)[,2]
  br_cells <- hge_cells[as.numeric(tpoly[!duplicated(tpoly)])]
  
  to_rm <- c(812, 813, 841)
  
  #remove unneeded objects
  rm(hge)
  rm(nrng)
  rm(nrng_sp)
  
  c_rm <- which(br_cells %in% to_rm)
  if (length(c_rm) > 0)
  {
    overlap_cells <- br_cells[-c_rm]  
  } else {
    overlap_cells <- br_cells
  }
  
  #remove unneeded objects
  rm(br_cells)
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
  cells <- cc_df$cell
  
  #retain rows that match selected cells
  m_mf <- data[which(data$cell %in% cells),]
  
  #remove unneeded objects
  rm(cc_df)
  rm(data)
  
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
                           # n1W = NA,
                           n0 = NA,
                           n0i = NA,
                           njd = NA,
                           njd1 = NA,
                           njd0 = NA,
                           njd0i = NA,
                           t_mat)
  
  #save to rds object
  setwd(paste0(dir, 'Bird_Phenology/Data/Processed/halfmax_juv_', RUN_DATE))
  saveRDS(halfmax_df, file = paste0('halfmax_juv_', args, '.rds'))
  stop('Range not suitable for modeling!')
}

if (NROW(m_mf) == 0)
{
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
                           # n1W = NA,
                           n0 = NA,
                           n0i = NA,
                           njd = NA,
                           njd1 = NA,
                           njd0 = NA,
                           njd0i = NA,
                           t_mat)
  
  #save to rds object
  setwd(paste0(dir, 'Bird_Phenology/Data/Processed/halfmax_juv_', RUN_DATE))
  saveRDS(halfmax_df, file = paste0('halfmax_juv_', args, '.rds'))
  stop('No data in breeding range!')
}


# Add juv info ------------------------------------------------------------

#add col for juveniles (1 = juv, 0 = adult)
m_mf$juv <- NA
m_mf$juv[which(m_mf$age == 2)] <- 1 #hatching year bird
m_mf$juv[-which(m_mf$age == 2)] <- 0

#exclude young bird incapable of flight or missing band_id
to.na <- which(m_mf$age == 4 | is.na(m_mf$band_id))
if (length(to.na) > 0)
{
  m_mf2 <- m_mf[-to.na,]  
} else {
  m_mf2 <- m_mf
}

#all band_ids numeric
m_mf2$band_id <- as.numeric(m_mf2$band_id)

years <- sort(unique(m_mf2$year))
cells <- sort(unique(m_mf2$cell))
ncell <- length(cells)
nyr <- length(years)


# setup data object -------------------------------------------------------

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
                         # n1W = NA,
                         n0 = NA,
                         n0i = NA,
                         njd = NA,
                         njd1 = NA,
                         njd0 = NA,
                         njd0i = NA,
                         t_mat)


#create dir for figs if doesn't exist
ifelse(!dir.exists(paste0(dir, 'Bird_Phenology/Figures/halfmax/juv_', RUN_DATE)),
       dir.create(paste0(dir, 'Bird_Phenology/Figures/halfmax/juv_', RUN_DATE)),
       FALSE)

setwd(paste0(dir, 'Bird_Phenology/Figures/halfmax/juv_', RUN_DATE))


#loop through each species, year, cell and extract half-max parameter
counter <- 1
for (j in 1:nyr)
{
  #j <- 3
  yspdata <- dplyr::filter(m_mf2, year == years[j])
  
  for (k in 1:ncell)
  {
    #k <- 24
    print(paste0('species: ', args, ', year: ', j, ', cell: ', k))
    
    #filter by cell and breeding status for capture station
    cyspdata <- dplyr::filter(yspdata, cell == cells[k], yr_br_status == 'B')
    
    #adults
    cyspdata_0 <- dplyr::filter(cyspdata, juv == 0)
    #juv
    cyspdata_1 <- dplyr::filter(cyspdata, juv == 1)
    
    #only first capture for each juv
    ujuv <- unique(as.numeric(cyspdata_1$band_id))
    jidx <- rep(NA, length(ujuv))
    tna <- c()
    if (length(ujuv) > 0)
    {
      for (b in 1:length(ujuv))
      {
        #b <- 43
        tidx <- which(cyspdata_1$band_id == ujuv[b])
        if (length(tidx) > 1)
        {
          #first day obs idx
          jidx[b] <- tidx[which.min(cyspdata_1$day[tidx])]
        } else {
          jidx[b] <- tidx
        }
      }
    }
    
    #bind adult and juv data together
    cyspdata_j <- rbind(cyspdata_0, cyspdata_1[jidx,])
    
    #number of surveys where juv were captured
    n1 <- sum(cyspdata_j$juv)
    #number of surveys where juv were not captured (but adult was)
    n0 <- sum(cyspdata_j$juv == 0)
    # #number of detections that came before jday 60
    # n1W <- sum(cyspdata_j$juv * as.numeric(cyspdata_j$day < 60))
    #number of unique days with juv captures
    njd1 <- length(unique(cyspdata_j$day[which(cyspdata_j$juv == 1)]))
    #number of unique days with no juv captures (but adults captures)
    njd0 <- length(unique(cyspdata_j$day[which(cyspdata_j$juv == 0)]))
    #number of total unique days
    njd <- length(unique(cyspdata_j$day))
    
    if (n1 > 0)
    {
      #number of unique days of adult captures before first juv capture
      njd0i <- length(unique(cyspdata_j$day[which(cyspdata_j$juv == 0 & cyspdata_j$day < 
                                                    min(cyspdata_j$day[which(cyspdata_j$juv == 1)]))]))
      #number of adult captures before first juv capture
      n0i <- length(which(cyspdata_j$juv == 0 & 
                            cyspdata_j$day < min(cyspdata_j$day[which(cyspdata_j$juv == 1)])))
    } else {
      njd0i <- 0
      n0i <- 0
    }
    
    halfmax_df$n1[counter] <- n1
    # halfmax_df$n1W[counter] <- n1W
    halfmax_df$n0[counter] <- n0
    halfmax_df$n0i[counter] <- n0i
    halfmax_df$njd[counter] <- njd
    halfmax_df$njd1[counter] <- njd1
    halfmax_df$njd0[counter] <- njd0
    halfmax_df$njd0i[counter] <- njd0i
    
    #defaults for rstanarm are 0.95 and 15
    DELTA <- 0.95
    TREE_DEPTH <- 15
    
    if (n1 > 4 & n0 > 4 & njd0i > 2 & njd1 > 2 & njd > 19)
    {
      fit2 <- rstanarm::stan_gamm4(juv ~ s(day, k = 20),
                                   data = cyspdata_j,
                                   family = binomial(link = "logit"),
                                   algorithm = 'sampling',
                                   iter = ITER,
                                   chains = CHAINS,
                                   cores = CHAINS,
                                   adapt_delta = DELTA,
                                   control = list(max_treedepth = TREE_DEPTH))
      
      # fit2 <- rstanarm::stan_glm(juv ~ day,
      #                              data = cyspdata_j,
      #                              family = binomial(link = "logit"),
      #                              algorithm = 'sampling',
      #                              iter = ITER,
      #                              chains = CHAINS,
      #                              cores = CHAINS,
      #                              adapt_delta = DELTA,
      #                              control = list(max_treedepth = TREE_DEPTH))
      
      #calculate diagnostics
      num_diverge <- rstan::get_num_divergent(fit2$stanfit)
      num_tree <- rstan::get_num_max_treedepth(fit2$stanfit)
      num_BFMI <- length(rstan::get_low_bfmi_chains(fit2$stanfit))
      
      #rerun model if divergences occurred
      while (sum(c(num_diverge, num_BFMI)) > 0 & DELTA <= 0.98)
      {
        DELTA <- DELTA + 0.01
        
        fit2 <- rstanarm::stan_gamm4(juv ~ s(day, k = 20),
                                     data = cyspdata_j,
                                     family = binomial(link = "logit"),
                                     algorithm = 'sampling',
                                     iter = ITER,
                                     chains = CHAINS,
                                     cores = CHAINS,
                                     adapt_delta = DELTA,
                                     control = list(max_treedepth = TREE_DEPTH))
        
        # fit2 <- rstanarm::stan_glm(juv ~ day,
        #                            data = cyspdata_j,
        #                            family = binomial(link = "logit"),
        #                            algorithm = 'sampling',
        #                            iter = ITER,
        #                            chains = CHAINS,
        #                            cores = CHAINS,
        #                            adapt_delta = DELTA,
        #                            control = list(max_treedepth = TREE_DEPTH))
        
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
      predictDays <- range(cyspdata_j$day)[1]:range(cyspdata_j$day)[2]
      newdata <- data.frame(day = predictDays)
      
      #predict response
      dfit <- rstanarm::posterior_linpred(fit2, newdata = newdata, transform = T)
      halfmax_fit <- rep(NA, ((ITER/2)*CHAINS))
      tlmax <- rep(NA, ((ITER/2)*CHAINS))
      #day at which probability of occurence is half local maximum value
      for (L in 1:((ITER/2)*CHAINS))
      {
        rowL <- as.vector(dfit[L,])
        #first detection
        fd <- min(cyspdata_j$day[which(cyspdata_j$juv == 1)])
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
        #local mins before max (global and local mins)
        lmin_idx <- c(which.min(rowL[1:lmax2_idx]), 
                      which(diff(sign(diff(rowL[1:lmax2_idx]))) == 2) + 1)
        lmin <- predictDays[lmin_idx]
        #local min nearest to local max
        lmin2_idx <- lmin_idx[which.min(lmax2 - lmin)]
        lmin2 <- predictDays[lmin2_idx]
        
        #value at local max - value at min (typically 0)
        dmm <- rowL[lmax2_idx] - rowL[lmin2_idx]
        #all positions less than or equal to half diff between max and min + value min
        tlm <- which(rowL <= ((dmm/2) + rowL[lmin2_idx]))
        #which of these come before max and after or at min
        
        vgm <- tlm[which(tlm < lmax2_idx & tlm >= lmin2_idx)]
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
      
      pdf(paste0(args, '_', years[j], '_', cells[k], '_juv.pdf'))
      plot(predictDays, UCI_dfit, type = 'l', col = 'red', lty = 2, lwd = 2,
           ylim = c(0, max(UCI_dfit)),
           main = paste0(args, ' - ', years[j], ' - ', cells[k]),
           xlab = 'Julian Day', ylab = 'Probability of capture')
      lines(predictDays, LCI_dfit, col = 'red', lty = 2, lwd = 2)
      lines(predictDays, mn_dfit, lwd = 2)
      dd <- cyspdata_j$juv
      dd[which(dd == 1)] <- max(UCI_dfit)
      points(cyspdata_j$day, dd, col = rgb(0,0,0,0.25))
      abline(v = mn_hm, col = rgb(0,0,1,0.5), lwd = 2)
      abline(v = LCI_hm, col = rgb(0,0,1,0.5), lwd = 2, lty = 2)
      abline(v = UCI_hm, col = rgb(0,0,1,0.5), lwd = 2, lty = 2)
      legend('topleft',
             legend = c('Model fit', 'CI fit', 'Half max', 'CI HM'),
             col = c('black', 'red', rgb(0,0,1,0.5), rgb(0,0,1,0.5)),
             lty = c(1,2,1,2), lwd = c(2,2,2,2), cex = 1.3)
      dev.off()
      
      # #alternative visualization
      pdf(paste0(args, '_', years[j], '_', cells[k], '_juv_realizations.pdf'))
      plot(NA, xlim = c(range(cyspdata_j$day)[1], range(cyspdata_j$day)[2]), 
           ylim = c(0, quantile(dfit, 0.999)),
           xlab = 'Julian Day', ylab = 'Probability of capture')
      for (L in 1:((ITER/2)*CHAINS))
      {
        lines(range(cyspdata_j$day)[1]:range(cyspdata_j$day)[2], as.vector(dfit[L,]), type = 'l', col = rgb(0,0,0,0.025))
      }
      for (L in 1:((ITER/2)*CHAINS))
      {
        abline(v = halfmax_fit[L], col = rgb(1,0,0,0.05))
      }
      dev.off()
      ########################
    }
    counter <- counter + 1
  } #k
} #j

#order by year then cell
OUT <- halfmax_df[order(halfmax_df[,'year'], halfmax_df[,'cell']),]

#save to rds object
setwd(paste0(dir, '/Bird_Phenology/Data/Processed/halfmax_juv_', RUN_DATE))
saveRDS(OUT, file = paste0('halfmax_juv_', args, '.rds'))


# runtime -----------------------------------------------------------------

time <- proc.time() - tt
rtime <- round(time[3]/60, 2)
paste0('Runtime (minutes): ', rtime)


print('I completed!')
