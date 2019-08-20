######################
# 5 - bern model juvs MAPS
#
######################



# runtime -----------------------------------------------------------------

tt <- proc.time()


# Top-level dir -----------------------------------------------------------

#desktop/laptop
#dir <- '~/Google_Drive/R/'

#Xanadu
dir <- '/UCHC/LABS/Tingley/phenomismatch/'


# other dir ------------------------------------------------------------

#IAR data
arr_master_dir <- 'arrival_master_2019-05-26'

#run date
RUN_DATE <- '2019-08-20'



# load packages -----------------------------------------------------------

library(dplyr)
library(plyr)
library(dggridR)
library(rstanarm)
library(rstan)



# Get args fed to script --------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
#args <- 'Vireo_olivaceus'
#args <- 'Agelaius_phoeniceus'
#args <- 'Dumetella_carolinensis'



# model settings ----------------------------------------------------------

#model settings
ITER <- 3000
CHAINS <- 4

# ITER <- 10
# CHAINS <- 1



# read in data ------------------------------------------------------------

#read in - RDS create with 1-query-db.R in wing_chord_changes project
setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))

data <- readRDS('MAPS-age-filled.rds')
colnames(data)[grep('sci_name', colnames(data))] <- 'species'
#add underscore to species naems
data$species <- gsub(' ', '_', data$species)

#master arrival data (from IAR output)
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', arr_master_dir))

arr_master <- readRDS(paste0(arr_master_dir, '.rds'))



# create grid -------------------------------------------------------------

hexgrid6 <- dggridR::dgconstruct(res = 6)
data$cell <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                         in_lon_deg = data$lng, 
                                         in_lat_deg = data$lat)[[1]]




# filter data -------------------------------------------------------------
head(arr_master)
#filter by species, breeding range, and pre-IAR data
arr_br_range <- dplyr::filter(arr_master, species == args, breed_cell == TRUE, !is.na(mean_pre_IAR))

#keep only MAPS obs where there are IAR arrival estimates
m_mrg <- dplyr::inner_join(data, arr_br_range, by = c('species', 'cell', 'year'))
m_mf <- dplyr::select(m_mrg, common_name, species, cell, 
                     year, day, age, true_age, mean_pre_IAR, 
                     sd_pre_IAR, mean_post_IAR, sd_post_IAR)

#add col for juveniles (1 = juv, 0 = adult)
m_mf$juv <- NA
m_mf$juv[which(m_mf$true_age == 0)] <- 1
m_mf$juv[-which(m_mf$true_age == 0)] <- 0

years <- sort(unique(m_mf$year))
cells <- sort(unique(m_mf$cell))
ncell <- length(cells)
nyrs <- length(years)



# setup data object -------------------------------------------------------

t_mat <- matrix(data = NA, nrow = ncell*nyrs, ncol = ((ITER/2)*CHAINS))
colnames(t_mat) <- paste0('iter_', 1:((ITER/2)*CHAINS))
halfmax_df <- data.frame(species = args, 
                         year = rep(years, each = ncell), 
                         cell = rep(cells, nyrs), 
                         max_Rhat = NA,
                         min_neff = NA,
                         num_diverge = NA,
                         num_tree = NA,
                         num_BFMI = NA,
                         delta = NA,
                         tree_depth = NA,
                         n1 = NA,
                         # n1W = NA,
                         n0 = NA,
                         n0i = NA,
                         njd1 = NA,
                         njd0 = NA,
                         njd0i = NA,
                         t_mat)


#create dir for figs if doesn't exist
ifelse(!dir.exists(paste0(dir, 'Bird_Phenology/Figures/halfmax/juvs_', RUN_DATE)),
       dir.create(paste0(dir, 'Bird_Phenology/Figures/halfmax/juvs_', RUN_DATE)),
       FALSE)

setwd(paste0(dir, 'Bird_Phenology/Figures/halfmax/juvs_', RUN_DATE))


#loop through each species, year, cell and extract half-max parameter
counter <- 1
for (j in 1:nyrs)
{
  #j <- 1
  ydata <- dplyr::filter(m_mf, year == years[j])
  
  for (k in 1:ncell)
  {
    #k <- 1
    print(paste0('species: ', args, ', year: ', j, ', cell: ', k))
    
    cydata <- dplyr::filter(m_mf, cell == cells[k])
    
    #number of surveys where brood patch was detected
    n1 <- sum(cydata$juv)
    #number of surveys where brood patch was not detected
    n0 <- sum(cydata$juv == 0)
    # #number of detections that came before jday 60
    # n1W <- sum(cydata$juv * as.numeric(cydata$day < 60))
    #number of unique days with detections
    njd1 <- length(unique(cydata$day[which(cydata$juv == 1)]))
    #number of unique days with non-detection
    njd0 <- length(unique(cydata$day[which(cydata$juv == 0)]))
    
    
    if (n1 > 0)
    {
      #number of unique days of non-detections before first detection
      njd0i <- length(unique(cydata$day[which(cydata$juv == 0 & cydata$day < 
                                                  min(cydata$day[which(cydata$juv == 1)]))]))
      #number of non-detections before first detection
      n0i <- length(which(cydata$juv == 0 & 
                            cydata$day < min(cydata$day[which(cydata$juv == 1)])))
    } else {
      njd0i <- 0
      n0i <- 0
    }
    
    halfmax_df$n1[counter] <- n1
    # halfmax_df$n1W[counter] <- n1W
    halfmax_df$n0[counter] <- n0
    halfmax_df$n0i[counter] <- n0i
    halfmax_df$njd1[counter] <- njd1
    halfmax_df$njd0[counter] <- njd0
    halfmax_df$njd0i[counter] <- njd0i
    
    #defaults for rstanarm are 0.95 and 15
    DELTA <- 0.95
    TREE_DEPTH <- 15
    
    #br thresholds
    #if (n1 > 29 & n1W < (n1/50) & n0 > 29 & njd0i > 29 & njd1 > 19)
    
    if (n1 > 15 & n0 > 30 & njd0i > 10 & njd1 > 10)
    {
      fit2 <- rstanarm::stan_gamm4(juv ~ s(day), 
                                   data = cydata,
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
      
      # #rerun model if things didn't go well
      # while (sum(c(num_diverge, num_tree, num_BFMI)) > 0 & DELTA <= 0.98)
      # {
      #   DELTA <- DELTA + 0.01
      #   TREE_DEPTH <- TREE_DEPTH + 1
      #   
      #   fit2 <- rstanarm::stan_gamm4(bp ~ s(day),
      #                                data = cydata,
      #                                family = binomial(link = "logit"),
      #                                algorithm = 'sampling',
      #                                iter = ITER,
      #                                chains = CHAINS,
      #                                cores = CHAINS,
      #                                adapt_delta = DELTA,
      #                                control = list(max_treedepth = TREE_DEPTH))
      #   
      #   num_diverge <- rstan::get_num_divergent(fit2$stanfit)
      #   num_tree <- rstan::get_num_max_treedepth(fit2$stanfit)
      #   num_BFMI <- length(rstan::get_low_bfmi_chains(fit2$stanfit))
      # }
      
      halfmax_df$num_diverge[counter] <- num_diverge
      halfmax_df$num_tree[counter] <- num_tree
      halfmax_df$num_BFMI[counter] <- num_BFMI
      halfmax_df$delta[counter] <- DELTA
      halfmax_df$tree_depth[counter] <- TREE_DEPTH
      
      #generate predict data
      predictDays <- range(cydata$day)[1]:range(cydata$day)[2]
      newdata <- data.frame(day = predictDays)
      
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
      
      cydata2 <- cydata
      #summary(fit2)
      mn_dfit <- apply(dfit, 2, mean)
      LCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.025))
      UCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.975))
      mn_hm <- mean(halfmax_fit)
      LCI_hm <- quantile(halfmax_fit, probs = 0.025)
      UCI_hm <- quantile(halfmax_fit, probs = 0.975)
      
      pdf(paste0(args, '_', years[j], '_', cells[k], '_juvs.pdf'))
      plot(predictDays, UCI_dfit, type = 'l', col = 'red', lty = 2, lwd = 2,
           ylim = c(0, max(UCI_dfit)),
           main = paste0(args, ' - ', years[j], ' - ', cells[k]),
           xlab = 'Julian Day', ylab = 'Brood Patch Probability')
      lines(predictDays, LCI_dfit, col = 'red', lty = 2, lwd = 2)
      lines(predictDays, mn_dfit, lwd = 2)
      cydata2$juv[which(cydata2$juv == 1)] <- max(UCI_dfit)
      points(cydata2$day, cydata2$juv, col = rgb(0,0,0,0.25))
      abline(v = mn_hm, col = rgb(0,0,1,0.5), lwd = 2)
      abline(v = LCI_hm, col = rgb(0,0,1,0.5), lwd = 2, lty = 2)
      abline(v = UCI_hm, col = rgb(0,0,1,0.5), lwd = 2, lty = 2)
      # legend('topleft',
      #        legend = c('Cubic fit', 'CI fit', 'Half max', 'CI HM'),
      #        col = c('black', 'red', rgb(0,0,1,0.5), rgb(0,0,1,0.5)),
      #        lty = c(1,2,1,2), lwd = c(2,2,2,2), cex = 1.3)
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
setwd(paste0(dir, '/Bird_Phenology/Data/Processed/halfmax_juvs_', RUN_DATE))
saveRDS(halfmax_df, file = paste0('halfmax_juvs_', args, '.rds'))



# runtime -----------------------------------------------------------------

time <- proc.time() - tt
rtime <- round(time[3]/60, 2)
paste0('Runtime (minutes): ', rtime)


print('I completed!')


