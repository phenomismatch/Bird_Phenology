####################
# 8 - logit-cubic-breeding
#
# *logit cubic for ebird breeding codes
#
####################


# load packages -----------------------------------------------------------

library(rstanarm)
library(dplyr)


# directory ---------------------------------------------------------------

#dir <- '~/Google_Drive/R/'

#Xanadu
dir <- '/home/CAM/cyoungflesh/phenomismatch/'



# get args ----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
#args <- 'Vireo_olivaceus'

RUN_DATE <- '2019-01-16'



# read in data ------------------------------------------------------------


#IAR input data (to get relevant cells)
DATE_MA <- '2018-11-12'

setwd(paste0(dir, 'Bird_phenology/Data/Processed/IAR_input_', DATE_MA))
df_master <- readRDS(paste0('IAR_input-', DATE_MA, '.rds'))



#read in ebird breeding code data
DATE_BC <- '2019-01-15'

setwd(paste0(dir, 'Bird_phenology/Data/Processed/breeding_cat_query_', DATE_BC))
temp_bc <- readRDS(paste0('ebird_NA_breeding_cat_', args, '.rds'))
temp_master <- dplyr::filter(df_master, species == args)



# model -------------------------------------------------------------------

#Breeding categories - http://www.ctbirdatlas.org/Surveys-Breeding-codes.htm
#C1 - bird observed, not breeding (though there are many surveys without a breeding category which I believe would fulfill C1)
#C2 - breeding possible
#C3 - breeding probable
#C4 - breeding confirmed


#model settings
#ITER <- 2500
ITER <- 10
CHAINS <- 4


#only cells that are in IAR input
kp_cells <- unique(temp_master$cell)
temp_bc_f <- dplyr::filter(temp_bc, cell %in% kp_cells)

yrs <- unique(temp_bc_f$year)
nyr <- length(yrs)

cells <- unique(temp_bc_f$cell)
ncell <- length(cells)

t_mat <- matrix(data = NA, nrow = ncell*nyr, ncol = ((ITER/2)*CHAINS))
colnames(t_mat) <- paste0('iter_', 1:((ITER/2)*CHAINS))
halfmax_df <- data.frame(species = args, 
                         year = rep(yrs, each = ncell), 
                         cell = rep(cells, nyr), 
                         max_Rhat = NA,
                         min_neff = NA,
                         sh = NA,
                         t_mat)


#new data for plotting model fit
predictDays <- scale(c(1:200), scale = FALSE)
predictDays2 <- scale(c(1:200)^2, scale = FALSE)
predictDays3 <- scale(c(1:200)^3, scale = FALSE)
newdata <- data.frame(sjday = predictDays, sjday2 = predictDays2, sjday3 = predictDays3, shr = 0)


counter <- 1
for (j in 1:nyr)
{
  #j <- 13
  t_yr <- dplyr::filter(temp_bc_f, year == yrs[j])
  
  for (k in 1:ncell)
  {
    #k <- 44
    t_cell <- dplyr::filter(t_yr, cell == cells[k])
    
    #new column with 'probable' or 'confirmed' breeding
    t_cell$br <- as.numeric(t_cell$bba_category == 'C3' | t_cell$bba_category == 'C4')
    
    #remove instances where bird was not observed
    # to.rm <- which(t_cell$br == 0)
    # t_cell2 <- t_cell[-to.rm,]
    
    #sample only some breeding obs to look at feasibility of lower threshold (half original)
    # breeds <- which(t_cell$br == 1)
    # stm <- sample(which(t_cell$br == 1), 10)
    # t_cell2 <- t_cell[-stm,]
    
    t_cell2 <- t_cell
    
    #bird not seen - fill with zeros
    na.ind <- which(is.na(t_cell2$br))
    t_cell2$br[na.ind] <- 0
    
    #number of surveys where breeding was detected (confirmed or probable)
    n1 <- sum(t_cell2$br)
    #number of surveys where breeding was not detected (bird not seen breeding or not seen)
    n0 <- sum(t_cell2$br == 0)
    
    #check
    #n0+n1 == NROW(t_cell2)
    
    #number of detections that came before jday 60
    n1W <- sum(t_cell2$br * as.numeric(t_cell2$day < 60))
    
    #number of unique days with detections
    njd1 <- length(unique(t_cell2$day[which(t_cell2$br == 1)]))
    #number of unique days of non-detections before first detection
    njd0i <- length(unique(t_cell2$day[which(t_cell2$br == 0 & t_cell2$day < 
                                               min(t_cell2$day[which(t_cell2$br == 1)]))]))
    
    print(paste0('species: ', args, ', year: ', j, ', cell: ', k, ', br obs: ', n1))
    
    #different thresholds from arrival models
    if (n1 > 20 & n1W < (n1/50) & n0 > 29 & njd1 > 15 & njd0i > 29)
    {
      fit2 <- rstanarm::stan_glm(br ~ sjday + sjday2 + sjday3 + shr,
                                 data = t_cell2,
                                 family = binomial(link = "logit"),
                                 algorithm = 'sampling',
                                 iter = ITER,
                                 chains = CHAINS,
                                 cores = CHAINS)
      
      dfit <- rstanarm::posterior_linpred(fit2, newdata = newdata, transform = T)
      halfmax_fit <- rep(NA, ((ITER/2)*CHAINS))
      
      for (L in 1:((ITER/2)*CHAINS))
      {
        #L <- 3000
        rowL <- as.vector(dfit[L,])
        halfmax_fit[L] <- min(which(rowL > (max(rowL)/2)))
      }
      
      ########################
      #PLOT MODEL FIT AND DATA
      
      #summary(fit2)
      setwd(paste0(dir, 'Bird_Phenology/Results/Plots'))
      mn_dfit <- apply(dfit, 2, mean)
      LCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.025))
      UCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.975))
      mn_hm <- mean(halfmax_fit)
      LCI_hm <- quantile(halfmax_fit, probs = 0.025)
      UCI_hm <- quantile(halfmax_fit, probs = 0.975)
      
      pdf(paste0(args, '_', yrs[j], '_', cells[k], '_breeding.pdf'))
      plot(UCI_dfit, type = 'l', col = 'red', lty = 2, lwd = 2,
           ylim = c(0, max(UCI_dfit)),
           main = paste0(args, ' - ', yrs[j], ' - ', cells[k]),
           xlab = 'Julian Day', ylab = 'Detection Probability')
      lines(LCI_dfit, col = 'red', lty = 2, lwd = 2)
      lines(mn_dfit, lwd = 2)
      t_cell2$br[which(t_cell2$br == 1)] <- max(UCI_dfit)
      points(t_cell2$day, t_cell2$br, col = rgb(0,0,0,0.25))
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
      halfmax_df$sh[counter] <- round(shapiro.test(halfmax_fit)$p.value, 2)
    }
    counter <- counter + 1
  } #k
} #j


setwd(paste0(dir, '/Bird_Phenology/Data/Processed/halfmax_breeding_', RUN_DATE))
saveRDS(halfmax_df, file = paste0('halfmax_df_breeding_', args, '.rds'))



