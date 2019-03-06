library(rstanarm)
library(dplyr)

#import data for species
spdata <- readRDS('/Users/Jacob/Downloads/Phenomismatch_data/ebird_NA_phen_proc_Vireo_olivaceus.rds')
spdata$ssjday <- scale(spdata$sjday)
spdata$ssjday2 <- scale(spdata$sjday^2)
spdata$ssjday3 <- scale(spdata$sjday^3)
spdata$ssjday4 <- scale(spdata$sjday^4)
spdata$ssjday5 <- scale(spdata$sjday^5)
spdata$ssjday6 <- scale(spdata$sjday^6)
spdata$ssjday7 <- scale(spdata$sjday^7)

cells <- unique(spdata$cell)
ncel <- length(cells)

years <- min(spdata$year):max(spdata$year)
nyr <- length(years)

formulas <- list(f3 = 'detect ~ ssjday + ssjday2 + ssjday3 + shr',
                 f4 = 'detect ~ ssjday + ssjday2 + ssjday3 + ssjday4 + shr',
                 f5 = 'detect ~ ssjday + ssjday2 + ssjday3 + ssjday4 + ssjday5 + shr',
                 f6 = 'detect ~ ssjday + ssjday2 + ssjday3 + ssjday4 + ssjday5 + ssjday6 + shr',
                 f7 = 'detect ~ ssjday + ssjday2 + ssjday3 + ssjday4 + ssjday5 + ssjday6 + ssjday7 + shr')


# Create newdata ---------------------------------------------------

#calculate polynomial, then center data
predictDays <- scale(c(1:200))
predictDays2 <- scale(c(1:200)^2)
predictDays3 <- scale(c(1:200)^3)
predictDays4 <- scale(c(1:200)^4)
predictDays5 <- scale(c(1:200)^5)
predictDays6 <- scale(c(1:200)^6)
predictDays7 <- scale(c(1:200)^7)

newdata <- data.frame(ssjday = predictDays, ssjday2 = predictDays2, ssjday3 = predictDays3, 
                      ssjday4 = predictDays4, ssjday5 = predictDays5, ssjday6 = predictDays6,
                      ssjday7 = predictDays7, shr = 0)

fit_diag <- halfmax_matrix_list <- list()
ppcs <- list()



# fit logit cubic ---------------------------------------------------------

#number of iterations each model should be run
ITER <- 2000


#loop through each species, year, cell and extract half-max parameter
for(f in 1:length(formulas)){
  halfmax_matrix_list[[f]] <- list()
  fit_diag[[f]] <- list()
  ppccs[[f]] <- list()
for (j in 1:nyr)
{
  #j <- 10
  halfmax_matrix_list[[f]][[j]] <- matrix(data = NA, nrow = ncel, ncol = (ITER*2))
  fit_diag[[f]][[j]] <- list()
  ppccs[[f]][[j]] <- list()
  yspdata <- spdata[which(spdata$year == years[j]), ]
  
  for (k in 1:ncel)
  {
    #k <- 1
    print(paste(j,k))
    fit_diag[[f]][[j]][[k]] <- list(maxRhat = NA, mineffsSize = NA)
    ppcs[[f]][[j]][[k]] <- list(bin1 = NA, bin2 = NA, bin3 = NA, bin4 = NA, bin5 = NA, bin6 = NA,
                           bin7 = NA, bin8 = NA, bin9 = NA, bin10 = NA)
#    cyspdata <- yspdata[which(yspdata$cell6 == cells[k]), ]
    cyspdata <- yspdata[which(yspdata$cell == cells[k]), ]
    
    #number of surveys where species was detected
    n1 <- sum(cyspdata$detect)
    #number of surveys where species was not detected
    n0 <- sum(cyspdata$detect == 0)
    #number of detections that came before jday 60
    n1W <- sum(cyspdata$detect * as.numeric(cyspdata$day < 60))
    
    if (n1 > 29 & n1W < (n1/50) & n0 > 29)
    {
      fit2 <- stan_glm(formulas[[f]],
                       data = cyspdata,
                       family = binomial(link = "logit"),
                       algorithm = 'sampling',
                       iter = ITER,
                       chains = 4,
                       cores = 4)
      
      dfit <- posterior_linpred(fit2, newdata = newdata, transform = T)
      halfmax_fit <- rep(NA, ITER*2)
      
      for (L in 1:(ITER*2))
      {
        #L <- 1
        rowL <- as.vector(dfit[L,])
        halfmax_fit[L] <- min(which(rowL > (max(rowL)/2)))
      }
      
      pp <- rstanarm::posterior_predict(fit2)
      bin_number <- 1 + cyspdata$day %/% 20
      
      for(bin in 1:10){
        binpredict <- rowSums(pp[, bin_number == bin])
        truenumber <- sum(cyspdata$detect[bin_number == bin])
        ppcs[[f]][[j]][[k]][[bin]] <- (sum(binpredict < truenumber) + 0.5*sum(binpredict == truenumber))/4000
      }
      
      halfmax_matrix_list[[f]][[j]][k,] <- halfmax_fit
      fit_diag[[f]][[j]][[k]]$maxRhat <- max(summary(fit2)[, "Rhat"])
      fit_diag[[f]][[j]][[k]]$mineffsSize <- min(summary(fit2)[, "n_eff"])
    }
  }
}
}

