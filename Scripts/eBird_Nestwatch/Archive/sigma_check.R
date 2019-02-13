#checks and compares IAR models using different parameteriations for sigma (scaling parameter for combined spatial and non-spatial random effects) - DEPRECATED AS THIS IAR MODEL NO LONGER BEING USED



# packages -----------------------------------------------------------------

library(MCMCvis)


# Read in and proces - original -------------------------------------------

setwd("~/Google_Drive/R/Bird_Phenology/Data/Processed/IAR_output_2018-11-12")

fls <- list.files()[grep('rds', list.files())]

C_minimus <- readRDS(fls[1])
E_virescens <- readRDS(fls[2])
V_olivaceus <- readRDS(fls[3])

med_C_minimus <- MCMCpstr(C_minimus, params = 'mu', func = median)[[1]]
sd_C_minimus <- MCMCpstr(C_minimus, params = 'mu', func = sd)[[1]]


med_E_virescens <- MCMCpstr(E_virescens, params = 'mu', func = median)[[1]]
sd_E_virescens <- MCMCpstr(E_virescens, params = 'mu', func = sd)[[1]]


med_V_olivaceus <- MCMCpstr(V_olivaceus, params = 'mu', func = median)[[1]]
sd_V_olivaceus <- MCMCpstr(V_olivaceus, params = 'mu', func = sd)[[1]]



# Read in and proces - sigma -------------------------------------------

setwd("~/Google_Drive/R/Bird_Phenology/Data/Processed/Sigma")

fls_s <- list.files()

C_minimus_sigma <- readRDS(fls_s[1])
E_virescens_sigma <- readRDS(fls_s[2])
V_olivaceus_sigma <- readRDS(fls_s[3])

med_C_minimus_sigma <- MCMCpstr(C_minimus_sigma, params = 'mu', func = median)[[1]]
sd_C_minimus_sigma <- MCMCpstr(C_minimus_sigma, params = 'mu', func = sd)[[1]]


med_E_virescens_sigma <- MCMCpstr(E_virescens_sigma, params = 'mu', func = median)[[1]]
sd_E_virescens_sigma <- MCMCpstr(E_virescens_sigma, params = 'mu', func = sd)[[1]]


med_V_olivaceus_sigma <- MCMCpstr(V_olivaceus_sigma, params = 'mu', func = median)[[1]]
sd_V_olivaceus_sigma <- MCMCpstr(V_olivaceus_sigma, params = 'mu', func = sd)[[1]]




# Read in and proces - sigma sigma -------------------------------------------

setwd("~/Google_Drive/R/Bird_Phenology/Data/Processed/Sigma_sigma")

fls_ss <- list.files()

C_minimus_sigma2 <- readRDS(fls_ss[1])
E_virescens_sigma2 <- readRDS(fls_ss[2])
V_olivaceus_sigma2 <- readRDS(fls_ss[3])

med_C_minimus_sigma2 <- MCMCpstr(C_minimus_sigma2, params = 'mu', func = median)[[1]]
sd_C_minimus_sigma2 <- MCMCpstr(C_minimus_sigma2, params = 'mu', func = sd)[[1]]


med_E_virescens_sigma2 <- MCMCpstr(E_virescens_sigma2, params = 'mu', func = median)[[1]]
sd_E_virescens_sigma2 <- MCMCpstr(E_virescens_sigma2, params = 'mu', func = sd)[[1]]


med_V_olivaceus_sigma2 <- MCMCpstr(V_olivaceus_sigma2, params = 'mu', func = median)[[1]]
sd_V_olivaceus_sigma2 <- MCMCpstr(V_olivaceus_sigma2, params = 'mu', func = sd)[[1]]





# C_minimus diff ----------------------------------------------------------

par(mfrow = c(3,2))

#difference between original mean and sigma mean
hist(med_C_minimus - med_C_minimus_sigma, 
     main = 'C_minimus - MEAN: orig - sigma',
     xlab = 'Difference')
hist(sd_C_minimus - sd_C_minimus_sigma, 
     main = 'C_minimus - SD: orig - sigma',
     xlab = 'Difference')

#difference between original mean and sigma2 mean
hist(med_C_minimus - med_C_minimus_sigma2, main = 'C_minimus - MEAN: orig - sigma2')
hist(sd_C_minimus - sd_C_minimus_sigma2, main = 'C_minimus - SD: orig - sigma2')

# E_virescens diff ----------------------------------------------------------

#difference between original mean and sigma mean
hist(med_E_virescens - med_E_virescens_sigma, 
     main = 'E_virescens - MEAN: orig - sigma',
     xlab = 'Difference')
hist(sd_E_virescens - sd_E_virescens_sigma, 
     main = 'E_virescens - SD: orig - sigma',
     xlab = 'Difference')

#difference between original mean and sigma2 mean
hist(med_E_virescens - med_E_virescens_sigma2, main = 'E_virescens - MEAN: orig - sigma2')
hist(sd_E_virescens - sd_E_virescens_sigma2, main = 'E_virescens - SD: orig - sigma2')

# V_olivaceus diff ----------------------------------------------------------

#difference between original mean and sigma mean
hist(med_V_olivaceus - med_V_olivaceus_sigma, 
     main = 'V_olivaceus - MEAN: orig - sigma',
     xlab = 'Difference')
hist(sd_V_olivaceus - sd_V_olivaceus_sigma, 
     main = 'V_olivaceus - SD: orig - sigma',
     xlab = 'Difference')

#difference between original mean and sigma2 mean
hist(med_V_olivaceus - med_V_olivaceus_sigma2, main = 'V_olivaceus - MEAN: orig - sigma2')
hist(sd_V_olivaceus - sd_V_olivaceus_sigma2, main = 'V_olivaceus - SD: orig - sigma2')




# PPO - C_minimus - original ---------------------------------------------------------

#beta0
PR <- rnorm(10000, 120, 10)
MCMCtrace(C_minimus, 
          params = 'beta0',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#rho
PR <- rbeta(10000, 0.5, 0.5)
MCMCtrace(C_minimus, 
          params = 'rho',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#sigma
PR <- rnorm(10000, 0, 3)
PR <- PR[which(PR > 0)]
MCMCtrace(C_minimus, 
          params = 'sigma',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)



# PPO - C_minimus - sigma ---------------------------------------------------------

#beta0
PR <- rnorm(10000, 120, 10)
MCMCtrace(C_minimus_sigma, 
          params = 'beta0',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#rho
PR <- rbeta(10000, 0.5, 0.5)
MCMCtrace(C_minimus_sigma, 
          params = 'rho',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#mu_sigma
PR <- rnorm(10000, 0, 3)
PR <- PR[which(PR > 0)]
MCMCtrace(C_minimus_sigma, 
          params = 'mu_sigma',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)


# PPO - C_minimus - sigma2 ---------------------------------------------------------

#beta0
PR <- rnorm(10000, 120, 10)
MCMCtrace(C_minimus_sigma2, 
          params = 'beta0',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#rho
PR <- rbeta(10000, 0.5, 0.5)
MCMCtrace(C_minimus_sigma2, 
          params = 'rho',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#mu_sigma
PR <- rnorm(10000, 0, 3)
PR <- PR[which(PR > 0)]
MCMCtrace(C_minimus_sigma2, 
          params = 'mu_sigma',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#sigma_sigma
PR <- rnorm(10000, 0, 3)
PR <- PR[which(PR > 0)]
MCMCtrace(C_minimus_sigma2, 
          params = 'sigma_sigma',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

# PPO - E_virescens - original ---------------------------------------------------------

#beta0
PR <- rnorm(10000, 120, 10)
MCMCtrace(E_virescens, 
          params = 'beta0',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#rho
PR <- rbeta(10000, 0.5, 0.5)
MCMCtrace(E_virescens, 
          params = 'rho',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#sigma
PR <- rnorm(10000, 0, 3)
PR <- PR[which(PR > 0)]
MCMCtrace(E_virescens, 
          params = 'sigma',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)



# PPO - E_virescens - sigma ---------------------------------------------------------

#beta0
PR <- rnorm(10000, 120, 10)
MCMCtrace(E_virescens_sigma, 
          params = 'beta0',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#rho
PR <- rbeta(10000, 0.5, 0.5)
MCMCtrace(E_virescens_sigma, 
          params = 'rho',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#mu_sigma
PR <- rnorm(10000, 0, 3)
PR <- PR[which(PR > 0)]
MCMCtrace(E_virescens_sigma, 
          params = 'mu_sigma',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)


# PPO - E_virescens - sigma2 ---------------------------------------------------------

#beta0
PR <- rnorm(10000, 120, 10)
MCMCtrace(E_virescens_sigma2, 
          params = 'beta0',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#rho
PR <- rbeta(10000, 0.5, 0.5)
MCMCtrace(E_virescens_sigma2, 
          params = 'rho',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#mu_sigma
PR <- rnorm(10000, 0, 3)
PR <- PR[which(PR > 0)]
MCMCtrace(E_virescens_sigma2, 
          params = 'mu_sigma',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#sigma_sigma
PR <- rnorm(10000, 0, 3)
PR <- PR[which(PR > 0)]
MCMCtrace(E_virescens_sigma2, 
          params = 'sigma_sigma',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)


# PPO - V_olivaceus - original ---------------------------------------------------------

#beta0
PR <- rnorm(10000, 120, 10)
MCMCtrace(V_olivaceus, 
          params = 'beta0',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#rho
PR <- rbeta(10000, 0.5, 0.5)
MCMCtrace(V_olivaceus, 
          params = 'rho',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#sigma
PR <- rnorm(10000, 0, 3)
PR <- PR[which(PR > 0)]
MCMCtrace(V_olivaceus, 
          params = 'sigma',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)



# PPO - V_olivaceus - sigma ---------------------------------------------------------

#beta0
PR <- rnorm(10000, 120, 10)
MCMCtrace(V_olivaceus_sigma, 
          params = 'beta0',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#rho
PR <- rbeta(10000, 0.5, 0.5)
MCMCtrace(V_olivaceus_sigma, 
          params = 'rho',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#mu_sigma
PR <- rnorm(10000, 0, 3)
PR <- PR[which(PR > 0)]
MCMCtrace(V_olivaceus_sigma, 
          params = 'mu_sigma',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)


# PPO - V_olivaceus - sigma2 ---------------------------------------------------------

#beta0
PR <- rnorm(10000, 120, 10)
MCMCtrace(V_olivaceus_sigma2, 
          params = 'beta0',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#rho
PR <- rbeta(10000, 0.5, 0.5)
MCMCtrace(V_olivaceus_sigma2, 
          params = 'rho',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#mu_sigma
PR <- rnorm(10000, 0, 3)
PR <- PR[which(PR > 0)]
MCMCtrace(V_olivaceus_sigma2, 
          params = 'mu_sigma',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#sigma_sigma
PR <- rnorm(10000, 0, 3)
PR <- PR[which(PR > 0)]
MCMCtrace(V_olivaceus_sigma2, 
          params = 'sigma_sigma',
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)
