#look at pairs
pairs_stan <- function(chain, stan_model, pars)
{
  energy <- as.matrix(sapply(get_sampler_params(stan_model, inc_warmup = F), 
                             function(x) x[,"energy__"]))
  pars <- extract(stan_model, pars = pars, permuted = F)
  df <- data.frame(energy[,chain], pars[,chain,])
  names(df)[1] <- "energy"
  GGally::ggpairs(df, title = paste0("Chain", chain), 
                  lower = list(continuous = GGally::wrap("points", alpha = 0.2)))                    
}
pairs_stan(3, fit, c('sigma_beta0', 'sigma_phi', 'sigma_y_true', 'sigma_gamma', 
                     'alpha_gamma', 'beta_gamma'))


#There is a trade-off between phi[i,j] and y_true[i,j] - this is the cause of the tree depth exceeds

#Results from normal IAR and and BYM2-style models are nearly identical



setwd('~/Desktop/')
fit <- readRDS('Vireo_olivaceus-bym-test-2019-11-14.rds')

yt_ch <- MCMCvis::MCMCchains(fit, params = 'y_true')
phi_ch <- MCMCvis::MCMCchains(fit, params = 'phi')
#plot(yt_ch, phi_ch, pch = '.')

#DATA$ii_obs[11:12,1:20]
#data not avail
par(mfrow = c(2,2))
grep('y_true[12,7]', colnames(yt_ch), fixed = TRUE)
plot(yt_ch[,975], phi_ch[,975], pch = '.', main = 'BYM nd')

#data avail
grep('y_true[12,8]', colnames(yt_ch), fixed = TRUE)
plot(yt_ch[,976], phi_ch[,976], pch = '.', main = 'BYM d')

#BYM2 model
setwd("~/Google_Drive/R/Bird_Phenology/Data/Processed/IAR_output_2019-05-26")
fit2 <- readRDS('Vireo_olivaceus-2019-05-26-iar-stan_output.rds')

yt_ch <- MCMCvis::MCMCchains(fit2, params = 'y_true')
phi_ch <- MCMCvis::MCMCchains(fit2, params = 'phi')
#plot(yt_ch, phi_ch, pch = '.')

#data not avail
grep('y_true[7,12]', colnames(yt_ch), fixed = TRUE)
plot(yt_ch[,114], phi_ch[,114], pch = '.', main = 'BYM2 nd')

#data avail
grep('y_true[8,12]', colnames(yt_ch), fixed = TRUE)
plot(yt_ch[,131], phi_ch[,131], pch = '.', main = 'BYM2 d')



#POSTERIOR SUMMARIES
#BYM
#data not avail
MCMCvis::MCMCsummary(fit, params = 'y_true\\[12,7\\]', 
                     ISB = FALSE, round = 2)
#data avail
MCMCvis::MCMCsummary(fit, params = 'y_true\\[12,8\\]', 
                     ISB = FALSE, round = 2)

#BYM2
#data not avail
MCMCvis::MCMCsummary(fit2, params = 'y_true\\[7,12\\]', 
                     ISB = FALSE, round = 2)
#data avail
MCMCvis::MCMCsummary(fit2, params = 'y_true\\[8,12\\]', 
                     ISB = FALSE, round = 2)




#ADDITIONAL CHECKS
#IAR
yt_f <- MCMCvis::MCMCpstr(fit, params = 'y_true', fun = mean)[[1]]
#BYM2-style
yt_f2 <- MCMCvis::MCMCpstr(fit2, params = 'y_true', fun = mean)[[1]]
#transpose because put years in rows in later models
yt_f22 <- t(yt_f2)
#difference is very slight
hist(yt_f - yt_f22)

#IAR
yt_f <- MCMCvis::MCMCpstr(fit, params = 'y_true', fun = var)[[1]]
#BYM2-style
yt_f2 <- MCMCvis::MCMCpstr(fit2, params = 'y_true', fun = var)[[1]]
#transpose because put years in rows in later models
yt_f22 <- t(yt_f2)
#difference is very slight
hist(yt_f - yt_f22)

#IAR
yt_f <- MCMCvis::MCMCpstr(fit, params = 'beta0', fun = mean)[[1]]
#BYM2-style
yt_f2 <- MCMCvis::MCMCpstr(fit2, params = 'beta0', fun = mean)[[1]]
#transpose because put years in rows in later models
yt_f22 <- t(yt_f2)
#difference is very slight
hist(yt_f - yt_f22)
MCMCvis::MCMCsummary(fit, params = 'beta0', round = 2)
MCMCvis::MCMCsummary(fit2, params = 'beta0', round = 2)


#IAR
yt_f <- MCMCvis::MCMCpstr(fit, params = 'beta_gamma', fun = mean)[[1]]
#BYM2-style
yt_f2 <- MCMCvis::MCMCpstr(fit2, params = 'beta_gamma', fun = mean)[[1]]
#transpose because put years in rows in later models
yt_f22 <- t(yt_f2)
#difference is very slight
hist(yt_f - yt_f22)
# MCMCvis::MCMCsummary(fit, params = 'beta_gamma', round = 2)
# MCMCvis::MCMCsummary(fit2, params = 'beta_gamma', round = 2)

