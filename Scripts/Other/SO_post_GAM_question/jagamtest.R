library(mgcv)
library(MASS)
library(rjags)

b <- gam(accel~s(times,k=40), data=mcycle, method = 'REML') 
jb <- jagam(accel~s(times,k=40), data=mcycle, file = 'test.jags') 


#fit JAGS model
jm <- jags.model("test.jags", data = jb$jags.data,
                   inits = jb$jags.ini, n.adapt = 2000, n.chains = 4)
update(jm, n.iter = 2000)
out <- coda.samples(jm, variable.names = c("b", "rho", "scale", "mu"), 
                    n.iter = 10000, thin = 10)

#check model
MCMCvis::MCMCsummary(out, round = 3)

#mean
mu_ch <- MCMCvis::MCMCchains(out, params = 'mu')
#95% CI
mu_mn <- MCMCvis::MCMCpstr(out, params = 'mu')[[1]]
mu_sd <- MCMCvis::MCMCpstr(out, params = 'mu', func = sd)[[1]]


#plot fit
plot(mcycle$times, mcycle$accel, pch = 19, 
     col = rgb(0,0,0,0.3))
#new data to predict at
pd <- data.frame(times = seq(min(mcycle$times),
                             max(mcycle$times),
                             length = 500))
lines(mcycle$times, mu_mn, col = 'green', lwd = 2)
lines(mcycle$times, mu_mn + 1.96 * mu_sd, col = 'green', lwd = 2, lty = 2)
lines(mcycle$times, mu_mn - 1.96 * mu_sd, col = 'green', lwd = 2, lty = 2)
for (i in 1:100)
{
  lines(mcycle$times, mu_ch[i,], col = rgb(1,0,0,0.1))
}
