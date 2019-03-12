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
pairs_stan(3, fit, c('beta0', 'sigma_beta', 'mu_beta', 'sigma_nu', 'rho', 'mu_sn'))

# pars = c('alpha', 'mu_alpha', 'sigma_alpha', 'beta',
#          'sigma_beta', 'alpha2', 'beta2',
#          'sigma_x_true', 'x_true'),