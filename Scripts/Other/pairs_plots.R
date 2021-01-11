#look at pairs
library(rstan)

pairs_stan <- function(chain, stan_model, pars)
{
  # stan_model <- tt
  # pars <- c('sigma_beta0', 'sigma_phi', 'sigma_y_true', 'sigma_gamma', 
  #           'alpha_gamma', 'beta_gamma')
  # chain <- 3
  
  energy <- as.matrix(sapply(get_sampler_params(stan_model, inc_warmup = F), 
                             function(x) x[,"energy__"]))
  tpars <- extract(stan_model, pars = pars, permuted = F)
  df <- data.frame(energy[,chain], tpars[,chain,])
  names(df)[1] <- "energy"
  divergent <- get_sampler_params(stan_model, 
                                  inc_warmup=FALSE)[[1]][,'divergent__']
  df$divergent <- divergent
  
  nd_df <- dplyr::filter(df, divergent == 0)
  d_df <- dplyr::filter(df, divergent == 1)
  cn <- which(colnames(nd_df) == 'divergent')
  tplt <- rbind(nd_df[,-cn], d_df[,-cn])
  
  pairs(tplt, 
        col = c(rep(rgb(0,0,0,0.1), NROW(nd_df)),
        rep(rgb(1,0,0,0.1), NROW(d_df))), 
        pch = 19)
}
pairs_stan(3, tt, c('sigma_beta0', 'sigma_phi', 'sigma_y_true', 'sigma_gamma', 
                     'alpha_gamma', 'beta_gamma'))





#need to mod to display log for x when specified
plt_fun <- function(fit, p1, p2, log_p1 = FALSE, log_p2 = FALSE)
{
  #extract chains
  x_ch <- MCMCvis::MCMCchains(fit, params = p1)
  y_ch <- MCMCvis::MCMCchains(fit, params = p2)
  
  #get divergent transitions
  divergent <- get_sampler_params(fit, inc_warmup = FALSE)[[1]][,'divergent__']
  div_idx <- which(divergent == 1)
  
  if (length(div_idx) > 0)
  {
    x_ch2 <- x_ch[-div_idx, ]
    y_ch2 <- y_ch[-div_idx, ]
    x_ch_div <- x_ch[div_idx, ]
    y_ch_div <- y_ch[div_idx, ]
  } else {
    x_ch2 <- x_ch
    y_ch2 <- y_ch
    x_ch_div <- NULL
    y_ch_div <- NULL
  }
  
  if (log_p1 == TRUE)
  {
    x_ch_plt <- log(x_ch2)
    x_ch_div_plt <- log(x_ch_div)
  } else {
    x_ch_plt <- x_ch2
    x_ch_div_plt <- x_ch_div
  }
  
  if (log_p2 == TRUE)
  {
    y_ch_plt <- log(y_ch2)
    y_ch_div_plt <- log(y_ch_div)
    YLAB <- paste0('log ', p2)
  } else {
    y_ch_plt <- y_ch2
    y_ch_div_plt <- y_ch_div
    YLAB <- paste0(p2)
  }
  
  par(mfrow = c(3,3))
  for (i in 1:NCOL(x_ch2))
  {
    #i <- 1
    plot(x_ch2[,i], y_ch_plt, 
         col = rgb(0,0,0,0.1), pch = 19, 
         xlab = colnames(x_ch2)[i], ylab = YLAB)
    if (length(div_idx) > 0)
    {
      points(x_ch_div[,i], y_ch_div_plt, 
             col = rgb(0,1,0,0.3), pch = 19)
    }
  }
}


setwd("~/Desktop/cnc_Pipilo_erythrophthalmus_test")
tt <- readRDS('fit6.rds')
plt_fun(fit = tt, p1 = 'gamma', p2 = 'sigma_gamma', log_p2 = TRUE)
plt_fun(fit = tt, p1 = 'beta0', p2 = 'sigma_beta0', log_p2 = TRUE)
plt_fun(fit = tt, p1 = 'phi', p2 = 'sigma_phi', log_p2 = TRUE)
plt_fun(fit = tt, p1 = 'y_true', p2 = 'sigma_y_true', log_p2 = TRUE)

