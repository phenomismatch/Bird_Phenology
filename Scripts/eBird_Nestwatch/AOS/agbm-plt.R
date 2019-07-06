
INPUT <- fit5

RNG <- rng_beta_gamma

mu_rep_mn <- MCMCvis::MCMCpstr(INPUT, params = 'mu_rep', 
                               func = mean)[[1]]
mu_rep_LCI <- MCMCvis::MCMCpstr(INPUT, params = 'mu_rep', 
                                func = function(x) quantile(x, probs = 0.025))[[1]]
mu_rep_UCI <- MCMCvis::MCMCpstr(INPUT, params = 'mu_rep', 
                                func = function(x) quantile(x, probs = 0.975))[[1]]

# #data to plot model fit
# FIT_PLOT <- data.frame(x_sim = seq(from = RNG[1], 
#                                    to = RNG[2], length.out = 100), 
#                        mu_rep_mn, mu_rep_LCI, mu_rep_UCI)

#mean and +- 1 sd for data
y_true_mn <- MCMCvis::MCMCpstr(INPUT, params = 'y_', func = mean)[[1]]
y_true_sd <- MCMCvis::MCMCpstr(INPUT, params = 'y_true', func = sd)[[1]]
y_true_LCI <- y_true_mn - y_true_sd
y_true_UCI <- y_true_mn + y_true_sd

x_true_mn <- MCMCvis::MCMCpstr(INPUT, params = 'x_true', func = mean)[[1]]
x_true_sd <- MCMCvis::MCMCpstr(INPUT, params = 'x_true', func = sd)[[1]]
x_true_LCI <- x_true_mn - x_true_sd
x_true_UCI <- x_true_mn + x_true_sd



#mean and +- 1 sd for data
y_obs_mn <- DATA$y_obs
y_obs_sd <- DATA$sd_y
y_obs_LCI <- y_obs_mn - y_obs_sd
y_obs_UCI <- y_obs_mn + y_obs_sd

x_obs_mn <- DATA$x_obs
x_obs_sd <- DATA$sd_x
x_obs_LCI <- x_obs_mn - x_obs_sd
x_obs_UCI <- x_obs_mn + x_obs_sd


DATA_PLOT <- data.frame(y_true_mn, y_true_LCI, y_true_UCI,
                        x_true_mn, x_true_LCI, x_true_UCI)
DATA_PLOT <- data.frame(y_obs_mn, y_obs_LCI, y_obs_UCI,
                        x_obs_mn, x_obs_LCI, x_obs_UCI)

plt <- ggplot(data = DATA_PLOT, aes(x_obs_mn, y_obs_mn), color = 'black', alpha = 0.6) +
  # geom_ribbon(data = FIT_PLOT, 
  #             aes(x = x_sim, ymin = mu_rep_LCI, ymax = mu_rep_UCI),
  #             fill = 'grey', alpha = 0.6,
  #             inherit.aes = FALSE) +
  # geom_line(data = FIT_PLOT, aes(x_sim, mu_rep_mn), color = 'red',
  #           alpha = 0.9,
  #           inherit.aes = FALSE,
  #           size = 1.4) +
  geom_point(data = DATA_PLOT,
             aes(x_obs_mn, y_obs_mn), color = 'black',
             inherit.aes = FALSE, size = 5, alpha = 0.4) +
  geom_errorbar(data = DATA_PLOT,
                aes(ymin = y_obs_LCI, ymax = y_obs_UCI), #width = 0.05,
                color = 'black', alpha = 0.2) +
  geom_errorbarh(data = DATA_PLOT,
                 aes(xmin = x_obs_LCI, xmax = x_obs_UCI), #height = 0.05,
                 color = 'black', alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
    axis.ticks.length= unit(0.2, 'cm')) #length of axis tick
setwd('~/Desktop')
ggsave(paste0('ag_bg.pdf'), plt)
