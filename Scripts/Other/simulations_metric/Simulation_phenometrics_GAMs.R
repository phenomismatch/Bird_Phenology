# Code to simulate data-varying temporal trends and look at phenological metrics using hierarchical and independent-year GAMs

# Load packages -----------------------------------------------------------
library(mgcv)
library(gratia)
library(dplyr)
library(tictoc)

# set number of sims (100 = 2 hours) ---------------------------------------------------
n.sim <- 100

# import processed data ---------------------------------------------------

#import eBird query data for species
data <- readRDS(paste0('ebird_arrival_query_Catharus_ustulatus.rds'))

#filter to 2 key cells
#617 is migratory-type functional form
#567 is breeding-type functional form
simdata1 <- dplyr::filter(data, cell == 567, year >= 2014, year <= 2018)
simdata1 <- transform(simdata1, year = factor(year, ordered = FALSE))
simdata2 <- dplyr::filter(data, cell == 617, year >= 2014, year <= 2018)
simdata2 <- transform(simdata2, year = factor(year, ordered = FALSE))

# function ----------------------------------------------------------------

hm_fun <- function(fit_res = gi_plt, data = cyspdata2)
{
  years <- unique(data$year)
  tdf <- data.frame(year = factor(years),
                  hm_jday = NA,
                  max_jday = NA,
                  max_val = NA,
                  first = NA,
                  median = NA)
  for (i in 1:length(years))
  {
    #i <- 1
    #row for specified year for data
    t_idx <- which(data$year == years[i])
    tdata <- data[t_idx,]
    #first detection
    fd <- min(tdata$jday[which(tdata$detect == 1)])
    #median detection
    md <- median(tdata$jday[which(tdata$detect == 1)])
    #row for specified year for predictions
    t_idx2 <- which(fit_res$year == years[i])
    tfit <- fit_res[t_idx2,]
    #local maximum(s)
    #from: stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
    lmax_idx <- which(diff(sign(diff(tfit$mfit))) == -2) + 1
    lmax <- tfit$jday[lmax_idx]
    
    #first local max to come after first detection
    flm <- which(lmax > fd)
    if (length(flm) > 0)
    {
      #first local max to come after first detection
      lmax2_idx <- lmax_idx[min(flm)]
      lmax2 <- lmax[min(flm)]
    } else {
      #no local max
      lmax2_idx <- which.max(tfit$mfit)
      lmax2 <- tfit$jday[which.max(tfit$mfit)]
    }
    #local mins before max (global and local mins)
    lmin_idx <- c(which.min(tfit$mfit[1:lmax2_idx]), 
                  which(diff(sign(diff(tfit$mfit[1:lmax2_idx]))) == 2) + 1)
    lmin <- tfit$jday[lmin_idx]
    #local min nearest to local max
    lmin2_idx <- lmin_idx[which.min(lmax2 - lmin)]
    lmin2 <- tfit$jday[lmin2_idx]
    
    #value at local max - value at min (typically 0)
    dmm <- tfit$mfit[lmax2_idx] - tfit$mfit[lmin2_idx]
    #all positions less than or equal to half diff between max and min + value min
    tlm <- which(tfit$mfit <= ((dmm/2) + tfit$mfit[lmin2_idx]))
    #which of these come before max and after or at min
    vgm <- tlm[which(tlm < lmax2_idx & tlm >= lmin2_idx)]
    #insert halfmax (first day for situations where max is a jday = 1)
    if (length(vgm) > 0)
    {
      tdf$hm_jday[i] <- tfit$jday[max(vgm)]
    } else {
      tdf$hm_jday[i] <- tfit$jday[1]
    }

    tdf$first[i] <- fd
    tdf$median[i] <- md
    tdf$max_jday[i] <- lmax2
    tdf$max_val[i] <- tfit$mfit[lmax2]
  }
  return(tdf)
}


# Fit the "truth" model in two cells ------------------------------------------------------------

#fit hierarchical GAMs just cause
fit_sim1 <- mgcv::gam(detect ~ s(jday, k = 30) +
                        s(jday, year, k = 30, bs = 'fs', m = 1) +
                        s(year, bs = 're', k = 30),
                      data = simdata1,
                      family = binomial(link = "logit"),
                      method = 'REML')
fit_sim2 <- mgcv::gam(detect ~ s(jday, k = 30),
                      data = simdata2,
                      family = binomial(link = "logit"),
                      method = 'REML')

#predict data
nd <- data.frame(jday = rep(1:200, length(unique(simdata1$year))), 
                 year = rep(unique(simdata1$year), each = 200))

pred_sim1 <- predict(fit_sim1, type = 'response', 
                     newdata = nd,
                     se.fit=TRUE)
pred_sim2 <- predict(fit_sim2, type = 'response', 
                     newdata = nd,
                     se.fit=TRUE)

gi_plt1 <- data.frame(nd, mfit = pred_sim1$fit, mfit_se = pred_sim1$se.fit)
hm_gi1 <- hm_fun(fit = gi_plt1, data = simdata1)
plot(probs_sim1$jday, probs_sim1$prob, type = "l", xlab = "jday", ylab = "Pr", main = "Local Breeder cell (truth)")
stats_sim1 <- hm_gi1[hm_gi1$year == max(as.character(unique(simdata1$year))), ]

# sim 2 stats
gi_plt2 <- data.frame(nd, mfit = pred_sim2$fit, mfit_se = pred_sim2$se.fit)
hm_gi2 <- hm_fun(fit = gi_plt2, data = simdata2)
probs_sim2 <- data.frame(jday = 1:200, prob = pred_sim2$fit[nd$year == max(as.character(unique(simdata2$year)))])
plot(probs_sim2$jday, probs_sim2$prob, type = "l", xlab = "jday", ylab = "Pr", main = "Passage Migrant cell (truth)")
stats_sim2 <- hm_gi2[hm_gi2$year == max(as.character(unique(simdata2$year))), ]


# Simulation loop: varying the trend and number of points for a local breeder cell ------------------------------------------------------------

# get the seasonal sampling frame for birder activity (surveys aren't equally likely across the year)
pr.survey <- table(data$jday) / dim(data)[1]

# define the number of surveys in each sim (approximately covering the range of observed variation from 2005 to 2017)
n.surveys <- seq(from = 500, to = 10000, by = 750)

#storage devices
sim.summary1 <- data.frame(model = NA,  
                           sim = NA,
                           mean_diff_hm = NA,
                           mean_diff_max = NA,
                           trend_hm = NA,
                           trend_max = NA)
sim.metrics1 <- data.frame(model = NA,
                            sim = NA,
                             year = NA,  
                             hm_jday = NA,
                             max_jday = NA,
                             max_val = NA,
                             first = NA,
                             median = NA,
                             diff_hm = NA,
                             diff_max = NA,
                             n.surveys = NA)
sim.raw1 <- data.frame(model = NA,
                       sim = NA,
                       year = NA,
                       jday = NA,
                       mfit = NA,
                       mfit_se = NA)

# loop objects
nd.i <- data.frame(jday = 1:200)
nd.h <- data.frame(jday = rep(1:200, length(n.surveys)), 
                        year = rep(1:length(n.surveys), each = 200))
# start loop
for(j in 1:n.sim) {
  tic(msg = j)
  sim.pred.i <- data.frame(year = NA, jday = NA, mfit = NA, mfit_se = NA)
  sim.data.h <- data.frame(year = NA, jday = NA, detect = NA)
  for(i in 1:length(n.surveys)) {
    surveys.i <- sort((sample(1:199, size = rev(n.surveys)[i], replace = T, prob = pr.survey)))
    surveys.shift <- surveys.i + (i - 1) # this translates the pheno-curve over time
    surveys.shift[surveys.shift > 200] <- 200
    occ.i <- rbinom(n = length(surveys.i), size = 1, prob = probs_sim1$prob[surveys.shift])
    x.i <- data.frame(year = rep(i, length(surveys.i)), jday = surveys.i, detect = occ.i)
    sim.data.h <- rbind(sim.data.h, x.i)
    # run single-year GAM
    fit_sim.i <- mgcv::gam(detect ~ s(jday, k = 30),
                                data = x.i,
                                family = binomial(link = "logit"),
                                method = 'REML')
    pred_sim.i <- predict(fit_sim.i, type = 'response', 
                               newdata = nd.i,
                               se.fit=TRUE)
    y.i <- data.frame(year = i,
                           jday = nd.i$jday,
                           mfit = pred_sim.i$fit,
                           mfit_se = pred_sim.i$se.fit)
    # store predictions in object
    sim.pred.i <- rbind(sim.pred.i, y.i)
  }
  sim.pred.i <- sim.pred.i[-1, ]
  sim.raw1 <- rbind(sim.raw1, data.frame(model = "independent", sim = j, sim.pred.i))
  sim.data.h <- sim.data.h[-1, ]
  fit_sim.h <- mgcv::gam(detect ~ s(jday, k = 30) +
                                s(jday, year, k = 30, bs = 'fs', m = 1) +
                                s(year, bs = 're', k = 30),
                              data = sim.data.h,
                              family = binomial(link = "logit"),
                              method = 'REML')
  
  pred_sim.h <- predict(fit_sim.h, type = 'response', 
                             newdata = nd.h,
                             se.fit=TRUE)
  
  gi_plt.h <- data.frame(nd.h, mfit = pred_sim.h$fit, mfit_se = pred_sim.h$se.fit)
  sim.raw1 <- rbind(sim.raw1, data.frame(model = "hierarchical", sim = j, sim.pred.i))
  hm_gi.h <- hm_fun(fit_res = gi_plt.h, data = sim.data.h)
  hm_gi.h$diff_hm = hm_gi.h$hm_jday - stats_sim1$hm_jday + 0:12
  hm_gi.h$diff_max = hm_gi.h$max_jday - stats_sim1$max_jday + 0:12
  hm_gi.h$n.surveys <- rev(n.surveys)
  
  hm_gi.i <- hm_fun(fit_res = sim.pred.i, data = sim.data.h)
  hm_gi.i$diff_hm = hm_gi.i$hm_jday - stats_sim1$hm_jday + 0:12
  hm_gi.i$diff_max = hm_gi.i$max_jday - stats_sim1$max_jday + 0:12
  hm_gi.i$n.surveys <- rev(n.surveys)
  sim.metrics1 <- rbind(sim.metrics1, data.frame(model = "independent", sim = j, hm_gi.i), data.frame(model = "hierarchical", sim = j, hm_gi.h))
  
  summary.j <- data.frame(model = c("independent", "hierarchical"), 
                          sim = j,
                          mean_diff_hm = c(mean(hm_gi.i$diff_hm), mean(hm_gi.h$diff_hm)),
                          mean_diff_max = c(mean(hm_gi.i$diff_max), mean(hm_gi.h$diff_max)),
                          trend_hm = c(coef(lm(hm_jday ~ as.numeric(year), data = hm_gi.i))[2], coef(lm(hm_jday ~ as.numeric(year), data = hm_gi.h))[2]),
                          trend_max = c(coef(lm(max_jday ~ as.numeric(year), data = hm_gi.i))[2], coef(lm(max_jday ~ as.numeric(year), data = hm_gi.h))[2]))
  sim.summary1 <- rbind(sim.summary1, summary.j)
  toc(quiet = F)
}
sim.raw1 <- sim.raw1[-1, ]
sim.metrics1 <- sim.metrics1[-1, ]
sim.summary1 <- sim.summary1[-1, ]

# Conclusion: 
write.csv(sim.raw1, file = "sim_raw_lb.csv")
write.csv(sim.metrics1, file = "sim_metrics_lb.csv")
write.csv(sim.summary1, file = "sim_summary_lb.csv")

# Quick statistics
tapply(sim.summary1$mean_diff_hm, INDEX = sim.summary1$model, mean)
tapply(sim.summary1$mean_diff_max, INDEX = sim.summary1$model, mean)
tapply(sim.metrics1$hm_jday, INDEX = sim.metrics1$model, FUN = sd)
tapply(sim.metrics1$max_jday, INDEX = sim.metrics1$model, FUN = sd)
tapply(sim.summary1$trend_hm, INDEX = sim.summary1$model, mean)
tapply(sim.summary1$trend_max, INDEX = sim.summary1$model, mean)




# Simulation loop: repeat for a passage migrant cell ------------------------------------------------------------

#storage devices
sim.summary2 <- data.frame(model = NA,  
                           sim = NA,
                           mean_diff_hm = NA,
                           mean_diff_max = NA,
                           trend_hm = NA,
                           trend_max = NA)
sim.metrics2 <- data.frame(model = NA,
                           sim = NA,
                           year = NA,  
                           hm_jday = NA,
                           max_jday = NA,
                           max_val = NA,
                           first = NA,
                           median = NA,
                           diff_hm = NA,
                           diff_max = NA,
                           n.surveys = NA)
sim.raw2 <- data.frame(model = NA,
                       sim = NA,
                       year = NA,
                       jday = NA,
                       mfit = NA,
                       mfit_se = NA)

# start loop
for(j in 1:n.sim) {
  tic(msg = j)
  sim.pred.i <- data.frame(year = NA, jday = NA, mfit = NA, mfit_se = NA)
  sim.data.h <- data.frame(year = NA, jday = NA, detect = NA)
  for(i in 1:length(n.surveys)) {
    surveys.i <- sort((sample(1:199, size = rev(n.surveys)[i], replace = T, prob = pr.survey)))
    surveys.shift <- surveys.i + (i - 1) # this translates the pheno-curve over time
    surveys.shift[surveys.shift > 200] <- 200
    occ.i <- rbinom(n = length(surveys.i), size = 1, prob = probs_sim2$prob[surveys.shift])
    x.i <- data.frame(year = rep(i, length(surveys.i)), jday = surveys.i, detect = occ.i)
    sim.data.h <- rbind(sim.data.h, x.i)
    # run single-year GAM
    fit_sim.i <- mgcv::gam(detect ~ s(jday, k = 30),
                           data = x.i,
                           family = binomial(link = "logit"),
                           method = 'REML')
    pred_sim.i <- predict(fit_sim.i, type = 'response', 
                          newdata = nd.i,
                          se.fit=TRUE)
    y.i <- data.frame(year = i,
                      jday = nd.i$jday,
                      mfit = pred_sim.i$fit,
                      mfit_se = pred_sim.i$se.fit)
    # store predictions in object
    sim.pred.i <- rbind(sim.pred.i, y.i)
  }
  sim.pred.i <- sim.pred.i[-1, ]
  sim.raw2 <- rbind(sim.raw1, data.frame(model = "independent", sim = j, sim.pred.i))
  sim.data.h <- sim.data.h[-1, ]
  fit_sim.h <- mgcv::gam(detect ~ s(jday, k = 30) +
                           s(jday, year, k = 30, bs = 'fs', m = 1) +
                           s(year, bs = 're', k = 30),
                         data = sim.data.h,
                         family = binomial(link = "logit"),
                         method = 'REML')
  
  pred_sim.h <- predict(fit_sim.h, type = 'response', 
                        newdata = nd.h,
                        se.fit=TRUE)
  
  gi_plt.h <- data.frame(nd.h, mfit = pred_sim.h$fit, mfit_se = pred_sim.h$se.fit)
  sim.raw2 <- rbind(sim.raw2, data.frame(model = "hierarchical", sim = j, sim.pred.i))
  hm_gi.h <- hm_fun(fit_res = gi_plt.h, data = sim.data.h)
  hm_gi.h$diff_hm = hm_gi.h$hm_jday - stats_sim2$hm_jday + 0:12
  hm_gi.h$diff_max = hm_gi.h$max_jday - stats_sim2$max_jday + 0:12
  hm_gi.h$n.surveys <- rev(n.surveys)
  
  hm_gi.i <- hm_fun(fit_res = sim.pred.i, data = sim.data.h)
  hm_gi.i$diff_hm = hm_gi.i$hm_jday - stats_sim2$hm_jday + 0:12
  hm_gi.i$diff_max = hm_gi.i$max_jday - stats_sim2$max_jday + 0:12
  hm_gi.i$n.surveys <- rev(n.surveys)
  sim.metrics2 <- rbind(sim.metrics2, data.frame(model = "independent", sim = j, hm_gi.i), data.frame(model = "hierarchical", sim = j, hm_gi.h))
  
  summary.j <- data.frame(model = c("independent", "hierarchical"), 
                          sim = j,
                          mean_diff_hm = c(mean(hm_gi.i$diff_hm), mean(hm_gi.h$diff_hm)),
                          mean_diff_max = c(mean(hm_gi.i$diff_max), mean(hm_gi.h$diff_max)),
                          trend_hm = c(coef(lm(hm_jday ~ as.numeric(year), data = hm_gi.i))[2], coef(lm(hm_jday ~ as.numeric(year), data = hm_gi.h))[2]),
                          trend_max = c(coef(lm(max_jday ~ as.numeric(year), data = hm_gi.i))[2], coef(lm(max_jday ~ as.numeric(year), data = hm_gi.h))[2]))
  sim.summary2 <- rbind(sim.summary2, summary.j)
  toc(quiet = F)
}
sim.raw2 <- sim.raw2[-1, ]
sim.metrics2 <- sim.metrics2[-1, ]
sim.summary2 <- sim.summary2[-1, ]

# Write results files 
write.csv(sim.raw2, file = "sim_raw_pm.csv")
write.csv(sim.metrics2, file = "sim_metrics_pm.csv")
write.csv(sim.summary2, file = "sim_summary_pm.csv")

# Quick statistics
tapply(sim.summary2$mean_diff_hm, INDEX = sim.summary1$model, mean)
tapply(sim.summary2$mean_diff_max, INDEX = sim.summary1$model, mean)
tapply(sim.metrics2$hm_jday, INDEX = sim.metrics1$model, FUN = sd)
tapply(sim.metrics2$max_jday, INDEX = sim.metrics1$model, FUN = sd)
tapply(sim.summary2$trend_hm, INDEX = sim.summary1$model, mean)
tapply(sim.summary2$trend_max, INDEX = sim.summary1$model, mean)



# Results plot: to summarize 95%, 50%, and median quantiles ------------------------------------------------------------

pdf(file = "Simulation_results.pdf", width = 5, height = 10)
par(tcl = -0.2, mgp = c(2, 0.4, 0), mfrow = c(2, 1), mar = c(4, 4, 2, 1))
plot(1,1, xlim = c(-3, 3), ylim = c(0.5,8.5), yaxt = "n", ylab = "", xlab = "Bias (days)", type = "n")
title(main = "Mean bias in estimators", adj = 0, line = 0.2)
axis(side = 2, at = 1:8, labels = rep(c("LB", "PM"), 4), las = 1)
abline(v = 0, lty = 2, col = "gray")
mtext("Individual", side = 2, line = 3, adj = 0.8)
mtext("Hierarchcial", side = 2, line = 3, adj = 0.2)
mtext("half-max", side = 2, line = 2, adj = .92)
mtext("max", side = 2, line = 2, adj = 0.63)
mtext("half-max", side = 2, line = 2, adj = 0.35)
mtext("max", side = 2, line = 2, adj = 0.1)
arrows(x0 = unlist(tapply(
  sim.summary2$mean_diff_hm,
  INDEX = sim.summary2$model,
  quantile,
  probs = c(0.025, 0.25, 0.5))),
  x1 = unlist(tapply(
    sim.summary2$mean_diff_hm,
    INDEX = sim.summary2$model,
    quantile,
    probs = c(0.975, 0.75, 0.5))),
  y0 = c(4, 4, 4, 8, 8, 8),
  lwd = c(1, 4, 10),
  length = 0)
arrows(x0 = unlist(tapply(
  sim.summary1$mean_diff_hm,
  INDEX = sim.summary1$model,
  quantile,
  probs = c(0.025, 0.25, 0.5))),
  x1 = unlist(tapply(
    sim.summary1$mean_diff_hm,
    INDEX = sim.summary1$model,
    quantile,
    probs = c(0.975, 0.75, 0.5))),
  y0 = c(3, 3, 3, 7, 7, 7),
  lwd = c(1, 4, 10),
  length = 0)
arrows(x0 = unlist(tapply(
  sim.summary2$mean_diff_max,
  INDEX = sim.summary2$model,
  quantile,
  probs = c(0.025, 0.25, 0.5))),
  x1 = unlist(tapply(
    sim.summary2$mean_diff_max,
    INDEX = sim.summary2$model,
    quantile,
    probs = c(0.975, 0.75, 0.5))),
  y0 = c(2, 2, 2, 6, 6, 6),
  lwd = c(1, 4, 10),
  length = 0)
arrows(x0 = unlist(tapply(
  sim.summary1$mean_diff_max,
  INDEX = sim.summary1$model,
  quantile,
  probs = c(0.025, 0.25, 0.5))),
  x1 = unlist(tapply(
    sim.summary1$mean_diff_max,
    INDEX = sim.summary1$model,
    quantile,
    probs = c(0.975, 0.75, 0.5))),
  y0 = c(1, 1, 1, 5, 5, 5),
  lwd = c(1, 4, 10),
  length = 0)

plot(1,1, xlim = c(-2, 0), ylim = c(0.5,8.5), yaxt = "n", ylab = "", xlab = "Bias (days per year)", type = "n")
title(main = "Bias in temporal trend", adj = 0, line = 0.2)
axis(side = 2, at = 1:8, labels = rep(c("LB", "PM"), 4), las = 1)
abline(v = -1, lty = 2, col = "gray")
mtext("Individual", side = 2, line = 3, adj = 0.8)
mtext("Hierarchcial", side = 2, line = 3, adj = 0.2)
mtext("half-max", side = 2, line = 2, adj = .92)
mtext("max", side = 2, line = 2, adj = 0.63)
mtext("half-max", side = 2, line = 2, adj = 0.35)
mtext("max", side = 2, line = 2, adj = 0.1)
arrows(x0 = unlist(tapply(
  sim.summary2$trend_hm,
  INDEX = sim.summary2$model,
  quantile,
  probs = c(0.025, 0.25, 0.5))),
  x1 = unlist(tapply(
    sim.summary2$trend_hm,
    INDEX = sim.summary2$model,
    quantile,
    probs = c(0.975, 0.75, 0.5))),
  y0 = c(4, 4, 4, 8, 8, 8),
  lwd = c(1, 4, 10),
  length = 0)
arrows(x0 = unlist(tapply(
  sim.summary1$trend_hm,
  INDEX = sim.summary1$model,
  quantile,
  probs = c(0.025, 0.25, 0.5))),
  x1 = unlist(tapply(
    sim.summary1$trend_hm,
    INDEX = sim.summary1$model,
    quantile,
    probs = c(0.975, 0.75, 0.5))),
  y0 = c(3, 3, 3, 7, 7, 7),
  lwd = c(1, 4, 10),
  length = 0)
arrows(x0 = unlist(tapply(
  sim.summary2$trend_max,
  INDEX = sim.summary2$model,
  quantile,
  probs = c(0.025, 0.25, 0.5))),
  x1 = unlist(tapply(
    sim.summary2$trend_max,
    INDEX = sim.summary2$model,
    quantile,
    probs = c(0.975, 0.75, 0.5))),
  y0 = c(2, 2, 2, 6, 6, 6),
  lwd = c(1, 4, 10),
  length = 0)
arrows(x0 = unlist(tapply(
  sim.summary1$trend_max,
  INDEX = sim.summary1$model,
  quantile,
  probs = c(0.025, 0.25, 0.5))),
  x1 = unlist(tapply(
    sim.summary1$trend_max,
    INDEX = sim.summary1$model,
    quantile,
    probs = c(0.975, 0.75, 0.5))),
  y0 = c(1, 1, 1, 5, 5, 5),
  lwd = c(1, 4, 10),
  length = 0)
dev.off()