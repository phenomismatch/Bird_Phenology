##########################################
# look at estimates for cells with MAPS vs. eBird data
#
# MAPS estimates seem to be much later than eBird estimates
# seems to be due to lack of clear 'hump' in data
##########################################

#mrg_f6 from bj model (all rows)
#filter data
bj_data <- dplyr::filter(mrg_f6, !is.na(bj_IAR_mean))
br_data <- dplyr::filter(mrg_f6, !is.na(br_GAM_mean))

#fit models
bj_fit <- lm(bj_IAR_mean ~ gr_mn, data = bj_data)
bj_data$resid <- residuals(bj_fit)

br_fit <- lm(br_GAM_mean ~ gr_mn, data = br_data)
br_data$resid <- residuals(br_fit)


#plots
#BJ data
#GREEN = cells with MAPS data
bj_IAR <- dplyr::filter(bj_data, VALID_juv_GAM != TRUE)
plot(bj_IAR$gr_mn, bj_IAR$bj_IAR_mean, 
     xlab = 'Greenup', ylab = 'IAR-derived fledge date',
     main = 'IAR-derived data',
     pch = 19, col = rgb(0,0,0,0.2),
     xlim = c(90, 150), ylim = c(115, 205))
abline(bj_fit, col = 'red')
bj_juv <- dplyr::filter(bj_data, VALID_juv_GAM == TRUE)
points(bj_juv$gr_mn, bj_juv$bj_IAR_mean, col = rgb(0,1,0,0.4), pch = 19)

#BR data
plot(br_data$gr_mn, br_data$br_GAM_mean, 
     xlab = 'Greenup', ylab = 'GAM-derived fledge date (from eBird data)',
     main = 'GAM-derived data',
     pch = 19, col = rgb(0,0,0,0.2),
     xlim = c(90, 150), ylim = c(115, 205))
abline(br_fit, col = 'red')

dplyr::filter(mrg_f6, VALID_br_GAM == TRUE, VALID_juv_GAM == TRUE)
plot(bj_data$br_GAM_mean, bj_data$juv_GAM_mean, 
     pch = 19)




#read in juv GAM data (to get plmax)
tt <- readRDS('juv_master_2020-12-04.rds')
tt2 <- dplyr::select(tt, species, year, cell, plmax)
bj_juv2 <- dplyr::left_join(bj_juv, tt2)
#look at residuals in relation to plmax
plot(bj_juv2$plmax, bj_juv2$resid, 
     xlab = 'Proportion of GAM iterations with local maximum',
     ylab = 'Residual from breeding ~ greenup model above',
     pch = 19, col = rgb(0,0,0,0.5))

#only three MAPS with plmax > 0.99
dplyr::filter(bj_juv2, plmax > 0.99)
