#############################
#MAPS exploration
#
#
#############################



# maps_data <- readRDS('MAPS_age_filled.rds')


#What characterizes breeding? Brood patch?



#--------------------------------#
#distirbution of hatch year birds caught over the season

maps_f <- dplyr::filter(maps_data, true_age == 0)
plyr::count(maps_f, 'Sci_name')

oc_0 <- filter(maps_f, sci_name == 'Oreothlypis celata')
plyr::count(oc_0, c('Cell', 'Year'))
tt <- dplyr::filter(oc_0, Cell == 444, Year == 2011)

plot(density(tt$Jday))

#--------------------------------#


#--------------------------------#
#logistic model for presence of brood patch

maps_adults <- dplyr::filter(maps_data, age %in% c('1', '5', '6', '7', '8'))

plyr::count(maps_adults, 'sci_name')

C_ustulatus <- dplyr::filter(maps_adults, sci_name == 'Catharus ustulatus')

#breeding defined as:
#BP 2-4
#CP 1-3

cu_f <- dplyr::filter(C_ustulatus, sex == 'F')
cu_m <- dplyr::filter(C_ustulatus, sex == 'M')

#add breeder col
C_ustulatus$br <- NA
cu_f$br <- NA

#male female
plot(density(C_ustulatus[which(C_ustulatus$Sex == 'M'),]$day))
lines(density(C_ustulatus[which(C_ustulatus$Sex == 'F'),]$day), col = 'red')

#cloacal protuberance probably not a good metric
plot(density(cu_m[which(cu_m$cloacal_pro <= 2),]$day))
lines(density(cu_m[which(cu_m$cloacal_pro > 2),]$day), col = 'red')

#brood patch maybe better
plot(density(cu_f[which(cu_f$brood_patch <= 2 | cu_f$brood_patch == 5),]$day))
lines(density(cu_f[which(cu_f$brood_patch > 2 & cu_f$brood_patch < 5),]$day), 
      col = 'red')

cu_f[which(cu_f$brood_patch > 2 & cu_f$brood_patch < 5),]$br <- 1
cu_f[which(cu_f$brood_patch <= 2 | cu_f$brood_patch == 5),]$br <- 0

#DATA <- dplyr::filter(cu_f, Jday < DAY, Jday > 50)
DATA <- cu_f

DATA$jday <- DATA$day
DATA$jday2 <- DATA$day^2
DATA$jday3 <- DATA$day^3

plyr::count(cu_f, c('Cell', 'Year'))


library(rstanarm)
ITER <- 1000
#ITER <- 10
CHAINS <- 3

#defaults for rstanarm are 0.95 and 15
DELTA <- 0.95
TREE_DEPTH <- 15

fit2 <- rstanarm::stan_gamm4(br ~ s(jday), 
                             data = cyspdata,
                             family = binomial(link = "logit"),
                             algorithm = 'sampling',
                             iter = ITER,
                             chains = CHAINS,
                             cores = CHAINS,
                             adapt_delta = DELTA,
                             control = list(max_treedepth = TREE_DEPTH))


#calculate diagnostics
num_diverge <- rstan::get_num_divergent(fit2$stanfit)
num_tree <- rstan::get_num_max_treedepth(fit2$stanfit)
num_BFMI <- length(rstan::get_low_bfmi_chains(fit2$stanfit))

#rerun model if things didn't go well
while (sum(c(num_diverge, num_tree, num_BFMI)) > 0 & DELTA <= 0.98)
{
  DELTA <- DELTA + 0.01
  TREE_DEPTH <- TREE_DEPTH + 1
  
  fit2 <- rstanarm::stan_gamm4(detect ~ s(jday) + shr,
                               data = cyspdata,
                               family = binomial(link = "logit"),
                               algorithm = 'sampling',
                               iter = ITER,
                               chains = CHAINS,
                               cores = CHAINS,
                               adapt_delta = DELTA,
                               control = list(max_treedepth = TREE_DEPTH))
  
  num_diverge <- rstan::get_num_divergent(fit2$stanfit)
  num_tree <- rstan::get_num_max_treedepth(fit2$stanfit)
  num_BFMI <- length(rstan::get_low_bfmi_chains(fit2$stanfit))
}

halfmax_df$num_diverge[counter] <- num_diverge
halfmax_df$num_tree[counter] <- num_tree
halfmax_df$num_BFMI[counter] <- num_BFMI
halfmax_df$delta[counter] <- DELTA
halfmax_df$tree_depth[counter] <- TREE_DEPTH

#generate predict data
predictDays <- range(cyspdata$jday)[1]:range(cyspdata$jday)[2]
newdata <- data.frame(jday = predictDays, shr = 0)

#predict response
dfit <- rstanarm::posterior_linpred(fit2, newdata = newdata, transform = T)
halfmax_fit <- rep(NA, ((ITER/2)*CHAINS))

for (L in 1:((ITER/2)*CHAINS))
{
  #L <- 1
  rowL <- as.vector(dfit[L,])
  halfmax_fit[L] <- predictDays[min(which(rowL > (max(rowL)/2)))]
}

########################
#PLOT MODEL FIT AND DATA

#summary(fit2)
mn_dfit <- apply(dfit, 2, mean)
LCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.025))
UCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.975))
mn_hm <- mean(halfmax_fit)
LCI_hm <- quantile(halfmax_fit, probs = 0.025)
UCI_hm <- quantile(halfmax_fit, probs = 0.975)

pdf(paste0(args, '_', years[j], '_', cells[k], '_arrival.pdf'))
plot(predictDays, UCI_dfit, type = 'l', col = 'red', lty = 2, lwd = 2,
     ylim = c(0, max(UCI_dfit)),
     main = paste0(args, ' - ', years[j], ' - ', cells[k]),
     xlab = 'Julian Day', ylab = 'Detection Probability')
lines(predictDays, LCI_dfit, col = 'red', lty = 2, lwd = 2)
lines(predictDays, mn_dfit, lwd = 2)
cyspdata$detect[which(cyspdata$detect == 1)] <- max(UCI_dfit)
points(cyspdata$day, cyspdata$detect, col = rgb(0,0,0,0.25))
abline(v = mn_hm, col = rgb(0,0,1,0.5), lwd = 2)
abline(v = LCI_hm, col = rgb(0,0,1,0.5), lwd = 2, lty = 2)
abline(v = UCI_hm, col = rgb(0,0,1,0.5), lwd = 2, lty = 2)
legend('topleft',
       legend = c('Cubic fit', 'CI fit', 'Half max', 'CI HM'),
       col = c('black', 'red', rgb(0,0,1,0.5), rgb(0,0,1,0.5)),
       lty = c(1,2,1,2), lwd = c(2,2,2,2), cex = 1.3)
dev.off()






# phenology as it related to age ------------------------------------------

# jfit <- lm(MAPS_age$Jday ~ MAPS_age$True_age)
# summary(jfit)
# 
# plot(MAPS_age$True_age, MAPS_age$Jday, pch = '.', col = rgb(0,0,0, 0.5))
# abline(jfit, col = 'red')
# a1 <- dplyr::filter(MAPS_age, True_age == 1)
# a2 <- dplyr::filter(MAPS_age, True_age == 2)
# a3 <- dplyr::filter(MAPS_age, True_age == 3)
# a4 <- dplyr::filter(MAPS_age, True_age == 4)
# a5 <- dplyr::filter(MAPS_age, True_age == 5)
# plot(density(a1$Jday))
# lines(density(a2$Jday))
# lines(density(a3$Jday))
# lines(density(a4$Jday))
# lines(density(a5$Jday))



#*breeding timing
#Each species/cell/year:
#maybe each capture different row, 1 for brood patch/CP, 0 if not -> fit a logistic to this for each cell/year?
#y[i] ~ bern(p[i])
#logit(p[i]) = alpha + beta * DAY[i]

#*phenology and age
#could fit logistic regression with age as random effect -> each age class would have different beta param
#for each cell/year: predict start of breeding for each age class (age 1, age 2, age 3+)
#calculate derived qty - difference in halfmax estimates
#y[i] ~ bern(p[i])
#logit(p[i]) = alpha[id[i]] + beta[id[i]] * DAY[i]

#*difference in arrival between males and females - how this maps onto phlogeny/traits




# other metric changes with age -------------------------------------------


#maybe standardized by wing chord

#*fat content and age - Fat[i] ~ alpha[id[i]] + beta1[id[i]] * Age[i] + beta2[id[i]] * Day[i] + beta3[id[i]] * Day[i]^2
plot(MAPS_age$True_age, MAPS_age$Fat_content, col = rgb(0,0,0, 0.5))
plot(MAPS_age$Jday, MAPS_age$Fat_content, col = rgb(0,0,0, 0.5))
summary(lm(Fat_content ~ True_age + poly(Jday, 2, raw = TRUE), data = MAPS_age))

#*Weight and age - Weight[i] ~ alpha[id[i]] + beta1[id[i]] * Age[i] + beta2[id[i]] * Day[i] + beta3[id[i]] * Day[i]^2
plot(MAPS_age$True_age, MAPS_age$Weight, col = rgb(0,0,0, 0.5))
plot(MAPS_age$Jday, MAPS_age$Weight, col = rgb(0,0,0, 0.5))
summary(lm(Weight ~ True_age + poly(Jday, 2, raw = TRUE), data = MAPS_age))




# has age distribution of birds changed over time? ------------------------

#more older birds over time? more younger birds over time?
