#####
#Look at alternatives for logit cubic - SCRATCH
#
#####



# zero inflated beta -----------------------------------------------------

zeros <- cyspdata$day[which(cyspdata$detect == 0)]
ones <- cyspdata$day[which(cyspdata$detect == 1)]

ds <- sort(unique(cyspdata$day))
prop_d <- rep(NA, length(ds))
for (p in 1:length(ds))
{
  #p <- 2
  td <- filter(cyspdata, day == ds[p])
  prop_d[p] <- mean(td$detect)
}

plot(ds, prop_d, type = 'l')

ma <- function(x, n = 5)
{
  stats::filter(x, rep(1/n, n), sides = 2)
}
ma_prop_d <- as.numeric(ma(prop_d))

plot(ds, ma_prop_d, type = 'l')

hist(prop_d)



din <- data.frame(ds = ds, ds2 = ds^2, ds3 = ds^3, pd = ma_prop_d)
din <- data.frame(ds = scale(ds, scale = FALSE), ds2 = scale(ds^2, scale = FALSE), 
                  ds3 = scale(ds^3, scale = FALSE), pd = ma_prop_d)
din2 <- din[-c(1,2,170,171),]

zib <- brm(bf(pd ~ ds, zi ~ ds),
           data = din2, family = zero_inflated_beta())
summary(zib)


newdata <- data.frame(ds = 1:200)

pv <- predict(zib, newdata)

plot(pv[,1], ylim = c(0, 1))
lines(din2$ds, din2$pd)

str(din2)
str(zinb)

range(din2$ma_prop_d)

library(brms)

zinb <- read.csv("https://paul-buerkner.github.io/data/fish.csv")
head(zinb)


fit_zinb1 <- brm(count ~ persons + child + camper, 
                 data = zinb, family = zero_inflated_poisson())

summary(fit_zinb1)
plot(marginal_effects(fit_zinb1), ask = FALSE)

fit_zinb2 <- brm(bf(count ~ persons + child + camper, zi ~ child), 
                 data = zinb, family = zero_inflated_poisson())

summary(fit_zinb2)
plot(marginal_effects(fit_zinb2), ask = FALSE)


zinb$count <- rbeta(NROW(zinb), shape1 = 1, shape2 = 1)
fit_zinb1 <- brm(count ~ persons + child + camper, 
                 data = zinb, family = zero_inflated_beta())




# ARMA model --------------------------------------------------------------

ds <- sort(unique(cyspdata$day))
prop_d <- rep(NA, length(ds))
for (p in 1:length(ds))
{
  #p <- 2
  td <- filter(cyspdata, day == ds[p])
  prop_d[p] <- mean(td$detect)
}

plot(ds, prop_d, type = 'l')

DATA <- data.frame(x = ds, y = prop_d)


library(brms)
hist(DATA$pd)

#ARMA
fit.ARMA <- brms::brm(y  ~ 1, 
                      autocor = cor_arma(~ day, p = 1, q = 1),
                      chains = 3,
                      cores = 3,
                      DATA)

summary(fit.ARMA)


#BSTS - use other package
fit <- brm(pd ~ day, 
           data = DATA, 
           autocor = cor_bsts(),
           chains = 3,
           cores = 3)

summary(fit)


#GP
fit.GP <- brms::brm(bf(y  ~ 1 + 
                         s(x, k = 100, bs = "gp")),
                    family = Beta(link = 'logit'),
                    chains = 3,
                    cores = 3,
                    DATA)

summary(fit.GP)
plot(fit.GP)


#PLOT results

#FIT <- fit.GP
FIT <- fit.ARMA
newdata <- data.frame(pd = rep(NA, NROW(DATA)), day = DATA$day)
pred <- as.data.frame(predict(FIT, newdata, nsamples = 500))

names(pred) <- c("est", "se", "lower", "upper")
pred$pd <- DATA$pd
pred$day <- DATA$day
ggplot(pred, aes(day, est, ymin = lower, ymax = upper)) + 
  geom_smooth(stat = 'identity') +
  geom_line(aes(day, pd), col = 'red', inherit.aes = FALSE)
#geom_line(aes(day, pd), col = 'red')


# data 
set.seed(250)
t_drift <- arima.sim(list(order = c(1,0,0), ar = 0.8), n = 50) + 0.50 * seq(1,50)
t_drift_df <- data.frame(y = as.matrix(t_drift))
t_drift_df$x <- seq(1, 50)


bayes_fit <- brm(
  y ~ 1,
  autocor = cor_ma(~ x, q = 5),
  family = Beta(link = 'logit'),
  data = DATA,
  chains = 3,
  cores = 3)

summary(bayes_fit)
newdata <- rbind(bayes_fit$data, data.frame(y = rep(NA, 6), x = 200:205))
pred <- as.data.frame(predict(bayes_fit, newdata, nsamples = 100))

names(pred) <- c("est", "se", "lower", "upper")
pred$y <- newdata$y
pred$x <- seq_len(nrow(newdata))
ggplot(pred, aes(x, est, ymin = lower, ymax = upper)) + 
  geom_smooth(stat = "identity") + 
  geom_line(aes(x, y), inherit.aes = FALSE, col = 'red')








# Process test data ----------------------------------------------------------------

#Agelaius phoeniceus - 2014 - 536
#Agelaius phoeniceus - 2014 - 538
#Agelaius phoeniceus - 2014 - 564
#Agelaius phoeniceus - 2014 - 566

# saveRDS(spdata2, 'ex_data.rds')
# spdata2 <- readRDS('ex_data.rds')

YEAR <- 2014
CELL <- 566

cyspdata <- dplyr::filter(spdata2, year == YEAR, cell == CELL)



# logistic ----------------------------------------------------------------

zeros <- cyspdata$day[which(cyspdata$detect == 0)]
ones <- cyspdata$day[which(cyspdata$detect == 1)]

ds <- sort(unique(cyspdata$day))
prop_d <- rep(NA, length(ds))
for (p in 1:length(ds))
{
  #p <- 2
  td <- dplyr::filter(cyspdata, day == ds[p])
  prop_d[p] <- mean(td$detect)
}

plot(ds, prop_d, type = 'l')


#moving average to smooth time series
ma <- function(x, n = 5)
{
  stats::filter(x, rep(1/n, n), sides = 2)
}
ma_prop_d <- as.numeric(ma(prop_d))
plot(ds, ma_prop_d, type = 'l')

DATA <- data.frame(prop_d, ds)

#fit model
#Asym/(1+exp((xmid-input)/scal))

fit <- nls(prop_d ~ SSlogis(ds, Asym, xmid, scal),
           control = nls.control(maxiter = 10000, minFactor = 1/5000000),
           data = DATA)


#extract params
CI_Asym <- nlstools::confint2(fit, 'Asym')
CI_xmid <- nlstools::confint2(fit, 'xmid')
CI_scal <- nlstools::confint2(fit, 'scal')
Asym <- coef(fit)[1]
xmid <- coef(fit)[2]
scal <- coef(fit)[3]
x <- ds
y <- Asym/(1+exp((xmid-x)/scal))
pdt <- data.frame(x,y)




# logit -------------------------------------------------------------------

ITER <- 1000
CHAINS <- 3
DELTA <- 0.95
TREE_DEPTH <- 17

tfit <- rstanarm::stan_glm(detect ~ sjday + sjday2 + sjday3 + sjday4 + sjday5 + sjday6 + shr, 
                           data = cyspdata,
                           family = binomial(link = "logit"),
                           algorithm = 'sampling',
                           iter = ITER,
                           chains = CHAINS,
                           cores = CHAINS,
                           adapt_delta = DELTA,
                           control = list(max_treedepth = TREE_DEPTH))

# plots logit -------------------------------------------------------------

predictDays <- range(cyspdata$sjday)[1]:range(cyspdata$sjday)[2]
predictDays2 <- predictDays^2
predictDays3 <- predictDays^3
predictDays4 <- predictDays^4
predictDays5 <- predictDays^5
predictDays6 <- predictDays^6
predictDays7 <- predictDays^7

newdata <- data.frame(sjday = predictDays, sjday2 = predictDays2,
                      sjday3 = predictDays3, sjday4 = predictDays4,
                      sjday5 = predictDays5, sjday6 = predictDays6,
                      sjday7 = predictDays7, shr = 0)

dfit <- rstanarm::posterior_linpred(tfit, newdata = newdata, transform = T)
halfmax_fit <- rep(NA, ((ITER/2)*CHAINS))

for (L in 1:((ITER/2)*CHAINS))
{
  #L <- 1
  rowL <- as.vector(dfit[L,])
  halfmax_fit[L] <- predictDays[min(which(rowL > (max(rowL)/2)))]
}

cyspdata2 <- cyspdata
mn_dfit <- apply(dfit, 2, mean)
LCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.025))
UCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.975))
mn_hm <- mean(halfmax_fit)
LCI_hm <- quantile(halfmax_fit, probs = 0.025)
UCI_hm <- quantile(halfmax_fit, probs = 0.975)

plot(predictDays, UCI_dfit, type = 'l', col = 'red', lty = 2, lwd = 3,
     ylim = c(-(max(UCI_dfit)/5), max(UCI_dfit)),
     main = paste0(args, ' - ', YEAR, ' - ', CELL),
     xlab = 'Julian Day', ylab = 'Detection Probability')
lines(predictDays, LCI_dfit, col = 'red', lty = 2, lwd = 3)
lines(predictDays, mn_dfit, lwd = 3)
cyspdata2$detect[which(cyspdata2$detect == 1)] <- max(UCI_dfit)
points(cyspdata2$day, cyspdata2$detect, col = rgb(0,0,0,0.25))

ds <- sort(unique(cyspdata$day))
prop_d <- rep(NA, length(ds))
for (p in 1:length(ds))
{
  #p <- 2
  td <- dplyr::filter(cyspdata, day == ds[p])
  prop_d[p] <- mean(td$detect)
}

ma <- function(x, n = 5)
{
  stats::filter(x, rep(1/n, n), sides = 2)
}
ma_prop_d <- ma(prop_d) * max(UCI_dfit)

lines(ds, ma_prop_d, col = rgb(0.5,0.5,0.5,0.9), lwd = 2)

segments(x0 = mn_hm, x1 = mn_hm, y0 = 0, y1 = 1, col = rgb(0,0,1,0.5), lwd = 3)
segments(x0 = LCI_hm, x1 = LCI_hm, y0 = 0, y1 = 1, col = rgb(0,0,1,0.5), lwd = 3, lty = 2)
segments(x0 = UCI_hm, x1 = UCI_hm, y0 = 0, y1 = 1, col = rgb(0,0,1,0.5), lwd = 3, lty = 2)



# GAMM --------------------------------------------------------------------

tfit_gamm <- rstanarm::stan_gamm4(detect ~ s(sjday) + shr, 
                                  data = cyspdata,
                                  family = binomial(link = "logit"),
                                  algorithm = 'sampling',
                                  iter = ITER,
                                  chains = CHAINS,
                                  cores = CHAINS,
                                  adapt_delta = DELTA,
                                  control = list(max_treedepth = TREE_DEPTH))



# plots GAM --------------------------------------------------------------

predictDays <- range(cyspdata$sjday)[1]:range(cyspdata$sjday)[2]
predictDays2 <- predictDays^2
predictDays3 <- predictDays^3
predictDays4 <- predictDays^4
predictDays5 <- predictDays^5
predictDays6 <- predictDays^6
predictDays7 <- predictDays^7

newdata <- data.frame(sjday = predictDays, sjday2 = predictDays2,
                      sjday3 = predictDays3, sjday4 = predictDays4,
                      sjday5 = predictDays5, sjday6 = predictDays6,
                      sjday7 = predictDays7, shr = 0)

dfit <- rstanarm::posterior_linpred(tfit_gamm, newdata = newdata, transform = T)
halfmax_fit <- rep(NA, ((ITER/2)*CHAINS))

for (L in 1:((ITER/2)*CHAINS))
{
  #L <- 1
  rowL <- as.vector(dfit[L,])
  halfmax_fit[L] <- predictDays[min(which(rowL > (max(rowL)/2)))]
}

cyspdata2 <- cyspdata
mn_dfit <- apply(dfit, 2, mean)
LCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.025))
UCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.975))
mn_hm <- mean(halfmax_fit)
LCI_hm <- quantile(halfmax_fit, probs = 0.025)
UCI_hm <- quantile(halfmax_fit, probs = 0.975)


plot(predictDays, UCI_dfit, type = 'l', col = 'red', lty = 2, lwd = 3,
     ylim = c(-(max(UCI_dfit)/5), max(UCI_dfit)),
     main = paste0(args, ' - ', YEAR, ' - ', CELL, ' - GAM'),
     xlab = 'Julian Day', ylab = 'Detection Probability')
lines(predictDays, LCI_dfit, col = 'red', lty = 2, lwd = 3)
lines(predictDays, mn_dfit, lwd = 3)
cyspdata2$detect[which(cyspdata2$detect == 1)] <- max(UCI_dfit)
points(cyspdata2$day, cyspdata2$detect, col = rgb(0,0,0,0.25))

ds <- sort(unique(cyspdata$day))
prop_d <- rep(NA, length(ds))
for (p in 1:length(ds))
{
  #p <- 2
  td <- dplyr::filter(cyspdata, day == ds[p])
  prop_d[p] <- mean(td$detect)
}

ma <- function(x, n = 5)
{
  stats::filter(x, rep(1/n, n), sides = 2)
}
ma_prop_d <- ma(prop_d) * max(UCI_dfit)

lines(ds, ma_prop_d, col = rgb(0.5,0.5,0.5,0.9), lwd = 2)

segments(x0 = mn_hm, x1 = mn_hm, y0 = 0, y1 = 1, col = rgb(0,0,1,0.5), lwd = 3)
segments(x0 = LCI_hm, x1 = LCI_hm, y0 = 0, y1 = 1, col = rgb(0,0,1,0.5), lwd = 3, lty = 2)
segments(x0 = UCI_hm, x1 = UCI_hm, y0 = 0, y1 = 1, col = rgb(0,0,1,0.5), lwd = 3, lty = 2)




# plot logistic function --------------------------------------------------

#95% CI around xmid from logistic function
rect(CI_xmid[1], 0, CI_xmid[2], max(UCI_dfit), col = rgb(1,0,0,0.2), border = rgb(1,0,0,0.2))

#plot model fit
lines(pdt$x, (pdt$y * max(UCI_dfit)), col = rgb(1,0,0,0.6), lwd = 3)