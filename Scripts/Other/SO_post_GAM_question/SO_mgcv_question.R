library(mgcv)

#Setophaga_tigrina
#2002
#619

#https://fromthebottomoftheheap.net/2014/06/16/simultaneous-confidence-intervals-for-derivatives/
#https://fromthebottomoftheheap.net/2011/06/12/additive-modelling-and-the-hadcrut3v-global-mean-temperature-series/

#response
y <- c(rep(0, 89), 1, 0, 1, 0, 0, 1, rep(0, 13), 1, 0, 0, 1, 
       rep(0, 10), 1, 0, 0, 1, 1, 0, 1, rep(0,4), 1, rep(0,3),  
       1, rep(0, 3), 1, rep(0, 10), 1, rep(0, 4), 1, 0, 1, 0, 0, 
       rep(1, 4), 0, rep(1, 5), rep(0, 4), 1, 1, rep(0, 46))

#predictor
x <- c(1, 1, 1, 2, 4, 5, 5, 6, 6, 6, 10, 10, 12, 12, 14, 16, 
       rep(19, 6), 24, 26, 27, 28, 28, rep(33, 3), 34, 35, 36, 
       37, 40, 40, 45, 47, 48, 50, 55, 55, 61, 62, 68, rep(75, 5), 
       76, rep(82, 3), 83, 87, 89, 90, 90, 94, 95, rep(96, 4), 
       rep(97, 3), 99, 102, 103, 104, 106, rep(110, 3), 111, 113, 
       113, 115, 117, 117, 118, 119, 120, 122, 122, 123, 123, 
       rep(124, 10), rep(125, 11), rep(126, 6), rep(127, 4), 
       rep(128, 5), rep(129, 3), rep(130, 5), rep(131, 10), 
       rep(132, 5), rep(133, 3), 134, 134, rep(135, 6), rep(136, 3),
       137, 137, rep(138, 9), rep(139, 4), 140, rep(141, 3), 142, 
       142, 143, 143, 144, rep(145, 3), rep(146, 3), 147, 147, 149,
       150, 152, 152, rep(157, 3), 158, 158, 159, 159, 163, 164, 
       167, 173, 173, 175, 179, 180, 180, 183, 185, 186, 186, 187, 
       188, 191, 194, 197, 197)


#fit model
fit <- mgcv::gam(y ~ s(x, bs = 'ad', k = 50),
                  method = 'REML',
                  select = TRUE,
                  family = binomial(link = "logit"))

#new data to predict at
x_pred <- data.frame(x = 0:200)

#predict
pp <- predict(fit, x_pred, type = 'response', se.fit = TRUE)

#extract Xp matrix
Xp <- predict(fit, x_pred, type = "lpmatrix")

#posterior realizations
n <- 100
br <- MASS::mvrnorm(n, coef(fit), fit$Vc)

#plot data
plot(x, y, col = rgb(0,0,0,0.25), ylim = c(0,1))
#plot model fit and 95% CI
lines(x_pred$x, pp$fit, col = 'green', lwd = 2)
lines(x_pred[,1], pp$fit - 1.97 * pp$se.fit, col = 'green', lwd = 2, lty = 2)
lines(x_pred[,1], pp$fit + 1.97 * pp$se.fit, col = 'green', lwd = 2, lty = 2)

#plot posterior realizations of model fit
for (i in 1:n)
{ 
#i <- 1
  fits <- Xp %*% br[i,]
  lines(x_pred$x, boot::inv.logit(fits), col = rgb(1,0,0,0.1))
}

#extract mean and 95% CI from posterior realizations
pm <- Xp %*% t(br)
xp_mn <- apply(boot::inv.logit(pm), 1, mean)
xp_LCI <- apply(boot::inv.logit(pm), 1, 
                function(x) quantile(x, probs = 0.025))
xp_UCI <- apply(boot::inv.logit(pm), 1, 
                function(x) quantile(x, probs = 0.975))

lines(x_pred$x, xp_mn, col = 'blue', lwd = 2)
lines(x_pred$x, xp_LCI, col = 'blue', lwd = 2, lty = 2)
lines(x_pred$x, xp_UCI, col = 'blue', lwd = 2, lty = 2)

