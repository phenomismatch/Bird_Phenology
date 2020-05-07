
#new data to predict at
x_pred <- data.frame(jday = 0:200,
                     shr = 0)

#predict
pp <- predict(fit2, x_pred, type = 'response', se.fit = TRUE)

#extract Xp matrix
Xp <- predict(fit2, x_pred, type = "lpmatrix")

#posterior realizations
N <- 1e3
br <- rmvn(N, coef(fit2), fit2$Vc)
br <- MASS::mvrnorm(N, coef(fit2), fit2$Vc)

#extract mean and 95% CI from posterior realizations
pm <- Xp %*% t(br)

#plot data
plot(cyspdata$jday, cyspdata$detect, col = rgb(0,0,0,0.25), ylim = c(0,1))
#plot model fit and 95% CI
lines(x_pred$jday, pp$fit, col = 'green', lwd = 2)
lines(x_pred$jday, pp$fit - 1.97 * pp$se.fit, col = 'green', lwd = 2, lty = 2)
lines(x_pred$jday, pp$fit + 1.97 * pp$se.fit, col = 'green', lwd = 2, lty = 2)

dfit <- t(pm)
predictDays <- x_pred$jday
halfmax_fit <- rep(NA, N)
tlmax <- rep(NA, N)
for (L in 1:NROW(dfit))
{
  #L <- 842
  rowL <- as.vector(boot::inv.logit(dfit[L,]))
  lines(predictDays, rowL, col = rgb(0,0,0,0.1))
  #first detection
  fd <- min(cyspdata$jday[which(cyspdata$detect == 1)])
  #abline(v = fd, col = 'red')
  #local maximum(s)
  #from: stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
  lmax_idx <- which(diff(sign(diff(rowL))) == -2) + 1
  lmax <- predictDays[lmax_idx]
  #first local max to come after first detection
  flm <- which(lmax > fd)
  if (length(flm) > 0)
  {
    #first local max to come after first detection
    lmax2_idx <- lmax_idx[min(flm)]
    lmax2 <- lmax[min(flm)]
    tlmax[L] <- TRUE
  } else {
    #no local max
    lmax2_idx <- which.max(rowL)
    lmax2 <- predictDays[which.max(rowL)]
    tlmax[L] <- FALSE
  }
  #VV
  #local mins before max (global min + local min; could be the same)
  lmin_idx <- c(which.min(rowL[1:lmax2_idx]), 
                which(diff(sign(diff(rowL[1:lmax2_idx]))) == 2) + 1)
  lmin <- predictDays[lmin_idx]
  #local min nearest to local max
  lmin2_idx <- lmin_idx[which.min(lmax2 - lmin)]
  lmin2 <- predictDays[lmin2_idx]
  #^^
  
  #position of min value before max - typically, where 0 is
  #lmin_idx <- which.min(rowL[1:lmax2_idx])
  #value at local max - value at min (typically 0)
  dmm <- rowL[lmax2_idx] - rowL[lmin2_idx]
  #all positions less than or equal to half diff between max and min + value min
  tlm <- which(rowL <= ((dmm/2) + rowL[lmin2_idx]))
  #which of these come before max and after or at min
  vgm <- tlm[which(tlm < lmax2_idx & tlm >= lmin2_idx)]
  #insert halfmax (first day for situations where max is a jday = 1)
  if (length(vgm) > 0)
  {
    halfmax_fit[L] <- predictDays[max(vgm)]
  } else {
    halfmax_fit[L] <- predictDays[1]
  }
  abline(v = halfmax_fit[L], col = rgb(1,0,0,0.1))
}

abline(v = mean(halfmax_fit), col = 'lightblue', lwd = 3)
abline(v = quantile(halfmax_fit, 0.025), col = 'lightblue', lty = 2)
abline(v = quantile(halfmax_fit, 0.975), col = 'lightblue', lty = 2)
hist(halfmax_fit)
which(halfmax_fit > 150)

mean(halfmax_fit)
sd(halfmax_fit)
quantile(halfmax_fit, 0.025)
quantile(halfmax_fit, 0.975)
