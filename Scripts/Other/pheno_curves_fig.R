x <- seq(-8, 8, length=100)
hx <- dnorm(x, 0 , 1.2)
h2 <- dnorm(x, 0 , 1)
h3 <- dnorm(x, -4 , 1)
to.rm <- which(x < -3.5 | x > 3.5)
to.rm2 <- which(x < -7.5 | x > -0.5)

#color from upload to: https://imagecolorpicker.com/
plot(x[-to.rm], hx[-to.rm], type="l", lty = 1, 
     lwd = 6, xlab = '', col = '#B4C7E7', ylab="", main="", 
     ylim = c(0, 0.5), xlim = c(-8, 4))
lines(x[-to.rm], h2[-to.rm] * 0.75, lty = 1, 
      lwd = 6, col = '#F8CBAE')

plot(x[-to.rm], hx[-to.rm], type="l", lty = 1, 
     lwd = 6, xlab = '', xlim = c(-8, 4), 
     col = '#B4C7E7', ylab="", 
     main="", ylim = c(0, 0.5))
lines(x[-to.rm2], h3[-to.rm2] * 0.75, lty = 1, lwd = 6, col = '#F8CBAE')




legend("topright", inset=.05, title="Distributions",
       labels, lwd=2, lty=c(1, 1, 1, 1, 2), col=colors)
