## DATA GENERATION
# reproducibility
library(bayesplot)
set.seed(1839)

# matches the stan function
inv_logit <- function(x)
{
  exp(x) / (1 + exp(x))
}

n <- 10000 # sample size
day <- runif(n, 0, 200) # predictor value
day2 <- day^2
day3 <- day^3

# theta
a_zo <- -10
b1_zo <- 0.05
b2_zo <- 0.001
b3_zo <- -0.0000048

theta <- inv_logit(a_zo + b1_zo * day + b2_zo * day2 + b3_zo * day3)
#hist(theta)
plot(day, theta)


# mu
a_b <- -1
b1_b <- 0.05
b2_b <- 0.0001
b3_b <- -0.0000015

mu <- inv_logit(a_b + b1_b * day + b2_b * day2 + b3_b * day3)
# hist(mu)
plot(day, mu)

# calculating phi
phi <- 1

# calculating shape parameters for beta distribution
p <- mu * phi
q <- phi - mu * phi


# simulate y
y <- rep(NA, n) # initialize empty vector of outcomes
for (i in 1:n) 
{
  #bernouli
  is.discrete <- rbinom(1, 1, theta[i])
  
  #draw from beta if 1, else 0
  if (is.discrete == 1)
  {
    y[i] <- rbeta(1, p[i], q[i]) 
  } else {
    y[i] <- 0
  }
}


# hist(y, breaks = 100)


tt <- data.frame(day, y)
stt <- tt[order(tt$day),]
plot(stt$day, stt$y, type = 'l')

ma <- function(x, n = 10)
{
  stats::filter(x, rep(1/n, n), sides = 2)
}
y_ma <- as.numeric(ma(stt$y))

plot(stt$day, y_ma, type = 'l')







library(rstan)

DATA <- list(y = y, n = length(y), day = day, day2 = day2, day3 = day3)

model <- stan_model('zoib.stan')
m_fit <- rstan::sampling(model, 
                         data = DATA, 
                         pars = c('alpha_zo', 'beta1_zo', 'beta2_zo', 'beta3_zo', 
                                  'alpha_b', 'beta1_b', 'beta2_b', 'beta3_b', 'mu', 'phi'), 
                         cores = 3)
m_fit

draws <- as.matrix(m_fit)
draws <- draws[, !(colnames(draws) == 'lp__')]

true <- c(a_zo, b1_zo, b2_zo, b3_zo,
          a_b, b1_b, b2_b, b3_b)

mcmc_recover_intervals(draws, true)
