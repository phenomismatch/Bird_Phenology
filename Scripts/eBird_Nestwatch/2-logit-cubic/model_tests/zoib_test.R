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
x <- rnorm(n) # predictor value

# coefficients for alpha
a_b0 <- 0
a_b1 <- 0.3

# calculating alpha, adding error
alpha <- inv_logit(a_b0 + a_b1 * x)

# coefficients for gamma
g_b0 <- 0
g_b1 <- 0.3

# calculating gamma, adding error
gamma <- inv_logit(g_b0 + g_b1 * x)

# coefficients for mu
m_b0 <- -0.6
m_b1 <- 0.5

# calculating mu, adding error 
mu <- inv_logit(m_b0 + m_b1 * x)

# coefficients for phi
p_b0 <- -0.5
p_b1 <- 0.7

# calculating phi, adding error
phi <- exp(p_b0 + p_b1 * x)

# calculating shape parameters for beta distribution
p <- mu * phi
q <- phi - mu * phi

# calculate, using alpha, if their score comes from bernoulli trial
y_is_discrete <- rbinom(n, 1, alpha)

# calculate y scores
y <- rep(NA, n) # initialize empty vector of outcomes
for (i in 1:n) {
  if (y_is_discrete[i] == 1) { # if y came from bernoulli...
    y[i] <- rbinom(1, 1, gamma[i]) # being 1 is determined by its gamma
  } else {
    y[i] <- rbeta(1, p[i], q[i]) # if it did not, draw from beta
  }
}

hist(y, breaks = 100)

library(rstan)

stan_d <- list(y = y, n = n, x = x)

m_init <- stan_model('zoib.stan')
m_fit <- rstan::sampling(m_init, data = stan_d, 
                  pars = c('coef_a', 'coef_g', 'coef_m', 'coef_p'), 
                  cores = 4)
m_fit

draws <- as.matrix(m_fit)
draws <- draws[, !(colnames(draws) == 'lp__')]

true <- c(a_b0, a_b1, 
          g_b0, g_b1, 
          m_b0, m_b1, 
          p_b0, p_b1)

mcmc_recover_intervals(draws, true)
