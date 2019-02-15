#look at the effect of inserting arbitrary values for sigma_y_obs
#that value makes no different as to the estimate of y_true, as expected

library(rstan)

#simulate data
sim_y <- rnorm(100, 50, 5)
#simulate obs error
obs_err <- rnorm(100, 0, 3)
#full y
full_y <- sim_y + obs_err

#remove some values
to.rm <- seq(1, 100, by = 5)
y_obs <- full_y[-to.rm]
#get indices
ind <- 1:100
full_sigma <- rep(3, 100)
full_sigma[to.rm] <- 0.1

#create data list for Stan
DATA <- list(N = 100,
             y_obs = y_obs,
             sigma_y = full_sigma,
             ii_obs = ind[-to.rm],
             lobs = length(ind[-to.rm]),
             ii_mis = to.rm,
             lmis = length(to.rm))


# Stan model --------------------------------------------------------------

test <- '
data {
int<lower = 0> N;
int<lower = 0> lobs;
int<lower = 0> lmis;
real y_obs[lobs];
real<lower = 0> sigma_y[N];
int<lower = 0> ii_obs[lobs];
int<lower = 0> ii_mis[lmis];
}

parameters {
real y_mis[lmis];
real mu;
real<lower = 0> sigma_yt;
real y_true[N];
}

transformed parameters {
real y[N];

y[ii_obs] = y_obs;
y[ii_mis] = y_mis;
}

model {
y ~ normal(y_true, sigma_y);
y_true ~ normal(mu, sigma_yt);
}'



# Run model ---------------------------------------------------------------

DELTA <- 0.95
TREE_DEPTH <- 15
STEP_SIZE <- 0.1
CHAINS <- 3
ITER <- 5000

fit <- rstan::stan(model_code = test,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('y', 'y_true', 'sigma_yt'),
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))

MCMCvis::MCMCsummary(fit, n.eff = TRUE)
MCMCvis::MCMCplot(fit, params = 'y_true', 
         rank = TRUE,
         sz_labels = 0.5)


