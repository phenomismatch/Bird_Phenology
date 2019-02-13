
library(rstan)

#sim data
alpha <- 
beta <- 
alpha1 <- 
beta1 <- 
sigma <- 
sigma1 <- 
y_true <- 

#feed to model  
y_obs <-
y_sd <- 
x <- 
x2 <-


#create data list for Stan
DATA <- list(y_obs = ,
             y_sd = ,
             x = ,
             x2 = 
             N = length(Y_OBS))

#model
test <- '
data {
  int<lower = 0> N;
  vector[N] y_obs;
  vector[N] y_sd;
  vector[N] x;
  vector[N] x2;
}

parameters {
  vector[N] y_true;  
  real alpha;
  real alpha1;
  real beta;
  real beta1;
  real sigma;
  real sigma1;
}

model {
  y_obs ~ normal(y_true, sigma_y);              \\ observation model
  y_true ~ normal(alpha1 + beta1 * x1, sigma1);
  y_true ~ normal(alpha2 + beta2 * x2, sigma2);
}
'

#stan forums relevant page
#https://discourse.mc-stan.org/t/cv-indexing-and-multiple-tilde-assignment/2123/4

#[y_{true}, \alpha_{1}, \beta_{1}, \sigma_{1}, \alpha_{2}, \beta_{2}, \sigma_{2} \mid y_{obs}, \sigma_{y}, x_{1}, x_{2}] \propto \prod_{i=1}^{I} N(y_{obs[i]} \mid y_{true[i]}, \sigma_{y}) \times N(y_{true[i]} \mid \alpha_{1}, \beta_{1}, \sigma_{1}, x_{1[i]}) \times N(y_{true[i]} \mid  \alpha_{2}, \beta_{2}, \sigma_{2}, x_{2[i]}) \times priors



# Run model ---------------------------------------------------------------

fit <- rstan::stan(model_code = test,
                   data = DATA,
                   chains = 3,
                   iter = 500,
                   cores = 1,
                   pars = 'y')

MCMCvis::MCMCsummary(fit, round = 2, n.eff = TRUE)

