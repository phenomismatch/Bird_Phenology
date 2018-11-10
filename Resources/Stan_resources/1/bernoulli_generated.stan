data {
  int<lower = 0> N;
  int<lower = 0, upper = 1> y[N];
}
parameters {
  real<lower = 0, upper = 1> theta;
}
model {
  y ~ bernoulli(theta);  
}
generated quantities {
  int n_rep;
  {
    int y_rep[N];
    for (n in 1:N)
      y_rep[n] <- bernoulli_rng(theta);
    n_rep <- sum(y_rep);
  }
}
