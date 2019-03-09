data {
  int n;
  vector[n] day;
  vector[n] day2; 
  vector<lower=0, upper=1>[n] y;
}

parameters {
  real alpha_zo;
  real beta1_zo;
  real beta2_zo;
  real alpha_b;
  real beta1_b;
  real beta2_b;
}

transformed parameters {
  vector<lower=0, upper=1>[n] alpha;
  vector[n] mu;
  vector<lower=0>[n] phi;
  vector<lower=0>[n] p;
  vector<lower=0>[n] q;
  
  theta = inv_logit(alpha_zo + beta1_zo * day + beta2_zo * day2);
  mu = inv_logit(alpha_b + beta1_b * day + beta2_b * day2);
  p = mu .* phi;
  q = phi - mu .* phi;
}

model {
  // zero-inflated beta likelihood
  for (i in 1:n) 
  {
    if (y[i] == 0)
    {
      target += log(theta[i]);
    } else {
      target += log1m(theta[i]) + beta_lpdf(y[i] | p[i], q[i]);
    }
  }
}

