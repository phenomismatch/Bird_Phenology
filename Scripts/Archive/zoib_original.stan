data {
  int n;
  vector[n] x; 
  vector<lower=0, upper=1>[n] y;
}

parameters {
  vector[2] coef_a;
  vector[2] coef_g;
  vector[2] coef_m;
  vector[2] coef_p;
}

transformed parameters {
  vector<lower=0, upper=1>[n] alpha;
  vector<lower=0, upper=1>[n] gamma;
  vector[n] mu;
  vector<lower=0>[n] phi;
  vector<lower=0>[n] p;
  vector<lower=0>[n] q;

  alpha = inv_logit(coef_a[1] + coef_a[2] * x);
  gamma = inv_logit(coef_g[1] + coef_g[2] * x);
  mu = inv_logit(coef_m[1] + coef_m[2] * x);
  phi = exp(coef_p[1] + coef_p[2] * x);
  p = mu .* phi;
  q = phi - mu .* phi;
}

model {
  coef_a ~ normal(0, 1);
  coef_g ~ normal(0, 1);
  coef_m ~ normal(0, 1);
  coef_p ~ normal(0, 1);
  
  // zero one inflated beta likelihood
  for (i in 1:n) {
    if (y[i] == 0) {
      target += log(alpha[i]) + log1m(gamma[i]);
    } else if (y[i] == 1) {
      target += log(alpha[i]) + log(gamma[i]);
    } else {
      target += log1m(alpha[i]) + beta_lpdf(y[i] | p[i], q[i]);
    }
  }
}
