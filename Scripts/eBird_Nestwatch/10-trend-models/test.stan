data {
int<lower = 0> N;                                     // number of obs
int<lower = 0> US;                                    // number of species
real<lower = 0, upper = 200> x_obs[N];                // mean halfmax IAR
real<lower = 0> sigma_x[N];                           // sd halfmax IAR
int<lower = 1, upper = US> cn_id[N];                  // species ids
real<lower = 1, upper = 17> year[N];
}

parameters {
real<lower = 0, upper = 200> x_true[N];                           //true arrival
real mu_alpha;
real mu_beta;
vector[US] alpha;
vector[US] beta;
real<lower = 0> sigma_x_true;
vector<lower = 0>[2] sigma_ab;
corr_matrix[2] rho;
}

transformed parameters {
vector[2] v_mu_ab = [mu_alpha, mu_beta]';             // vector with mu_alpha and mu_beta
vector[2] v_ab[US];                           // vector with alpha and beta
cov_matrix[2] SRS_sigma_ab_rho;               // covariance matrix

for (j in 1:US)
{
  v_ab[j, 1:2] = [alpha[j], beta[j]]';
}

SRS_sigma_ab_rho = quad_form_diag(rho, sigma_ab);
}

model {

vector[N] mu;

target += normal_lpdf(mu_alpha | 0, 10);
target += normal_lpdf(mu_beta | 0, 3);
target += normal_lpdf(sigma_x_true | 0, 5);
target += normal_lpdf(sigma_ab | 0, 3);
target += lkj_corr_lpdf(rho | 2);

for (i in 1:N)
{
  mu[i] = alpha[cn_id[i]] + beta[cn_id[i]] * year[i];
}

// observation model - modeling true state as a function of some observed state
target += normal_lpdf(x_obs | x_true, sigma_x);
target += normal_lpdf(x_true | mu, sigma_x_true);

target += multi_normal_lpdf(v_ab | v_mu_ab, SRS_sigma_ab_rho);
}