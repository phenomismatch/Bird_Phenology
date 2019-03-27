#FROM: https://groups.google.com/forum/#!topic/stan-users/HaPPJ5dFUKw
#WRITTEN BY: Richard McElreath
#See also: Page 39 Stan Users Guide
#See also: https://github.com/ssp3nc3r/rethinking/blob/master/chapter13.Rmd



# example Cholesky mixed-effects model (non-centered parameterization)

# sim data

library(MASS)
N <- 1000
N_groups <- 100
group <- rep(1:N_groups,each=N/N_groups)
B <- mvrnorm(N_groups,c(-1,0.5),matrix(c(1,-0.5,-0.5,1),2,2))
x <- rnorm(N)
y <- rnorm(N,B[group,1]+B[group,2]*x,2)

# implied model: y ~ (1+x|group) + x

# first ordinary centered case, for comparison
# this is only partly centered, because means are in linear model, not in prior

mc1 <- '
data{
    int N;
    int N_groups;
    real x[N];
    real y[N];
    int group[N];
}
transformed data{
    vector[2] zeros;
    zeros[1] <- 0;
    zeros[2] <- 0;
}
parameters{
    real a;
    real b;
    real<lower=0> sigma;
    vector<lower=0>[2] sigma_group;
    corr_matrix[2] Rho;
    vector[2] v[N_groups];
}
model{
    matrix[2,2] SIGMA;
    vector[N] mu;
    SIGMA <- quad_form_diag(Rho,sigma_group);
    v ~ multi_normal( zeros , SIGMA );
    a ~ normal(0,10);
    b ~ normal(0,10);
    Rho ~ lkj_corr(2);
    sigma ~ cauchy(0,2);
    sigma_group ~ cauchy(0,2);
    for ( i in 1:N )
        mu[i] <- a + v[group[i],1] + 
                (b + v[group[i],2])*x[i];
    y ~ normal(mu,sigma);
}'

# now Cholesky non-centered version

mc2 <- "
data{
    int N;
    int N_groups;
    real x[N];
    real y[N];
    int group[N];
}
parameters{
    real a;
    real b;
    real<lower=0> sigma;
    vector<lower=0>[2] sigma_group;
    cholesky_factor_corr[2] L_Rho;
    matrix[2,N_groups] z;       // z-scores for constructing varying effects
}
transformed parameters {
    matrix[N_groups,2] v;       // same indexing as: vector[2] v[N-groups]
    matrix[2,2] Rho;
    v <- (diag_pre_multiply(sigma_group,L_Rho) * z)';   // note transpose
    Rho <- L_Rho * L_Rho';
}
model{
    vector[N] mu;
    to_vector(z) ~ normal(0,1); // sample z-scores for varying effects
    a ~ normal(0,10);
    b ~ normal(0,10);
    L_Rho ~ lkj_corr_cholesky(2);
    sigma ~ cauchy(0,2);
    sigma_group ~ cauchy(0,2);
    for ( i in 1:N )
        mu[i] <- a + v[group[i],1] + 
                (b + v[group[i],2])*x[i];
    y ~ normal(mu,sigma);
}"

# fit both

dat <- list(
    N = N,
    N_groups = N_groups,
    y = y,
    x = x,
    group = group
)

m1 <- stan( model_code=mc1 , data=dat , chains=2 , warmup=1000 , iter=3000 , refresh=500 )

m2 <- stan( model_code=mc2 , data=dat , chains=2 , warmup=1000 , iter=3000 , refresh=500 )

# compare inference
library(rethinking)
coeftab(m1,m2)

