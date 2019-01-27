data {
  int<lower=1> N;                  // number of observations
  real S[N];                       // saccade amplitude
  int<lower=6,upper=13> E[N];    
  int<lower=8,upper=12> meanE[N];    
  int<lower=1,upper=3> posu[N];    
  int<lower=1> J;                  // number of subjects
  int<lower=1,upper=J> id[N];      // subject id
}

parameters {
  vector[6] beta;               // fixed-effects parameters
  vector[3] alpha;              // compression weight
  real<lower=0> sigma_e;        // residual std
  vector<lower=0>[9] sigma_u;   // random effects standard deviations
  cholesky_factor_corr[9] L_u;  // L_u is the Choleski factor of the correlation matrix
  matrix[9,J] z_u;              // random effect matrix
}

transformed parameters {
  matrix[9,J] u;
  u = diag_pre_multiply(sigma_u, L_u) * z_u; // use Cholesky to set correlation
}

model {
  real mu; // conditional mean of the dependent variable

  //priors
  L_u ~ lkj_corr_cholesky(2);   // LKJ prior for the correlation matrix
  to_vector(z_u) ~ normal(0,1); // before Cholesky, random effects are normal variates with SD=1
  sigma_u ~ cauchy(0, 1);       // SD of random effects (vectorized)
  sigma_e ~ cauchy(0, 2);       // prior for residual standard deviation
  beta[1:3] ~ normal(0, 2);     // intercepts
  beta[4:6] ~ normal(1, 2);     // slopes
  alpha ~ normal(0, 1);

  //likelihood
  for (i in 1:N){
    mu = beta[posu[i]] + beta[3+posu[i]]*(alpha[posu[i]] * meanE[i] + (1-alpha[posu[i]])*E[i]);
    S[i] ~ normal(mu, sigma_e);
  }
}

