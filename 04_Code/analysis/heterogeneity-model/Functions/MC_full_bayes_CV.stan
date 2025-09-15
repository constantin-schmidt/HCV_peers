// 1) Inputs to the program
data {
  int<lower=0> m;           // Number of time points
  int<lower=0> k;           // Number of latent factors
  int<lower=0> n;           // Number of units
  int<lower=0> n_treat;     // Number of treated observations
  int<lower=0> v;           // Number of spline parameters
  int<lower=0> y[n, m];     // Y matrix with units in rows and times in columns
  int<lower=0> g[n];        // Vector of treatment start time periods
  int<lower=0> p[n, m];     // Current exposure matrix with units in rows and times in columns
  int<lower=0> ld[n, m];    // Lockdown / treatment interaction
  real cps[n, m, v];        // Cumulative exposures matrix with units in rows and times in columns
  int<lower=0> atp;         // First time period ODNs could employ peers
  int<lower=0> o[n];        // Observations left out for cross-validation
}
  
// 2) Transformations of inputs
transformed data {
  int delta_p[n, m - atp + 1];

  for (i in 1:n) {
    for (j in atp:m) {
      delta_p[i, j - atp + 1] = p[i, j] - p[i, j - 1];
    }
  }
}

// 3) Definition of parameters
parameters {
  // Outcome model
  vector[n] d0;
  row_vector[m-1] c0;
  matrix[n, k] L;
  matrix[m, k] FS;
  vector[v] theta1;
  real theta2;
  
  real<lower=0> sqrt_disp_0;
  real<lower=0> sqrt_disp_a;
  
  // Intervention assignment model
  real gamma0;
  real gamma1;
  vector[k] gamma2;
}

// 4) Transformations of parameters
transformed parameters {
  // Outcome model
  matrix[n, m] Ups;
  matrix[n, m] lin_pred_count;
  matrix[n, m] lin_pred;
  matrix<lower=0>[n, m] Mu;
  
  real<lower=0> phi_0 = 1/sqrt_disp_0^2;
  real<lower=0> phi_a = 1/sqrt_disp_a^2;
  
  // Intervention assignment model
  matrix[n, m - atp + 1] lin_pred_assgn;
  matrix<lower=0>[n, m - atp + 1] Mu_assgn;

  // Outcome model
  Ups = L * FS';
  lin_pred_count = rep_matrix(d0, m);
  
  for (j in 2:m) {
    lin_pred_count[, j] += c0[j-1];
  }

  lin_pred_count += Ups;
  lin_pred = lin_pred_count + theta2 * to_matrix(ld);
  
  for (l in 1:v) {
    lin_pred += theta1[l] * to_matrix(cps[, , l]);
  }

  Mu = exp(lin_pred);
  
  // Intervention assignment model
  lin_pred_assgn = rep_matrix(gamma0 + gamma1 * d0 + L * gamma2, m - atp + 1);
  Mu_assgn = exp(lin_pred_assgn);
}

// 5) Model
model {
  // Priors
  
  // Outcome model
  d0 ~ normal(0, 50);
  c0 ~ normal(0, 10);
  to_vector(FS) ~ normal(0, 10);
  to_vector(L) ~ normal(0, 50);
  theta1 ~ normal(0, 10);
  theta2 ~ normal(0, 10);
  sqrt_disp_0 ~ normal(0, .11);
  sqrt_disp_a ~ normal(0, .1165);
  
  // Intervention assignment model
  gamma0 ~ normal(0, 10);
  gamma1 ~ normal(0, 10);
  gamma2 ~ normal(0, 10);
  
  // Likelihood
  for (i in 1:n) {
    // Outcome model with cross-validation
    for (j in 1:m) {
      if (j != o[i]) {
        if (j < g[i]) {
          y[i, j] ~ neg_binomial_2(Mu[i, j], phi_0);
        } else {
          y[i, j] ~ neg_binomial_2(Mu[i, j], phi_a);
        }
      }
    }

    // Intervention assignment model
    for (j in 1:m - atp + 1) {
      if (j != o[i]) {
        delta_p[i,j] ~ poisson(Mu_assgn[i,j]);
      }
    }
  }
}

// 6) Generated quantities
generated quantities {
  // No generated quantities
}
