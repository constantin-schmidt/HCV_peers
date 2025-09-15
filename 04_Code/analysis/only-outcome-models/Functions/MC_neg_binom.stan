// 1) Inputs to the program
data {
  int<lower=0> m;           // Number of time points
  int<lower=0> k;           // Number of latent factors
  int<lower=0> n;           // Number of units
  int<lower=0> v;           // Number of spline parameters
  int<lower=0> y[n, m];     // Y matrix with units in rows and times in columns
  int<lower=0> g[n];        // Vector of treatment start time periods
  int<lower=0> ld[n, m];    // Lockdown / treatment interaction
  real cps[n, m, v];        // Cumulative exposures matrix with units in rows and times in columns
}

// 2) Transformations of inputs
transformed data {
  // No data transformations in this model
}

// 3) Definition of parameters
parameters {
  vector[n] d0;              // Unit-specific intercepts
  row_vector[m-1] c0;        // Time-specific effects
  matrix[n, k] L;            // Loadings
  matrix[m, k] FS;           // Factor scores
  vector[v] theta1;          // Vector of spline parameters for cumulative exposure
  real theta2;               // Effect of lockdown/treatment interaction
  real<lower=0> sqrt_disp_0; // Dispersion for observations with j < g[n]
  real<lower=0> sqrt_disp_a; // Dispersion for observations with j >= g[n]
}

// 4) Transformations of parameters
transformed parameters {
  matrix[n, m] Ups;         // Intermediate predictor
  matrix[n, m] lin_pred_count;      // Linear predictor without treatment effects
  matrix[n, m] lin_pred;    // Linear predictor
  matrix<lower=0>[n, m] Mu; // Mean parameter for negative binomial
  
  // Dispersion parameters
  real<lower=0> phi_0 = 1/sqrt_disp_0^2;
  real<lower=0> phi_a = 1/sqrt_disp_a^2;

  // Compute Ups as the product of L and the transpose of FS
  Ups = L * FS';

  // Initialize lin_pred_count with the intercepts
  lin_pred_count = rep_matrix(d0, m);
  
  // Add time-specific effects to lin_pred_count
  for (j in 2:m) {
    lin_pred_count[, j] += c0[j-1];
  }

  // Add the latent factor contributions (Ups)
  lin_pred_count += Ups;

  // Compute the expectation of the potential outcome
  lin_pred = lin_pred_count
             + theta2 * to_matrix(ld);
             
  // Add the effects of cumulative exposures
  for (l in 1:v) {
    lin_pred += theta1[l] * to_matrix(cps[, , l]);
  }

  // Exponentiate to get the mean parameter for the negative binomial
  Mu = exp(lin_pred);
}

// 5) Model
model {
  // Priors
  d0 ~ normal(0, 50);
  c0 ~ normal(0, 10);
  to_vector(FS) ~ normal(0, 10);
  to_vector(L) ~ normal(0, 50);
  theta1 ~ normal(0, 10);
  theta2 ~ normal(0, 10);
  sqrt_disp_0 ~ normal(0, .111);
  sqrt_disp_a ~ normal(0, .111);
  
  // Likelihood
  for (i in 1:n) {
    for (j in 1:m) {
      if (j < g[i]) {
        y[i, j] ~ neg_binomial_2(Mu[i, j], phi_0);
      } else {
        y[i, j] ~ neg_binomial_2(Mu[i, j], phi_a);
      }
    }
  }
}

// 6) Generated quantities
generated quantities {
  matrix[n, m] Count;                     // Baseline prediction

  // Compute the baseline prediction
  Count = exp(lin_pred_count);
}
