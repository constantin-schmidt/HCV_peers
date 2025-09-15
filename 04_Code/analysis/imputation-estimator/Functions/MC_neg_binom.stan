// 1) Inputs to the program
data {
  int<lower = 0> m;    // Number of time points
  int<lower = 0> k;    // Number of latent factors
  int<lower = 0> n;    // Number of units
  int<lower = 2> g[n]; // Treatment start times
  int<lower = 0> y[n,m]; // Y(0) matrix with units in rows and times in columns
}

// 2) Transformations of inputs
transformed data{
  
}

// 3) Definition of parameters
parameters {
  vector[n] d0;           // Time-specific specific intercepts
  row_vector[m-1] c0;     // ODN effects (time-specific)
  matrix[n,k] L;          // Loadings (latent factor loadings)
  matrix[m,k] FS;         // Factor scores (latent factor contributions across time)
  real<lower=0> sqrt_disp;  // Setting appropriate prior for dispersion parameter
}

// 4) Transformations of parameters
transformed parameters {
  matrix[n,m] Ups;               // Intermediate predictor
  row_vector<lower=0>[m] Mu[n];  // Negative binomial mean
  
  // Dispersion parameter
  real<lower=0> phi = 1/sqrt_disp^2; // Dispersion parameter for the negative binomial

  // Latent predictors
  Ups = L * FS';

  // Mean structure
  for (i in 1:n) {
    Mu[i,1] = exp(d0[i] + Ups[i,1]);
    for (t in 2:m) {
      Mu[i,t] = exp(d0[i] + c0[t-1] + Ups[i,t]);
    }
  }
}

// 5) Model
model {
  to_vector(L) ~ normal(0, 50);  
  to_vector(FS) ~ normal(0, 10); 
  d0 ~ normal(0, 50);   
  c0 ~ normal(0, 10);
  sqrt_disp ~ normal(0, .11);

  // Likelihood terms
  for (i in 1:n) {
    for (t in 1:g[i]-1) {
      y[i,t] ~ neg_binomial_2(Mu[i,t], phi);
    }
  }
}

// 6) Generated quantities
generated quantities {
  // Can add generated quantities if needed
}
