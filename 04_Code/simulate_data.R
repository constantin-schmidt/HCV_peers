##  Simulates some synthetic data ##
####################################
library(tidyr)
library(dplyr)
library(splines)
set.seed(42)

##  Set parameter values  ##
# --- parameters for A_it ---
delta0        <- -2.5
# --- parameters for Y_it ---
phi0 <- 16                 # NB dispersion under no exposure (size in rnbinom)
phi1 <- 18                 # NB dispersion under exposure
theta1 <- 0.3              # lockdown effect in ψ(·)
# f(U_it) = κ_i + β_t + λ_i V_t
beta_t_amp <- 0.10          # seasonal amplitude for β_t
# V_t: linear trend (scaled)
# spline coefficients (p=3 cubic, b*=3 internal knots -> 6 basis coefs)
wb <- c( 0.009, 0.05, .28, 0.13, 0.18, 0.05)
# --- parameters for Y_it and A_it ---
delta_kappa   <- -.4
delta_lambda  <- 0.15

##  ODN ##
ODN <- c(
  "North 1","North 2","North 3","North 4",
  "South 1","South 2","South 3","South 4",
  "East 1","East 2","East 3","East 4",
  "West 1","West 2","West 3","West 4",
  "North-East 1","North-East 2",
  "North-West 1","North-West 2",
  "South-East 1","South-West 1"
)
n_odns <- length(ODN)

##  Month  ##
Month <- format(
  seq(as.Date("2016-06-01"), as.Date("2021-05-01"), by = "month"),
  "%Y-%m")

##  Expand  ##
HEP.dta <- expand_grid(
  ODN = ODN,
  Month  = Month
)

##  Number of Peers ##
t_min <- as.Date("2018-01-01") 

# 5 ODNs with zero peers forever
no_peer_odns <- sample(unique(HEP.dta$ODN), 5)

# prep dates and order
HEP.dta <- HEP.dta |> 
  mutate(MonthDate = as.Date(paste0(Month, "-01"))) |> 
  arrange(ODN, MonthDate) 
# unit-level latent traits: U_i = (kappa_i, lambda_i) with scalar lambda
unit_latent <- tibble( ODN = unique(HEP.dta$ODN),
                       kappa_i = rnorm(n_odns, 3, .5),
                       lambda_i = rnorm(n_odns, 0, .5) )

# log(mu_it) = δ0 + g(U_i) with g(U_i)=δκ κ_i + δλ λ_i  (constant over t for unit i)
HEP.dta <- HEP.dta |>
  left_join(unit_latent, by = "ODN") |>
  group_by(ODN) |>
  arrange(MonthDate, .by_group = TRUE) |>
  mutate(
    # Unit mu
    log_mu = delta0 + delta_kappa * kappa_i + delta_lambda * lambda_i,
    mu     = exp(log_mu),
    
    # Poisson draw per row; 0 before t_min
    Mit_raw = if_else(MonthDate < t_min, 0L, rpois(n(), mu)),
    
    # Accumulate from 0 at t_min automatically
    Ait = cumsum(Mit_raw),
    
    # Enforce always-zero ODNs and cap at 3
    NumberOfPeers = if_else(ODN %in% no_peer_odns, 0L, Ait),
    
    Dit = as.integer(NumberOfPeers > 0),
    
    # first month with peers > 0 (recycled across rows)
    InterventionStart = if (any(NumberOfPeers > 0))
      Month[which(NumberOfPeers > 0)[1]] else NA_character_
  ) |>
  ungroup()


##  Lockdown indicator  ##
HEP.dta <- HEP.dta |>
  mutate(first_lockdown = if_else(Month %in% c("2020-03", "2020-04", "2020-05"),
                                  1, 0))


##  Outcome ##
# --- build unit/time pieces used in f(U_it) ---
HEP.dta <- HEP.dta |>
  arrange(ODN, MonthDate) |>
  group_by(ODN) |>
  mutate(
    # cumulative exposure z_it = sum_{j=1}^t a_ij
    z = cumsum(NumberOfPeers)
  ) |>
  ungroup() |>
  mutate(
    mon = as.integer(format(MonthDate, "%m")),
    beta_t = beta_t_amp * sin(2*pi*mon/12),                      # seasonality
    Vt = as.numeric(MonthDate - min(MonthDate))/30.44,           # months since start
    Vt = scale(Vt)[,1]                                           # scaled
  )

# --- spline s(z) with p=3, b*=3 ---
# internal knots at quartiles of positive exposure (3 knots)
posz <- HEP.dta$z[HEP.dta$z > 0]
knots <- if (length(posz) >= 4) quantile(posz, probs = c(.25,.5,.75)) else c(1,2,3)
Bs <- bs(HEP.dta$z,
         degree = 3,
         knots = knots,
         Boundary.knots = c(0, max(HEP.dta$z)),
         intercept = FALSE)

s_z <- as.numeric(Bs %*% wb)

# --- f(U_it) = κ_i + β_t + λ_i V_t ---
HEP.dta <- HEP.dta |>
  mutate(
    fUit = kappa_i + beta_t + lambda_i * Vt,
    # ψ(ā_it, X_it; θ) = s(z_it) + θ1 * 1{first lockdown}; (no X term)
    psi  = s_z + theta1 * first_lockdown,
    
    # log q^0_it and log q^1_it
    log_q0 = fUit,
    log_q1 = fUit + psi,
    
    # choose regime by whether cumulative exposure up to t is zero
    use_exposed = z > 0,
    q = ifelse(use_exposed, exp(log_q1), exp(log_q0)),
    phi = ifelse(use_exposed, phi1, phi0),
    
    # NegBin draw: size = phi, mean = q
    received_blueteq = rnbinom(n(), size = phi, mu = q)
  )

##  Save  ##
HEP.dta <- HEP.dta |>
  select(ODN, Month, NumberOfPeers, Dit, InterventionStart, first_lockdown, 
         received_blueteq)
save(HEP.dta, file = "./03_Data/HEP_dta.RData")
