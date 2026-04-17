
# simplified version of debug_working that uses fixed z_init

# Load libraries
library(tidyverse)
library(Matrix)
library(here)
library(sf)
library(devtools)
library(dplyr)
library(fs)
library(here)
library(lwgeom) 
library(spacetime)
library(gstat)

# Complete package
load_all()

# Set seed
set.seed(3)
# Turn spherical geometry off
sf_use_s2(FALSE)

# Clock time
t0 <- Sys.time()

# --- Settings -----------------------------------------------------------------
# model parameters
# ell_km <- 120          # RFF length-scale in **kilometres**
# tau2   <- 0.5         # RW1 variance in feature space
p_init_true <- 0.05
z_init_true <- qlogis(p_init_true)

# inference parameters
D        <- 50         # number of random frequencies (try 200–500)
max_iter <- 5           # EM iterations
z_eps    <- 1e-6        # threshold for |z| to use n/4 in E[ω]
omega_floor <- 1e-10    # floor on ω to avoid 1/ω explosions
jitter_S <- 1e-10       # diagonal jitter for innovation covariance

# prediction parameters
nx <- 50
ny <- 50
t_vec <- 2000:2025
t_num <- length(t_vec)

# posterior draws
num_post_draws <- 100

# --------------------------- Load & filter data ----------------------------
# Read in prevalence data
# dat <- read.csv("R_ignore/R_scripts/data/all_mutations_get_prevalence.csv") |>
#   mutate(collection_day = as.Date(collection_day)) |>
#   select(survey_id, study_id, longitude, latitude, year, country_name, numerator, denominator, prevalence, mutation) 

# Read Africa shape files
africa_shp_admin1 <- readRDS(file = "R_ignore/R_scripts/data/sf_admin1_africa.rds")

# --- Spatial grid ---------------------------------------------------------
nx_sim <- 20
ny_sim <- 20
xs_sim <- seq(28, 48, length.out = nx_sim)
ys_sim <- seq(-4.6, 18, length.out = ny_sim)
grid_sim <- expand.grid(lon = xs_sim, lat = ys_sim)
S <- as.matrix(grid_sim[, c("lon", "lat")])
M <- nrow(S)

# --- Spatial covariance for latent fields --------------------------------
dist_mat <- as.matrix(dist(S))
ell_sim <- 3
sigma2_rw   <- 0.3   # innovation variance for year-to-year changes

K_rw   <- sigma2_rw   * exp(-(dist_mat^2) / (2 * ell_sim^2))
diag(K_rw)   <- diag(K_rw) + 1e-8
L_rw   <- t(chol(K_rw))

# --- Generate latent field over time -------------------------------------
z_true_array <- array(NA_real_, dim = c(nx_sim, ny_sim, t_num))
p_true_array <- array(NA_real_, dim = c(nx_sim, ny_sim, t_num))

# Initial year
z_prev <- rep(z_init_true, M)

for (tt in 1:t_num) {
  if (tt > 1) {
    innovation <- as.numeric(L_rw %*% rnorm(M))
    z_prev <- z_prev + innovation   # random walk over time
  }
  
  p_prev <- plogis(z_prev)
  
  z_true_array[, , tt] <- matrix(z_prev, nrow = nx_sim, ncol = ny_sim, byrow = FALSE)
  p_true_array[, , tt] <- matrix(p_prev, nrow = nx_sim, ncol = ny_sim, byrow = FALSE)
}

# --- Sample observations for each year -----------------------------------
samp_years <- 2018:2023
n_samp_years <- length(samp_years)
n_obs_per_year <- 50

dat_list <- vector("list", n_samp_years)

for (tt in 1:n_samp_years) {
  obs_idx <- sort(sample(seq_len(M), n_obs_per_year))
  denom <- sample(c(20, 50, 100), n_obs_per_year, replace = TRUE)
  
  p_tt <- as.vector(p_true_array[, , match(samp_years, t_vec)[tt]])
  y_obs <- rbinom(n_obs_per_year, size = denom, prob = p_tt[obs_idx])
  
  dat_list[[tt]] <- data.frame(
    longitude = S[obs_idx, 1],
    latitude = S[obs_idx, 2],
    year = samp_years[tt],
    numerator = y_obs,
    denominator = denom,
    mutation = "fake"
  )
}

dat <- dplyr::bind_rows(dat_list)

# --------------------------- Fliter admin regions ----------------------------
target_crs <- 4326
africa_shp_admin1 <- st_transform(africa_shp_admin1, crs = target_crs)

# Define East Africa box and crop Admin1 shapefile
bbox_east_africa <- get_east_africa_bbox(target_crs)

# Using the cropped Admin 1 data for East Africa
admin1_regions <- st_make_valid(africa_shp_admin1) |> st_crop(bbox_east_africa)

# # --- Add K13 overall data -----------------------------------------------------
# dat_with_k13 <- add_combined_k13(dat)
# 
# # --- Define K13 mutants -------------------------------------------------------
# args <- commandArgs(trailingOnly = TRUE)
# mut <- "k13:comb"

# Per-mutation clock start
mut_t0 <- Sys.time()

# --- Filter data based on mutation --------------------------------------------
# dat_sub <- dat_with_k13 |>
#   filter(mutation == mut) |>
#   filter(is.finite(numerator), is.finite(denominator), denominator > 0) |>
#   mutate(
#     p_hat = (numerator + 0.5) / (denominator + 1),
#     z     = qlogis(p_hat)
#   )

dat_sub <- dat |>
  filter(is.finite(numerator), is.finite(denominator), denominator > 0) |>
  mutate(
    p_hat = (numerator + 0.5) / (denominator + 1),
    z     = qlogis(p_hat)
  )

# --- Project to local Cartesian (km) ------------------------------------------
# Center projection on study area for good local properties
lon0 <- median(dat_sub$longitude, na.rm = TRUE)
lat0 <- median(dat_sub$latitude,  na.rm = TRUE)

# Local Lambert Azimuthal Equal-Area (units: meters)
crs_laea <- paste0(
  "+proj=laea +lat_0=", lat0,
  " +lon_0=", lon0,
  " +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
)

# Make sf points and transform
pts_ll <- st_as_sf(dat_sub, coords = c("longitude", "latitude"), crs = 4326)
pts_xy <- st_transform(pts_ll, crs_laea)

# Extract projected coordinates in **km**
xy_m <- st_coordinates(pts_xy)
xy_km <- xy_m / 1000
dat_sub$Xkm <- xy_km[, 1]
dat_sub$Ykm <- xy_km[, 2]

# ==============================================================================
# Load spatial and temporal parameters
# ==============================================================================

ell_km <- 3
tau2   <- 0.3
#z_init <- qlogis(0.5)
z_init <- z_init_true

# --- Random Fourier Features --------------------------------------------------
# Draw D frequencies ω_j ~ N(0, ℓ⁻² I₂), ℓ in **km**
Omega <- matrix(rnorm(D * 2, mean = 0, sd = 1 / ell_km), ncol = 2)  # [D x 2]

# Helper to build feature matrix Φ(S) for coords S = [x_km, y_km] in km
phi_of_coords <- function(S_km) {
  OS <- as.matrix(S_km) %*% t(Omega)                 # (n x D)
  cbind(cos(OS), sin(OS)) * (1 / sqrt(D))            # (n x 2D)
}

# --- Create grid for prediction ---------------------------------------------------
# Define xlim and ylim for each mutation
xs <- xs_sim
ys <- ys_sim
nx <- nx_sim
ny <- ny_sim
grid_ll <- grid_sim
grid_km <- as.matrix(grid_ll[, c("lon", "lat")])
lon_min <- min(xs)
lon_max <- max(xs)
lat_min <- min(ys)
lat_max <- max(ys)

grid_xy <- as.matrix(grid_ll[, c("lon", "lat")])
dat_sub$Xkm <- dat_sub$longitude
dat_sub$Ykm <- dat_sub$latitude

# Φ on projected grid (km)
OS_full   <- grid_km %*% t(Omega)
Phi_full  <- cbind(cos(OS_full), sin(OS_full)) * (1 / sqrt(D))

# --- Organize by time ---------------------------------------------------------
# For each time index, collect projected coords and Binomial stats

obs_by_t <- vector("list", t_num)
for (ti in 1:t_num) {
  dat_t <- filter(dat_sub, year == t_vec[ti])
  if (nrow(dat_t) > 0) {
    S_t    <- cbind(dat_t$Xkm, dat_t$Ykm)          # (m_t x 2) projected coords (km)
    Phi_t  <- phi_of_coords(S_t)                   # (m_t x 2D)
    y_t    <- as.numeric(dat_t$numerator)          # positives
    n_t    <- as.numeric(dat_t$denominator)        # totals
    kappa  <- y_t - 0.5 * n_t                      # centered counts
    mu_t   <- rep(z_init, length(y_t))             # constant logit mean
    obs_by_t[[ti]] <- list(Phi = Phi_t, y = y_t, n = n_t, kappa = kappa, mu = mu_t)
  }
}

# --------------------------- State-space model -----------------------------
# Random walk in feature space: w_t = w_{t-1} + ξ_t,  ξ_t ~ N(0, τ² I)
Qw   <- Diagonal(2 * D, tau2)
m_w  <- matrix(0, 2 * D, t_num)          # smoothed means (initially 0)
m0_w <- rep(0, 2 * D)                    # prior mean
C0_w <- Matrix(0, 2 * D, 2 * D, sparse = FALSE)  # prior covariance (zero variance)

# --------------------------- EM loop ---------------------------------------
m_pred <- vector("list", t_num)
C_pred <- vector("list", t_num)
m_filt <- matrix(NA, 2 * D, t_num)
C_filt <- vector("list", t_num)

for (i in 1:max_iter) {
  message(sprintf("\nEM iteration %d/%d", i, max_iter))
  
  # ===== E-step: Expected PG weights =======================================
  Omega_list <- vector("list", t_num)
  r_list     <- vector("list", t_num)
  
  for (ti in 1:t_num) {
    ob <- obs_by_t[[ti]]
    if (!is.null(ob)) {
      Phi_t  <- ob$Phi
      mu_t   <- ob$mu
      n_t    <- ob$n
      kappa  <- ob$kappa
      
      # Current linear predictor z_t(s) = μ + Φ w_t
      z_t <- as.vector(mu_t + Phi_t %*% m_w[, ti])
      
      # Expected PG weights: E[ω] = (n / (2|z|)) * tanh(|z|/2)
      abs_z <- pmax(abs(z_t), z_eps)
      omega <- 0.5 * n_t / abs_z * tanh(0.5 * abs_z)
      near0 <- (abs(z_t) < z_eps) | !is.finite(omega)
      omega[near0] <- n_t[near0] / 4
      omega <- pmax(omega, omega_floor)
      
      Omega_list[[ti]] <- omega
      r_list[[ti]]     <- kappa / omega
    }
  }
  
  # ===== M-step: KF + RTS in feature space =================================
  # Observation model: r_t = μ_t + Φ_t w_t + ε_t,  ε_t ~ N(0, Ω_t^{-1})
  m_prev <- m0_w
  C_prev <- C0_w
  ll_sum <- 0
  
  # ---------- Forward (Kalman filter) ----------
  for (ti in 1:t_num) {
    m_t_pred <- m_prev
    C_t_pred <- C_prev + Qw
    
    ob <- obs_by_t[[ti]]
    if (is.null(ob)) {
      m_t_filt <- m_t_pred
      C_t_filt <- C_t_pred
    } else {
      Phi_t <- ob$Phi
      mu_t  <- ob$mu
      omega <- Omega_list[[ti]]
      r_t   <- r_list[[ti]]
      
      m_t <- length(omega)
      R_t <- Diagonal(m_t, 1 / omega)               # Var(ε_t) = Ω_t⁻¹
      y_c <- r_t - mu_t                             # centered pseudo-obs
      
      # Innovation covariance S_t = Φ P Φ^T + R
      A   <- Phi_t %*% C_t_pred
      S_t <- forceSymmetric(A %*% t(Phi_t) + R_t)
      diag(S_t) <- diag(S_t) + jitter_S             # ensure PD
      
      # Innovation ν_t
      v_t <- y_c - as.vector(Phi_t %*% m_t_pred)
      
      # Log-likelihood increment (pseudo-Gaussian)
      Rch <- chol(S_t)
      logdetS <- 2 * sum(log(diag(Rch)))
      quad <- drop(t(solve(Rch, solve(t(Rch), v_t))) %*% v_t)
      ll_sum <- ll_sum - 0.5 * (logdetS + quad + m_t * log(2*pi))
      
      # Kalman gain + update
      B   <- C_t_pred %*% t(Phi_t)
      K_t <- t(solve(S_t, t(B)))
      m_t_filt <- as.vector(m_t_pred + K_t %*% v_t)
      C_t_filt <- forceSymmetric(C_t_pred - K_t %*% A)
    }
    
    # Store + propagate
    m_pred[[ti]] <- m_t_pred
    C_pred[[ti]] <- C_t_pred
    m_filt[, ti] <- m_t_filt
    C_filt[[ti]] <- C_t_filt
    m_prev <- m_t_filt
    C_prev <- C_t_filt
  }
  
  # ---------- Backward (RTS smoother) ----------
  m_smooth <- m_filt
  C_smooth <- C_filt
  
  for (ti in (t_num - 1):1) {
    Ptt  <- C_filt[[ti]]
    Pnpt <- C_pred[[ti + 1]]
    mnpt <- m_pred[[ti + 1]]
    
    Rch <- chol(forceSymmetric(Pnpt))
    J_t <- Ptt %*% solve(Rch, solve(t(Rch), Diagonal(nrow(Pnpt))))
    
    m_smooth[, ti] <- as.vector(
      m_filt[, ti] + J_t %*% (m_smooth[, ti + 1] - mnpt)
    )
    
    C_smooth[[ti]] <- forceSymmetric(
      Ptt + J_t %*% (C_smooth[[ti + 1]] - Pnpt) %*% t(J_t)
    )
  }
  
  # Update EM parameters
  m_w_old <- m_w
  m_w     <- m_smooth
  num <- sqrt(sum((m_w - m_w_old)^2))
  den <- sqrt(sum(m_w_old^2) + 1e-12)
  rel_change <- num / den
  
  message(sprintf("  pseudo-Gaussian log-lik: %.3f,  rel_change(w): %.3e", ll_sum, rel_change))
  
}

# ================================================
# Joint posterior sampling of β_t and surfaces
# ================================================

# Storage:
# - beta_draws:  (2D x t_num x B)      feature-space trajectories
# - p_draws:     (nx x ny x t_num x B) prevalence surfaces
beta_draws <- array(NA_real_, dim = c(2 * D, t_num, num_post_draws))
p_draws    <- array(NA_real_, dim = c(nx, ny, t_num, num_post_draws))

for (b in 1:num_post_draws) {
  message(sprintf("Posterior draw %s of %s", b, num_post_draws))
  
  # --------------------------------------------
  # 1. Backward sampling in feature space (β_t)
  # --------------------------------------------
  beta_b <- matrix(NA_real_, nrow = 2 * D, ncol = t_num)
  
  # Sample final time T from N(m_filt[,T], C_filt[[T]])
  T_idx <- t_num
  P_T   <- forceSymmetric(C_filt[[T_idx]])
  Rch_T <- chol(P_T)
  eta_T <- rnorm(2 * D)
  beta_b[, T_idx] <- as.numeric(m_filt[, T_idx] + Rch_T %*% eta_T)
  
  # Backward recursion for t = T-1,...,1
  for (ti in (t_num - 1):1) {
    # Filtered mean/cov at time t
    a_t <- m_filt[, ti]                 # β_{t|t}
    P_t <- forceSymmetric(C_filt[[ti]]) # P_{t|t}
    
    # One-step-ahead prediction at time t+1
    a_tp1 <- m_pred[[ti + 1]]           # β_{t+1|t}
    R_tp1 <- forceSymmetric(C_pred[[ti + 1]])  # P_{t+1|t}
    
    # Smoothing gain J_t = P_{t|t} (P_{t+1|t})^{-1}
    Rch <- chol(R_tp1)
    J_t <- P_t %*% solve(Rch, solve(t(Rch), Diagonal(nrow(R_tp1))))
    
    # Conditional mean and covariance of β_t | β_{t+1}, y
    mean_t <- as.numeric(a_t + J_t %*% (beta_b[, ti + 1] - a_tp1))
    Var_t  <- forceSymmetric(P_t - J_t %*% R_tp1 %*% t(J_t))
    
    # Draw β_t ~ N(mean_t, Var_t)
    Rch_t <- chol(Var_t)
    eta_t <- rnorm(2 * D)
    beta_b[, ti] <- as.numeric(mean_t + Rch_t %*% eta_t)
  }
  
  # Store β path
  beta_draws[, , b] <- beta_b
  
  # --------------------------------------------
  # 2. Map β_t path to full surfaces z(s,t), p(s,t)
  # --------------------------------------------
  for (ti in 1:t_num) {
    z_vec <- as.numeric(z_init + Phi_full %*% beta_b[, ti])   # length M = nx*ny
    z_mat <- matrix(z_vec, nrow = nx, ncol = ny, byrow = FALSE)
    p_mat <- plogis(z_mat)
    
    p_draws[, , ti, b] <- p_mat
  }
}

# ============================================================
# Calculate mean, median and difference in CI from draws
# ============================================================
p_post_mean = apply(p_draws, MARGIN = c(1, 2, 3), FUN = mean)
p_post_median = apply(p_draws, MARGIN = c(1, 2, 3), FUN = median)
p_post_CI <- apply(p_draws, MARGIN = c(1, 2, 3), FUN = function(x) {
  diff(quantile(x, probs = c(0.025, 0.975)))
})
p_post_lower <- apply(p_draws, MARGIN = c(1, 2, 3), FUN = function(x) {
  quantile(x, probs = 0.025, na.rm = TRUE)
})
p_post_upper <- apply(p_draws, MARGIN = c(1, 2, 3), FUN = function(x) {
  quantile(x, probs = 0.975, na.rm = TRUE)
})

model_output <- list(
  p_post_mean = p_post_mean,
  p_post_median = p_post_median,
  p_post_CI = p_post_CI,
  
  # Grid information
  xs = xs,
  ys = ys,
  lon_min = lon_min,
  lon_max = lon_max,
  lat_min = lat_min,
  lat_max = lat_max,
  t_vec = t_vec,
  
  # Original data
  data_subset = dat_sub,
  
  # Posterior draws
  num_posteriors = num_post_draws,
  
  # Model settings
  settings = list(ell_km = ell_km, tau2 = tau2, D = D)
)

plot_times <- seq(2000, 2025, by = 5)

time_idx <- match(plot_times, t_vec)
l <- list()
for (i in seq_along(plot_times)) {
  ti <- time_idx[i]
  l[[i]] <- data.frame(
    x = grid_ll$lon,
    y = grid_ll$lat,
    year = plot_times[i],
    p_true = as.vector(p_true_array[, , ti]),
    z_true = as.vector(z_true_array[, , ti]),
    p_mean = as.vector(p_post_mean[, , ti]),
    p_lo = as.vector(p_post_lower[, , ti]),
    p_hi = as.vector(p_post_upper[, , ti])
  )
}
plot_df <- bind_rows(l) |>
  mutate(p_width  = p_hi - p_lo)


plot1 <- ggplot(plot_df, aes(x = x, y = y, fill = p_true)) +
  geom_raster() +
  geom_point(data = filter(dat, year %in% plot_times),
             aes(x = longitude, y = latitude, fill = numerator / denominator),
             colour = "grey", pch = 21, size = 1.0) +
  scale_fill_viridis_c(name = "", limits = c(0, 1)) +
  ggtitle("True prevalence") +
  facet_wrap(~year, nrow = 1)

plot2 <- ggplot(plot_df, aes(x = x, y = y, fill = p_mean)) +
  geom_raster() +
  geom_point(data = filter(dat, year %in% plot_times),
             aes(x = longitude, y = latitude, fill = numerator / denominator),
             colour = "grey", pch = 21, size = 1.5) +
  scale_fill_viridis_c(name = "", limits = c(0, 1)) +
  ggtitle("Posterior mean prevalence") +
  facet_wrap(~year, nrow = 1)

plot3 <- ggplot(plot_df, aes(x = x, y = y, fill = p_width)) +
  geom_raster() +
  scale_fill_viridis_c(name = "") +
  ggtitle("95% interval width") +
  facet_wrap(~year, nrow = 1)

cowplot::plot_grid(plot1, plot2, plot3, ncol = 1)
