# ============================================================
# RFF + PG-EM with Kalman filter/RTS in projected (Cartesian) coords
# ============================================================
# - Reads WHO GET prevalence data
# - Filters to a mutation and spatial window
# - Projects lon/lat -> local LAEA (km), runs RFF + PG-EM in km
# - Reconstructs prevalence on a lon/lat grid for plotting
# - Overlays observations with an in-window/out-of-window styling
# ============================================================

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
set.seed(1)
# Turn spherical geometry off
sf_use_s2(FALSE)

# Clock time
t0 <- Sys.time()

# --- Helper Function ----------------------------------------------------------
# Build consistent filenames per (mutation, length_space, length_time, n_features)
cache_path <- function(mut, what, ext = c("rds", "parquet", "rds.gz")) {
  ext <- match.arg(ext)
  fn  <- sprintf("%s_%s.%s",
                 gsub("[:/ ]", "_", mut),
                 what,
                 ext)
  file.path(CACHE_DIR, fn)
}

# --- Settings -----------------------------------------------------------------
# model parameters
# ell_km <- 120          # RFF length-scale in **kilometres**
# tau2   <- 0.5         # RW1 variance in feature space
p_init <- 0.001
z_init <- qlogis(p_init)

# inference parameters
D        <- 500         # number of random frequencies (try 200–500)
max_iter <- 5           # EM iterations
z_eps    <- 1e-6        # threshold for |z| to use n/4 in E[ω]
omega_floor <- 1e-10    # floor on ω to avoid 1/ω explosions
jitter_S <- 1e-10       # diagonal jitter for innovation covariance

# prediction parameters
nx <- 200
ny <- 200
t_vec <- 1995:2025
t_num <- length(t_vec)

plot_times <- seq(2012, 2025, by = 1)

# posterior draws
num_post_draws <- 100 

# parameters for exceedance plot
cell_size_km <- 5
half_m <- (cell_size_km * 1000) / 2
q_thr <- 0.80   # probability threshold for exceedance probability
bootstrap    <- 500    # bootstrap resamples

# --------------------------- Load & filter data ----------------------------
CACHE_DIR <- "R_ignore/R_scripts/outputs/model_outputs/GRFF_model_output_key_WHO_mutations"
dir_create(CACHE_DIR)

# Read in prevalence data
dat <- read.csv("R_ignore/R_scripts/data/all_who_get_prevalence.csv") |>
  mutate(collection_day = as.Date(collection_day)) |>
  select(longitude, latitude, year, numerator, denominator, prevalence, mutation) 

# Read Africa shape files
africa_shp_admin1 <- readRDS(file = "R_ignore/R_scripts/data/sf_admin1_africa.rds")

# --------------------------- Fliter admin regions ----------------------------
target_crs <- 4326
africa_shp_admin1 <- st_transform(africa_shp_admin1, crs = target_crs)

# Define East Africa box and crop Admin1 shapefile
bbox_east_africa <- st_bbox(c(
  xmin = 28.48,
  xmax = 48.43,
  ymin = -4.6,
  ymax = 15.29
), crs = target_crs)

# Using the cropped Admin 1 data for East Africa
admin1_regions <- st_make_valid(africa_shp_admin1) |> st_crop(bbox_east_africa)

# --- Add K13 overall data -----------------------------------------------------
k13_any_site_year <- dat %>%
  # keep only K13 mutations; adjust this filter to your naming convention
  filter(grepl("^k13", mutation, ignore.case = TRUE)) %>%
  group_by(year, longitude, latitude) %>%
  summarise(
    numerator = sum(numerator, na.rm = TRUE),
    denominator = max(denominator, na.rm = TRUE),
    prevalence = numerator / denominator *100,
    .groups = "drop"
  ) %>%
  transmute(
    year,
    longitude = longitude,
    latitude  = latitude,
    mutation  = "k13:comb", 
    numerator = numerator,
    denominator = denominator,
    prevalence = prevalence
  )

dat_with_k13 <- dat %>%
  bind_rows(k13_any_site_year) %>%
  arrange(year, longitude, latitude, mutation)

# --- Define K13 mutants -------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
# mut <- args[[1]]
mut <- "k13:675:V"

# Per-mutation clock start
mut_t0 <- Sys.time()

# --- Filter data based on mutation --------------------------------------------
dat_sub <- dat_with_k13 |>
  filter(mutation == mut) |>
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
# Variogram analysis to identify optimal length space parameter (ell_km) 
# and temporal variance (tau2)
# ==============================================================================

# --- Spatial variogram to estimate the spatial length-scale ell_km ------------
# Construct a spatiotemporal object (STIDF) from the data subset
# Spatial coordinates in kilometres (projected coordinates)
sp_pts <- sp::SpatialPoints(dat_sub[, c("Xkm", "Ykm")])

# Create time index where time is treated as numeric across all years in dat_sub
time_idx <- as.POSIXct(paste0(dat_sub$year, "-07-01"), tz = "UTC")

# Combine spatial points + time index + observed values into a spacetime object
st_obj <- spacetime::STIDF(sp_pts, time_idx, data = data.frame(z = dat_sub$z))

# Compute spatiotemporal variogram, then isolate the ~1-year temporal slice
# Small tolerance around 365 days to catch "approximately one year apart" pairs
eps_days <- 0.5

# Compute spatiotemporal variogram
vg_st_1year <- variogramST(
  z ~ 1,
  data   = st_obj,
  tlags  = c(365 - eps_days, 365 + eps_days, 3 * 365),  # boundaries for ≤2 years
  tunit  = "days",
  cutoff = 1500,
  width  = 25
)

# Extract only the ~1-year temporal lag bin
vg_1y <- subset(vg_st_1year, timelag == 365)
# Keep only spatial variogram columns
vg_1y_sp <- vg_1y[, c("dist", "gamma", "np")]

# Defensive casting to numeric (variogramST output can sometimes be factor-like)
vg_1y_sp$dist  <- as.numeric(vg_1y_sp$dist)
vg_1y_sp$gamma <- as.numeric(vg_1y_sp$gamma)
vg_1y_sp$np    <- as.numeric(vg_1y_sp$np)

# Assign gstat variogram class so gstat plotting/fitting methods work
class(vg_1y_sp) <- c("gstatVariogram", "data.frame")
# Quick diagnostic plot of the empirical spatial variogram at ~1-year lag
plot(vg_1y_sp)

# Fit a parametric variogram model to estimate the spatial range parameter
# Reasonable starting values
nugget0 <- min(vg_1y_sp$gamma)
psill0  <- max(vg_1y_sp$gamma) - nugget0
range0  <- 500

# Initialize a Gaussian variogram model
vgm_init <- vgm(
  psill  = psill0,
  model  = "Gau",
  range  = range0,
  nugget = nugget0
)

# Fit the model to the empirical variogram points
vgm_fit <- fit.variogram(vg_1y_sp, vgm_init)

# Extract fitted Gaussian range parameter (gstat’s "range" for model "Gau")
a_hat  <- vgm_fit$range[vgm_fit$model == "Gau"]

# Convert gstat Gaussian range parameter to an RFF-style length-scale ell
ell_km_hat <- a_hat / sqrt(2)

# Plot empirical variogram points with fitted model curve overlay
plot(
  gamma ~ dist,
  data = vg_1y_sp,
  xlab = "Distance (km)",
  ylab = "Semivariance",
  main = "Empirical variogram (≈1-year pairs) with exponential fit",
  pch = 1
)

# Create a smooth model-implied variogram curve across the observed distances
vg_line <- variogramLine(
  vgm_fit,
  maxdist = max(vg_1y_sp$dist),
  n = 200
)

# Overlay fitted variogram curve
lines(vg_line$dist, vg_line$gamma, lwd = 2, col = "black")

# --- Temporal variogram to estimate RW1 variance (tau^2) ----------------------

# Extract time index from spatiotemporal object
t_idx <- spacetime::index(st_obj@time)
# Maximum temporal separation across all observation pairs (in days)
max_lag_days <- as.numeric(diff(range(t_idx)), units = "days")

# Define temporal lag bins from 0 to maximum lag; here we bin by 180 days = 6 months
tlags_days <- seq(0, max_lag_days, by = 180)

# Compute spatiotemporal variogram
vg_st_time <- variogramST(
  z ~ 1,
  data   = st_obj,
  tlags  = tlags_days,  # covers ALL time lags in your data
  tunit  = "days",
  cutoff = 25,          # only pairs with distance ≤ 50 km
  width  = 25           # 1 spatial bin: [0, 50) km
)

# Sanity check: should be a single spatial bin (~25 km midpoint)
unique(vg_st_time$dist)

# Extract temporal variogram (collapse spatial dimension)
vg_t <- vg_st_time[, c("timelag", "gamma", "np")]
names(vg_t)[names(vg_t) == "timelag"] <- "dist"   # rename for gstat

# Remove bins with undefined or non-finite estimates
vg_t <- vg_t[is.finite(vg_t$dist) & is.finite(vg_t$gamma), ]

# Ensure numeric columns (defensive programming)
vg_t$dist  <- as.numeric(vg_t$dist)   # time lag in days
vg_t$gamma <- as.numeric(vg_t$gamma)
vg_t$np    <- as.numeric(vg_t$np)

# Set class so gstat methods (plot, fit.variogram, etc.) work correctly
class(vg_t) <- c("gstatVariogram", "data.frame")

# Visual inspection of temporal variogram
plot(gamma ~ dist, data = vg_t,
     xlab = "Time lag (days)", ylab = "Semivariance",
     main = "Temporal variogram (pairs ≤ 25 km)")

# Estimate RW variance via weighted linear regression (model: gamma(h) = 0.5 * tau^2 * h)
lm_rw0 <- lm(gamma ~ 0 + dist, data = vg_t, weights = np)
summary(lm_rw0)

# Estimated slope
b_hat <- coef(lm_rw0)[["dist"]]   # slope

# Convert slope to RW variance parameter
tau2_hat_day <- 2 * b_hat
tau2_hat <- tau2_hat_day * 365

# Final plot with fitted line
plot(gamma ~ dist, data = vg_t,
     xlab = "Time lag (days)", ylab = "Semivariance",
     main = "Temporal variogram with line anchored at 0")
abline(a = 0, b = b_hat, lwd = 2)

# ==============================================================================
# Spatiotemporal model
# ==============================================================================
tau2 <- tau2_hat
ell_km <- ell_km_hat

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
if (mut == "k13:675:V"){
  xlim = c(28, 37)
  ylim = c(-4, 5)
  # xlim = c(25, 40)
  # ylim = c(-7, 18)
} else if (mut == "k13:comb"){
  xlim = c(28, 48)
  ylim = c(-4.60, 18)
} else if (mut == "k13:622:I"){
  xlim = c(32, 45)
  ylim = c(4, 17)
} else if (mut == "k13:561:H"){
  xlim = c(28, 32)
  ylim = c(-4, 0)
}

lon_min <- xlim[1] # bb["xmin"]
lon_max <- xlim[2] # bb["xmax"]
lat_min <- ylim[1] # bb["ymin"]
lat_max <- ylim[2] # bb["ymax"]

# Grid
xs <- seq(lon_min, lon_max, length.out = nx)
ys <- seq(lat_min, lat_max, length.out = ny)
M  <- nx * ny

# Build lon/lat grid, then project to km for features
grid_ll <- expand.grid(lon = xs, lat = ys)
grid_sf <- st_as_sf(grid_ll, coords = c("lon", "lat"), crs = 4326)
grid_xy <- st_transform(grid_sf, crs_laea)
grid_m  <- st_coordinates(grid_xy)
grid_km <- grid_m / 1000

# Φ on projected grid (km)
OS_full   <- grid_km %*% t(Omega)
Phi_full  <- cbind(cos(OS_full), sin(OS_full)) * (1 / sqrt(D))

# Storage for means (and later, CIs) at t_vec
z_post_mean <- array(NA, dim = c(nx, ny, length(t_vec)))
p_post_mean <- array(NA, dim = c(nx, ny, length(t_vec)))

# uncertainty containers
z_post_sd   <- array(NA, dim = c(nx, ny, length(t_vec)))
p_post_lo   <- array(NA, dim = c(nx, ny, length(t_vec)))
p_post_hi   <- array(NA, dim = c(nx, ny, length(t_vec)))

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
C0_w <- Matrix(0, 2 * D, 2 * D, sparse = FALSE)  # prior covariance (zero variance, i.e. initial frequency known exactly)

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
    m_smooth[, ti] <- as.vector(m_filt[, ti] + J_t %*% (m_smooth[, ti + 1] - mnpt))
  }
  
  # Update EM parameters
  m_w_old <- m_w
  m_w     <- m_smooth
  num <- sqrt(sum((m_w - m_w_old)^2))
  den <- sqrt(sum(m_w_old^2) + 1e-12)
  rel_change <- num / den
  
  message(sprintf("  pseudo-Gaussian log-lik: %.3f,  rel_change(w): %.3e", ll_sum, rel_change))
  
  # ============================================================
  # Plot on-the-fly for checking
  # ============================================================
  
  if (FALSE) {
    
    for (k in seq_along(plot_times)) {
      t_idx <- which(t_vec == plot_times[k])
      
      # Mean z on grid
      z_mean_vec <- as.numeric(z_init + Phi_full %*% m_w[, t_idx])
      
      # Var[z] on grid: diag(Phi_full P_t Phi_full')
      # Efficient diagonal via row-wise quadratic form:
      #  V = rowSums( (Phi_full %*% P_t) * Phi_full )
      P_t <- C_smooth[[t_idx]]
      PhiP <- Phi_full %*% P_t               # [M x 2D]
      z_var_vec <- rowSums(PhiP * Phi_full)  # elementwise product, row-sum
      z_sd_vec  <- sqrt(pmax(z_var_vec, 0))
      
      # Transform to prevalence and 95% pointwise intervals
      p_mean_vec <- plogis(z_mean_vec)
      p_lo_vec   <- plogis(z_mean_vec - 1.96 * z_sd_vec)
      p_hi_vec   <- plogis(z_mean_vec + 1.96 * z_sd_vec)
      
      # Store as [nx x ny]
      z_post_mean[, , k] <- matrix(z_mean_vec, nrow = nx, ncol = ny, byrow = FALSE)
      p_post_mean[, , k] <- matrix(p_mean_vec, nrow = nx, ncol = ny, byrow = FALSE)
      z_post_sd[,   , k] <- matrix(z_sd_vec,   nrow = nx, ncol = ny, byrow = FALSE)
      p_post_lo[,   , k] <- matrix(p_lo_vec,   nrow = nx, ncol = ny, byrow = FALSE)
      p_post_hi[,   , k] <- matrix(p_hi_vec,   nrow = nx, ncol = ny, byrow = FALSE)
    }
    
    # Long df for raster plotting (lon/lat axes)
    p_long <- do.call(rbind, lapply(seq_along(plot_times), function(k) {
      data.frame(
        x = rep(xs, times = ny),
        y = rep(ys, each  = nx),
        p = as.vector(p_post_mean[, , k]),
        t = factor(plot_times[k], levels = plot_times)
      )
    }))
    
    # Base facet raster
    plot1 <- ggplot() + theme_bw() +
      geom_raster(aes(x = x, y = y, fill = p), data = p_long) +
      coord_fixed(expand = FALSE) +
      scale_fill_viridis_c(option = "magma", name = "Prevalence") +
      facet_wrap(~ t, nrow = 2) +
      labs(x = "Longitude", y = "Latitude") + 
      ggtitle(sprintf("EM step %s", i))
    
    print(plot1)
    
  }
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
# Exceedance probability surface: Pr{ p(s,t) > p_thresh }
# ============================================================
p_thresh_list <- c("1" = 0.01, "5" = 0.05, "10" = 0.10)

exceed_post_p_thresh <- array(
  NA_real_, dim = c(nx, ny, length(t_vec), length(p_thresh_list))
)

count <- 1
for (p_thresh in p_thresh_list) {
  
  exceed_post <- apply(p_draws, MARGIN = c(1, 2, 3), FUN = function(x) {
    mean(x > p_thresh)
  })
  exceed_post <- array(exceed_post, dim = c(nx, ny, length(t_vec)))
  
  exceed_post_p_thresh[, , , count] <- exceed_post
  count <- count + 1
}

# ============================================================
# Calculate mean, median and difference in CI from draws
# ============================================================
p_post_mean = apply(p_draws, MARGIN = c(1, 2, 3), FUN = mean)
p_post_median = apply(p_draws, MARGIN = c(1, 2, 3), FUN = median)
p_post_CI <- apply(p_draws, MARGIN = c(1, 2, 3), FUN = function(x) {
  diff(quantile(x, probs = c(0.025, 0.975)))
})

# ============================================================
# Map grid cells to admin units (unchanged)
# ============================================================
grid_sf <- st_as_sf(
  expand.grid(x = xs, y = ys),
  coords = c("x", "y"),
  crs = st_crs(admin1_regions)
) %>%
  mutate(pixel_index = dplyr::row_number())   # index BEFORE join / filtering

points_to_regions <- st_join(
  grid_sf,
  admin1_regions %>% dplyr::select(id_1, name_0, name_1)
) %>%
  st_drop_geometry() %>%
  drop_na(id_1)  

# ============================================================
# Time-series from 4D array: p_draws[nx, ny, t, draws]
# ============================================================

T <- t_num

region_draws <- vector("list", length(plot_times))

for (k in seq_along(plot_times)) {
  # index of this year in the 3rd dimension of p_draws
  t_idx <- which(t_vec == plot_times[k])
  if (length(t_idx) != 1L) stop("plot_times not aligned with t_vec")
  
  # slice: [nx, ny, draws] -> [M, B] with locations as rows, draws as columns
  arr_k   <- p_draws[ , , t_idx, ]              # [nx, ny, B]
  draws_k <- matrix(arr_k, nrow = M, ncol = num_post_draws)  # [M, B]
  
  # keep only pixels that fall inside admin1 regions
  draws_k <- draws_k[points_to_regions$pixel_index, , drop = FALSE]
  
  # make a data.frame with one column per draw
  df <- as.data.frame(draws_k)
  colnames(df) <- paste0("draw_", seq_len(num_post_draws))
  df$id_1 <- points_to_regions$id_1
  
  # long format + region-wise mean per draw
  df_long <- df %>%
    tidyr::pivot_longer(
      cols      = starts_with("draw_"),
      names_to  = "draw",
      values_to = "p"
    ) %>%
    dplyr::group_by(id_1, draw) %>%
    dplyr::summarise(
      mean_p = mean(p, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(t = plot_times[k])
  
  region_draws[[k]] <- df_long
}

df_posterior_timeseries <- dplyr::bind_rows(region_draws)

# ============================================================
# Create a list of results to save
# ============================================================
model_output <- list(
  p_post_mean = p_post_mean,
  p_post_median = p_post_median,
  p_post_CI = p_post_CI,

  time_series = df_posterior_timeseries,
  
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
  
  #Exceedance prob
  exceed_post_list = exceed_post_p_thresh,
  
  # Model settings
  settings = list(ell_km = ell_km, tau2 = tau2, D = D)
)

# Use your cache_path function to create a filename
output_filename <- cache_path(
  mut = mut, 
  what = "optimized_ell_tau2_full_model_output", 
  ext = "rds"
)

# Save the single list object to an .rds file
saveRDS(model_output, file = output_filename)

message(sprintf("Saved full model output for %s", mut)) 

mut_elapsed <- Sys.time() - mut_t0
message(sprintf("Elapsed for %s: %.2f min",
                mut, as.numeric(mut_elapsed, units = "mins")))

total_elapsed <- Sys.time() - t0
message(sprintf("Total elapsed: %.2f min",
                as.numeric(total_elapsed, units = "mins")))
