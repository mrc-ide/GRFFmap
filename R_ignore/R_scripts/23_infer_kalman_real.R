# ============================================================
# RFF + PG-EM with Kalman filter/RTS in projected (Cartesian) coords
# ============================================================
# - Reads WHO GET prevalence data
# - Filters to a mutation and spatial window
# - Projects lon/lat -> local LAEA (km), runs RFF + PG-EM in km
# - Reconstructs prevalence on a lon/lat grid for plotting
# - Overlays observations with an in-window/out-of-window styling
# ============================================================

library(tidyverse)
library(Matrix)
library(here)
library(sf)
library(devtools)
library(dplyr)
library(fs)

load_all()

set.seed(1)

# build consistent filenames per (mutation, length_space, length_time, n_features)
cache_path <- function(mut, lenS, lenT, what, ext = c("rds","parquet","rds.gz")) {
  ext <- match.arg(ext)
  fn <- sprintf("%s_lenS%s_lenT%s_%s.%s",
                gsub("[:/ ]", "_", mut), lenS, lenT, what, ext)
  file.path(CACHE_DIR, fn)
}

# --------------------------- Settings --------------------------------------
# model parameters
ell_km <- 80          # RFF length-scale in **kilometres**
tau2   <- 0.1         # RW1 variance in feature space
p_init <- 0.01
z_init <- qlogis(p_init)

# inference parameters
D        <- 300         # number of random frequencies (try 200–500)
max_iter <- 5           # EM iterations
z_eps    <- 1e-6        # threshold for |z| to use n/4 in E[ω]
omega_floor <- 1e-10    # floor on ω to avoid 1/ω explosions
jitter_S <- 1e-10       # diagonal jitter for innovation covariance

# prediction parameters
nx <- 400
ny <- 400
t_vec <- 2014:2025
t_num <- length(t_vec)

# plotting parameters
plot_times <- seq(2014, 2025, by = 1) # should be within t_vec

# --------------------------- Load & filter data ----------------------------

CACHE_DIR <- "R_ignore/R_scripts/outputs/GRFF_kalman_cache_annual"
dir_create(CACHE_DIR)

# read in prevalence data
dat <- read.csv("R_ignore/R_scripts/data/all_who_get_prevalence_africa.csv") |>
  mutate(collection_day = as.Date(collection_day),
         year_group = year(collection_day))

# which mutation we are focusing on
all_who_mutations <- c("k13:675:V", "k13:446:I", "k13:458:Y", "k13:476:I",   "k13:493:H",   "k13:539:T",
                       "k13:543:T",  "k13:553:L",   "k13:561:H",   "k13:574:L",  "k13:580:Y",
                       "k13:441:L", "k13:449:A",   "k13:469:F",   "k13:481:V",
                       "k13:515:K", "k13:527:H",  "k13:537:I", "k13:537:D", "k13:538:V",  "k13:568:G")

# --------------------------- Run spatiotemporal ----------------------------
for (mut in all_who_mutations){
  dat_sub <- dat |>
    filter(mutation == mut) |>
    filter((latitude  > -5) & (latitude  < 20)) |>
    filter((longitude > 25) & (longitude < 40)) |>
    filter(is.finite(numerator), is.finite(denominator), denominator > 0)
  
  # --------------------------- Project to local Cartesian (km) ---------------
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
  
  # --------------------------- Random Fourier Features -----------------------
  # Draw D frequencies ω_j ~ N(0, ℓ⁻² I₂), ℓ in **km**
  Omega <- matrix(rnorm(D * 2, mean = 0, sd = 1 / ell_km), ncol = 2)  # [D x 2]
  
  # Helper to build feature matrix Φ(S) for coords S = [x_km, y_km] in km
  phi_of_coords <- function(S_km) {
    OS <- as.matrix(S_km) %*% t(Omega)                 # (n x D)
    cbind(cos(OS), sin(OS)) * (1 / sqrt(D))            # (n x 2D)
  }
  
  # --------------------------- Create grid for prediction -----------------------
  
  # Bounding box in lon/lat (pad so raster isn't clipped)
  # pad <- 1
  # lon_min <- min(dat_sub$longitude, na.rm = TRUE) - pad
  # lon_max <- max(dat_sub$longitude, na.rm = TRUE) + pad
  # lat_min <- min(dat_sub$latitude,  na.rm = TRUE) - pad
  # lat_max <- max(dat_sub$latitude,  na.rm = TRUE) + pad
  lon_min <- -20
  lon_max <- 55
  lat_min <- -35
  lat_max <- 38
  
  # Grid
  # xs <- seq(lon_min, lon_max, length.out = nx)
  # ys <- seq(lat_min, lat_max, length.out = ny)
  xs <- seq(lon_min, lon_max, length.out = nx)
  ys <- seq(lat_min, lat_max, length.out = ny)
  M <- nx * ny
  
  # Build lon/lat grid, then project to km for features
  grid_ll <- expand.grid(lon = xs, lat = ys)
  grid_sf <- st_as_sf(grid_ll, coords = c("lon", "lat"), crs = 4326)
  grid_xy <- st_transform(grid_sf, crs_laea)
  grid_m  <- st_coordinates(grid_xy)
  grid_km <- grid_m / 1000
  
  # Φ on projected grid (km)
  OS_full   <- grid_km %*% t(Omega) # Dot product of a location's coordinates and a random frequency vector (output: M x D matrix)
  Phi_full  <- cbind(cos(OS_full), sin(OS_full)) * (1 / sqrt(D)) # Feature Matrix: approximate a covariance function
  
  # Storage for means (and later, CIs) at plot_times
  z_post_mean <- array(NA, dim = c(nx, ny, length(plot_times)))
  p_post_mean <- array(NA, dim = c(nx, ny, length(plot_times)))
  
  # uncertainty containers
  z_post_sd   <- array(NA, dim = c(nx, ny, length(plot_times)))
  p_post_lo   <- array(NA, dim = c(nx, ny, length(plot_times)))
  p_post_hi   <- array(NA, dim = c(nx, ny, length(plot_times)))
  
  # --------------------------- Organize by time ------------------------------
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
    #
    Omega_list <- vector("list", t_num)
    r_list     <- vector("list", t_num)
    for (ti in 1:t_num) {
      ob <- obs_by_t[[ti]]
      if (!is.null(ob)) {
        Phi_t  <- ob$Phi
        mu_t   <- ob$mu
        n_t    <- ob$n
        kappa  <- ob$kappa
        
        # Current linear predictor z_t(s) = μ + Φ w_t (current prev map)
        # Phi_t is the weight at time t
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
    # Reconstruct z_t(s) and prevalence on lon/lat grid (for plotting)
    # ============================================================
    
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
  }
  # ============================================================
  # Exceedance probability surface: Pr{ p(s,t) > p_thresh }
  # (Drop-in block; no changes to earlier code required)
  # ============================================================
    
  # ---- Choose threshold(s) in prevalence units ---
  p_thresh_list <- c(0.10, 0.05, 0.01)
  exceed_post_list = list()
  count = 1
  for (p_thresh in p_thresh_list){
    z_thresh <- qlogis(p_thresh)
    
    # ---- Compute exceedance on the lon/lat grid for plot_times ----
    exceed_post <- array(NA, dim = c(nx, ny, length(plot_times)))
    
    for (k in seq_along(plot_times)) {
      t_idx <- which(t_vec == plot_times[k])
      
      # Mean z on grid from smoothed feature weights
      z_mean_vec <- as.numeric(z_init + Phi_full %*% m_w[, t_idx])
      
      # Var[z] on grid: diag(Phi_full %*% P_t %*% t(Phi_full))
      P_t   <- C_smooth[[t_idx]]                 # smoothed state covariance at time t
      PhiP  <- Phi_full %*% P_t                  # [M x 2D]
      z_var_vec <- rowSums(PhiP * Phi_full)      # fast diagonal
      z_var_vec <- pmax(z_var_vec, 0)            # guard tiny negatives
      z_sd_vec  <- sqrt(z_var_vec)
      
      # Exceedance under Gaussian approx on logit scale:
      # Pr{p > p_thresh} = Pr{z > z_thresh} = 1 - Phi((z_thresh - mean)/sd)
      z_std <- ifelse(z_sd_vec > 0, (z_thresh - z_mean_vec) / z_sd_vec, NA_real_)
      exc_vec <- ifelse(
        is.na(z_std),
        as.numeric(z_mean_vec > z_thresh),       # if sd==0, degenerate
        1 - pnorm(z_std)
      )
      
      exceed_post[, , k] <- matrix(exc_vec, nrow = nx, ncol = ny, byrow = FALSE)
    }
    exceed_post_list[[count]] = exceed_post
    count = count + 1
  }
  
  # Create a list of results to save
  model_output <- list(
    # Essential predictions
    p_mean = p_post_mean,
    p_lower = p_post_lo,
    p_upper = p_post_hi,
    
    # Grid information
    xs = xs,
    ys = ys,
    plot_times = plot_times,
    
    # Original data
    data_subset = dat_sub,
    
    # Core model components (optional, for flexibility)
    weights_mean = m_w,
    weights_cov = C_smooth,
    random_frequencies = Omega,
    
    #Exceedance prob
    exceed_prob_10_perc = exceed_post_list[[1]],
    exceed_prob_5_perc = exceed_post_list[[2]],
    exceed_prob_1_perc = exceed_post_list[[3]],
    
    # Model settings
    settings = list(ell_km = ell_km, tau2 = tau2, D = D)
  )
  
  # Use your cache_path function to create a filename
  output_filename <- cache_path(
    mut = mut, 
    lenS = ell_km, 
    lenT = tau2, 
    what = "full_model_output", 
    ext = "rds"
  )
  
  # Save the single list object to an .rds file
  saveRDS(model_output, file = output_filename)
  
  message(sprintf("Saved full model output for %s", mut))
}
