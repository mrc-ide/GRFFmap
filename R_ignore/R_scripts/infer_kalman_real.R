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

set.seed(1)

# --------------------------- Settings --------------------------------------
# model parameters
ell_km <- 80           # RFF length-scale in **kilometres**
tau2   <- 0.1        # RW1 variance in feature space
p_init <- 0.01
z_init <- qlogis(p_init)

# inference parameters
D        <- 200         # number of random frequencies (try 200–500)
max_iter <- 5           # EM iterations
z_eps    <- 1e-6        # threshold for |z| to use n/4 in E[ω]
omega_floor <- 1e-10    # floor on ω to avoid 1/ω explosions
jitter_S <- 1e-10       # diagonal jitter for innovation covariance

# prediction parameters
nx <- 120
ny <- 120
t_vec <- 1995:2030
t_num <- length(t_vec)

# plotting parameters
plot_times <- seq(1995, 2030, by = 1) # should be within t_vec
t_win <- 2           # window around each facet time (years). Points within this window are displayed

# --------------------------- Load & filter data ----------------------------
dat <- read.csv("R_ignore/R_scripts/data/all_who_get_prevalence_africa.csv")

dat_sub <- dat |>
  filter(mutation == "k13:622:I") |>
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
pad <- 1
lon_min <- min(dat_sub$longitude, na.rm = TRUE) - pad
lon_max <- max(dat_sub$longitude, na.rm = TRUE) + pad
lat_min <- min(dat_sub$latitude,  na.rm = TRUE) - pad
lat_max <- max(dat_sub$latitude,  na.rm = TRUE) + pad

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

# ============================================================
# Final plotting: lower / mean / upper prevalence surfaces
# ============================================================

# --- Build dataframe of observed points for overlay ------------------------
points_df <- tidyr::crossing(
  t = plot_times,
  dat_sub
) |>
  mutate(
    t = factor(t, levels = plot_times),
    in_window = abs(year - as.numeric(as.character(t))) <= t_win,
    p_obs = numerator / denominator,
    t = year(collection_day)
  )

# --- Helper to convert [nx x ny x K] array -> long data frame --------------
make_raster_long <- function(arr3, xs, ys, times) {
  do.call(rbind, lapply(seq_along(times), function(k) {
    data.frame(
      x = rep(xs, times = ny),
      y = rep(ys, each  = nx),
      p = as.vector(arr3[, , k]),
      t = factor(times[k], levels = times)
    )
  }))
}

# --- Build long data for lower / mean / upper ------------------------------
p_long_lo   <- make_raster_long(p_post_lo,   xs, ys, plot_times)
p_long_mean <- make_raster_long(p_post_mean, xs, ys, plot_times)
p_long_hi   <- make_raster_long(p_post_hi,   xs, ys, plot_times)

shape_Africa <- readRDS("R_ignore/R_scripts/data/shapefiles/sf_admin0_africa.rds")
shape_water  <- sf::st_read("R_ignore/R_scripts/data/shapefiles/africa_water_bodies.shp", quiet = TRUE)

# Use the same local LAEA you already defined for modeling, or a generic projected CRS (e.g., 3857)
crs_laea <- paste0(
  "+proj=laea +lat_0=", lat0,
  " +lon_0=", lon0,
  " +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
)

fix_geom_planar <- function(g, to_crs = crs_laea, snap_size = 10) {
  # 1) project to planar meters
  g_p <- st_transform(g, to_crs)
  # 2) make valid first (repairs bow ties/self-intersections)
  g_p <- st_make_valid(g_p)
  # 3) set precision in meters (grid size you want to keep)
  g_p <- st_set_precision(g_p, snap_size)
  # 4) snap to that grid to kill duplicate/near-duplicate vertices
  g_p <- st_snap_to_grid(g_p, snap_size)
  # 5) make valid again after snapping (paranoid but helpful)
  g_p <- st_make_valid(g_p)
  # 6) back to lon/lat
  st_transform(g_p, 4326)
}

# Example usage on your shapes (currently in 4326):
shape_Africa_fix <- fix_geom_planar(shape_Africa, crs_laea, snap_size = 10)  # 10 m grid
shape_water_fix  <- fix_geom_planar(shape_water,  crs_laea, snap_size = 10)

pp_cols <- c(
  "#eafff5",  # 0% (very light mint, not gray)
  "#c9fae7",  # ~0.2
  "#a6f0d7",  # ~0.5
  "#8fe3c7",  # ~1
  "#66CDAA",  # 2
  "#bfe6a1",  # 5
  "khaki2",   # 10
  "#f3c86b",  # 20
  "orange",   # 40
  "red"       # 100
)
pp_vals <- rescale(c(0, 0.2, 0.5, 1, 2, 5, 10, 20, 40, 100))  # emphasize <10%
pts_pos <- dat_sub |> 
  dplyr::filter(prevalence > 5) |>
  sf::st_as_sf(coords = c("longitude","latitude"), crs = 4326, remove = FALSE)

buf_km   <- 100
bbox_poly <- pts_pos |>
  sf::st_transform(3857) |>
  sf::st_bbox() |> sf::st_as_sfc() |>
  sf::st_buffer(buf_km * 1000) |>
  sf::st_transform(4326)

bb   <- sf::st_bbox(bbox_poly)
xlim <- c(bb["xmin"], bb["xmax"])
ylim <- c(bb["ymin"], bb["ymax"])

# --- Helper to generate a consistent ggplot with point overlays ------------
plot_layer <- function(p_long_df, title_text, shp, shp_water) {
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = p*100), data = p_long_df) +
    geom_sf(data = shp,       linewidth = 0.2, fill = NA, color = "white") +
    geom_sf(data = shp_water, fill = "white", colour = NA) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    # In-window observations (filled circles, viridis fill)
    geom_point(
      aes(x = longitude, y = latitude, fill = p_obs*100),
      data = points_df %>% filter(t>2013&t<2026) %>% arrange(p_obs),
      shape = 21, colour = "grey60", size = 2
    ) +
    scale_fill_gradientn(colours = pp_cols,
                         values  = pp_vals,
                         limits = c(0, 100), name = "Prevalence (%)") +
    facet_wrap(~ t, nrow = 2) +
    labs(x = "Longitude", y = "Latitude", title = title_text)
}

# --- Create and print the three separate plots -----------------------------
p_mean <- p_long_mean %>%
  mutate(t = as.numeric(as.character(t))) %>%
  filter(t > 2013&t<2026)
plot_lower <- plot_layer(p_long_lo,   "Lower (2.5%) prevalence", shape_Africa, shape_water)
plot_mean  <- plot_layer(p_mean, "Mean prevalence", shape_Africa, shape_water)
plot_upper <- plot_layer(p_long_hi,   "Upper (97.5%) prevalence", shape_Africa, shape_water)

print(plot_lower)
print(plot_mean)
print(plot_upper)
