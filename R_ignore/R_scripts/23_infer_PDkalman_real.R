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
library(here)
library(lwgeom) 

load_all()

set.seed(1)
sf_use_s2(FALSE)

# clock time
t0 <- Sys.time()

# build consistent filenames per (mutation, length_space, length_time, n_features)
cache_path <- function(mut, lenS, lenT, what, ext = c("rds","parquet","rds.gz")) {
  ext <- match.arg(ext)
  fn <- sprintf("%s_lenS%s_lenT%s_%s.%s",
                gsub("[:/ ]", "_", mut), lenS, lenT, what, ext)
  file.path(CACHE_DIR, fn)
}

# NA-robust row mean for 0/1 matrices
.row_mean01 <- function(X01) {
  Dk   <- rowSums(!is.na(X01))
  succ <- rowSums(X01, na.rm = TRUE)
  ifelse(Dk > 0L, succ / Dk, NA_real_)
}

# Per-time exceedance from draws: Pk is [M x D] (M = nx*ny)
exceed_prob_from_draws <- function(Pk, p_thr) {
  .row_mean01(Pk > p_thr)  # length M
}

# --------------------------- Load & filter data ----------------------------

CACHE_DIR <- "R_ignore/R_scripts/outputs/GRFF_kalman_partner_drug_cache_annual_1995_2024/"
dir_create(CACHE_DIR)

# read in prevalence data
dat <- read.csv("R_ignore/R_scripts/data/partner_drug_calc_get_prevalence.csv") |>
  mutate(collection_day = as.Date(collection_day)) |>
  select(longitude, latitude, year, numerator, denominator, prevalence, mutation)

# read Africa shape files
shape_Africa <- readRDS("R_ignore/R_scripts/data/shapefiles/sf_admin0_africa.rds")
shape_water  <- sf::st_read("R_ignore/R_scripts/data/shapefiles/africa_water_bodies.shp", quiet = TRUE)
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

# --------------------------- Fliter Africa shape ----------------------------
grey_countries <- c(
  "Egypt", "Morocco", "Libya", "Tunisia", "Algeria",
  "Cabo Verde", "Lesotho", "Mauritius", "Seychelles"
)

shape_malaria_africa <- shape_Africa %>%
  filter(!name_0 %in% grey_countries)

# --- Create white raster for outside Africa ------------------------------
bbox_info <- build_africa_bbox_and_crop(shape_malaria_africa, shape_water, buf_km = 200)
xlim        <- bbox_info$xlim
ylim        <- bbox_info$ylim

lon_min <- bbox_info$xlim[1]
lon_max <- bbox_info$xlim[2]
lat_min <- bbox_info$ylim[1]
lat_max <- bbox_info$ylim[2]

# --------------------------- Settings --------------------------------------
# inference parameters
D        <- 500         # number of random frequencies (try 200–500)
max_iter <- 5           # EM iterations
z_eps    <- 1e-6        # threshold for |z| to use n/4 in E[ω]
omega_floor <- 1e-10    # floor on ω to avoid 1/ω explosions
jitter_S <- 1e-10       # diagonal jitter for innovation covariance
ell_km <- 500           # RFF length-scale in **kilometres** #80 for k13
tau2   <- 0.1           # RW1 variance in feature space
p_init <- 0.001         
z_init <- qlogis(p_init)

# prediction parameters
nx <- 200
ny <- 200
t_vec <- 1995:2024
t_num <- length(t_vec)

# plotting parameters
plot_times <- seq(1995, 2024, by = 1) # should be within t_vec

# posterior draws
n_post_draws <- 1000

# --------------------------- Define PD mutants ----------------------------
# which mutation we are focusing on
all_who_mutations <- c("mdr1C:86:N", "crt:76:T")
all_who_mutations <- "crt:76:T"
for (mut in all_who_mutations){
  mut_t0 <- Sys.time()  # per-mutation clock stars
  
  # --------------------------- Run spatiotemporal ----------------------------
  dat_sub <- dat |>
    filter(mutation == mut) |>
    filter(is.finite(numerator), is.finite(denominator), denominator > 0)
  
  # --- Project to local Cartesian (km) ----------------------------------------
  proj <- project_to_local_laea(dat_sub,
                                lon_col = "longitude",
                                lat_col = "latitude",
                                crs_ll  = 4326)
  
  xy_km        <- proj$xy_km
  dat_sub$Xkm  <- xy_km[, 1]
  dat_sub$Ykm  <- xy_km[, 2]
  crs_laea     <- proj$crs_laea 
  lon0         <- proj$lon0 
  lat0         <- proj$lat0
  
  # --- Random Fourier Features ------------------------------------------------
  # Draw D frequencies ω_j ~ N(0, ℓ⁻² I₂), ℓ in **km**
  Omega <- matrix(rnorm(D * 2, mean = 0, sd = 1 / ell_km), ncol = 2)  # [D x 2]
  
  # Helper to build feature matrix Φ(S) for coords S = [x_km, y_km] in km
  phi_of_coords <- function(S_km) {
    OS <- as.matrix(S_km) %*% t(Omega)                 # (n x D)
    cbind(cos(OS), sin(OS)) * (1 / sqrt(D))            # (n x 2D)
  }
  
  # --------------------------- Create grid for prediction -----------------------
  # Grid
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
  p_post_mean_draw <- array(NA, dim = c(nx, ny, length(plot_times)))
  
  # uncertainty containers
  z_post_sd   <- array(NA, dim = c(nx, ny, length(plot_times)))
  p_post_lo   <- array(NA, dim = c(nx, ny, length(plot_times)))
  p_post_hi   <- array(NA, dim = c(nx, ny, length(plot_times)))
  p_post_lo_draw   <- array(NA, dim = c(nx, ny, length(plot_times)))
  p_post_hi_draw  <- array(NA, dim = c(nx, ny, length(plot_times)))
  
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
  }
  
  # ============================================================
  # Reconstruct z_t(s) and prevalence on lon/lat grid (for plotting)
  # ============================================================
  posterior_draws_list <- list()
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
    # p_mean_vec <- plogis(z_mean_vec)
    # p_lo_vec   <- plogis(z_mean_vec - 1.96 * z_sd_vec)
    # p_hi_vec   <- plogis(z_mean_vec + 1.96 * z_sd_vec)
    # 
    # Posterior draws
    message(sprintf("Time Series Draws: %d ", n_post_draws))
    L_t <- chol(P_t)
    w_draws <- matrix(m_w[, t_idx], ncol = 1) + L_t %*% matrix(rnorm(ncol(P_t) * n_post_draws), ncol = n_post_draws)
    z_draws <-as.matrix( z_init + Phi_full %*% w_draws)
    p_draws <- plogis(z_draws)
    posterior_draws_list[[k]] <- p_draws 
    
    # Calculate mean with CrI
    p_mean_exact <- rowMeans(p_draws)
    p_lo_exact   <- apply(p_draws, 1, quantile, 0.025, na.rm=TRUE)
    p_hi_exact   <- apply(p_draws, 1, quantile, 0.975, na.rm=TRUE)
    
    # Store as [nx x ny]
    z_post_mean[, , k] <- matrix(z_mean_vec, nrow = nx, ncol = ny, byrow = FALSE)
    #p_post_mean[, , k] <- matrix(p_mean_vec, nrow = nx, ncol = ny, byrow = FALSE)
    z_post_sd[,   , k] <- matrix(z_sd_vec,   nrow = nx, ncol = ny, byrow = FALSE)
    #p_post_lo[,   , k] <- matrix(p_lo_vec,   nrow = nx, ncol = ny, byrow = FALSE)
    #p_post_hi[,   , k] <- matrix(p_hi_vec,   nrow = nx, ncol = ny, byrow = FALSE)
    p_post_mean_draw[, , k] <- matrix(p_mean_exact, nx, ny, byrow = FALSE)
    p_post_lo_draw[,   , k] <- matrix(p_lo_exact,   nx, ny, byrow = FALSE)
    p_post_hi_draw[,   , k] <- matrix(p_hi_exact,   nx, ny, byrow = FALSE)
  }
  
  # ============================================================
  # Time-series
  # ============================================================
  # Map grid to admin units
  grid_sf <- st_as_sf(expand.grid(x = xs, y = ys), coords = c("x", "y"), crs = st_crs(admin1_regions))
  points_to_regions <- st_join(grid_sf, admin1_regions %>% dplyr::select(id_1, name_0, name_1)) %>%
    st_drop_geometry() %>%
    mutate(pixel_index = 1:n()) %>%
    drop_na(id_1)
  
  # Combine posterior draws across time
  region_draws <- list()
  for (k in seq_along(plot_times)) {
    draws_k <- posterior_draws_list[[k]][points_to_regions$pixel_index, ]  # [pixels x draws]
    
    df <- as.data.frame(draws_k)
    df$id_1 <- points_to_regions$id_1
    
    df_long <- df %>%
      pivot_longer(cols = starts_with("V"), names_to = "draw", values_to = "p") %>%
      group_by(id_1, draw) %>%
      summarise(mean_p = mean(p, na.rm = TRUE), .groups = "drop") %>%
      mutate(t = plot_times[k])
    
    region_draws[[k]] <- df_long
  }
  
  df_posterior_timeseries <- bind_rows(region_draws)
  
  # ============================================================
  # Exceedance probability surface: Pr{ p(s,t) > p_thresh }
  # ============================================================
  # ---- Choose threshold(s) in prevalence units ---
  p_thresh_list <- c("50" = 0.5, "20" = 0.2 ,"10" = 0.10, "5" = 0.05, "1" = 0.01)
  
  # Calculate area km2 per cell
  lat_mid <- (lat_min + lat_max) / 2
  # 1 degree lat ≈ 111.32 km and 1 degree lon ≈ 111.32 * cos(latitude)
  km_per_deg_lat <- 111.32
  km_per_deg_lon <- 111.32 * cos(lat_mid * pi / 180)
  dx_deg <- (lon_max - lon_min) / (nx - 1)  # cell width in degrees lon
  dy_deg <- (lat_max - lat_min) / (ny - 1)  # cell height in degrees lat
  dx_km <- dx_deg * km_per_deg_lon
  dy_km <- dy_deg * km_per_deg_lat
  cell_area_km2 <- dx_km * dy_km
  
  # initalize list for exceedance prob based on Gaussian approx
  # exceed_post_gaussian_approx_list <- vector("list", length(p_thresh_list))
  # names(exceed_post_gaussian_approx_list) <- names(p_thresh_list)
  
  # initalize list for exceedance prob based on draws
  exceed_post_draws_list <- vector("list", length(p_thresh_list))
  names(exceed_post_draws_list) <- names(p_thresh_list)
  
  count = 1
  for (p_thresh in p_thresh_list){
    z_thresh <- qlogis(p_thresh)
    
    # ---- Compute exceedance on the lon/lat grid for plot_times ----
    #exceed_post_gaussian_approx <- array(NA, dim = c(nx, ny, length(plot_times)))
    exceed_post_draw <- array(NA_real_, dim = c(nx, ny, length(plot_times)))
    
    for (k in seq_along(plot_times)) {
      # ============================================================
      # Exceedance probability: via Gaussian approximation
      # ============================================================
      #t_idx <- which(t_vec == plot_times[k])
      
      # Mean z on grid from smoothed feature weights
      #z_mean_vec <- as.numeric(z_init + Phi_full %*% m_w[, t_idx])
      
      # Var[z] on grid: diag(Phi_full %*% P_t %*% t(Phi_full))
      # P_t   <- C_smooth[[t_idx]]                 # smoothed state covariance at time t
      # PhiP  <- Phi_full %*% P_t                  # [M x 2D]
      # z_var_vec <- rowSums(PhiP * Phi_full)      # fast diagonal
      # z_var_vec <- pmax(z_var_vec, 0)            # guard tiny negatives
      # z_sd_vec  <- sqrt(z_var_vec)
      # 
      # Exceedance under Gaussian approx on logit scale:
      # Pr{p > p_thresh} = Pr{z > z_thresh} = 1 - Phi((z_thresh - mean)/sd)
      # z_std <- ifelse(z_sd_vec > 0, (z_thresh - z_mean_vec) / z_sd_vec, NA_real_)
      # exc_vec <- ifelse(
      #   is.na(z_std),
      #   as.numeric(z_mean_vec > z_thresh),       # if sd==0, degenerate
      #   1 - pnorm(z_std)
      # )
      
      #exceed_post_gaussian_approx[, , k] <- matrix(exc_vec, nrow = nx, ncol = ny, byrow = FALSE)
      
      # ============================================================
      # Exceedance probability: via draws
      # ============================================================
      Pk <- posterior_draws_list[[k]]            # [M x D] of prevalences from draws
      theta_hat <- exceed_prob_from_draws(Pk, p_thresh)   # length M
      exceed_post_draw[ , , k] <- matrix(theta_hat, nrow = nx, ncol = ny, byrow = FALSE)
    }
    
    #exceed_post_gaussian_approx_list[[count]] = exceed_post_gaussian_approx
    exceed_post_draws_list[[count]] = exceed_post_draw
    count = count + 1
  }
  
  # Create a list of results to save
  model_output <- list(
    # Essential predictions
    #p_mean = p_post_mean,
    #p_lower = p_post_lo,
    #p_upper = p_post_hi,
    
    p_post_mean_draw = p_post_mean_draw,
    p_post_lower_draw = p_post_lo_draw,
    p_post_upper_draw = p_post_hi_draw,
    
    # Grid information
    xs = xs,
    ys = ys,
    lon_min = lon_min,
    lon_max = lon_max,
    lat_min = lat_min,
    lat_max = lat_max,
    plot_times = plot_times,
    
    # Original data
    data_subset = dat_sub,
    
    # Core model components (optional, for flexibility)
    weights_mean = m_w,
    weights_cov = C_smooth,
    random_frequencies = Omega,
    
    # Posterior draws
    num_posterior_draws = n_post_draws,
    posterior_draws = posterior_draws_list,
    
    # Time-series
    time_series = df_posterior_timeseries,
    
    #Exceedance prob
    # exceed_post_gaussian_approx_list = exceed_post_gaussian_approx_list,
    exceed_post_draw_list = exceed_post_draws_list,
    # 
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
  
  mut_elapsed <- Sys.time() - mut_t0
  message(sprintf("Elapsed for %s: %.2f min",
                  mut, as.numeric(mut_elapsed, units = "mins")))
}

total_elapsed <- Sys.time() - t0
message(sprintf("Total elapsed: %.2f min",
                as.numeric(total_elapsed, units = "mins")))
