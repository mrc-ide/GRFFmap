# ============================================================
# Calculating distance and semivariance
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

# --- Define K13 mutants -------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
mut <- args[1]
mut <- "k13:comb"
# "crt:76:T"
# "mdr1:86:Y"

# Per-mutation clock start
mut_t0 <- Sys.time()

# --------------------------- Load & filter data ----------------------------
CACHE_DIR <- "R_ignore/R_scripts/outputs/model_outputs/variogram_distances"
dir_create(CACHE_DIR)

VARIOGRAM_DIR <- "R_ignore/R_scripts/outputs/model_outputs/variogram_distances"

# Read in prevalence data
dat <- read.csv("R_ignore/R_scripts/data/all_mutations_get_prevalence.csv") |>
  mutate(collection_day = as.Date(collection_day)) |>
  select(survey_id, study_id, collection_day, longitude, latitude, year, country_name, numerator, denominator, prevalence, mutation) 

# --- Add K13 overall data -----------------------------------------------------
dat_with_k13 <- add_combined_k13(dat)

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
# Spatiotemporal analysis
# ==============================================================================

# --- Spatial variogram to estimate the spatial length-scale ell_km ------------
# Construct a spatiotemporal object (STIDF) from the data subset
# Spatial coordinates in kilometres (projected coordinates)
sp_pts <- sp::SpatialPoints(dat_sub[, c("Xkm", "Ykm")])

# Create time index where time is treated as numeric across all years in dat_sub
if (mut == "k13:comb"){
  time_idx <- as.POSIXct(paste0(dat_sub$year, "-07-01"), tz = "UTC")
} else {
  time_idx <- as.POSIXct(dat_sub$collection_day, tz = "UTC")
}
# Combine spatial points + time index + observed values into a spacetime object
st_obj <- spacetime::STIDF(sp_pts, time_idx, data = data.frame(z = dat_sub$z))

# Compute spatiotemporal variogram, then isolate the ~1-year temporal slice
# Small tolerance around 365 days to catch "approximately one year apart" pairs
eps_days <- 0.5

# Compute spatiotemporal variogram
coords <- as.matrix(dat_sub[, c("Xkm","Ykm")])
dmax <- max(sp::spDists(coords, longlat = FALSE), na.rm = TRUE)

vg_st_1year <- variogramST(
  z ~ 1,
  data   = st_obj,
  tlags  = c(365 - eps_days, 365 + eps_days, 3 * 365),  # boundaries for ≤2 years
  tunit  = "days",
  cutoff = dmax,
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
vg_1y_sp <- subset(vg_1y_sp, is.finite(dist) & is.finite(gamma) & is.finite(np) & np > 0)
class(vg_1y_sp) <- c("gstatVariogram","data.frame")

message("Spatial variogram (1-year lag) calculated.")

# Saving outputs
saveRDS(vg_1y_sp, file.path(CACHE_DIR, paste0(mut, "_variogram_distance_space.rds")))

total_elapsed <- Sys.time() - mut_t0
message(sprintf("Total elapsed: %.2f min",
                as.numeric(total_elapsed, units = "mins")))

# --- Spatial variogram to estimate the spatial length-scale ell_km ------------
# vg_1y_sp <- readRDS(file.path(VARIOGRAM_DIR, paste0(mut, "_variogram_distance_space.rds")))

df <- subset(vg_1y_sp, is.finite(dist) & is.finite(gamma) & dist > 0)
h  <- df$dist
g  <- df$gamma
w <- sqrt(df$np)

nugget_fix <- 0.8

gamma_gauss_nug <- function(h, psill, range) {
  nugget_fix + psill * (1 - exp(-(h / range)^2))
}

obj_wls_nug <- function(par_log) {
  psill <- exp(par_log[1])
  range <- exp(par_log[2])
  pred  <- gamma_gauss_nug(h, psill, range)
  sum(w * (g - pred)^2)
}

psill0 <- max(1e-8, median(g, na.rm = TRUE) - nugget_fix)
range0 <- median(h, na.rm = TRUE)

fit <- optim(
  par    = log(c(psill0, range0)),
  fn     = obj_wls_nug,
  method = "Nelder-Mead"
)

psill_hat  <- exp(fit$par[1])
range_hat  <- exp(fit$par[2])
nugget_hat <- nugget_fix
ell_km_hat <- range_hat / sqrt(2)

hgrid <- seq(0, max(h), length.out = 500)
pred  <- gamma_gauss_nug(hgrid, psill_hat, range_hat)

message("Estimated nugget = ", nugget_hat)
message("Estimated length-scale ell_km = ", ell_km_hat)

# Collect spatial variogram plot data
vg_df <- data.frame(
  dist   = hgrid,
  fitted = pred
)
obs_df <- data.frame(
  dist  = h,
  gamma = g,
  np    = df$np
)
spatial_vg_plot_data <- list(
  obs_df  = obs_df,    # observed semivariance points
  vg_df   = vg_df,     # fitted curve
  nugget  = nugget_hat,
  ell_km  = ell_km_hat
)

# ==============================================================================
# Temporal Variogram analysis
# ==============================================================================

# Build all 1-year-apart pairs
pairs_1yr <- dat_sub |>
  mutate(id = row_number()) |>
  inner_join(
    dat_sub |> mutate(id2 = row_number()),
    by = character()
  ) |>
  filter(year.y == year.x + 1) |>
  mutate(
    dx   = Xkm.y - Xkm.x,
    dy   = Ykm.y - Ykm.x,
    dist = sqrt(dx^2 + dy^2),
    dz   = z.y - z.x
  )

# Restrict to nearby pairs
dist_max <- 75   # km; tune this

pairs_local <- pairs_1yr |>
  filter(dist <= dist_max)

# Estimate tau^2
tau2_hat <- mean(pairs_local$dz^2, na.rm = TRUE)
message("Estimated RW1 variance tau2 = ", tau2_hat)

# Collect fitted hyperparameters
hyperparams <- list(
  mutation      = mut,
  ell_km        = ell_km_hat,
  tau2_year     = tau2_hat,
  tau2_day      = tau2_hat_day,
  psill         = psill_hat,
  range_km      = range_hat,
  fitted_at     = Sys.time()
)

# Save as RDS
saveRDS(hyperparams, file.path(VARIOGRAM_DIR, paste0(mut, "_variogram_hyperparams.rds")))
saveRDS(spatial_vg_plot_data,  file.path(VARIOGRAM_DIR, paste0(mut, "_spatial_vg_plot_data.rds")))
