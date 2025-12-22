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
# Distance calculation for variogram analysis 
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
# unique(vg_st_time$dist)

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

message("Temporal variogram calculated.")

# Saving outputs
saveRDS(vg_1y_sp, file.path(CACHE_DIR, paste0(mut, "_variogram_distance_space.rds")))
saveRDS(vg_t, file.path(CACHE_DIR, paste0(mut, "_variogram_distance_time.rds")))

total_elapsed <- Sys.time() - mut_t0
message(sprintf("Total elapsed: %.2f min",
                as.numeric(total_elapsed, units = "mins")))

# ==============================================================================
# Variogram analysis to identify optimal length space parameter (ell_km) 
# and temporal variance (tau2)
# ==============================================================================

# --- Spatial variogram to estimate the spatial length-scale ell_km ------------
# vg_1y_sp <- readRDS(file.path(VARIOGRAM_DIR, paste0(mut, "_variogram_distance_space.rds")))

df <- subset(vg_1y_sp, is.finite(dist) & is.finite(gamma) & dist > 0)
h  <- df$dist
g  <- df$gamma
w <- sqrt(df$np)

# nugget_fix <- 0.1

# gamma_gauss_fixnug <- function(h, psill, range) {
#   nugget_fix + psill * (1 - exp(-(h / range)^2))
# }

# obj_wls_fixnug <- function(par_log) {
#   psill <- exp(par_log[1])
#   range <- exp(par_log[2])
#   pred  <- gamma_gauss_fixnug(h, psill, range)
#   sum(w * (g - pred)^2)
# }

# # starts (partial sill should exclude nugget)
# psill0 <- max(1e-8, median(g, na.rm = TRUE) - nugget_fix)
# range0 <- median(h, na.rm = TRUE)

# fit <- optim(
#   par    = log(c(psill0, range0)),
#   fn     = obj_wls_fixnug,
#   method = "Nelder-Mead"
# )

# psill_hat  <- exp(fit$par[1])
# range_hat  <- exp(fit$par[2])
# nugget_hat <- nugget_fix

# ell_km_hat <- range_hat / sqrt(2)

# hgrid <- seq(0, max(h), length.out = 500)
# pred  <- nugget_fix + psill_hat * (1 - exp(-(hgrid / range_hat)^2))

gamma_gauss_nonug <- function(h, psill, range) {
  psill * (1 - exp(-(h / range)^2))   # no nugget term
}

obj_wls_nonug <- function(par_log) {
  psill <- exp(par_log[1])
  range <- exp(par_log[2])
  pred  <- gamma_gauss_nonug(h, psill, range)
  sum(w * (g - pred)^2)
}

# starts
psill0 <- max(1e-8, median(g, na.rm=TRUE))   # no need to subtract nugget
range0 <- median(h, na.rm=TRUE)

fit <- optim(log(c(psill0, range0)), obj_wls_nonug, method="Nelder-Mead")

psill_hat <- exp(fit$par[1])
range_hat <- exp(fit$par[2])
ell_km_hat <- range_hat / sqrt(2)

hgrid <- seq(0, max(h), length.out = 500)
pred  <- gamma_gauss_nonug(hgrid, psill_hat, range_hat)

png(paste0(CACHE_DIR, "/", mut, "_variogram_1year_spatial.png"), width = 6, height = 5, units = "in", res = 300)
plot(g ~ h, xlab="Distance (km)", ylab="Semivariance", pch=1)
lines(hgrid, pred, lwd=2)
dev.off()

message("Estimated length-scale ell_km = ", ell_km_hat)

# --- Temporal variogram to estimate RW1 variance (tau^2) ----------------------

# vg_t <- readRDS(file.path(VARIOGRAM_DIR, paste0(mut, "_variogram_distance_time.rds")))

# Estimate RW variance via weighted linear regression (model: gamma(h) = 0.5 * tau^2 * h)
lm_rw0 <- lm(gamma ~ 0 + dist, data = vg_t, weights = np)

# Estimated slope
b_hat <- coef(lm_rw0)[["dist"]]   # slope

# Convert slope to RW variance parameter
tau2_hat_day <- 2 * b_hat
tau2_hat <- tau2_hat_day * 365

# Final plot with fitted line
png(paste0(CACHE_DIR, "/", mut, "_temporal_variogram_rw_fit.png"), width = 6, height = 5, units = "in", res = 300)
plot(
  gamma ~ dist,
  data = vg_t,
  xlab = "Time lag (days)",
  ylab = "Semivariance",
  main = "Temporal variogram with RW fit (intercept fixed at 0)"
)
# Fitted RW line: gamma(h) = b_hat * h
abline(a = 0, b = b_hat, lwd = 2)
dev.off()

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
saveRDS(
  hyperparams,
  file.path(VARIOGRAM_DIR, paste0(mut, "_variogram_hyperparams.rds"))
)