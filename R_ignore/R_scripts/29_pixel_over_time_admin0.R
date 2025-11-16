# ── Packages ────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(sf)
  library(ggplot2)
  library(units)
  library(rlang)  # for .data pronoun
  library(fs)
  library(devtools)
  library(purrr)
})

load_all()
sf_use_s2(FALSE)

# --- Paths & constants -----------------------------------------------------------
CACHE_DIR      <- "R_ignore/R_scripts/outputs/GRFF_kalman_cache_annual_2012_2025"
OUT_BASE       <- "GRFF_updated_maps"
OUT_PLOT_DIR <- "pixel_over_time_admin0"
dir_create(c(paste0("R_ignore/R_scripts/outputs/plots/", OUT_PLOT_DIR)))

# --- Cache helpers -----------------------------------------------------------
cache_path <- function(mut, lenS, lenT, what, ext = c("rds","parquet")) {
  ext <- match.arg(ext)
  fn  <- sprintf("%s_lenS%s_lenT%s_%s.%s",
                 gsub("[:/ ]", "_", mut), lenS, lenT, what, ext)
  file.path(CACHE_DIR, fn)
}

make_or_load <- function(mut, ell_km, tau2) {
  cached_rds <- cache_path(mut, ell_km, tau2, "full_model_output", "rds") 
  stopifnot(file_exists(cached_rds))
  readRDS(cached_rds)
}

make_raster_long <- function(arr3, xs, ys, times) {
  nx <- length(xs); ny <- length(ys)
  
  df <- do.call(rbind, lapply(seq_along(times), function(k) {
    data.frame(
      x = rep(xs, times = ny),
      y = rep(ys, each  = nx),
      p = as.vector(arr3[ , , k]),
      t = factor(times[k], levels = times)
    )
  }))
  
  df %>%
    dplyr::mutate(t = as.numeric(as.character(t))) %>%
    dplyr::filter(t %in% plot_times)
}

# --- Data I/O --------------------------------------------------------
shape_Africa <- readRDS("R_ignore/R_scripts/data/shapefiles/sf_admin0_africa.rds")
shape_water  <- sf::st_read("R_ignore/R_scripts/data/shapefiles/africa_water_bodies.shp", quiet = TRUE)

dat <- read.csv("R_ignore/R_scripts/data/all_who_get_prevalence.csv") |>
  mutate(
    collection_day  = as.Date(collection_day),
    collection_year = lubridate::year(collection_day)
  )

# --------------------------- Add K13 overall data ----------------------------
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

# --- Parameters ---------------------------------------------------------------
# choose a single run to (re)plot from cache
ell_km <- 80          # RFF length-scale in **kilometres**
tau2   <- 0.1         # RW1 variance in feature space
plot_times <- seq(2012, 2023, by = 1)

# prediction parameters
nx <- 200
ny <- 200

# Area parameters
q_thr <- 0.80      # exceedance-probability threshold
D     <- 1000       # resamples for ribbon

EA_long_lim <- c(15, 52)
EA_lat_lim  <- c(-12, 18)
x_range     <- c(-20, 55)  # lon extent used to subset points
y_range     <- c(-35, 38)  # not used below, but kept for reference
plot_crs = 4326

all_who_mutations <- c(
  "k13:comb","k13:622:I", "k13:469:Y", "k13:675:V", "k13:446:I", "k13:458:Y", 
  "k13:476:I",   "k13:493:H", "k13:539:T", "k13:543:T",  "k13:553:L", 
  "k13:561:H",   "k13:574:L", "k13:580:Y", "k13:441:L", "k13:449:A", 
  "k13:469:F",   "k13:481:V", "k13:515:K", "k13:527:H",  "k13:537:I", 
  "k13:537:D", "k13:538:V",  "k13:568:G")

area_ts_list <- list()
for (mut in all_who_mutations){
  print(paste0("Processing ", mut, "........."))
  # subset data for this mutation and AOI
  dat_sub <- dat_with_k13 |>
    filter(mutation == mut,
           longitude > x_range[1], longitude < x_range[2],
           year > 2011)
  
  # quick skip if no positive prevalence
  if (!any(dat_sub$prevalence > 0, na.rm = TRUE)) next
  
  # load cached model output (skip cleanly if missing)
  cached <- try(make_or_load(mut, ell_km, tau2), silent = TRUE)
  
  p_post_mean_draw <-  cached$p_post_mean_draw
  xs <- cached$xs
  ys <- cached$ys
  p_post_mean_long <- make_raster_long(p_post_mean_draw, xs, ys, plot_times) %>%
    rename(post_mean = p)
  
  # make points in WGS84 (change crs if you used something else)
  p_post_mean_draw_sf <- p_post_mean_long %>%
    sf::st_as_sf(coords = c("x", "y"), crs = 4326, remove = FALSE)
  
  # make sure shapes are in same CRS
  shape_Africa_4326 <- st_transform(shape_Africa, 4326)
  
  # keep only the admin0 name / id columns you need
  shp_admin0 <- shape_Africa_4326 %>%
    select(id_0, name_0)   # or whatever your column is called
  
  # spatial join: each pixel gets the country it falls in
  pts_with_country <- st_join( p_post_mean_draw_sf, shp_admin0, left = FALSE) %>%
    mutate(pixel_id = interaction(x, y, drop = TRUE))
  
  country_to_plot <- "Uganda"   # replace with your admin0_name
  
  df_plot <- pts_with_country %>%
    filter(name_0 == country_to_plot) %>%
    st_drop_geometry()
  
  ggplot(df_plot, aes(x = t, y = post_mean, group = pixel_id)) +
    geom_line(alpha = 0.1) +
    labs(
      title = paste("Posterior mean trajectories per pixel in", country_to_plot, "for", mut),
      x = "Year",
      y = "Posterior mean prevalence"
    ) +
    theme_bw()
  
}
