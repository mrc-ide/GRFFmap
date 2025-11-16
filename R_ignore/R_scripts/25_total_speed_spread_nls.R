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
})

load_all()

# --- Paths & constants -----------------------------------------------------------
CACHE_DIR      <- "R_ignore/R_scripts/outputs/GRFF_kalman_cache_annual_2012_2025"
OUT_BASE       <- "GRFF_updated_maps"
OUT_SPEED   <- file.path("speed_estimate")
OUT_PLOT_DIR <- "total_speed_spread_nls"
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
  
  exceed_prob_draws <- cached$exceed_post_draw_list
  exceed_prob_draws_1 <- exceed_prob_draws$`1`
  exceed_prob_draws_5 <- exceed_prob_draws$`5`
  exceed_prob_draws_10 <- exceed_prob_draws$`10`
  
  posterior_draws = cached$posterior_draws
  
  exceed_prob_draws_long_10 <- make_raster_long(exceed_prob_draws_10, xs, ys, plot_times)
  exceed_prob_draws_long_5 <- make_raster_long(exceed_prob_draws_5, xs, ys, plot_times)
  exceed_prob_draws_long_1 <- make_raster_long(exceed_prob_draws_1, xs, ys, plot_times)
  
  lat_min <- cached$lat_min
  lat_max <- cached$lat_max
  lon_min <- cached$lon_min
  lon_max <- cached$lon_max
  
  # Combine into a single data frame
  df_pixel_exceed_prob <- exceed_prob_draws_long_1 %>%
    rename(P_ex_p1 = p) %>%
    left_join(exceed_prob_draws_long_5 %>% rename(P_ex_p5 = p),
              by = c("x", "y", "t")) %>%
    left_join(exceed_prob_draws_long_10 %>% rename(P_ex_p10 = p),
              by = c("x", "y", "t"))
  
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
  
  #--- Area calculation based on draws --------------------------------
  df_area_over_time_exceed <- df_pixel_exceed_prob %>%
    # Identify which pixels are "high-risk" at each time point
    mutate(p1 = P_ex_p1 > q_thr,
           p5 = P_ex_p5 > q_thr,
           p10 = P_ex_p10 > q_thr) %>%
    group_by(t) %>%
    # Calculate the total high-risk area for that time point
    summarise(
      # Count the number of high-risk pixels
      n_p1 = sum(p1),
      n_p5 = sum(p5),
      n_p10 = sum(p10),
      # Calculate the total area
      total_area_p1 = n_p1 * cell_area_km2,
      total_area_p5 = n_p5 * cell_area_km2,
      total_area_p10 = n_p10 * cell_area_km2,
      .groups = "drop"
    )

  # Fit Nonlinear Least Squares (NLS) regression
  model_fits <- list()
  for (p_thr in c(0.01, 0.05, 0.10)) {
    # Map threshold to column name
    p_thr_col <- switch(as.character(p_thr),
                        "0.01" = "total_area_p1",
                        "0.05" = "total_area_p5",
                        "0.1"  = "total_area_p10")
    
    # y vector for starting values
    y <- df_area_over_time_exceed[[p_thr_col]]
    A_max <- max(y, na.rm = TRUE)
    t_centered <- df_area_over_time_exceed$t - min(df_area_over_time_exceed$t, na.rm = TRUE)
    df_tmp <- data.frame(t = df_area_over_time_exceed$t, t_centered = t_centered, y = y)
    start_list <- list(
      a  = sqrt(A_max / pi) * 1.05,  # radius s.t. pi*a^2 ≈ max area
      b  = 0.1,
      t0 = 0
    )
    
    # Build formula: total_area_pX ~ (pi * a^2) * (1 - exp(-b * (t - t0)))^2
    fit_formula <- as.formula(
      "y ~ pi * a^2 * (1 - exp(-b * (t_centered - t0)))^2"
    )
    
    fit_model <- nls(
      fit_formula,
      data    = df_tmp,
      start   = start_list,
      algorithm = "port",
      lower   = c(a = 0,   b = 0,   t0 = -5),
      upper   = c(a = Inf, b = 10,  t0 = max(t_centered) + 5),
      control = nls.control(maxiter = 500, minFactor = 1e-10, warnOnly = TRUE)
    )
    
    model_fits[[p_thr_col]] <- fit_model
  }
  
  # Generate prediction points
  df_area_over_time_exceed$t_centered <- df_area_over_time_exceed$t - min(df_area_over_time_exceed$t)
  time_range <- seq(min(df_area_over_time_exceed$t_centered), max(df_area_over_time_exceed$t_centered), length.out = 100)
  predictions_exceed <- data.frame(t_centered = time_range)
  predictions_exceed$area_fit_p1 <- predict(model_fits$total_area_p1, newdata = predictions_exceed)
  predictions_exceed$area_fit_p5 <- predict(model_fits$total_area_p5, newdata = predictions_exceed)
  predictions_exceed$area_fit_p10 <- predict(model_fits$total_area_p10, newdata = predictions_exceed)
  
  df_area_exceed_long <- df_area_over_time_exceed %>%
    select(t, total_area_p1, total_area_p5, total_area_p10) %>%
    pivot_longer(
      cols      = starts_with("total_area_p"),
      names_to  = "threshold",
      values_to = "total_area"
    ) %>%
    mutate(
      threshold = recode(
        threshold,
        "total_area_p1"  = "p > 0.01",
        "total_area_p5"  = "p > 0.05",
        "total_area_p10" = "p > 0.10"
      )
    )
  
  # Fitted curves in long format
  predictions_exceed_long <- predictions_exceed %>%
    pivot_longer(
      cols      = starts_with("area_fit_p"),
      names_to  = "threshold",
      values_to = "area_fit"
    ) %>%
    mutate(
      t = t_centered + min(df_area_over_time_exceed$t),
      threshold = recode(
        threshold,
        "area_fit_p1"  = "p > 0.01",
        "area_fit_p5"  = "p > 0.05",
        "area_fit_p10" = "p > 0.10"
      )
    )
  
  model_fit_exceed_nls <- ggplot() +
    # observed points + lines
    geom_point(
      data = df_area_exceed_long,
      aes(x = t, y = total_area, color = threshold),
      size = 1.5
    ) +
    geom_line(
      data = df_area_exceed_long,
      aes(x = t, y = total_area, color = threshold),
      linewidth = 0.8
    ) +
    # fitted curves (dashed)
    geom_line(
      data = predictions_exceed_long,
      aes(x = t, y = area_fit, color = threshold),
      linewidth = 0.9,
      linetype = "dashed"
    ) +
    labs(
      x     = "Year",
      y     = "Total high-risk area (km²)",
      color = "Prevalence threshold"
    ) +
    theme_bw() +
    theme(
      panel.border    = element_rect(color = "grey80", fill = NA),
      legend.position = "bottom"
    )
  save_figs(file.path(OUT_PLOT_DIR, paste0(mut, "_speed_spread_time_total_area_exceed_nls")), model_fit_exceed_nls, width = 10)
  
  # --- Area over time calculated from the prevalence draws ------------------
  # Loop through all time points and draws to generate the uncertainty
  area_draws <- list()
  for (k in seq_along(plot_times)) {
    draws_k <- posterior_draws[[k]] # [pixels x draws]
    t_k <- plot_times[k]
    
    # For each draw (column), classify pixels and sum area
    for (d in 1:D) {
      # The 'p' for this specific draw d at time t_k (p_{i, t, d})
      p_for_draw_d <- draws_k[, d]
      
      # For this specific draw, calculate the total area where p > tau
      # The logic: 1. (p_for_draw_d > p_thr) converts to TRUE/FALSE (1/0)
      #           2. sum() counts the number of exceeding pixels
      area_for_draw_p1 <- sum(p_for_draw_d > 0.01) * cell_area_km2
      area_for_draw_p5 <- sum(p_for_draw_d > 0.05) * cell_area_km2
      area_for_draw_p10 <- sum(p_for_draw_d > 0.1) * cell_area_km2
      
      area_draws[[length(area_draws) + 1]] <- data.frame(
        t = t_k,
        draw = d,
        area_p1 = area_for_draw_p1,
        area_p5 = area_for_draw_p5,
        area_p10 = area_for_draw_p10
      )
    }
  }
  
  df_area_draws <- bind_rows(area_draws)
  # Compute Quantiles for the Ribbon Plot
  df_area_quantiles <- df_area_draws %>%
    group_by(t) %>%
    summarise(
      area_p1_mean = mean(area_p1),
      area_p1_lower = quantile(area_p1, 0.025), # 2.5th percentile
      area_p1_upper = quantile(area_p1, 0.975), # 97.5th percentile
      area_p5_mean = mean(area_p5),
      area_p5_lower = quantile(area_p5, 0.025), # 2.5th percentile
      area_p5_upper = quantile(area_p5, 0.975), # 97.5th percentile
      area_p10_mean = mean(area_p10),
      area_p10_lower = quantile(area_p10, 0.025), # 2.5th percentile
      area_p10_upper = quantile(area_p10 , 0.975), # 97.5th percentile
      .groups = "drop"
    )
  
  # Fit Nonlinear Least Squares (NLS) regression
  model_fits_draws <- list()
  for (p_thr in c(0.01, 0.05, 0.10)) {
    # Map threshold to column name
    p_thr_col <- switch(as.character(p_thr),
                        "0.01" = "area_p1_mean",
                        "0.05" = "area_p5_mean",
                        "0.1"  = "area_p10_mean")
    
    # y vector for starting values
    y <- df_area_quantiles[[p_thr_col]]
    A_max <- max(y, na.rm = TRUE)
    t_centered <- df_area_quantiles$t - min(df_area_quantiles$t, na.rm = TRUE)
    df_tmp <- data.frame(t = df_area_quantiles$t, t_centered = t_centered, y = y)
    start_list <- list(
      a  = sqrt(A_max / pi) * 1.05,  # radius s.t. pi*a^2 ≈ max area
      b  = 0.1,
      t0 = 0
    )
    
    # Build formula: total_area_pX ~ (pi * a^2) * (1 - exp(-b * (t - t0)))^2
    fit_formula <- as.formula(
      "y ~ pi * a^2 * (1 - exp(-b * (t_centered - t0)))^2"
    )
    
    fit_model <- nls(
      fit_formula,
      data    = df_tmp,
      start   = start_list
      # start   = start_list,
      # algorithm = "port",
      # lower   = c(a = 0,   b = 0,   t0 = -5),
      # upper   = c(a = Inf, b = 10,  t0 = max(t_centered) + 5),
      # control = nls.control(maxiter = 500, minFactor = 1e-10, warnOnly = TRUE)
    )
    
    model_fits_draws[[p_thr_col]] <- fit_model
  }
  
  # Generate prediction points
  df_area_draws$t_centered <- df_area_draws$t - min(df_area_draws$t)
  time_range <- seq(min(df_area_draws$t_centered), max(df_area_draws$t_centered), length.out = 100)
  predictions_draws <- data.frame(t_centered = time_range)
  predictions_draws$area_fit_p1 <- predict(model_fits_draws$area_p1_mean, newdata = predictions_draws)
  predictions_draws$area_fit_p5 <- predict(model_fits_draws$area_p5_mean, newdata = predictions_draws)
  predictions_draws$area_fit_p10 <- predict(model_fits_draws$area_p10_mean, newdata = predictions_draws)
  
  df_area_draws_long <- df_area_quantiles %>%
    # 1) Gather all area_* columns into long format
    pivot_longer(
      cols = starts_with("area_"),
      names_to = c("thr_raw", "stat"),
      names_pattern = "area_(p[0-9]+)_(mean|lower|upper)",
      values_to = "area"
    ) %>%
    # 2) Make nice threshold labels
    mutate(
      threshold = dplyr::recode(
        thr_raw,
        "p1"  = "p > 0.01",
        "p5"  = "p > 0.05",
        "p10" = "p > 0.10"
      )
    ) %>%
    select(t, threshold, stat, area) %>%
    # 3) Spread mean / lower / upper back out as columns
    pivot_wider(
      names_from  = stat,
      values_from = area,
      names_prefix = "area_"  # gives area_mean, area_lower, area_upper
    ) %>%
    arrange(threshold, t)
  
  # Fitted curves in long format
  predictions_draws_long <- predictions_exceed %>%
    pivot_longer(
      cols      = starts_with("area_fit_p"),
      names_to  = "threshold",
      values_to = "area_fit"
    ) %>%
    mutate(
      t = t_centered + min(df_area_over_time_exceed$t),
      threshold = recode(
        threshold,
        "area_fit_p1"  = "p > 0.01",
        "area_fit_p5"  = "p > 0.05",
        "area_fit_p10" = "p > 0.10"
      )
    )
  
  model_fit_draws_nls <- ggplot() +
    # ribbon from quantiles
    geom_ribbon(
      data = df_area_draws_long,
      aes(x = t, ymin = area_lower, ymax = area_upper, fill = threshold),
      alpha = 0.2
    ) +
    # mean line
    geom_line(
      data = df_area_draws_long,
      aes(x = t, y = area_mean, color = threshold)
    ) +
    # fitted curves
    geom_line(
      data = predictions_draws_long,
      aes(x = t, y = area_fit, color = threshold),
      linetype = "dashed"
    ) +
    theme_bw() +
    labs(x = "Year", y = "Total high-risk area (km²)",
         color = "Prevalence threshold", fill = "Prevalence threshold")
  save_figs(file.path(OUT_PLOT_DIR, paste0(mut, "_speed_spread_time_total_area_draws_nls")), model_fit_draws_nls, width = 10)
  
}
