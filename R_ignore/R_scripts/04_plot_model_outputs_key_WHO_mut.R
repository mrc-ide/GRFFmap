# --- Packages -----------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(fs)
  library(arrow)       # optional but nice for big cached tables
  library(devtools)
  library(dplyr)
  library(fs)
  library(lwgeom)
  library(viridis)
  library(scales)
})

devtools::load_all()

# --- Settings -----------------------------------------------------------------

# prediction parameters
nx <- 200
ny <- 200
t_vec <- 1995:2025
t_num <- length(t_vec)

# plotting parameters
plot_times <- seq(2012, 2023, by = 1) # should be within t_vec

# --- Load & filter data -------------------------------------------------------
CACHE_DIR <- "R_ignore/R_scripts/outputs/model_outputs/GRFF_model_output_key_WHO_mutations"

OUT_DIR <- "R_ignore/R_scripts/outputs/plots/key_k13_mutations"
dir_create(OUT_DIR)

OUT_MEAN_DIR   <- file.path("key_k13_mutations", "mean_prev_prediction_avg_2_year")
OUT_MEDIAN_DIR   <- file.path("key_k13_mutations", "median_prev_prediction_avg_2_year")
OUT_CI_DIR   <- file.path("key_k13_mutations", "CI_diff_prediction_avg_2_year")
OUT_EXCEED_DIR <- file.path("key_k13_mutations", "exceedance_prob_avg_2_year")

dir_create(c(paste0("R_ignore/R_scripts/outputs/plots/", OUT_MEAN_DIR),
             paste0("R_ignore/R_scripts/outputs/plots/", OUT_MEDIAN_DIR),
             paste0("R_ignore/R_scripts/outputs/plots/", OUT_CI_DIR),
             paste0("R_ignore/R_scripts/outputs/plots/", OUT_EXCEED_DIR)),
           recurse = TRUE)

dat <- read.csv("R_ignore/R_scripts/data/all_who_get_prevalence.csv") |>
  mutate(collection_day = as.Date(collection_day)) |>
  select(longitude, latitude, year, numerator, denominator, prevalence, mutation, country_name)

shape_Africa <- readRDS("R_ignore/R_scripts/data/sf_admin0_africa.rds")
shape_water  <- sf::st_read("R_ignore/R_scripts/data/shapefiles/africa_water_bodies.shp", quiet = TRUE)

# --- Add combined K13 data ----------------------------------------------------
dat_with_k13 <- add_combined_k13(dat)

# --------------------------- Define Mutations ----------------------------
all_who_mutations <- c("k13:comb", 
                       "k13:675:V", 
                       "k13:622:I", 
                       "k13:561:H")
# Note for Cecile: the following mutations have a wide grid where prev > 0:
# 441L, 449A, 469F, 473I, 539T, 553L, 561H, 568G, 574L, 622I

# --- Create white raster for outside Africa ------------------------------
CRS_LL     <- sf::st_crs(4326)  # WGS84 lon/lat for plotting
CRS_METRIC <- sf::st_crs(3857)  # planar for buffers/ops

shape_Africa_ll <- shape_Africa |> sf::st_transform(CRS_LL)
shape_water_ll  <- shape_water  |> sf::st_transform(CRS_LL)

# --- Define color scheme ------------------------------
pp_cols<- c(
  "#5E3A9B",  # dark purple (0)
  "#8cc4e0",  # medium-dark blue
  "#a8ecbf",  # mint-teal
  "palegreen3", # <-- ADDED GREEN
  "khaki2",
  "#f3c86b",
  "orange",
  "hotpink",   # <-- ADDED PINK
  "red",      
  "darkred"    
)

pp_vals <- rescale(c(0, 1, 5, 10, 20, 30, 50, 70, 90, 95 , 100))

# --- Loop over each mutation ------------------------------
mut <- "k13:675:V"
for (mut in all_who_mutations){
  message("Processing ", mut)
  clean_mut <- paste0("k13 ", gsub("^k13:(\\d+):([A-Za-z])$", "\\1\\2", mut))
  
  if (mut == "k13:comb"){
    clean_mut = ""
  }
  
  dat_sub <- dat_with_k13 |>
    filter(mutation == mut) |>
    filter(is.finite(numerator), is.finite(denominator), denominator > 0,
           ) |>
    mutate(
      longitude = as.numeric(longitude),
      latitude  = as.numeric(latitude)
    )
  
  # --- Define lat and lon bbox limits and axis breaks -------------------------
  if (mut == "k13:675:V"){
    mut_lims <- list(
      xlim = c(28, 37),
      ylim = c(-4, 5)
    )
    axis_break = 2
  } else if (mut == "k13:comb"){
    mut_lims <- list(
      xlim = c(28, 48),
      ylim = c(-4.60, 18)
    )
    axis_break = 5
  } else if (mut == "k13:622:I"){
    mut_lims <- list(
      xlim = c(32, 45),
      ylim = c(4, 17)
    )
    axis_break = 2
  } else if (mut == "k13:561:H"){
    mut_lims <- list(
      xlim = c(28, 32),
      ylim = c(-4, 0)
    )
    axis_break = 2
  }
  
  # --- Crop Africa polygons ---------------------------------------------------
  bbox_sfc_ll <- st_as_sfc(st_bbox(c(
    xmin = mut_lims[[1]][1], 
    xmax = mut_lims[[1]][2],
    ymin = mut_lims[[2]][1], 
    ymax = mut_lims[[2]][2]
  ), crs = CRS_LL))
  
  old_s2 <- sf_use_s2()
  sf_use_s2(FALSE)
  
  bbox_planar <- st_transform(bbox_sfc_ll, CRS_METRIC)
  
  shape_Africa_planar <- shape_Africa_ll |>
    st_transform(CRS_METRIC) |>
    st_make_valid() |>
    st_buffer(0) |>                # classic fix for minor self-intersections
    st_intersection(bbox_planar)   # or st_crop(bbox_planar)
  
  shape_water_planar <- shape_water_ll |>
    st_transform(CRS_METRIC) |>
    st_make_valid() |>
    st_buffer(0) |>
    st_intersection(bbox_planar)
  
  # Transform back to lon/lat for plotting
  shape_Africa_crop <- st_transform(shape_Africa_planar, 4326)
  shape_water_crop  <- st_transform(shape_water_planar,  4326)
  
  # Restore s2
  sf_use_s2(old_s2)
  
  africa_mask <- shape_Africa_crop %>%
    sf::st_union() %>%
    sf::st_make_valid()
  
  # --- Project to local Cartesian (km) ----------------------------------------
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
  
  # --- Build long data for lower / mean / upper -------------------------------
  cached <- make_or_load(mut)
  
  tau2 <- cached$settings$tau2
  ell_km <- cached$settings$ell_km
  
  message("Optimal ell_km = ", ell_km, " tau2 = ", tau2)

  dim <- cached$p_post_median
  p_post_mean <- cached$p_post_mean
  p_post_median <- cached$p_post_median
  p_post_CI <- cached$p_post_CI
  
  xs <- cached$xs
  ys <- cached$ys
  exceed_prob <- cached$exceed_post_list
  
  exceed_prob_1 <- exceed_prob[ , , , 1]
  exceed_prob_5 <- exceed_prob[ , , , 2]
  exceed_prob_10 <- exceed_prob[ , , , 3]
  
  # --- Build long data for lower / mean / upper -------------------------------
  plot_idx <- match(plot_times, t_vec)  # indices in 1:t_num

  p_long_median <- make_raster_long(p_post_median[ , , plot_idx, drop = FALSE], xs, ys, plot_times)
  p_long_mean <- make_raster_long(p_post_mean[ , , plot_idx, drop = FALSE], xs, ys, plot_times)
  p_long_CI <- make_raster_long(p_post_CI[ , , plot_idx, drop = FALSE], xs, ys, plot_times)
  
  exceed_prob_long_10 <- make_raster_long(exceed_prob_10[ , , plot_idx, drop = FALSE], xs, ys, plot_times)
  exceed_prob_long_5 <- make_raster_long(exceed_prob_5[ , , plot_idx, drop = FALSE], xs, ys, plot_times)
  exceed_prob_long_1 <- make_raster_long(exceed_prob_1[ , , plot_idx, drop = FALSE], xs, ys, plot_times)
  
  # --- Mask out sea as white background ---------------------------------------
  p_long_median <- mask_to_africa(p_long_median, africa_mask)
  p_long_mean <- mask_to_africa(p_long_mean, africa_mask)
  p_long_CI <- mask_to_africa(p_long_CI, africa_mask)
  
  exceed_probdraws_long_10 <- mask_to_africa(exceed_prob_long_10, africa_mask)
  exceed_prob_long_5  <- mask_to_africa(exceed_prob_long_5, africa_mask)
  exceed_prob_long_1  <- mask_to_africa(exceed_prob_long_1, africa_mask) 
  
  # --- Crop dataframes to bbox ------------------------------------------------
  p_long_median   <- crop_long_df(p_long_median, mut_lims[[1]], mut_lims[[2]])
  p_long_mean   <- crop_long_df(p_long_mean, mut_lims[[1]], mut_lims[[2]])
  p_long_CI   <- crop_long_df(p_long_CI, mut_lims[[1]], mut_lims[[2]])
  
  exceed_prob_long_10 <- crop_long_df(exceed_prob_long_10, mut_lims[[1]], mut_lims[[2]])
  exceed_prob_long_5  <- crop_long_df(exceed_prob_long_5,  mut_lims[[1]], mut_lims[[2]])
  exceed_prob_long_1  <- crop_long_df(exceed_prob_long_1,  mut_lims[[1]], mut_lims[[2]])
  
  # --- Create observed data points --------------------------------------------
  points_df <- dat_sub %>%
    filter(year >= min(plot_times), year <= max(plot_times)) %>%     # keep both years in each block
    mutate(
      p_obs = numerator / denominator,
      t     = year
    ) %>%
    filter(longitude >= mut_lims[[1]][1], longitude <= mut_lims[[1]][2],
           latitude  >= mut_lims[[2]][1],  latitude  <= mut_lims[[2]][2]) %>%
    transmute(
      x = longitude,
      y = latitude,
      t = t,
      p = p_obs * 100   # or prevalence, but keep consistent with your color scale
    )
  
  # --- Average exceed prob and pred prev every two years ----------------------
  p_long_median_avg_2y <- avg_2_year(p_long_median)
  p_long_mean_avg_2y <- avg_2_year(p_long_mean)
  p_long_CI_avg_2y <- avg_2_year(p_long_CI)
  
  exceed_prob_long_1_avg_2y <- avg_2_year(exceed_prob_long_1)
  exceed_prob_long_5_avg_2y <- avg_2_year(exceed_prob_long_5)
  exceed_prob_long_10_avg_2y <- avg_2_year(exceed_prob_long_10)
  
  points_df_2y <- points_df %>%
    dplyr::filter(t %in% plot_times) %>%
    add_year_group(year_col = t) %>%
    dplyr::mutate(t = year_group) %>%
    dplyr::select(-year_group)
  
  # --- Create and print the three separate plots ------------------------------
  sf::sf_use_s2(TRUE)
  plot_median <- plot_prev_layer(
    p_long_df   = p_long_median_avg_2y,
    title_text  = clean_mut,
    shp         = shape_Africa_crop,
    shp_water   = shape_water_crop,
    lims        = mut_lims,
    x_axis_break = axis_break,
    y_axis_break = axis_break,
    points_df   = points_df_2y,
    add_legend  = TRUE
  )
  
  plot_median_no_points  <- plot_prev_layer(
    p_long_df   = p_long_median_avg_2y,
    title_text  = clean_mut,
    shp         = shape_Africa_crop,
    shp_water   = shape_water_crop,
    lims        = mut_lims,
    x_axis_break = axis_break,
    y_axis_break = axis_break,
    add_legend  = TRUE
  )

  plot_mean <- plot_prev_layer(
    p_long_df   = p_long_mean_avg_2y, 
    title_text  = clean_mut,
    shp         = shape_Africa_crop,
    shp_water   = shape_water_crop,
    lims        = mut_lims,
    x_axis_break = axis_break,
    y_axis_break = axis_break,
    points_df   = points_df_2y,
    add_legend  = TRUE
  )
  
  plot_mean_no_points  <- plot_prev_layer(
    p_long_df   = p_long_mean_avg_2y, 
    title_text  = clean_mut,
    shp         = shape_Africa_crop,
    shp_water   = shape_water_crop,
    lims        = mut_lims,
    x_axis_break = axis_break,
    y_axis_break = axis_break,
    add_legend  = TRUE
  )
  
  plot_CI <- plot_prev_layer(
    p_long_df   = p_long_CI_avg_2y, 
    title_text  = clean_mut,
    shp         = shape_Africa_crop,
    shp_water   = shape_water_crop,
    lims        = mut_lims,
    x_axis_break = axis_break,
    y_axis_break = axis_break,
    points_df   = points_df_2y,
    add_legend = TRUE
    )
  
  plot_CI_no_points <- plot_prev_layer(
    p_long_df   = p_long_CI_avg_2y, 
    title_text  = clean_mut,
    shp         = shape_Africa_crop,
    shp_water   = shape_water_crop,
    lims        = mut_lims,
    x_axis_break = axis_break,
    y_axis_break = axis_break,
    add_legend = TRUE
  )
  
  # --- Plot exceedance probability --------------------------------------------
  plot_exceed_1 <- plot_exceedance(exceed_prob_long_1_avg_2y, 
                                   title_text   = clean_mut, 
                                   shp          = shape_Africa_crop,
                                   shp_water    = shape_water_crop,
                                   legend_title = expression(Pr(Prevalence >= 1*"%")),
                                   lims         = mut_lims,
                                   x_axis_break = axis_break,
                                   y_axis_break = axis_break
                                   )
  plot_exceed_5 <- plot_exceedance(exceed_prob_long_5_avg_2y, 
                                   title_text = clean_mut, 
                                   shp          = shape_Africa_crop,
                                   shp_water    = shape_water_crop,
                                   legend_title =expression(Pr(Prevalence >= 5*"%")),
                                   lims         = mut_lims,
                                   x_axis_break = axis_break,
                                   y_axis_break = axis_break)
  plot_exceed_10 <- plot_exceedance(exceed_prob_long_10_avg_2y, 
                                    title_text = clean_mut, 
                                    shp          = shape_Africa_crop,
                                    shp_water    = shape_water_crop,
                                    legend_title =expression(Pr(Prevalence >= 10*"%")),
                                    lims         = mut_lims,
                                    x_axis_break = axis_break,
                                    y_axis_break = axis_break)
  
  # --- Save plots -------------------------------------------------------------
  save_figs(file.path(OUT_MEAN_DIR, paste0(mut, "_mean_perc_prev")), plot_mean, width = 10)
  save_figs(file.path(OUT_MEAN_DIR, paste0(mut, "_mean_perc_prev_no_points")), plot_mean_no_points, width = 10)
  
  save_figs(file.path(OUT_MEDIAN_DIR, paste0(mut, "_median_perc_prev")), plot_median, width = 10)
  save_figs(file.path(OUT_MEDIAN_DIR, paste0(mut, "_median_perc_prev_no_points")), plot_median_no_points, width = 10)
  
  save_figs(file.path(OUT_CI_DIR, paste0(mut, "_CI_diff_prev")), plot_CI, width = 10)
  save_figs(file.path(OUT_CI_DIR, paste0(mut, "_CI_diff_prev_no_points")), plot_CI_no_points, width = 10)

  save_figs(file.path(OUT_EXCEED_DIR, paste0(mut, "_exceednace_prob_1")), plot_exceed_1, width = 10)
  save_figs(file.path(OUT_EXCEED_DIR, paste0(mut, "_exceednace_prob_5")), plot_exceed_5, width = 10)
  save_figs(file.path(OUT_EXCEED_DIR, paste0(mut, "_exceednace_prob_10")), plot_exceed_10, width = 10)
  
  print(paste0("Saved figures for ", mut)) 
}

