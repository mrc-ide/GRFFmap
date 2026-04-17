# ── Packages ────────────────────────────────────────────────────────────────────
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

load_all()

set.seed(1)

# --- Settings -----------------------------------------------------------------

# prediction parameters
nx <- 200
ny <- 200
t_vec <- 1995:2025
t_num <- length(t_vec)
plot_times <- seq(2001, 2024, by = 1) # should be within t_vec

# --------------------------- Load & filter data ----------------------------
CACHE_DIR <- "R_ignore/R_scripts/outputs/model_outputs/GRFF_model_output_parter_drug_mutations"

OUT_MEAN_DIR   <- file.path("key_partner_drug_mutations", "mean_prev_prediction_avg_2_year")
OUT_MEDIAN_DIR   <- file.path("key_partner_drug_mutations", "median_prev_prediction_avg_2_year")
OUT_CI_DIR   <- file.path("key_partner_drug_mutations", "CI_diff_prediction_avg_2_year")
OUT_SUPP <- file.path("supplemental", "key_partner_drug_mutations")

dir_create(c(paste0("R_ignore/R_scripts/outputs/plots/", OUT_MEAN_DIR),
             paste0("R_ignore/R_scripts/outputs/plots/", OUT_MEDIAN_DIR),
             paste0("R_ignore/R_scripts/outputs/plots/", OUT_CI_DIR),
             paste0("R_ignore/R_scripts/outputs/plots/", OUT_SUPP)),
           recurse = TRUE)

dat <- read.csv("R_ignore/R_scripts/data/all_mutations_get_prevalence.csv") |>
  mutate(collection_day = as.Date(collection_day)) |>
  select(study_id, survey_id, longitude, latitude, year, numerator, denominator, prevalence, mutation, country_name)

shape_Africa <- readRDS("R_ignore/R_scripts/data/shapefiles/sf_admin0_africa.rds")
shape_water  <- sf::st_read("R_ignore/R_scripts/data/shapefiles/africa_water_bodies.shp", quiet = TRUE)
africa_admin0 <- readRDS("R_ignore/R_scripts/data/sf_admin0_africa.rds")

# --------------------------- Define Mutations ----------------------------
all_pd_mutations <- c("mdr1:86:Y", "crt:76:T")

# --- Obtain Africa lat long ---------------------------------------------------
bb <- sf::st_bbox(shape_Africa)  # or shape_Africa_crop if that is full Africa

# after bb <- st_bbox(shape_Africa)
x_pad <- 2   # degrees of longitude
y_pad <- 2   # degrees of latitude

xlim_africa <- c(bb["xmin"] - x_pad, bb["xmax"] + x_pad)
ylim_africa <- c(bb["ymin"] - y_pad, bb["ymax"] + y_pad)

# --- Create white raster for outside Africa ------------------------------

# Construct spatial boundaries
bbox_info <- build_africa_bbox_and_crop(shape_Africa, shape_water, buf_km = 200)

# Raster-compatible land mask defining valid (on-land) spatial locations
africa_mask <- bbox_info$africa_mask

# Cropped vector geometries for Africa landmass and water bodies, restricted to the buffered bounding box
shape_Africa_crop <- bbox_info$africa_crop
shape_water_crop  <- bbox_info$water_crop

# Create a unified, valid land geometry for Africa
africa_mask_land <- shape_Africa_crop %>%
  sf::st_union() %>%
  sf::st_make_valid()

# Coordinate reference systems
CRS_LL     <- sf::st_crs(4326)
CRS_METRIC <- sf::st_crs(3857) 

# --- Grey out non-endemic malaria African countries ---------------------------
non_malaria_countries <- c(
  "Egypt", "Morocco", "Libya", "Tunisia", "Algeria",
  "Cabo Verde", "Lesotho", "Mauritius", "Seychelles"
)

shape_non_malaria <- africa_admin0 %>%
  dplyr::filter(name_0 %in% non_malaria_countries)

# --- Loop over each mutation ------------------------------
for (mut in all_pd_mutations){
  message("Processing ", mut)
  
  # Define clean version of mutation
  if (mut == "crt:76:T"){clean_mut <- gsub("^crt:(\\d+):([A-Za-z])$", "crt \\1\\2", mut)}
  if (mut == "mdr1:86:Y"){clean_mut <- gsub("^mdr1:(\\d+):([A-Za-z])$", "mdr1 \\1\\2", mut)}
  
  # Filter data
  dat_sub <- dat |>
    filter(mutation == mut) |>
    filter(is.finite(numerator), is.finite(denominator), denominator > 0) 
  
  # --- Build long data for lower / mean / upper ------------------------------
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
  
  xlim <- xlim_africa
  ylim <- ylim_africa
  
  # --- Build long data for lower / mean / upper -------------------------------
  plot_idx <- match(plot_times, t_vec)  # indices in 1:t_num
  
  p_long_median <- make_raster_long(p_post_median[ , , plot_idx, drop = FALSE], xs, ys, plot_times)
  p_long_mean <- make_raster_long(p_post_mean[ , , plot_idx, drop = FALSE], xs, ys, plot_times)
  p_long_CI <- make_raster_long(p_post_CI[ , , plot_idx, drop = FALSE], xs, ys, plot_times)
  
  # --- Mask out sea as white background ---------------------------------------
  p_long_median <- mask_to_africa(p_long_median, africa_mask)
  p_long_mean <- mask_to_africa(p_long_mean, africa_mask)
  p_long_CI <- mask_to_africa(p_long_CI, africa_mask)
  
  # --- Crop dataframes to bbox ------------------------------------------------
  p_long_median   <- crop_long_df(p_long_median, xlim, ylim)
  p_long_mean   <- crop_long_df(p_long_mean, xlim, ylim)
  p_long_CI   <- crop_long_df(p_long_CI, xlim, ylim)
  
  # --- Create observed data points --------------------------------------------
  points_df <- dat_sub %>%
    filter(year >= min(plot_times), year <= max(plot_times)) %>%     # keep both years in each block
    mutate(
      p_obs = numerator / denominator,
      t     = year
    ) %>%
    filter(longitude >= xlim[1], longitude <= xlim[2],
           latitude  >= ylim[1],  latitude  <= ylim[2]) %>%
    transmute(
      x = longitude,
      y = latitude,
      t = t,
      p = p_obs * 100   # or prevalence, but keep consistent with your color scale
    )
  
  # --- Average exceed prob and pred prev every two years ----------------------
  p_long_median_avg_2y <- avg_2_year_pd(p_long_median)
  p_long_mean_avg_2y <- avg_2_year_pd(p_long_mean)
  p_long_CI_avg_2y <- avg_2_year_pd(p_long_CI)
  
  points_df_2y <- points_df %>%
    dplyr::filter(t %in% plot_times) %>%
    add_year_group_pd(year_col = t) %>%
    dplyr::mutate(t = year_group) %>%
    dplyr::select(-year_group)
  
  mut_lims <- list()
  mut_lims[[1]] <- xlim
  mut_lims[[2]] <- ylim
  axis_break <- 20
  
  # --- Create and print the three separate plots ------------------------------
  sf::sf_use_s2(TRUE)
  
  plot_per_year_median <- plot_prev_layer(
    p_long_df   = p_long_median,
    title_text  = clean_mut,
    shp         = shape_Africa_crop,
    shp_water   = shape_water_crop,
    shp_non_malaria = shape_non_malaria,
    facet_n_row = 4,
    lims        = mut_lims,
    x_axis_break = axis_break,
    y_axis_break = axis_break,
    points_df   = points_df,
    add_legend  = TRUE
  )
  
  plot_binned_year_median <- plot_prev_layer(
    p_long_df   = p_long_median_avg_2y,
    title_text  = clean_mut,
    shp         = shape_Africa_crop,
    shp_water   = shape_water_crop,
    shp_non_malaria = shape_non_malaria,
    facet_n_row = 2,
    lims        = mut_lims,
    x_axis_break = axis_break,
    y_axis_break = axis_break,
    points_df   = points_df_2y,
    add_legend  = TRUE
  )
  
  plot_binned_year_median_no_points  <- plot_prev_layer(
    p_long_df   = p_long_median_avg_2y,
    title_text  = clean_mut,
    shp         = shape_Africa_crop,
    shp_water   = shape_water_crop,
    shp_non_malaria = shape_non_malaria,
    facet_n_row = 2,
    lims        = mut_lims,
    x_axis_break = axis_break,
    y_axis_break = axis_break,
    add_legend  = TRUE
  )
  
  plot_binned_year_mean <- plot_prev_layer(
    p_long_df   = p_long_mean_avg_2y, 
    title_text  = clean_mut,
    shp         = shape_Africa_crop,
    shp_water   = shape_water_crop,
    shp_non_malaria = shape_non_malaria,
    facet_n_row = 2,
    lims        = mut_lims,
    x_axis_break = axis_break,
    y_axis_break = axis_break,
    points_df   = points_df_2y,
    add_legend  = TRUE
  )
  
  plot_binned_year_mean_no_points  <- plot_prev_layer(
    p_long_df   = p_long_mean_avg_2y, 
    title_text  = clean_mut,
    shp         = shape_Africa_crop,
    shp_water   = shape_water_crop,
    shp_non_malaria = shape_non_malaria,
    facet_n_row = 2,
    lims        = mut_lims,
    x_axis_break = axis_break,
    y_axis_break = axis_break,
    add_legend  = TRUE
  )
  
  plot_binned_year_CI <- plot_CI_layer(
    p_long_df   = p_long_CI_avg_2y, 
    title_text  = clean_mut,
    shp         = shape_Africa_crop,
    shp_water   = shape_water_crop,
    shp_non_malaria = shape_non_malaria,
    facet_n_row = 2,
    lims        = mut_lims,
    x_axis_break = axis_break,
    y_axis_break = axis_break,
    points_df   = points_df_2y,
    viridis = TRUE,
    add_legend = TRUE
  )
  
  plot_binned_year_CI_no_points <- plot_CI_layer(
    p_long_df   = p_long_CI_avg_2y, 
    title_text  = clean_mut,
    shp         = shape_Africa_crop,
    shp_water   = shape_water_crop,
    shp_non_malaria = shape_non_malaria,
    facet_n_row = 2,
    lims        = mut_lims,
    x_axis_break = axis_break,
    y_axis_break = axis_break,
    viridis = TRUE,
    add_legend = TRUE
  )
  
  plot_height = 4.6
  # --- Save plots -----------------------------
  save_figs(file.path(OUT_MEAN_DIR, paste0(mut, "_mean_perc_prev")), plot_binned_year_mean, width = 10, height = plot_height)
  save_figs(file.path(OUT_MEAN_DIR, paste0(mut, "_mean_perc_prev_no_points")), plot_binned_year_mean_no_points, width = 10, height = plot_height)
  
  save_figs(file.path(OUT_MEDIAN_DIR, paste0(mut, "_median_perc_prev")), plot_binned_year_median, width = 10, height = plot_height)
  save_figs(file.path(OUT_MEDIAN_DIR, paste0(mut, "_median_perc_prev_no_points")), plot_binned_year_median_no_points, width = 10, height = plot_height)
  
  save_figs(file.path(OUT_SUPP, paste0(mut, "_median_perc_prev")), plot_per_year_median, width = 10, height = 10)
  
  save_figs(file.path(OUT_CI_DIR, paste0(mut, "_CI_diff_prev")), plot_binned_year_CI, width = 10, height = plot_height)
  save_figs(file.path(OUT_CI_DIR, paste0(mut, "_CI_diff_prev_no_points")), plot_binned_year_CI_no_points, width = 10, height = plot_height)
  
  
  print(paste0("Saved figures for ", mut)) 
}
