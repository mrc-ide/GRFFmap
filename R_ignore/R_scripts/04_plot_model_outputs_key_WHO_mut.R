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

cache_path <- function(mut, what, ext = c("rds","parquet")) {
  ext <- match.arg(ext)
  fn  <- sprintf("%s_lenS%s_lenT%s_%s.%s",
                 gsub("[:/ ]", "_", mut), what, ext)
  file.path(CACHE_DIR, fn)
}

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

make_or_load <- function(mut, ell_km, tau2) {
  cached_rds <- cache_path(mut, "optimized_ell_tau2_full_model_output", "rds") 
  stopifnot(file_exists(cached_rds))
  readRDS(cached_rds)
}

crop_long_df <- function(df, xlim, ylim) {
  dplyr::filter(df, x >= xlim[1], x <= xlim[2], y >= ylim[1], y <= ylim[2])
}

mask_to_africa <- function(df, africa_mask) {
  pts <- sf::st_as_sf(df, coords = c("x","y"), crs = 4326, remove = FALSE)
  inside <- lengths(sf::st_intersects(pts, africa_mask)) > 0  # sparse = TRUE by default
  df$p[!inside] <- NA_real_
  df
}

avg_2_year <- function(df){
  # Create a “year block” variable
  df <- df %>%
    mutate(
      year_block = case_when(
        t %in% c(2012, 2013) ~ "2012-2013",
        t %in% c(2014, 2015) ~ "2014-2015",
        t %in% c(2016, 2017) ~ "2016-2017",
        t %in% c(2018, 2019) ~ "2018-2019",
        t %in% c(2020, 2021) ~ "2020-2021",
        t %in% c(2022, 2023) ~ "2022-2023"
      )
    )
  
  # Average prevalence per pixel inside each 2-year block
  avg_two_year <- df %>%
    group_by(x, y, year_block) %>%
    summarise(mean_p = mean(p, na.rm = TRUE), .groups = "drop")
  
  return(avg_two_year)
}

add_year_block <- function(df) {
  df %>%
    mutate(
      year_block = case_when(
        t %in% c(2012, 2013) ~ "2012-2013",
        t %in% c(2014, 2015) ~ "2014-2015",
        t %in% c(2016, 2017) ~ "2016-2017",
        t %in% c(2018, 2019) ~ "2018-2019",
        t %in% c(2020, 2021) ~ "2020-2021",
        t %in% c(2022, 2023) ~ "2022-2023",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(year_block))
}

# --------------------------- Settings --------------------------------------
# model parameters
ell_km <- 120         # RFF length-scale in **kilometres**
tau2   <- 0.5

# inference parameters
# D        <- 300         # number of random frequencies (try 200–500)
# max_iter <- 5           # EM iterations
# z_eps    <- 1e-6        # threshold for |z| to use n/4 in E[ω]
# omega_floor <- 1e-10    # floor on ω to avoid 1/ω explosions
# jitter_S <- 1e-10       # diagonal jitter for innovation covariance

# prediction parameters
nx <- 100
ny <- 100
t_vec <- 1995:2025
t_num <- length(t_vec)

# plotting parameters
plot_times <- seq(2012, 2023, by = 1) # should be within t_vec
t_win <- 2           # window around each facet time (years). Points within this window are displayed

# --------------------------- Load & filter data ----------------------------
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

# --------------------------- Combine K13 data ----------------------------
k13_any_site_year <- dat %>%
  # keep only K13 mutations; adjust this filter to your naming convention
  filter(grepl("^k13", mutation, ignore.case = TRUE)) %>%
  group_by(year, longitude, latitude, country_name) %>%
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
    prevalence = prevalence,
    country_name = country_name
  )

dat_with_k13 <- dat %>%
  bind_rows(k13_any_site_year) %>%
  arrange(year, longitude, latitude, mutation)

# --------------------------- Define Mutations ----------------------------
all_who_mutations <- c("k13:comb", 
                       "k13:675:V", 
                       "k13:622:I", 
                       "k13:561:H", 
                       "k13:441:L")
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
    
  if (mut == "k13:675:V"){
    xlim = c(28, 37)
    ylim = c(-4, 5)
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
  
  # --- Crop Africa polygons ---------------
  bbox_sfc_ll <- st_as_sfc(st_bbox(c(
    xmin = xlim[1], xmax = xlim[2],
    ymin = ylim[1], ymax = ylim[2]
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
  
  # --- Build long data for lower / mean / upper ------------------------------
  cached <- make_or_load(mut, ell_km, tau2)
  cached <- readRDS("R_ignore/R_scripts/outputs/model_outputs/GRFF_model_output_key_WHO_mutations/k13_675_V_lenS120_lenT0.5_full_model_output.rds")
  
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
  
  # --- Build long data for lower / mean / upper ------------------------------
  plot_idx <- match(plot_times, t_vec)  # indices in 1:t_num

  p_long_median <- make_raster_long(p_post_median[ , , plot_idx, drop = FALSE], xs, ys, plot_times)
  p_long_mean <- make_raster_long(p_post_mean[ , , plot_idx, drop = FALSE], xs, ys, plot_times)
  p_long_CI <- make_raster_long(p_post_CI[ , , plot_idx, drop = FALSE], xs, ys, plot_times)
  
  exceed_prob_long_10 <- make_raster_long(exceed_prob_10[ , , plot_idx, drop = FALSE], xs, ys, plot_times)
  exceed_prob_long_5 <- make_raster_long(exceed_prob_5[ , , plot_idx, drop = FALSE], xs, ys, plot_times)
  exceed_prob_long_1 <- make_raster_long(exceed_prob_1[ , , plot_idx, drop = FALSE], xs, ys, plot_times)
  
  # # --- Mask out sea as white background ------------------------------
  # p_long_median <- mask_to_africa(p_long_median, africa_mask)
  # p_long_mean <- mask_to_africa(p_long_mean, africa_mask)
  # p_long_CI <- mask_to_africa(p_long_CI, africa_mask)
  # 
  # exceed_probdraws_long_10 <- mask_to_africa(exceed_prob_long_10, africa_mask)
  # exceed_prob_long_5  <- mask_to_africa(exceed_prob_long_5, africa_mask)
  # exceed_prob_long_1  <- mask_to_africa(exceed_prob_long_1, africa_mask) 
  # 
  # # --- Crop dataframes to bbox ------------------------------
  # p_long_median   <- crop_long_df(p_long_median, xlim, ylim)
  # p_long_mean   <- crop_long_df(p_long_mean, xlim, ylim)
  # p_long_CI   <- crop_long_df(p_long_CI, xlim, ylim)
  # 
  # exceed_prob_long_10 <- crop_long_df(exceed_prob_long_10, xlim, ylim)
  # exceed_prob_long_5  <- crop_long_df(exceed_prob_long_5,  xlim, ylim)
  # exceed_prob_long_1  <- crop_long_df(exceed_prob_long_1,  xlim, ylim)
  
  # --- Create observed data points ------------------------------
  points_df <- dat_sub %>%
    filter(year >= 2012, year <= 2023) %>%     # keep both years in each block
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
  
  # --- Average exceed prob and pred prev every two years --------------------
  p_long_median_avg_2y <- avg_2_year(p_long_median)
  p_long_mean_avg_2y <- avg_2_year(p_long_mean)
  p_long_CI_avg_2y <- avg_2_year(p_long_CI)
  
  exceed_prob_long_1_avg_2y <- avg_2_year(exceed_prob_long_1)
  exceed_prob_long_5_avg_2y <- avg_2_year(exceed_prob_long_5)
  exceed_prob_long_10_avg_2y <- avg_2_year(exceed_prob_long_10)
  
  points_df_2y <- add_year_block(points_df %>% filter(t %in% plot_times))
  
  # --- Helper to generate a consistent ggplot with point overlays ------------
  plot_theme <- function(p) {
    p + theme(
      strip.background = element_rect(fill = "white", color = NA),
      strip.text = element_text(color = "black", face = "plain"),
      panel.border = element_rect(color = "grey80", fill = NA),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      title = element_text(size = 12),
      axis.text.x = element_text(size = 8, angle = 30, hjust = 1),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(hjust = 0),
      legend.justification = "center"
    )
  }
  
  plot_layer <- function(p_long_df, title_text, shp = shape_Africa, shp_water = shape_water, add_points = FALSE, add_legend = FALSE) {
    p <- ggplot() +
      geom_raster(aes(x = x, y = y, fill = mean_p*100), data = p_long_df) +
      geom_sf(data = shp, linewidth = 0.2, fill = NA, color = "white") +
      geom_sf(data = shp_water, fill = "white", colour = NA) +
      coord_sf(
        xlim = unname(as.numeric(xlim)),
        ylim = unname(as.numeric(ylim)),
        expand = FALSE,
        crs = sf::st_crs(4326),          # reproject sf layers to WGS84
        default_crs = sf::st_crs(4326),  # interpret x/y (raster/points) as WGS84
        clip = "on"
      )
    if (add_points) {
      p <- p + geom_point(
        aes(x = x, y = y, fill = p),          # <- raw prevalence per point
        data = points_df_2y %>% arrange(p),   # <- NOT averaged
        shape = 21, colour = "darkgrey", size = 1.5, stroke = 0.2
      )
    }
    if (!add_legend){
      p <- p + theme(legend.position = "none")
    }
    else{
      p <- p + theme(legend.position = "bottom",
                     legend.key.width = unit(2, "cm"),
                     legend.key.height = unit(0.2, "cm"))
    }
    if (title_text != ""){
      p <- p + labs(title = title_text)
    }
    p <- p + scale_fill_gradientn(colours = pp_cols,
                                  values  = pp_vals,
                                  limits = c(0, 100), 
                                  name = "Prevalence (%)",
                                  breaks = seq(0, 100, by = 10),
                                  na.value = "white") +
      facet_wrap(~year_block, nrow = 1) +
      labs(x = "Longitude", y = "Latitude")
    if (mut == "k13:561:H"){
      p <- p + scale_x_continuous(breaks = seq(min(p_long_df$x)+1, max(p_long_df$x), by = 2))
    }
    plot_theme(p)
  }
  
  plot_exceedance <-function(exceed_prob, title_text, legend_title, shp = shape_Africa, shp_water = shape_water, add_points = FALSE, add_legend = FALSE) {
    p <- ggplot() +
      theme_bw() +
      annotate(
        "rect",
        xmin = xlim[1], xmax = xlim[2],
        ymin = ylim[1], ymax = ylim[2],
        fill = "white", colour = NA
      ) +
      geom_raster(aes(x = x, y = y, fill = mean_p*100), data = exceed_prob) +
      geom_sf(data = shape_Africa, linewidth = 0.2, fill = NA, color = "white") +
      geom_sf(data = shp_water, fill = "white", colour = NA) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = st_crs(4326)) +
      scale_fill_viridis_c(
        option = "mako",
        direction = 1,
        limits = c(0, 100),
        breaks = seq(0, 100, by = 10),
        name = legend_title,
        na.value = "white"
      ) +
      facet_wrap(~year_block, nrow = 1) +
      labs(
        x = "Longitude",
        y = "Latitude"
      ) +
      theme(legend.position = "bottom",
            legend.key.width = unit(2, "cm"),
            legend.key.height = unit(0.2, "cm"))
    if (title_text != ""){
      p <- p + labs(title = title_text)
    }
    if (mut == "k13:561:H"){
      p <- p + scale_x_continuous(breaks = seq(min(exceed_prob$x)+1, max(exceed_prob$x), by = 2))
    }
    plot_theme(p)
  }
  
  # --- Create and print the three separate plots -----------------------------
  
  plot_median <- plot_layer(p_long_median_avg_2y, clean_mut, shape_Africa_crop, shape_water_crop, add_points = TRUE, add_legend = TRUE)
  plot_median_no_points  <- plot_layer(p_long_median_avg_2y, clean_mut, shape_Africa, shape_water, add_points = FALSE, add_legend = TRUE)

  plot_mean <- plot_layer(p_long_mean_avg_2y, clean_mut, shape_Africa_crop, shape_water_crop, add_points = TRUE, add_legend = TRUE)
  plot_mean_no_points  <- plot_layer(p_long_mean_avg_2y, clean_mut, shape_Africa, shape_water, add_points = FALSE, add_legend = TRUE)
  
  plot_CI <- plot_layer(p_long_CI_avg_2y, clean_mut, shape_Africa_crop, shape_water_crop, add_points = TRUE, add_legend = TRUE)
  plot_CI_no_points <- plot_layer(p_long_CI_avg_2y, clean_mut, shape_Africa_crop, shape_water_crop, add_points = FALSE, add_legend = TRUE)
  
  plot_exceed_1 <- plot_exceedance(exceed_prob_long_1_avg_2y, title_text = clean_mut, legend_title =expression(Pr(Prevalence >= 1*"%")))
  plot_exceed_5 <- plot_exceedance(exceed_prob_long_5_avg_2y, title_text = clean_mut, legend_title =expression(Pr(Prevalence >= 5*"%")))
  plot_exceed_10 <- plot_exceedance(exceed_prob_long_10_avg_2y, title_text = clean_mut, legend_title =expression(Pr(Prevalence >= 10*"%")))
  
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

