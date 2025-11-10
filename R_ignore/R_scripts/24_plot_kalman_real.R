# ── Packages ────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(fs)
  library(arrow)       # optional but nice for big cached tables
  library(devtools)
  library(dplyr)
  library(fs)
  library(viridis)
  library(scales)
})

load_all()

set.seed(1)

cache_path <- function(mut, lenS, lenT, what, ext = c("rds","parquet")) {
  ext <- match.arg(ext)
  fn  <- sprintf("%s_lenS%s_lenT%s_%s.%s",
                 gsub("[:/ ]", "_", mut), lenS, lenT, what, ext)
  file.path(CACHE_DIR, fn)
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
    dplyr::filter(t > 2013 & t < 2026)
}

make_or_load <- function(mut, ell_km, tau2) {
  cached_rds <- cache_path(mut, ell_km, tau2, "full_model_output", "rds") 
  stopifnot(file_exists(cached_rds))
  readRDS(cached_rds)
}

# --------------------------- Settings --------------------------------------
# model parameters
ell_km <- 80           # RFF length-scale in **kilometres**
tau2   <- 0.1        # RW1 variance in feature space
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
t_win <- 2           # window around each facet time (years). Points within this window are displayed

# --------------------------- Load & filter data ----------------------------
CACHE_DIR <- "R_ignore/R_scripts/outputs/GRFF_kalman_cache_annual"
dir_create(CACHE_DIR)
OUT_PREV_DIR   <- file.path("prev_prediction_grouped_year_kalman")
OUT_EXCEED_DIR    <- file.path("exceedance_prob_grouped_year_kalman")
dir_create(c(paste0("R_ignore/R_scripts/outputs/plots/", OUT_PREV_DIR), paste0("R_ignore/R_scripts/outputs/plots/", OUT_EXCEED_DIR)), recurse = TRUE)

dat <- read.csv("R_ignore/R_scripts/data/all_who_get_prevalence.csv")
all_who_mutations <- c("k13:622:I", "k13:469:Y", "k13:675:V", "k13:446:I", "k13:458:Y", "k13:476:I",   "k13:493:H",   "k13:539:T",
                       "k13:543:T",  "k13:553:L",   "k13:561:H",   "k13:574:L",  "k13:580:Y",
                       "k13:441:L", "k13:449:A",   "k13:469:F",   "k13:481:V",
                       "k13:515:K", "k13:527:H",  "k13:537:I", "k13:537:D", "k13:538:V",  "k13:568:G")


shape_Africa <- readRDS("R_ignore/R_scripts/data/shapefiles/sf_admin0_africa.rds")
shape_water  <- sf::st_read("R_ignore/R_scripts/data/shapefiles/africa_water_bodies.shp", quiet = TRUE)

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

for (mut in all_who_mutations){
  dat_sub <- dat |>
    filter(mutation == mut) |>
    filter((latitude  > -5) & (latitude  < 20)) |>
    filter((longitude > 25) & (longitude < 40)) |>
    filter(is.finite(numerator), is.finite(denominator), denominator > 0)
  
  pts_pos <- dat_sub |> 
    dplyr::filter(prevalence > 5) 
  
  if (dim(pts_pos)[1] != 0){
    pts_pos <- pts_pos |>
      sf::st_as_sf(coords = c("longitude","latitude"), crs = 4326, remove = FALSE)
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
    p_post_lo <- cached$p_lower
    p_post_hi <- cached$p_upper
    p_post_mean <- cached$p_mean
    xs <- cached$xs
    ys <- cached$ys
    exceed_prob_10 <- cached$exceed_prob_10_perc
    exceed_prob_5 <- cached$exceed_prob_5_perc
    exceed_prob_1 <- cached$exceed_prob_1_perc
    
    # --- Build long data for lower / mean / upper ------------------------------
    p_long_lo   <- make_raster_long(p_post_lo,   xs, ys, plot_times)
    p_long_mean <- make_raster_long(p_post_mean, xs, ys, plot_times)
    p_long_hi   <- make_raster_long(p_post_hi,   xs, ys, plot_times)
    exceed_prob_long_10 <- make_raster_long(exceed_prob_10, xs, ys, plot_times)
    exceed_prob_long_5 <- make_raster_long(exceed_prob_5, xs, ys, plot_times)
    exceed_prob_long_1 <- make_raster_long(exceed_prob_1, xs, ys, plot_times)
    
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
    plot_theme <- function(p) {
      p + theme(
        strip.background = element_rect(fill = "white", color = NA),
        strip.text = element_text(color = "black", face = "plain"),
        panel.border = element_rect(color = "grey80", fill = NA),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        title = element_text(size = 14),
        axis.text.x = element_text(size = 9, angle = 30, hjust = 1),  # <-- key tweak
        axis.text.y = element_text(size = 9)
      )
    }
    
    plot_layer <- function(p_long_df, title_text, shp = shape_Africa, shp_water = shape_water, add_points = FALSE, add_legend = FALSE) {
      p <- ggplot() +
        geom_raster(aes(x = x, y = y, fill = p*100), data = p_long_df) +
        geom_sf(data = shp,       linewidth = 0.2, fill = NA, color = "white") +
        geom_sf(data = shp_water, fill = "white", colour = NA) +
        coord_sf(xlim = xlim, ylim = ylim, expand = FALSE)
      if (add_points) {
        p <- p + geom_point(
          aes(x = longitude, y = latitude, fill = p_obs*100),
          data = points_df %>% filter(t>2013&t<2026) %>% arrange(p_obs),
          shape = 21, colour = "grey60", size = 1
        )
      }
      if (!add_legend){
        p <- p + theme(legend.position = "none")
      }
      p <- p + scale_fill_gradientn(colours = pp_cols,
                                    values  = pp_vals,
                                    limits = c(0, 100), 
                                    name = "Prevalence (%)") +
        facet_wrap(~ t, nrow = 2) +
        labs(x = "Longitude", y = "Latitude", title = title_text)
      plot_theme(p)
    }
    
    plot_exceedance <-function(exceed_prob, title_text, legend_title, shp = shape_Africa, shp_water = shape_water, add_points = FALSE, add_legend = FALSE) {
      p <- ggplot() +
        theme_bw() +
        geom_raster(aes(x = x, y = y, fill = p), data = exceed_prob) +
        geom_sf(data = shape_Africa, linewidth = 0.2, fill = NA, color = "white") +
        geom_sf(data = shape_water, fill = "white", colour = NA) +
        coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
        geom_point(
          aes(x = longitude, y = latitude, fill = p_obs * 100),
          data = points_df %>% filter(t > 2013 & t < 2026) %>% arrange(p_obs),
          shape = 21, colour = "grey60", size = 1
        ) +
        scale_fill_viridis_c(
          option = "magma",
          name   = legend_title,
          limits = c(0, 1)
        ) +
        facet_wrap(~t, nrow = 2) +
        labs(
          x = "Longitude",
          y = "Latitude",
          title = title_text
        ) 
      plot_theme(p)
    }
    
    # --- Create and print the three separate plots -----------------------------
    
    plot_lower <- plot_layer(p_long_lo, paste0("Lower (2.5%) prevalence ", mut), shape_Africa, shape_water, add_points = TRUE,  add_legend = TRUE)
    plot_mean  <- plot_layer(p_long_mean, paste0("Mean prevalence ", mut), shape_Africa, shape_water, add_points = TRUE, add_legend = TRUE)
    plot_upper <- plot_layer(p_long_hi,   paste0("Upper (97.5%) prevalence ", mut), shape_Africa, shape_water, add_points = TRUE, add_legend = TRUE)
    
    plot_exceed_1 <- plot_exceedance(exceed_prob_long_1, title_text = "Exceedance Probability (Prevalence > 1%)", legend_title ="Pr(prevalence > 1%)")
    plot_exceed_5 <- plot_exceedance(exceed_prob_long_5, title_text = "Exceedance Probability (Prevalence > 5%)", legend_title ="Pr(prevalence > 5%)")
    plot_exceed_10 <- plot_exceedance(exceed_prob_long_10, title_text = "Exceedance Probability (Prevalence > 10%)", legend_title ="Pr(prevalence > 10%)")
    
    save_figs(file.path(OUT_PREV_DIR, paste0(mut, "_lower_perc_prev")), plot_lower, width = 10)
    save_figs(file.path(OUT_PREV_DIR, paste0(mut, "_mean_perc_prev")), plot_mean, width = 10)
    save_figs(file.path(OUT_PREV_DIR, paste0(mut, "_upper_perc_prev")), plot_upper, width = 10)
    
    save_figs(file.path(OUT_EXCEED_DIR, paste0(mut, "_exceednace_prob_1")), plot_exceed_1, width = 10)
    save_figs(file.path(OUT_EXCEED_DIR, paste0(mut, "_exceednace_prob_5")), plot_exceed_5, width = 10)
    save_figs(file.path(OUT_EXCEED_DIR, paste0(mut, "_exceednace_prob_10")), plot_exceed_10, width = 10)
    
    print(paste0("Saved figures for ", mut)) 
  }
  else{
    print(paste0(mut, " has no pos mutations")) 
  }
}

