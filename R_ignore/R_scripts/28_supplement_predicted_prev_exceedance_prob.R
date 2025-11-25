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
  library(cowplot)
  library(scales)
  library(matrixStats)
})

load_all()

set.seed(1)

cache_path <- function(mut, lenS, lenT, what, ext = c("rds","parquet")) {
  ext <- match.arg(ext)
  fn  <- sprintf("%s_lenS%s_lenT%s_%s.%s",
                 gsub("[:/ ]", "_", mut), lenS, lenT, what, ext)
  file.path(CACHE_DIR, fn)
}

make_raster_long <- function(arr3, xs, ys, t_vec, plot_times) {
  nx <- length(xs); ny <- length(ys)
  df <- do.call(rbind, lapply(seq_along(t_vec), function(k) {
    data.frame(
      x = rep(xs, times = ny),
      y = rep(ys, each  = nx),
      p = as.vector(arr3[ , , k]),
      t = t_vec[k]
    )
  }))
  df %>%
    dplyr::filter(t %in% plot_times) %>%
    dplyr::mutate(t = factor(t, levels = plot_times))
}

make_or_load <- function(mut, ell_km, tau2) {
  cached_rds <- cache_path(mut, ell_km, tau2, "full_model_output", "rds") 
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

posterior_se_fast <- function(draw_list) {
  # draw_list: list of identical-dimension arrays/matrices
  
  # 1. Stack draws along last dimension
  arr <- simplify2array(draw_list)        # dims: (d1, d2, ..., dK, D)
  dims <- dim(arr)
  D <- dims[length(dims)]
  
  # 2. Reshape into a matrix: pixels × draws
  mat <- matrix(arr, 
                nrow = prod(dims[-length(dims)]),
                ncol = D)
  
  # 3. Compute SD row-wise very efficiently
  sd_vec <- rowSds(mat)
  
  # 4. Standard error = SD / sqrt(D)
  se_vec <- sd_vec / sqrt(D)
  
  # 5. Reshape back to original spatial shape
  array(se_vec, dim = dims[-length(dims)])
}

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
    legend.justification = "center",
    panel.spacing = unit(0, "lines"),
    plot.margin = margin(t = 1, r = 2, b = 1, l = 2)
  )
}

plot_layer_combined <- function(p_long_df, title_text, shp = shape_Africa, shp_water = shape_water, add_points = FALSE, add_legend = TRUE) {
  p <- ggplot() +
    geom_raster(aes(x = x, y = y, fill = p*100), data = p_long_df) +
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
  p <- plot_theme(p)
  if (add_points) {
    p <- p + geom_point(
      aes(x = longitude, y = latitude, fill = p),  # p is already % units
      data  = points_df %>% arrange(p),
      shape = 21, colour = "darkgrey", size = 1.5, stroke = 0.2,
      inherit.aes = FALSE
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
    facet_wrap(~ t, nrow = 2) +
    labs(x = "Longitude", y = "Latitude")
  p
}

plot_exceedance_combined <-function(exceed_prob, title_text, legend_title, shp = shape_Africa, shp_water = shape_water, add_points = FALSE, add_legend = FALSE) {
  p <- ggplot() +
    theme_bw() +
    annotate(
      "rect",
      xmin = xlim[1], xmax = xlim[2],
      ymin = ylim[1], ymax = ylim[2],
      fill = "white", colour = NA
    ) +
    geom_raster(aes(x = x, y = y, fill = p*100), data = exceed_prob) +
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
    facet_wrap(~t, nrow = 2) +
    labs(
      x = "Longitude",
      y = "Latitude"
    ) 
  p <- plot_theme(p)
  p <- p + theme(legend.position = "bottom",
                 legend.key.width = unit(2, "cm"),
                 legend.key.height = unit(0.2, "cm")
  )
  if (title_text != ""){
    p <- p + labs(title = title_text)
  }
  p
}

plot_layer_CI_combined <- function(p_long_df, title_text, shp = shape_Africa, shp_water = shape_water, add_legend = TRUE) {
  p <- ggplot() +
    geom_raster(aes(x = x, y = y, fill = p*100), data = p_long_df) +
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
  p <- plot_theme(p)
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
  p <- p + scale_fill_viridis_c(
    option = "magma",
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10),
    name = "Prevalence (%)",
    na.value = "white"
  ) +
    facet_wrap(~ t, nrow = 2) +
    labs(x = "Longitude", y = "Latitude")
  p
}

# --------------------------- Settings --------------------------------------
# model parameters
ell_km <- 120          # RFF length-scale in **kilometres**
tau2   <- 0.5        # RW1 variance in feature space

# prediction parameters
nx <- 200
ny <- 200
t_vec <- 1995:2025
t_num <- length(t_vec)

# plotting parameters
plot_times <- seq(2012, 2023, by = 1) # should be within t_vec

# --------------------------- Load & filter data ----------------------------
CACHE_DIR <- "R_ignore/R_scripts/outputs/supplemental/GRFF_kalman_cache_annual_2012_2025"

OUT_COMBINED <- file.path("supplemental", "combined_pred_prev_exceedance_prob_grouped_year_kalman")
dir_create(c(paste0("R_ignore/R_scripts/outputs/plots/", OUT_COMBINED)), 
           recurse = TRUE)

dat <- read.csv("R_ignore/R_scripts/data/all_who_get_prevalence.csv") |>
  mutate(collection_day = as.Date(collection_day)) |>
  select(longitude, latitude, year, numerator, denominator, prevalence, mutation, country_name)

shape_Africa <- readRDS("R_ignore/R_scripts/data/shapefiles/sf_admin0_africa.rds")
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

# --- Create white raster for outside Africa -----------------------------------
CRS_LL     <- sf::st_crs(4326)  # WGS84 lon/lat for plotting
CRS_METRIC <- sf::st_crs(3857)  # planar for buffers/ops

shape_Africa_ll <- shape_Africa |> sf::st_transform(CRS_LL)
shape_water_ll  <- shape_water  |> sf::st_transform(CRS_LL)

# --- Crop Africa polygons -----------------------------------------------------

# Define East Africa box and crop Admin1 shapefile
bbox_east_africa <- st_bbox(c(
  xmin = 28,
  xmax = 48,
  ymin = -4.6,
  ymax = 18
), crs = 4326)

xlim <- c(bbox_east_africa["xmin"], bbox_east_africa["xmax"])
ylim <- c(bbox_east_africa["ymin"], bbox_east_africa["ymax"])

lon_min <- xlim[1] # bbox_east_africa["xmin"]
lon_max <- xlim[2] # bbox_east_africa["xmax"]
lat_min <- ylim[1] # bbox_east_africa["ymin"]
lat_max <- ylim[2] # bbox_east_africa["ymax"]

bbox_sfc_ll <- st_as_sfc(bbox_east_africa)

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


# --------------------------- Define Mutations ----------------------------
all_who_mutations <- c("k13:comb", "k13:675:V", "k13:622:I", "k13:469:Y", "k13:446:I", "k13:458:Y", "k13:476:I",   "k13:493:H",   "k13:539:T",
                       "k13:543:T",  "k13:553:L",   "k13:561:H",   "k13:574:L",  "k13:580:Y",
                       "k13:441:L", "k13:449:A",   "k13:469:F",   "k13:481:V",
                       "k13:515:K", "k13:527:H",  "k13:537:I", "k13:537:D", "k13:538:V",  "k13:568:G")

# Note for Cecile: the following mutations have a wide grid where prev > 0:
# 441L, 449A, 469F, 473I, 539T, 553L, 561H, 568G, 574L, 622I

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
for (mut in all_who_mutations){
  clean_mut <- paste0("k13 ", gsub("^k13:(\\d+):([A-Za-z])$", "\\1\\2", mut))
  
  dat_sub <- dat_with_k13 |>
    filter(mutation == mut) |>
    filter(is.finite(numerator), is.finite(denominator), denominator > 0) |>
    mutate(
      longitude = as.numeric(longitude),
      latitude  = as.numeric(latitude)
    ) |>
  filter(
    between(longitude, lon_min, lon_max),
    between(latitude, lat_min, lat_max)
  )
  
  pts_pos <- dat_sub |> 
    dplyr::filter(prevalence > 0) 
  
  if (dim(pts_pos)[1] != 0){
    # --------------------------- Create plotting grid ---------------
    pts_pos <- sf::st_as_sf(pts_pos,
                            coords = c("longitude", "latitude"),
                            crs = CRS_LL,
                            remove = FALSE)
    xlim <- sort(xlim)
    ylim <- sort(ylim)
    
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
    p_long_median <- make_raster_long(p_post_median, xs, ys, t_vec, plot_times)
    p_long_mean <- make_raster_long(p_post_mean, xs, ys, t_vec, plot_times)
    p_long_CI <- make_raster_long(p_post_CI, xs, ys, t_vec, plot_times)
    
    exceed_prob_long_10 <- make_raster_long(exceed_prob_10, xs, ys, t_vec, plot_times)
    exceed_prob_long_5 <- make_raster_long(exceed_prob_5, xs, ys, t_vec, plot_times)
    exceed_prob_long_1 <- make_raster_long(exceed_prob_1, xs, ys, t_vec, plot_times)
    
    # --- Mask out sea as white background ------------------------------
    p_long_median <- mask_to_africa(p_long_median, africa_mask)
    p_long_mean <- mask_to_africa(p_long_mean, africa_mask)
    p_long_CI <- mask_to_africa(p_long_CI, africa_mask)
    
    exceed_probdraws_long_10 <- mask_to_africa(exceed_prob_long_10, africa_mask)
    exceed_prob_long_5  <- mask_to_africa(exceed_prob_long_5, africa_mask)
    exceed_prob_long_1  <- mask_to_africa(exceed_prob_long_1, africa_mask) 
    
    # --- Crop dataframes to bbox ------------------------------
    p_long_median   <- crop_long_df(p_long_median, xlim, ylim)
    p_long_mean   <- crop_long_df(p_long_mean, xlim, ylim)
    p_long_CI   <- crop_long_df(p_long_CI, xlim, ylim)
    
    exceed_prob_long_10 <- crop_long_df(exceed_prob_long_10, xlim, ylim)
    exceed_prob_long_5  <- crop_long_df(exceed_prob_long_5,  xlim, ylim)
    exceed_prob_long_1  <- crop_long_df(exceed_prob_long_1,  xlim, ylim)
    
    # Rename t to year
    # p_long_median <- p_long_median %>% rename(year = t)
    # p_long_mean <- p_long_mean %>% rename(year = t)
    # exceed_prob_long_1 <- exceed_prob_long_1 %>% rename(year = t)
    # exceed_prob_long_5 <- exceed_prob_long_5 %>% rename(year = t)
    # exceed_prob_long_10 <- exceed_prob_long_10 %>% rename(year = t)
    # 
    # p_long_median$year       <- factor(p_long_median$year, levels = plot_times)
    # p_long_mean$year         <- factor(p_long_mean$year, levels = plot_times)
    # exceed_prob_long_1$year  <- factor(exceed_prob_long_1$year, levels = plot_times)
    # exceed_prob_long_5$year  <- factor(exceed_prob_long_5$year, levels = plot_times)
    # exceed_prob_long_10$year <- factor(exceed_prob_long_10$year, levels = plot_times)
    
    # --- Create observed data points ------------------------------
    points_df <- dat_sub %>%
      filter(year %in% plot_times) %>%
      mutate(
        p_obs = numerator / denominator,
        year  = factor(year, levels = plot_times),
        p     = p_obs * 100
      ) %>%
      filter(between(longitude, xlim[1], xlim[2]),
             between(latitude,  ylim[1], ylim[2])) %>%
      transmute(
        longitude, latitude, year, p_obs, prevalence, p
      )
    
    # --- Create and print the three separate plots -----------------------------
    
    plot_median_no_title_comb  <- plot_layer_combined(p_long_median, "", shape_Africa_crop, shape_water_crop, add_points = FALSE, add_legend = TRUE)
    plot_exceed_5_no_title_comb <- plot_exceedance_combined(exceed_prob_long_5, "", legend_title =expression(Pr(Prevalence >= 5*"%")))
    plot_CI_no_title_comb <- plot_layer_CI_combined(p_long_CI, "", shape_Africa_crop, shape_water_crop, add_legend = TRUE)
      
    prev_pred_exceed_5perc <- 
      ggdraw() +
      draw_label(
        clean_mut,
        fontface = "plain",
        size = 18,
        x = 0.5,
        y = 0.98,
        hjust = 0.5
      ) +
      draw_plot(
        cowplot::plot_grid(
          plot_median_no_title_comb,
          plot_CI_no_title_comb, 
          plot_exceed_5_no_title_comb,
          ncol = 1,
          labels = c("A", "B", "C"),
          label_size = 18,
          label_x = 0.02,
          label_y = 0.98,
          rel_heights = c(1, 1, 1)
        ),
        x = 0,
        y = 0,
        width = 1,
        height = 0.95
      )
    save_figs(file.path(OUT_COMBINED, paste0(mut, "_pred_prev_exceednace_5perc_no_title")), prev_pred_exceed_5perc, height = 14, width = 7)
    
    print(paste0("Saved figures for ", mut)) 
  }
  else{
    print(paste0(mut, " has no pos mutations")) 
  }
}

