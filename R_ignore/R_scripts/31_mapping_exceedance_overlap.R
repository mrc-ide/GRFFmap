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

cache_path <- function(mut, lenS, lenT, what, ext = c("rds","parquet")) {
  ext <- match.arg(ext)
  fn  <- sprintf("%s_lenS%s_lenT%s_%s.%s",
                 gsub("[:/ ]", "_", mut), lenS, lenT, what, ext)
  file.path(CACHE_DIR, fn)
}

make_raster_long <- function(arr3, xs, ys, plot_times) {
  nx <- length(xs); ny <- length(ys)
  df <- do.call(rbind, lapply(seq_along(plot_times), function(k) {
    data.frame(
      x = rep(xs, times = ny),
      y = rep(ys, each  = nx),
      p = as.vector(arr3[ , , k]),
      t = plot_times[k]
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

avg_2_year <- function(df){
  # Create a “year block” variable
  df <- df %>%
    mutate(
      year_block = case_when(
        t %in% c(2002, 2003) ~ "2002-2003",
        t %in% c(2004, 2005) ~ "2004-2005",
        t %in% c(2006, 2007) ~ "2006-2007",
        t %in% c(2008, 2009) ~ "2008-2009",
        t %in% c(2010, 2011) ~ "2010-2011",
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
#ell_km and p_init below within for loop
# inference parameters
D        <- 500         # number of random frequencies (try 200–500)
max_iter <- 5           # EM iterations
z_eps    <- 1e-6        # threshold for |z| to use n/4 in E[ω]
omega_floor <- 1e-10    # floor on ω to avoid 1/ω explosions
jitter_S <- 1e-10       # diagonal jitter for innovation covariance

# prediction parameters
nx <- 200
ny <- 200
t_vec <- 2012:2023
t_num <- length(t_vec)

# plotting parameters
plot_times <- seq(2012, 2023, by = 1) # should be within t_vec
t_win <- 2           # window around each facet time (years). Points within this window are displayed

# --------------------------- Load & filter data ----------------------------
CACHE_DIR <- "R_ignore/R_scripts/outputs/GRFF_kalman_PDcache_annual_2012_2025"
dir_create(CACHE_DIR)
OUT_PREV_DIR   <- file.path("prev_prediction_grouped_year_kalman")
#OUT_PREV_DIR_LO_UP <- file.path("prev_prediction_grouped_year_kalman/upper_lower_prev")
OUT_EXCEED_DIR <- file.path("exceedance_prob_grouped_year_kalman")

dir_create(c(paste0("R_ignore/R_scripts/outputs/plots/", OUT_PREV_DIR), 
             paste0("R_ignore/R_scripts/outputs/plots/", OUT_EXCEED_DIR)),
           recurse = TRUE)

dat <- read.csv("R_ignore/R_scripts/data/partner_drug_calc_get_prevalence.csv") |>
  mutate(collection_day = as.Date(collection_day)) |>
  select(longitude, latitude, year, numerator, denominator, prevalence, mutation, country_name)

shape_Africa <- readRDS("R_ignore/R_scripts/data/shapefiles/sf_admin0_africa.rds")
shape_water  <- sf::st_read("R_ignore/R_scripts/data/shapefiles/africa_water_bodies.shp", quiet = TRUE)

# --------------------------- Define Mutations ----------------------------
all_who_mutations <- c("mdr1C:86:N", "crt:76:T")
# --- Create white raster for outside Africa ------------------------------
CRS_LL     <- sf::st_crs(4326)  # WGS84 lon/lat for plotting
CRS_METRIC <- sf::st_crs(3857)  # planar for buffers/ops

shape_Africa_ll <- shape_Africa |> sf::st_transform(CRS_LL)
shape_water_ll  <- shape_water  |> sf::st_transform(CRS_LL)

#setup grid - use crt 76 for grid (all of africa ) for all
mut = "crt:76:T"
dat_sub <- dat |>
  filter(mutation == mut) |>
  filter(is.finite(numerator), is.finite(denominator), denominator > 0) |>
  mutate(
    longitude = as.numeric(longitude),
    latitude  = as.numeric(latitude)
  )

pts_pos <- dat_sub |> 
  dplyr::filter(prevalence > prev_grid_threshold) 
pts_pos <- sf::st_as_sf(pts_pos,
                        coords = c("longitude", "latitude"),
                        crs = CRS_LL,
                        remove = FALSE)
buf_km <- 200
bbox_ll <- pts_pos |>
  sf::st_transform(CRS_METRIC) |>
  sf::st_bbox() |> sf::st_as_sfc() |>
  sf::st_buffer(buf_km * 1000) |>
  sf::st_transform(CRS_LL)

bb      <- sf::st_bbox(bbox_ll)              # named numeric xmin/xmax/ymin/ymax
stopifnot(all(!is.na(bb)))                   # guard against empties
xlim    <- as.numeric(bb[c("xmin","xmax")])  # plain numerics
ylim    <- as.numeric(bb[c("ymin","ymax")])
xlim <- sort(xlim)
ylim <- sort(ylim)

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
shape_Africa_crop <- st_transform(shape_Africa_planar, 4326)
shape_water_crop  <- st_transform(shape_water_planar,  4326)
# Restore s2
sf_use_s2(old_s2)

africa_mask <- shape_Africa_crop %>%
  sf::st_union() %>%
  sf::st_make_valid()

exceedance_list <- list()

# --- Loop over each mutation ------------------------------
for (mut in all_who_mutations){
  if (mut == "crt:76:T"){clean_mut <- gsub("^crt:(\\d+):([A-Za-z])$", "crt \\1\\2", mut)}
  if (mut == "mdr1C:86:N"){clean_mut <- gsub("^mdr1C:(\\d+):([A-Za-z])$", "mdr1 \\1\\2", mut)}
  
  # --------------------------- Settings (cont.) --------------------------------------
  # model parameters
  if(mut == "mdr1C:86:N"){
    ell_km <- 500       # RFF length-scale in **kilometres** #80 for k13
    tau2   <- 0.1         # RW1 variance in feature space
    p_init <- 0.3       #0.001 for k13
    z_init <- qlogis(p_init)
  }
  if(mut == "crt:76:T"){
    ell_km <- 500       # RFF length-scale in **kilometres** #80 for k13
    tau2   <- 0.1         # RW1 variance in feature space
    p_init <- 0.07       #0.001 for k13
    z_init <- qlogis(p_init)
  }
  if (mut == "k13:comb"){
    #@CECILE -> add K13 comb params to get correct file pull and update cached file location if different
    CACHE_DIR <- "R_ignore/R_scripts/outputs/GRFF_kalman_PDcache_annual_2012_2025"
    
  }
  
  if (dim(pts_pos)[1] != 0){
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
    exceed_prob <- cached$exceed_post_draw_list
    xs <- cached$xs
    ys <- cached$ys
    ##Only pull of exceedance prob 10 for k13 
    exceed_prob_10 <- exceed_prob$`10`
    
    ###@CECILE -> need to add if statments to pull out the exceedancd  50 for mdr and crt
    ##pull out exceedance prob of 50 for crt and mdr
    
    exceed_prob_long_10 <- make_raster_long(exceed_prob_10, xs, ys, plot_times)
    exceed_prob_long_10 <- mask_to_africa(exceed_prob_long_10, africa_mask)
    
    exceedance_list[[mut]] <- exceed_prob_long_10
  }  
}
## extract three exceedance lists:
k13_exceed <- exceedance_list[["k13:comb"]]
mdr_exceed <- exceedance_list[["mdr1C:86:N"]]
crt_exceed <- exceedance_list[["crt:76:T"]]

k13_exceed <- k13_exceed %>% filter( t == 2023)
mdr_exceed <- mdr_exceed %>% filter( t == 2023)
crt_exceed <- crt_exceed %>% filter(t == 2023)

plot_theme <- function(p) {
  p + theme(
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(color = "black", face = "plain"),
    panel.border = element_rect(color = "grey80", fill = NA),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    title = element_text(size = 14),
    axis.text.x = element_text(size = 8, angle = 30, hjust = 1),  # <-- key tweak
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )
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
    geom_raster(aes(x = x, y = y, fill = p*100), data = exceed_prob) +
    geom_sf(data = shape_Africa, linewidth = 0.2, fill = NA, color = "white") +
    geom_sf(data = shp_water, fill = "white", colour = NA) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = st_crs(4326)) +
    scale_fill_viridis_c(
      option = "mako",
      direction = 1,   # optional (flip palette)
      limits = c(0, 100),
      breaks = seq(0, 100, by = 10),
      name = legend_title,
      na.value = "white"
    ) +
    facet_wrap(~t, nrow = 2) +
    labs(
      x = "Longitude",
      y = "Latitude"
    ) +
    theme(legend.position = "bottom",
          legend.key.width = unit(2, "cm"))
  if (title_text != ""){
    p <- p + labs(title = title_text)
  }
  plot_theme(p)
}

plot_exceed_10 <- plot_exceedance(exceed_prob_long_10, mut, legend_title =expression(Pr(Prevalence >= 10*"%")))


##########ADM1 OVERLAP ATTEMPT################################

# read Africa shape files
africa_shp_admin1 <- readRDS(file = "R_ignore/R_scripts/data/sf_admin1_africa.rds")
target_crs <- 4326
africa_shp_admin1 <- st_transform(africa_shp_admin1, crs = target_crs)

# Using the cropped Admin 1 data for East Africa
admin1_regions <- st_make_valid(africa_shp_admin1) 

summarize_exceedance_adm1 <- function(exceed_df, adm1_shp, threshold_prob = NULL) {
  # exceed_df must contain x, y, t, p (p in 0–1)
  
  pts_sf <- st_as_sf(
    exceed_df,
    coords = c("x", "y"),
    crs = 4326,
    remove = FALSE
  )
  sf_use_s2(FALSE)
  #THIS STEP IS V SLOW if you try to do all time points 
  joined <- st_join(pts_sf, admin1_regions)  # spatial overlay

  summary <- joined |>
    st_drop_geometry() |>
    dplyr::group_by(id_0, id_1, t) |>
    dplyr::summarise(
      mean_exceed = mean(p, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(summary)
}

mdr_adm1 <- summarize_exceedance_adm1(mdr_exceed, admin1_regions)
crt_adm1 <- summarize_exceedance_adm1(crt_exceed, admin1_regions)

adm1_wide <- k13_adm1 |>
  rename(mean_k13 = mean_exceed) |>
  inner_join(mdr_adm1 |> rename(mean_mdr = mean_exceed),
             by = c("id_0","id_1","t")) |>
  inner_join(crt_adm1 |> rename(mean_crt = mean_exceed),
             by = c("id_0","id_1","t"))

adm1_overlap <- adm1_wide |>
  mutate(
    all_high = (mean_k13 > 0.90 &
                  mean_mdr > 0.90 &
                  mean_crt > 0.90)
  )

adm1_wide <- mdr_adm1 |>
  rename(mean_mdr = mean_exceed) |>
  inner_join(crt_adm1 |> rename(mean_crt = mean_exceed),
             by = c("id_0","id_1","t"))

adm1_overlap <- adm1_wide |>
  mutate(
    all_high = (mean_mdr > 0.60 &
                  mean_crt > 0.60)
  ) %>%drop_na()

plot_adm1_overlap <- function(adm1_df, adm1_shp) {
  
  shp_plot <- adm1_shp |>
    left_join(adm1_df, by = c("id_0","id_1")) %>% drop_na(t)
  
  p <- ggplot(shp_plot) +
    geom_sf(aes(fill = all_high), color = "white", size = 0.2) +
    geom_sf(data = shape_Africa_ll, linewidth = 0.2, fill = NA, color = "black") +
    geom_sf(data = shape_Africa_crop, linewidth = 0.2, fill = NA, color = "grey40") +
    
    #geom_sf(data = shape_water_crop, fill = "white", colour = NA) +
    scale_fill_manual(
      values = c("FALSE" = "white", "TRUE" = "grey20"),
      breaks = c("TRUE"),                              # Only show TRUE in the legend
      labels = c("ADM1 regions with\n>60% mean exceedance\n(all markers)"), # The label text
      name = NULL 
    ) +
    facet_wrap(~ t) +
    theme_bw() +
    labs(title = "High-Confidence Multi-Locus Exceedance at ADM1 Level",
         x = "Longitude", y = "Latitude") +
    theme(legend.position = "bottom")
  plot_theme(p)
}

plot_adm1_overlap(adm1_overlap, admin1_regions)



# You must install and load this library first
# install.packages("ggpattern")
library(ggpattern)

plot_adm1_overlap <- function(adm1_df, adm1_shp) {
  
  shp_plot <- adm1_shp |>
    left_join(adm1_df, by = c("id_0","id_1")) %>% 
    drop_na(t)
  
  p <- ggplot(shp_plot) +
    
    # Use geom_sf_pattern instead of geom_sf
    geom_sf_pattern(
      aes(pattern = all_high),     # Map the TRUE/FALSE column to the pattern aesthetic
      fill            = "white",   # Background color of the polygons (white for all)
      color           = "white",   # Border color of the polygons
      pattern_color   = "grey20",  # Color of the actual hash lines
      pattern_fill    = "white",   # distinct from fill, often needed for specific pattern types
      pattern_angle   = 45,        # Angle of the hash lines
      pattern_density = 0.1,       # How thick the lines are (0 to 1)
      pattern_spacing = 0.02,      # Spacing between lines (adjust based on map scale)
      size            = 0.2
    ) +
    
    # Overlay the continent outline
    geom_sf(data = shape_Africa_crop, linewidth = 0.2, fill = NA, color = "grey40") +
    
    # Custom scale to handle the specific legend request
    scale_pattern_manual(
      values = c("FALSE" = "none", "TRUE" = "stripe"), # FALSE gets no pattern, TRUE gets hash
      breaks = c("TRUE"),                              # Only show TRUE in the legend
      labels = c("ADM1 regions with\n>60% mean exceedance\n(all markers)"), # The label text
      name = NULL                                      # Removes the legend title
    ) +
    
    facet_wrap(~ t) +
    theme_bw() +
    labs(title = "High-Confidence Multi-Locus Exceedance at ADM1 Level",
         x = "Longitude", y = "Latitude") +
    
    # Ensure the legend box is large enough to see the pattern
    theme(
      legend.position = "bottom",
      legend.key.size = unit(1.5, "cm") 
    )
  
  plot_theme(p)
}


####STRICT POINT WISE OVERLAP ATTEMPT#######################################
compute_joint_exceedance <- function(k13_df, mdr_df, crt_df) {
  
  # Join all three by spatial/time coordinates
  joint <- k13_df |>
    dplyr::rename(p_k13 = p) |>
    dplyr::inner_join(mdr_df |> dplyr::rename(p_mdr = p),
                      by = c("x","y","t")) |>
    dplyr::inner_join(crt_df |> dplyr::rename(p_crt = p),
                      by = c("x","y","t"))
  
  # Identify full-overlap pixels
  joint <- joint |>
    dplyr::mutate(
      overlap = (p_k13 == 1 & p_mdr == 1 & p_crt == 1)
    )
  
  return(joint)
}


compute_joint_exceedance_2 <- function(mdr_df, crt_df) {
  # Join all three by spatial/time coordinates
  joint <- mdr_df |>
    dplyr::rename(p_mdr = p)|>
    dplyr::inner_join(crt_df |> dplyr::rename(p_crt = p),
                      by = c("x","y","t"))
  
  # Identify full-overlap pixels
  joint <- joint |>
    dplyr::mutate(
      overlap = (p_mdr == 1 & p_crt == 1)
    )
  
  return(joint)
}

overlap <- compute_joint_exceedance_2(mdr_exceed, crt_exceed)
plot_joint_overlap <- function(overlap_df, shp = shape_Africa, shp_water = shape_water) {
  
  ggplot() +
    theme_bw() +
    annotate(
      "rect", xmin = xlim[1], xmax = xlim[2],
      ymin = ylim[1], ymax = ylim[2],
      fill = "white", colour = NA
    ) +
    geom_raster(
      data = overlap_df,
      aes(x = x, y = y, fill = overlap)
    ) +
    geom_sf(data = shp, linewidth = 0.2, fill = NA, color = "white") +
    geom_sf(data = shp_water, fill = "white", colour = NA) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = st_crs(4326)) +
    scale_fill_manual(
      values = c("FALSE" = "white", "TRUE" = "grey20"),
      name = "Full Overlap\n(100%)"
    ) +
    facet_wrap(~t, nrow = 2) +
    labs(x = "Longitude", y = "Latitude",
         title = "Regions With 100% Exceedance Across K13, mdr1, and crt") +
    theme(legend.position = "bottom") |>
    plot_theme()
}

plot_joint_overlap(overlap)

