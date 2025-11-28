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
  
  # Temporarily turn off s2 and ensure geometry is valid
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  africa_mask_valid <- sf::st_make_valid(africa_mask)
  
  inside <- lengths(sf::st_intersects(pts, africa_mask_valid)) > 0
  df$p[!inside] <- NA_real_
  df
}

avg_2_year <- function(df){
  # Create a “year block” variable
  df %>%
    mutate(
      year_block = case_when(
        t %in% c(2001, 2002) ~ "2001-2002",
        t %in% c(2003, 2004) ~ "2003-2004",
        t %in% c(2005, 2006) ~ "2005-2006",
        t %in% c(2007, 2008) ~ "2007-2008",
        t %in% c(2009, 2010) ~ "2009-2010",
        t %in% c(2011, 2012) ~ "2011-2012",
        t %in% c(2013, 2014) ~ "2013-2014",
        t %in% c(2015, 2016) ~ "2015-2016",
        t %in% c(2017, 2018) ~ "2017-2018",
        t %in% c(2019, 2020) ~ "2019-2020",
        t %in% c(2021, 2022) ~ "2021-2022",
        t %in% c(2023, 2024) ~ "2023-2024",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(year_block)) %>%              # <- drop NA here
    group_by(x, y, year_block) %>%
    summarise(mean_p = mean(p, na.rm = TRUE), .groups = "drop")
}

add_year_block <- function(df) {
  df %>%
    mutate(
      year_block = case_when(
        t %in% c(2001, 2002) ~ "2001-2002",
        t %in% c(2003, 2004) ~ "2003-2004",
        t %in% c(2005, 2006) ~ "2005-2006",
        t %in% c(2007, 2008) ~ "2007-2008",
        t %in% c(2009, 2010) ~ "2009-2010",
        t %in% c(2011, 2012) ~ "2011-2012",
        t %in% c(2013, 2014) ~ "2013-2014",
        t %in% c(2015, 2016) ~ "2015-2016",
        t %in% c(2017, 2018) ~ "2017-2018",
        t %in% c(2019, 2020) ~ "2019-2020",
        t %in% c(2021, 2022) ~ "2021-2022",
        t %in% c(2023, 2024) ~ "2023-2024"
      )
    ) %>%
    filter(!is.na(year_block))
}

# --------------------------- Load & filter data ----------------------------
CACHE_DIR <- "R_ignore/R_scripts/outputs/GRFF_kalman_partner_drug_cache_annual_1995_2024/"
OUT_PREV_DIR   <- file.path("partner_drug_prev_prediction_grouped_year_kalman")
#OUT_PREV_DIR_LO_UP <- file.path("prev_prediction_grouped_year_kalman/upper_lower_prev")
OUT_EXCEED_DIR <- file.path("partner_drug_exceedance_prob_grouped_year_kalman")

OUT_MEAN_DIR   <- file.path("partner_drug_mutations", "mean_prev_prediction_avg_2_year")
OUT_MEDIAN_DIR   <- file.path("partner_drug_mutations", "median_prev_prediction_avg_2_year")
OUT_CI_DIR   <- file.path("partner_drug_mutations", "CI_diff_prediction_avg_2_year")
OUT_EXCEED_DIR <- file.path("partner_drug_mutations", "exceedance_prob_avg_2_year")

dir_create(c(paste0("R_ignore/R_scripts/outputs/plots/", OUT_MEAN_DIR),
             paste0("R_ignore/R_scripts/outputs/plots/", OUT_MEDIAN_DIR),
             paste0("R_ignore/R_scripts/outputs/plots/", OUT_CI_DIR),
             paste0("R_ignore/R_scripts/outputs/plots/", OUT_EXCEED_DIR)),
           recurse = TRUE)
dat <- read.csv("R_ignore/R_scripts/data/partner_drug_calc_get_prevalence.csv") |>
  mutate(collection_day = as.Date(collection_day)) |>
  select(longitude, latitude, year, numerator, denominator, prevalence, mutation, country_name)

shape_Africa <- readRDS("R_ignore/R_scripts/data/shapefiles/sf_admin0_africa.rds")
shape_water  <- sf::st_read("R_ignore/R_scripts/data/shapefiles/africa_water_bodies.shp", quiet = TRUE)

# --------------------------- Define Mutations ----------------------------
all_who_mutations <- c("mdr1C:86:N", "crt:76:T")

ell_km <- 500       # RFF length-scale in **kilometres** #80 for k13
tau2   <- 0.1

# --- Create white raster for outside Africa ------------------------------
bbox_info <- build_africa_bbox_and_crop(shape_Africa, shape_water, buf_km = 200)
xlim        <- bbox_info$xlim
ylim        <- bbox_info$ylim
africa_mask <- bbox_info$africa_mask
shape_Africa_crop <- bbox_info$africa_crop
shape_water_crop  <- bbox_info$water_crop
africa_mask_land <- shape_Africa_crop %>%
  sf::st_union() %>%
  sf::st_make_valid()

CRS_LL     <- sf::st_crs(4326)
CRS_METRIC <- sf::st_crs(3857) 

plot_times <- seq(2001, 2024) # should be within t_vec

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

# --- Countries that will be greyed out ------------------------------
grey_countries <- c(
  "Egypt", "Morocco", "Libya", "Tunisia", "Algeria",
  "Cabo Verde", "Lesotho", "Mauritius", "Seychelles"
)

shape_grey <- shape_Africa %>%
  dplyr::filter(name_0 %in% grey_countries)

# --- Loop over each mutation ------------------------------
for (mut in all_who_mutations){
  if (mut == "crt:76:T"){clean_mut <- gsub("^crt:(\\d+):([A-Za-z])$", "crt \\1\\2", mut)}
  if (mut == "mdr1C:86:N"){clean_mut <- gsub("^mdr1C:(\\d+):([A-Za-z])$", "mdr1 \\1\\2", mut)}
  
  # --- Filter data --------------------------------------------------------
  dat_sub <- dat |>
    filter(mutation == mut) |>
    filter(is.finite(numerator), is.finite(denominator), denominator > 0) |>
    mutate(
      longitude = as.numeric(longitude),
      latitude  = as.numeric(latitude)
    )
    
  # --- Build long data for lower / mean / upper ------------------------------
  cached <- make_or_load(mut, ell_km, tau2)
  
  p_post_mean <- cached$p_post_mean_draw
  
  xs <- cached$xs
  ys <- cached$ys
  exceed_prob <- cached$exceed_post_draw_list
  
  #exceed_prob_1 <- exceed_prob$`1`
  exceed_prob_5 <- exceed_prob$`5`
  exceed_prob_10 <- exceed_prob$`10`
  exceed_prob_20 <- exceed_prob$`20`
  exceed_prob_50 <- exceed_prob$`50`
  
  # --- Build long data for lower / mean / upper ------------------------------
  p_long_mean <- make_raster_long(p_post_mean, xs, ys, plot_times)
  
  exceed_prob_long_50 <- make_raster_long(exceed_prob_50, xs, ys, plot_times)
  exceed_prob_long_20 <- make_raster_long(exceed_prob_20, xs, ys, plot_times)
  exceed_prob_long_10 <- make_raster_long(exceed_prob_10, xs, ys, plot_times)
  exceed_prob_long_5 <- make_raster_long(exceed_prob_5, xs, ys, plot_times)
  
  # --- Mask out sea as white background ------------------------------
  p_long_mean <- mask_to_africa(p_long_mean, africa_mask_land)
  
  exceed_prob_long_50 <- mask_to_africa(exceed_prob_long_50, africa_mask_land)
  exceed_prob_long_20 <- mask_to_africa(exceed_prob_long_20, africa_mask_land)
  exceed_prob_long_10 <- mask_to_africa(exceed_prob_long_10, africa_mask_land)
  exceed_prob_long_5  <- mask_to_africa(exceed_prob_long_5,  africa_mask_land)

  # # --- Crop dataframes to bbox ------------------------------
  # p_long_mean   <- crop_long_df(p_long_mean,   xlim, ylim)
  # 
  # exceed_prob_long_50 <- crop_long_df(exceed_prob_long_50, xlim, ylim)
  # exceed_prob_long_20 <- crop_long_df(exceed_prob_long_20, xlim, ylim)
  # exceed_prob_long_10 <- crop_long_df(exceed_prob_long_10, xlim, ylim)
  # exceed_prob_long_5  <- crop_long_df(exceed_prob_long_5,  xlim, ylim)
  # #exceed_prob_long_1  <- crop_long_df(exceed_prob_long_1,  xlim, ylim)
  
  # --- Create color plot ------------------------------
  points_df <- dat_sub %>%
    filter(year >= 2000, year <= 2024) %>%     # keep both years in each block
    mutate(
      p_obs = numerator / denominator,
      t     = factor(year, levels = plot_times)
    ) %>%
    filter(longitude >= xlim[1], longitude <= xlim[2],
           latitude  >= ylim[1],  latitude  <= ylim[2]) %>%
    transmute(
      x = longitude,
      y = latitude, 
      t = year,
      p = prevalence
    )
  
  # --- Average exceed prob and pred prev every two years --------------------
  #exceed_prob_long_1_avg_2y <- avg_2_year(exceed_prob_long_1)
  exceed_prob_long_5_avg_2y <- avg_2_year(exceed_prob_long_5)
  exceed_prob_long_10_avg_2y <- avg_2_year(exceed_prob_long_10)
  exceed_prob_long_20_avg_2y <- avg_2_year(exceed_prob_long_20)
  exceed_prob_long_50_avg_2y <- avg_2_year(exceed_prob_long_50)
  points_df_2y <- add_year_block(points_df)
  p_long_mean_avg_2y <- avg_2_year(p_long_mean)
  
  # --- Helper to generate a consistent ggplot with point overlays ------------
  plot_theme <- function(p) {
    p + theme(
      # make everything white
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      strip.background = element_rect(fill = "white", color = NA),
      panel.grid       = element_blank(),
      strip.text = element_text(color = "black", face = "plain"),
      panel.border = element_rect(color = "grey80", fill = NA),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 8),
      title        = element_text(size = 12),
      axis.text.x  = element_text(size = 8, angle = 30, hjust = 1),
      axis.text.y  = element_text(size = 8),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.title   = element_text(hjust = 0),
      legend.justification = "center"
    )
  }
  
  plot_layer <- function(p_long_df, title_text, shp = shape_Africa, shp_grey = shape_grey, shp_water = shape_water, add_points = FALSE, add_legend = FALSE) {
    p <- ggplot() +
      # fill entire bbox with white first (ocean/background)
      annotate(
        "rect",
        xmin = xlim[1], xmax = xlim[2],
        ymin = ylim[1], ymax = ylim[2],
        fill = "white", colour = NA
      ) +
      geom_raster(aes(x = x, y = y, fill = mean_p * 100), data = p_long_df) +
      geom_sf(data = shp,        linewidth = 0.2, fill = NA, color = "white") +
      geom_sf(data = shp_grey,   fill = "grey80", colour = NA) +
      geom_sf(data = shp_water,  fill = "white",  colour = NA) +
      coord_sf(
        xlim = unname(as.numeric(xlim)),
        ylim = unname(as.numeric(ylim)),
        expand = FALSE,
        crs = sf::st_crs(4326),
        default_crs = sf::st_crs(4326),
        clip = "on"
      )
    if (add_points) {
      p <- p + geom_point(
        aes(x = x, y = y, fill = p),
        data = points_df_2y %>% arrange(p),
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
      facet_wrap(~year_block, nrow = 2) +
      labs(x = "Longitude", y = "Latitude")
    plot_theme(p)
  }
  
  plot_exceedance <-function(exceed_prob, title_text, legend_title, shp = shape_Africa, shp_grey = shape_grey, shp_water = shape_water, add_points = FALSE, add_legend = FALSE) {
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
      geom_sf(data = shape_grey, fill = "grey80", colour = NA) +
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
      facet_wrap(~year_block, nrow = 2) +
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
  plot_mean  <- plot_layer(p_long_mean_avg_2y, clean_mut, shape_Africa_crop, shape_grey, shape_water_crop, add_points = TRUE, add_legend = TRUE)
  plot_mean_no_points  <- plot_layer(p_long_mean_avg_2y, clean_mut, shape_Africa_crop, shape_grey, shape_water_crop, add_points = FALSE, add_legend = TRUE)
  
  #plot_exceed_1 <- plot_exceedance(exceed_prob_long_1_avg_2y, title_text = clean_mut, legend_title =expression(Pr(Prevalence >= 1*"%")))
  plot_exceed_5 <- plot_exceedance(exceed_prob_long_5_avg_2y, title_text = clean_mut, legend_title =expression(Pr(Prevalence >= 5*"%")))
  plot_exceed_10 <- plot_exceedance(exceed_prob_long_10_avg_2y, title_text = clean_mut, legend_title =expression(Pr(Prevalence >= 10*"%")))
  plot_exceed_20 <- plot_exceedance(exceed_prob_long_20_avg_2y, title_text = clean_mut, legend_title =expression(Pr(Prevalence >= 20*"%")))
  plot_exceed_50 <- plot_exceedance(exceed_prob_long_50_avg_2y, title_text = clean_mut, legend_title =expression(Pr(Prevalence >= 50*"%")))
  
  save_figs(file.path(OUT_PREV_DIR, paste0(mut, "_mean_perc_prev")), plot_mean, width = 10)
  save_figs(file.path(OUT_PREV_DIR, paste0(mut, "_mean_perc_prev_no_points")), plot_mean_no_points, width = 10)
  
  #save_figs(file.path(OUT_EXCEED_DIR, paste0(mut, "_exceedance_prob_1")), plot_exceed_1, width = 10)
  save_figs(file.path(OUT_EXCEED_DIR, paste0(mut, "_exceedance_prob_5")), plot_exceed_5, width = 10)
  save_figs(file.path(OUT_EXCEED_DIR, paste0(mut, "_exceedance_prob_10")), plot_exceed_10, width = 10)
  save_figs(file.path(OUT_EXCEED_DIR, paste0(mut, "_exceedance_prob_20")), plot_exceed_20, width = 10)
  save_figs(file.path(OUT_EXCEED_DIR, paste0(mut, "_exceedance_prob_50")), plot_exceed_50, width = 10)
  print(paste0("Saved figures for ", mut)) 
  }
  else{
    print(paste0(mut, " has no pos mutations")) 
  }
}
