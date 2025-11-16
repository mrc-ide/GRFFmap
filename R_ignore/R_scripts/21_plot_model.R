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

add_year_group <- function(df, year_col = year) {
  dplyr::mutate(
    df,
    # 1. Create the year_group_numeric column
    year_group_numeric = dplyr::case_when(
      {{ year_col }} %in% 2012:2014 ~ 0L, # Use 0L for integer type
      {{ year_col }} %in% 2015:2017 ~ 1L,
      {{ year_col }} %in% 2018:2020 ~ 2L,
      {{ year_col }} %in% 2021:2023 ~ 3L,
      TRUE ~ NA_integer_
    ),
    
    # 2. Keep the original character bin for reference (optional)
    year_group_char = dplyr::case_when(
      {{ year_col }} %in% 2012:2014 ~ "2012-2014",
      {{ year_col }} %in% 2015:2017 ~ "2015-2017",
      {{ year_col }} %in% 2018:2020 ~ "2018-2020",
      {{ year_col }} %in% 2021:2023 ~ "2021-2023",
      TRUE ~ NA_character_
    )
  )
}

bin_prevalence <- function(v) {
  labs <- PREV_LEVELS()
  cut(
    v,
    breaks = c(-Inf, 0, 1, 5, 10, 20, 30, 40, Inf),
    labels = labs,
    right = TRUE,
    ordered_result = TRUE
  )
}

PREV_LEVELS <- function() {
  c("0", "0-1", "1-5", "5-10", "10-20", "20-30", "30-40", "40+")
}

# ── Paths & constants ───────────────────────────────────────────────────────────
CACHE_DIR      <- "R_ignore/R_scripts/outputs/GRFF_cache_group_years"
OUT_BASE       <- "GRFF_updated_maps_grouped_year"
OUT_PREV_DIR   <- file.path("prev_prediction_grouped_year")
OUT_EXC_DIR    <- file.path("exceedance_prob_grouped_year")
dir_create(c(CACHE_DIR, OUT_BASE, OUT_PREV_DIR, OUT_EXC_DIR), recurse = TRUE)

EA_long_lim <- c(15, 52)
EA_lat_lim  <- c(-12, 18)
x_range     <- c(-20, 55)  # lon extent used to subset points
y_range     <- c(-35, 38)  # not used below, but kept for reference

# ── Cache helpers ───────────────────────────────────────────────────────────────
cache_path <- function(mut, lenS, lenT, what, ext = c("rds","parquet")) {
  ext <- match.arg(ext)
  fn  <- sprintf("%s_lenS%s_lenT%s_%s.%s",
                 gsub("[:/ ]", "_", mut), lenS, lenT, what, ext)
  file.path(CACHE_DIR, fn)
}

make_or_load <- function(mut, lenS, lenT) {
  pred_rds <- cache_path(mut, lenS, lenT, "df_pred", "rds")
  meta_rds <- cache_path(mut, lenS, lenT, "meta",   "rds")
  stopifnot(file_exists(pred_rds), file_exists(meta_rds))
  list(df_pred = readRDS(pred_rds), meta = readRDS(meta_rds))
}

# ── Data I/O ────────────────────────────────────────────────────────────────────
shape_Africa <- readRDS("R_ignore/R_scripts/data/shapefiles/sf_admin0_africa.rds")
shape_water  <- sf::st_read("R_ignore/R_scripts/data/shapefiles/africa_water_bodies.shp", quiet = TRUE)

# read in prevalence data
dat <- read.csv("R_ignore/R_scripts/data/all_who_get_prevalence_africa.csv") |>
  mutate(collection_day = as.Date(collection_day),
         year_group = year(collection_day))

# Aggregate site-level prevalence
k13_site <- dat %>%
  group_by(latitude, longitude, study_name, country_name,
           site_name, collection_day, year, denominator, mutation) %>%
  summarise(prevalence = sum(prevalence, na.rm = TRUE), .groups = "drop")

# Bin prevalences by years
k13_grouped <- k13_site %>%
  add_year_group(year) %>%
  filter(!is.na(year_group_numeric)) %>%
  mutate(
    prevalence_bin = bin_prevalence(prevalence),
    prevalence_bin = factor(prevalence_bin, levels = PREV_LEVELS()),
    numerator = as.integer(prevalence*denominator/100)
  ) %>%
  arrange(prevalence) %>%
  select(!collection_day, year)

# choose a single run to (re)plot from cache
length_space <- 75 / 110
length_time  <- 3

all_who_mutations <- c("k13:622:I", "k13:469:Y", "k13:675:V", "k13:446:I", "k13:458:Y", "k13:476:I",   "k13:493:H",   "k13:539:T",
                       "k13:543:T",  "k13:553:L",   "k13:561:H",   "k13:574:L",  "k13:580:Y",
                       "k13:441:L", "k13:449:A",   "k13:469:F",   "k13:481:V",
                       "k13:515:K", "k13:527:H",  "k13:537:I", "k13:537:D", "k13:538:V",  "k13:568:G")

for (focal_mut in all_who_mutations){
  # subset raw points shown on the plots
  dat_sub <- k13_grouped |>
    filter(mutation == focal_mut,
           longitude > x_range[1], longitude < x_range[2]) #collection_year > 2013
  
  # ── Load cached model outputs (df_pred + meta) ──────────────────────────────────
  # (If you still need to run the model, do it in a separate script and save RDS.)
  cached <- make_or_load(focal_mut, length_space, length_time)
  df_pred <- cached$df_pred      # columns: x, y, t, `2.5%`, `50%`, `97.5%`, MOE, p_pred, exceed_prob_{1,5,10}, …
  meta    <- cached$meta         # list with limits, years, etc.
  
  # ── Clip grid to only map relevant areas ──────────────────────────────────
  # points you care about
  pts_pos <- dat_sub |> 
    dplyr::filter(prevalence > 5) |>
    sf::st_as_sf(coords = c("longitude","latitude"), crs = 4326, remove = FALSE)
  
  if (dim(pts_pos)[1] != 0){
    # 300 km bbox in meters, then back to lon/lat
    buf_km   <- 100
    bbox_poly <- pts_pos |>
      sf::st_transform(3857) |>
      sf::st_bbox() |> sf::st_as_sfc() |>
      sf::st_buffer(buf_km * 1000) |>
      sf::st_transform(4326)
    
    bb   <- sf::st_bbox(bbox_poly)
    xlim <- c(bb["xmin"], bb["xmax"])
    ylim <- c(bb["ymin"], bb["ymax"])
    
    # ── Small plot builders to avoid repetition ─────────────────────────────────────
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
    
    
    facet_theme <- theme_void() +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 9),
        legend.text  = element_text(size = 8),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0)
      )
    
    add_base_layers <- function(p, shp, shp_water, xlim, ylim) {
      p +
        geom_sf(data = shp, linewidth = 0.2, fill = NA) +  # original polygons
        geom_sf(data = shp_water, fill = "white",  colour = NA, inherit.aes = FALSE) +
        coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
        xlab("Longitude") + ylab("Latitude")
    }
    
    plot_prev_surface <- function(dfp, title, add_points = FALSE, points_df = NULL, add_legend = FALSE, shp, shp_water, xlim, ylim) {
      p <- ggplot(dfp) + facet_theme +
        geom_raster(aes(x = x, y = y, fill = 100 * p_pred)) +
        scale_fill_gradientn(colours = pp_cols,
                             values  = pp_vals,
                             limits = c(0, 100), name = "Prevalence (%)") +
        facet_wrap(~collection_year, nrow = 2, ncol = 5)
      p <- add_base_layers(p, shp, shp_water, xlim, ylim)
      if (isTRUE(add_points) && !is.null(points_df)) {
        points_df <- points_df %>%
          mutate(collection_year = factor(
            year_group_numeric, levels = 0:3,
            labels = c("2012-2014","2015-2017","2018-2020","2021-2023")
          )) %>%
          filter(between(longitude, xlim[1], xlim[2]),
                 between(latitude,  ylim[1], ylim[2]))
        
        p <- p + geom_point(
          data = points_df |> dplyr::arrange(prevalence), #dplyr::filter(collection_year > 2013)
          aes(x = longitude, y = latitude, fill = prevalence),
          pch = 21, colour = grey(0.5), size = 1, inherit.aes = FALSE
        ) +
          scale_colour_gradientn(
            colours = pp_cols, values = pp_vals,
            limits  = c(0, 100),
            oob     = scales::squish,
            name    = "Observed prevalence (%)",
            guide   = guide_colorbar(barheight = unit(3, "mm"), barwidth = unit(40, "mm"))
          )
      }
      if (isFALSE(add_legend)){
        p <- p + theme(legend.position = "none")
      }
      else{
        p <- p + theme(legend.position = "bottom")
      }
      p
    }
    
    plot_exceed_surface <- function(dfp, prob_col, add_legend = FALSE, shp, shp_water, xlim, ylim) {
      p <- ggplot(dfp) + facet_theme +
        geom_raster(aes(x = x, y = y, fill = .data[[prob_col]])) +
        scale_fill_viridis(
          option = "magma",         
          name   = "Exceedance\nProbability",
          limits = c(0, 1),
          labels = scales::percent_format(accuracy = 1)
        ) +
        geom_sf(data = shp, color = "white", linewidth = 0.1, fill = NA, inherit.aes = FALSE) +
        geom_sf( data = shp_water, fill = "white",  colour = NA, inherit.aes = FALSE) +
        coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
        facet_wrap(~collection_year, nrow = 2, ncol = 5) +
        labs(x = "Longitude", y = "Latitude")
      
      if (isFALSE(add_legend)){
        p <- p + theme(legend.position = "none")
      } else {
        p <- p + theme(legend.position = "bottom")
      }
      p
    }
    
    # Masked prevalence (clipped view)
    plot_masked_prev <- function(dfp, thr_col, add_legend = FALSE, shp, shp_water, xlim, ylim) {
      p <- ggplot(dfp) + facet_theme +
        geom_sf(data = shp, linewidth = 0.1, fill = NA, inherit.aes = FALSE) +
        geom_sf( data = shp_water, fill = "white",  colour = NA, inherit.aes = FALSE) +
        geom_raster(aes(
          x = x, y = y,
          fill  = p_pred,
          alpha = ifelse(.data[[thr_col]] > 0.8, "≥80%", "<80%")
        )) +
        scale_fill_gradientn(
          colours = pp_cols,
          values  = pp_vals,
          name    = "Prevalence",
          limits  = c(0, 1)  # adjust if p_pred is 0–1; use c(0,100) if you multiply by 100
        ) +
        scale_alpha_manual(values = c("<80%" = 0.3, "≥80%" = 1.0), guide = "none") +
        coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
        facet_wrap(~collection_year, nrow = 2, ncol = 5) +
        labs(x = "Longitude", y = "Latitude")
      if (isFALSE(add_legend)){
        p <- p + theme(legend.position = "none")
      }
      else{
        p <- p + theme(legend.position = "bottom")
      }
      p
    }
    
    # ── Prepare data for faceting (add year column) ─────────────────────────────────
    dfp <- df_pred |>
      mutate(collection_year = factor(
               t, levels = 0:3,
               labels = c("2012-2014","2015-2017","2018-2020","2021-2023")
             )) |>
      # collection_year > 2013,
      filter(
             x >= xlim[1], 
             x <= xlim[2], 
             y >= ylim[1], 
             y <= ylim[2])
    
    # ── Build plots ─────────────────────────────────────────────────────────────────
    plot1 <- plot_prev_surface(
      dfp, title = focal_mut, add_points = TRUE, add_legend = FALSE,
      points_df = dat_sub,
      shp = shape_Africa, shp_water = shape_water, xlim = xlim, ylim = ylim
    )
    
    plot2 <- plot_prev_surface(
      dfp, title = focal_mut, add_points = FALSE, add_legend = FALSE,
      shp = shape_Africa, shp_water = shape_water, xlim = xlim, ylim = ylim
    )
    
    plot3 <- plot_prev_surface(
      dfp, title = focal_mut, add_points = FALSE, add_legend = TRUE,
      shp = shape_Africa, shp_water = shape_water, xlim = xlim, ylim = ylim
    )
    
    plot3_points <- plot_prev_surface(
      dfp, title = focal_mut, add_points = TRUE, add_legend = TRUE,
      points_df = dat_sub,
      shp = shape_Africa, shp_water = shape_water, xlim = xlim, ylim = ylim
    )
    
    plot_masked_exceed_10 <- plot_masked_prev(
      dfp, thr_col = "exceed_prob_10", add_legend = FALSE,
      shp = shape_Africa, shp_water = shape_water, xlim = xlim, ylim = ylim
    )
    
    plot_masked_exceed_5 <- plot_masked_prev(
      dfp, thr_col = "exceed_prob_5", add_legend = FALSE,
      shp = shape_Africa, shp_water = shape_water, xlim = xlim, ylim = ylim
    )
    
    plot_masked_exceed_1 <- plot_masked_prev(
      dfp, thr_col = "exceed_prob_1", add_legend = FALSE,
      shp = shape_Africa, shp_water = shape_water, xlim = xlim, ylim = ylim
    )
    
    plot_masked_exceed_1_legend <- plot_masked_prev(
      dfp, thr_col = "exceed_prob_1", add_legend = TRUE,
      shp = shape_Africa, shp_water = shape_water, xlim = xlim, ylim = ylim
    )
    
    plot_exceed_10 <- plot_exceed_surface(
      dfp, prob_col = "exceed_prob_10", add_legend = FALSE,
      shp = shape_Africa, shp_water = shape_water, xlim = xlim, ylim = ylim
    )
    
    plot_exceed_5 <- plot_exceed_surface(
      dfp, prob_col = "exceed_prob_5", add_legend = FALSE,
      shp = shape_Africa, shp_water = shape_water, xlim = xlim, ylim = ylim
    )
    
    plot_exceed_1 <- plot_exceed_surface(
      dfp, prob_col = "exceed_prob_1", add_legend = FALSE,
      shp = shape_Africa, shp_water = shape_water, xlim = xlim, ylim = ylim
    )
    
    plot_exceed_1_legend <- plot_exceed_surface(
      dfp, prob_col = "exceed_prob_1", add_legend = TRUE,
      shp = shape_Africa, shp_water = shape_water, xlim = xlim, ylim = ylim
    )
    
    # ── Save plots (uses your save_figs(name, fig, ...)) ────────────────────────────
    # NOTE: save_figs expects a *base path* (without extension) if your version appends .png/.pdf.
    # If your save_figs already expects full file stems, we pass a stem without extension.
    
    save_figs(file.path(OUT_PREV_DIR, paste0(focal_mut, "_Africa_points_", length_space, "_", length_time)), plot1)
    
    save_figs(file.path(OUT_PREV_DIR, paste0(focal_mut, "_Africa_", length_space, "_", length_time)),
              plot2)
    
    save_figs(file.path(OUT_PREV_DIR, paste0(focal_mut, "_Africa_legend_", length_space, "_", length_time)),
              plot3)
    
    save_figs(file.path(OUT_PREV_DIR, paste0(focal_mut, "_Africa_legend_points_", length_space, "_", length_time)),
              plot3_points)
    
    save_figs(file.path(OUT_EXC_DIR, paste0(focal_mut, "_Africa_", length_space, "_", length_time, "_exceedance_prob_10_masked")),
              plot_masked_exceed_10)
    
    save_figs(file.path(OUT_EXC_DIR, paste0(focal_mut, "_Africa_", length_space, "_", length_time, "_exceedance_prob_5_masked")),
              plot_masked_exceed_5)
    
    save_figs(file.path(OUT_EXC_DIR, paste0(focal_mut, "_Africa_", length_space, "_", length_time, "_exceedance_prob_1_masked")),
              plot_masked_exceed_1)
    
    save_figs(file.path(OUT_EXC_DIR, paste0(focal_mut, "_Africa_legend_", length_space, "_", length_time, "_exceedance_prob_1_masked")),
              plot_masked_exceed_1_legend)
    
    save_figs(file.path(OUT_EXC_DIR, paste0(focal_mut, "_Africa_", length_space, "_", length_time, "_exceedance_prob_10")),
              plot_exceed_10)
    
    save_figs(file.path(OUT_EXC_DIR, paste0(focal_mut, "_Africa_", length_space, "_", length_time, "_exceedance_prob_5")),
              plot_exceed_5)
    
    save_figs(file.path(OUT_EXC_DIR, paste0(focal_mut, "_Africa_", length_space, "_", length_time, "_exceedance_prob_1")),
              plot_exceed_1)
    
    save_figs(file.path(OUT_EXC_DIR, paste0(focal_mut, "_Africa_legend_", length_space, "_", length_time, "_exceedance_prob_1")),
              plot_exceed_1_legend)
    
  } 
}

