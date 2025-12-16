# ===============================================================
# Time-Series Plot of Admin-1 Prevalence from GRFF Kalman Output
# ===============================================================

library(tidyverse)
library(Matrix)
library(here)
library(sf)
library(devtools)
library(fs)
library(ggrepel)
library(cowplot)

load_all()
sf_use_s2(FALSE)
set.seed(1)

# --- Define country based palette ----------------------------------------------------------

# base color palette for countries (will be repeated/trimmed as needed)
country_base_cols <- c(
  "#88CCEE","#DDCC77",
  "#CC6677","#882255","#AA4499","#DDDDDD","#000000",
  "#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
  "#D55E00","#CC79A7","#9A6324","#800000"
)

# --- Define cache and output directory ----------------------------------------
CACHE_DIR <- "R_ignore/R_scripts/outputs/model_outputs/GRFF_model_output_key_WHO_mutations"
OUT_PLOT_DIR <- "time_series"

dir_create(c(
  file.path("R_ignore/R_scripts/outputs/plots", OUT_PLOT_DIR)
))

# --- Load data ----------------------------------------------------------------
dat <- read.csv("R_ignore/R_scripts/data/all_who_get_prevalence.csv") |>
  mutate(collection_day = as.Date(collection_day))

# --- Filter data --------------------------------------------------------------
# Load Africa Shape file
africa_shp_admin1 <- readRDS("R_ignore/R_scripts/data/sf_admin1_africa.rds")

# --- Settings -----------------------------------------------------------------
plot_times    <- seq(2012, 2025, by = 1)
prev_cutoff   <- 0.05   # threshold for time series

# --- Fix shapes ---------------------------------------------------------------
target_crs <- 4326
africa_shp_admin1 <- st_transform(africa_shp_admin1, crs = target_crs)

# Define East Africa box and crop Admin1 shapefile
bbox_east_africa <- get_east_africa_bbox(target_crs)

# Using the cropped Admin 1 data for East Africa
admin1_regions <- africa_shp_admin1 |>
  st_make_valid() |>
  st_crop(bbox_east_africa)

# --- K13 Mutations of interests -----------------------------------------------
all_who_mutations <- c("k13:622:I", "k13:675:V", "k13:561:H")
fig5_plots <- list()

# --- Main loop over mutations -------------------------------------------------
for (mut in all_who_mutations){
  print(paste0("Processing ", mut, "........."))
  clean_mut <- gsub("^k13:(\\d+):([A-Za-z])$", "k13 \\1\\2", mut)
  
  # subset data for this mutation and AOI over plot_times
  dat_sub <- dat |>
    filter(mutation == mut,
           year >= min(plot_times), year <= max(plot_times)+1)
  
  # quick skip if no positive prevalence
  if (!any(dat_sub$prevalence > 0, na.rm = TRUE)) {
    message("  Skipping ", mut, " (no positive prevalence)")
    next
  }
  
  # load cached model output
  cached <- make_or_load(mut)
  xs <- cached$xs
  ys <- cached$ys
  df_posterior_timeseries <- cached$time_series
  p_median <- cached$p_post_median
    
  # --- Summarise posterior by (t, id_1) into median + 95% CI ------------------
  admin1_time_series_final <- df_posterior_timeseries %>%
    group_by(t, id_1) %>%
    summarise(
      lower_ci_bound = quantile(mean_p, probs = 0.025),
      mean_prevalence = quantile(mean_p, probs = 0.5),
      upper_ci_bound = quantile(mean_p, probs = 0.975),
      .groups = 'drop'
    ) %>%
    left_join(admin1_regions %>% st_drop_geometry() %>% distinct(id_1, name_0, name_1),
              by = "id_1") %>%
    filter(t >= min(plot_times), t <= max(plot_times)+1)

  # --- Summarise posterior by (t, id_1) into median + 95% CI ------------------
  annual_prev_plot_country <- ggplot(
    admin1_time_series_final,
    aes(
      x     = t,
      y     = mean_prevalence*100,
      group = id_1,
      color = name_0,
      fill  = name_0
    )
  ) +
    geom_ribbon(
      aes(ymin = lower_ci_bound*100, ymax = upper_ci_bound),
      alpha    = 0.10,
      linetype = "blank"
    ) +
    geom_line(linewidth = 0.7, alpha = 1) +
    facet_wrap(~ name_0) +
    geom_text_repel(
      data = admin1_time_series_final |>
        group_by(name_0, name_1) |>
        slice_max(order_by = t, n = 1) |>
        ungroup(),
      aes(label = name_1),
      color         = "black",
      hjust         = 1,
      nudge_x       = 0.5,
      size          = 2.5,
      show.legend   = FALSE,
      segment.color = "grey60",
      max.overlaps  = 10
    ) +
    labs(
      x        = "Year",
      y        = "Mean Prevalence",
      title    = paste("Annual Prevalence Trends by Country for", clean_mut),
      subtitle = "Each line represents an Admin 1 region"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # save_figs(
  #   file.path(OUT_PLOT_DIR, paste0(mut, "_country_prev_draws")),
  #   annual_prev_plot_country,
  #   width = 10
  # )
  
  # --- Plot 2: Plot only countries where prevalence ever exceeds cutoff ------------
  # which countries ever exceed prev_cutoff at region level
  select_countries <- admin1_time_series_final |>
    filter(mean_prevalence > prev_cutoff) |>
    distinct(name_0)
  
  adm1_timeseries_bycountry_data <- admin1_time_series_final |>
    semi_join(select_countries, by = "name_0") |>
    filter(
      !(name_1 == "Lake Kivu"      & mut %in% c("k13:561:H", "k13:675:V")),
      !(name_1 == "Lake Victoria"  & mut %in% c("k13:561:H", "k13:675:V"))
    )
  
  countries <- sort(unique(adm1_timeseries_bycountry_data$name_0))
  
  # repeat/trim base palette to match number of countries
  my_colors <- rep_len(country_base_cols, length(countries))
  names(my_colors) <- countries
  
  annual_prev_larger_plot_country <- adm1_timeseries_bycountry_data |>
    ggplot(
      aes(
        x     = t,
        y     = mean_prevalence,
        group = id_1,
        color = name_0,
        fill  = name_0
      )
    ) +
    geom_ribbon(
      aes(ymin = lower_ci_bound, ymax = upper_ci_bound),
      alpha    = 0.10,
      linetype = "blank"
    ) +
    geom_line(linewidth = 0.7, alpha = 1) +
    geom_text_repel(
      data = adm1_timeseries_bycountry_data |>
        group_by(name_0, name_1) |>
        slice_max(order_by = t, n = 1) |>
        ungroup(),
      aes(label = name_1),
      color         = "black",
      hjust         = 0,
      nudge_x       = 0.5,
      size          = 2.5,
      show.legend   = FALSE,
      segment.color = "grey60",
      max.overlaps  = 10
    ) +
    facet_wrap(~ name_0) +
    scale_color_manual(values = my_colors) +
    scale_fill_manual(values  = my_colors,
                      name = "Country") +
    scale_x_continuous(
      breaks = seq(
        from = min(adm1_timeseries_bycountry_data$t),
        to   = max(adm1_timeseries_bycountry_data$t),
        by   = 2
      )
    ) +
    labs(
      x     = "Year",
      y     = "Prevalence",
      title = clean_mut
    ) +
    theme_bw() + 
    theme(legend.position = "None")
  
  annual_prev_larger_plot_country <- plot_theme(annual_prev_larger_plot_country)
  
  
  save_figs(
    file.path(OUT_PLOT_DIR, paste0(mut, "_country_prev_larger_", prev_cutoff, "_draws")),
    annual_prev_larger_plot_country,
    width  = 8,
    height = 3
  )
  
  fig5_plots[[mut]] <- annual_prev_larger_plot_country
}
