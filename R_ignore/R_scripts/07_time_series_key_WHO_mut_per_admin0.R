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
country_colors <- c(
  "Uganda"    = "#0072B2",  # deep blue
  "Rwanda"    = "#D55E00",  # orange
  "Ethiopia"  = "#009E73",  # green
  "Eritrea"   = "#CC79A7",  # pink / magenta
  "Sudan"     = "#E69F00",  # golden yellow
  "Tanzania"  = "#882255",  # dark maroon
  "Democratic Republic of the Congo" = "#7F7F7F"  # neutral grey
)

# --- Define cache and output directory ----------------------------------------
CACHE_DIR <- "R_ignore/R_scripts/outputs/model_outputs/GRFF_model_output_key_WHO_mutations"
OUT_PLOT_DIR <- "time_series"

dir_create(c(
  file.path("R_ignore/R_scripts/outputs/plots", OUT_PLOT_DIR)
))

# --- Load data ----------------------------------------------------------------
dat <- read.csv("R_ignore/R_scripts/data/all_mutations_get_prevalence.csv") |>
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
key_mutations <- c("k13:622:I", "k13:675:V", "k13:561:H")
fig5_plots <- list()

# --- Main loop over mutations -------------------------------------------------
for (mut in key_mutations){
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
  shape_lookup <- admin1_regions %>%
    st_drop_geometry() %>%
    distinct(id_1, name_0, name_1)
  
  summ_df <- df_posterior_timeseries %>%
    group_by(t, id_1) %>%
    summarise(
      lower_ci_bound  = quantile(mean_p, probs = 0.025),
      mean_prevalence = quantile(mean_p, probs = 0.5),
      upper_ci_bound  = quantile(mean_p, probs = 0.975),
      .groups = "drop"
    )
  
  admin1_time_series_final <- summ_df %>%
    left_join(shape_lookup, by = "id_1") %>%
    filter(t >= min(plot_times), t <= max(plot_times) + 1)
  
  stopifnot(nrow(admin1_time_series_final) == nrow(summ_df))
  
  n_missing <- sum(is.na(admin1_time_series_final$name_0) | is.na(admin1_time_series_final$name_1))
  if (n_missing > 0) {
    message("WARNING: ", n_missing, " rows did not match an admin1 region (missing name_0/name_1).")
    
    # which id_1 values are missing?
    missing_ids <- admin1_time_series_final %>%
      filter(is.na(name_0) | is.na(name_1)) %>%
      distinct(id_1) %>%
      arrange(id_1)
    
    print(missing_ids)
    
    # optional: show a few (t, id_1) examples
    print(
      admin1_time_series_final %>%
        filter(is.na(name_0) | is.na(name_1)) %>%
        select(t, id_1) %>%
        arrange(t, id_1) %>%
        head(20)
    )
    
    # hard fail if you want strictness:
    stop("Join to admin1_regions failed for some id_1 values.")
  }
  
  # Also check the other direction: shape IDs not present in posterior summary
  ids_in_post  <- unique(summ_df$id_1)
  ids_in_shape <- unique(shape_lookup$id_1)
  
  missing_in_shape <- setdiff(ids_in_post, ids_in_shape)
  
  if (length(missing_in_shape) > 0) {
    message("Posterior has ", length(missing_in_shape), " id_1 values not present in shape_lookup.")
    print(head(missing_in_shape, 20))
  }
  
  # --- Plot only countries where prevalence ever exceeds cutoff ------------
  # which countries ever exceed prev_cutoff at region level
  select_countries <- admin1_time_series_final |>
    filter(mean_prevalence > prev_cutoff) |>
    distinct(name_0)
  
  adm1_timeseries_bycountry_data <- admin1_time_series_final |>
    semi_join(select_countries, by = "name_0") |>
    filter(
      !(name_1 == "Lake Kivu"      & mut %in% c("k13:561:H", "k13:675:V")),
      !(name_1 == "Lake Victoria"  & mut %in% c("k13:561:H", "k13:675:V")),
      !(name_1 == "Sud-Kivu"  & mut %in% c("k13:561:H", "k13:675:V")),
      !(name_1 == "Nord-Kivu"  & mut %in% c("k13:561:H", "k13:675:V"))
    )
  
  if (mut == "k13:561:H") {
    desired <- c("Tanzania", "Rwanda", "Uganda")
    
    adm1_timeseries_bycountry_data <- adm1_timeseries_bycountry_data |>
      mutate(
        name_0 = factor(
          name_0,
          levels = c(desired, setdiff(sort(unique(name_0)), desired))
        )
      )
  } else {
    adm1_timeseries_bycountry_data <- adm1_timeseries_bycountry_data |>
      mutate(name_0 = factor(name_0, levels = sort(unique(name_0))))
  }
  
  countries <- levels(adm1_timeseries_bycountry_data$name_0)
  my_colors <- country_colors[countries]
  
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
      nudge_x       = 1.5,
      size          = 2.5,
      show.legend   = FALSE,
      segment.color = "grey60",
      max.overlaps  = 10
    ) +
    facet_wrap(~ name_0) +
    coord_cartesian(ylim = c(0, 1)) +
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
    theme(legend.position = "None",
          panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
          panel.grid.minor = element_blank()
  )
  annual_prev_larger_plot_country <- plot_theme(annual_prev_larger_plot_country)
  
  
  save_figs(
    file.path(OUT_PLOT_DIR, paste0(mut, "_country_prev_larger_", prev_cutoff, "_draws")),
    annual_prev_larger_plot_country,
    width  = 9,
    height = 3.5
  )
  
  fig5_plots[[mut]] <- annual_prev_larger_plot_country
}
