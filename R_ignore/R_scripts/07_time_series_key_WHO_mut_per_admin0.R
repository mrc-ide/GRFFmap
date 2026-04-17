# ===============================================================
# Time-Series Plot of Admin-1 Prevalence from GRFF Kalman Output
# ===============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(here)
  library(sf)
  library(devtools)
  library(fs)
  library(ggrepel)
  library(cowplot)
})

load_all()
sf_use_s2(FALSE)
set.seed(1)

# --- Define country based palette ----------------------------------------------------------
country_colors <- c(
  "Burundi" = "#CC79A7",
  "Democratic Republic of the Congo" = "#7F7F7F",
  "Kenya" = "#56B4E9",
  "Uganda" = "#0072B2",
  "Sudan" = "#E69F00",
  "South Sudan" = "#009E73",
  "Rwanda" = "#D55E00",
  "Ethiopia" = "#332288",
  "Eritrea" = "#882255",
  "Tanzania" = "#44AA99"
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

# --- Add combined K13 data ----------------------------------------------------
dat <- add_combined_k13(dat)

# --- Filter data --------------------------------------------------------------
# Load Africa Shape file
africa_shp_admin1 <- readRDS("R_ignore/R_scripts/data/sf_admin1_africa.rds")

# --- Settings -----------------------------------------------------------------
plot_times    <- seq(2012, 2025, by = 1)
prev_cutoff   <- 0.1   # threshold for time series

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
key_mutations <- c("k13:622:I", "k13:675:V", "k13:561:H", "k13:comb")

all_mut_timeseries <- list()
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
    group_by(name_0, id_1, name_1) |>
    filter(any(mean_prevalence > prev_cutoff)) |>
    ungroup() |>
    filter(
      !(name_1 == "Lake Kivu"      & mut %in% c("k13:561:H", "k13:675:V", "k13:comb")),
      !(name_1 == "Lake Victoria"  & mut %in% c("k13:561:H", "k13:675:V", "k13:comb")),
      !(name_1 == "Sud-Kivu"       & mut %in% c("k13:561:H", "k13:675:V", "k13:comb")),
      !(name_1 == "Nord-Kivu"      & mut %in% c("k13:561:H", "k13:675:V", "k13:comb"))
    )
  # 
  # if (mut == "k13:561:H") {
  #   desired <- c("Tanzania", "Rwanda", "Uganda")
  #   
  #   adm1_timeseries_bycountry_data <- adm1_timeseries_bycountry_data |>
  #     mutate(
  #       name_0 = factor(
  #         name_0,
  #         levels = c(desired, setdiff(sort(unique(name_0)), desired))
  #       )
  #     )
  # } else {
  #   adm1_timeseries_bycountry_data <- adm1_timeseries_bycountry_data |>
  #     mutate(name_0 = factor(name_0, levels = sort(unique(name_0))))
  # }
  
  adm1_timeseries_bycountry_data <- adm1_timeseries_bycountry_data |>
    mutate(name_0 = factor(name_0, levels = sort(unique(name_0))))
  
  countries <- levels(adm1_timeseries_bycountry_data$name_0)
  my_colors <- country_colors[countries]
  
  adm1_timeseries_bycountry_data <- adm1_timeseries_bycountry_data %>%
    mutate(
      mean_prevalence = 100 * mean_prevalence,
      lower_ci_bound  = 100 * lower_ci_bound,
      upper_ci_bound  = 100 * upper_ci_bound
    )
  
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
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, by = 25),
      labels = function(x) paste0(x, "%")
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
  
  if (mut == "k13:561:H"){
    plot_height = 2.4
  } else{
    plot_height = 4
  }
  
  save_figs(
    file.path(OUT_PLOT_DIR, paste0(mut, "_country_prev_larger_", prev_cutoff, "_draws")),
    annual_prev_larger_plot_country,
    width  = 9,
    height = plot_height
  )
  
  # Save dataframe
  adm1_timeseries_bycountry_data <- adm1_timeseries_bycountry_data %>%
    mutate(mutation = mut)
  all_mut_timeseries[[mut]] <- adm1_timeseries_bycountry_data
}
combined_adm1_timeseries <- dplyr::bind_rows(all_mut_timeseries)
readr::write_csv(
  combined_adm1_timeseries,
  file.path("R_ignore/R_scripts/outputs/plots/", OUT_PLOT_DIR, "adm1_timeseries_bycountry_all_mutations.csv")
)

# Extract highest prev with corresponding admin1 and admin0 in 20212 and 2024
max_regions_2012_2024 <- combined_adm1_timeseries %>%
  filter(t %in% c(2012, 2024)) %>%
  group_by(mutation, t) %>%
  slice_max(
    order_by = mean_prevalence,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  select(mutation, peak_year = t, admin0 = name_0, admin1 = name_1)

max_regions_with_both_years <- max_regions_2012_2024 %>%
  left_join(
    combined_adm1_timeseries %>%
      filter(t %in% c(2012, 2024)) %>%
      select(
        mutation, t, name_0, name_1,
        mean_prevalence, lower_ci_bound, upper_ci_bound
      ),
    by = c("mutation", "admin0" = "name_0", "admin1" = "name_1")
  ) %>%
  mutate(year_label = paste0("y_", t)) %>%
  select(-t) %>%
  pivot_wider(
    names_from = year_label,
    values_from = c(mean_prevalence, lower_ci_bound, upper_ci_bound)
  ) %>%
  mutate(
    rate_change_per_year = (mean_prevalence_y_2024 - mean_prevalence_y_2012) / (2024 - 2012)
  ) %>%
  arrange(mutation, peak_year)

save_csv(name="time_series/max_prev_2012_2024_per_key_mut" , df = max_regions_with_both_years)
