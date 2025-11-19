# ============================================================
# RFF + PG-EM with Kalman filter/RTS in projected (Cartesian) coords
# ============================================================
# - Reads WHO GET prevalence data
# - Filters to a mutation and spatial window
# - Projects lon/lat -> local LAEA (km), runs RFF + PG-EM in km
# - Reconstructs prevalence on a lon/lat grid for plotting
# - Overlays observations with an in-window/out-of-window styling
# ============================================================

library(tidyverse)
library(Matrix)
library(here)
library(sf)
library(devtools)
library(dplyr)
library(fs)
library(here)
library(ggrepel)
library(cowplot)

load_all()
sf_use_s2(FALSE)
set.seed(1)

# build consistent filenames per (mutation, length_space, length_time, n_features)
cache_path <- function(mut, lenS, lenT, what, ext = c("rds","parquet","rds.gz")) {
  ext <- match.arg(ext)
  fn <- sprintf("%s_lenS%s_lenT%s_%s.%s",
                gsub("[:/ ]", "_", mut), lenS, lenT, what, ext)
  file.path(CACHE_DIR, fn)
}

make_or_load <- function(mut, ell_km, tau2) {
  cached_rds <- cache_path(mut, ell_km, tau2, "full_model_output", "rds") 
  stopifnot(file_exists(cached_rds))
  readRDS(cached_rds)
}

# --------------------------- Settings --------------------------------------
# plotting parameters
plot_times <- seq(2012, 2023, by = 1) # should be within t_vec

ell_km <- 80
tau2   <- 0.1
n_post_draws <- 1000

# --- Paths & constants -----------------------------------------------------------
CACHE_DIR      <- "R_ignore/R_scripts/outputs/GRFF_kalman_cache_annual_2012_2025"
OUT_BASE       <- "GRFF_updated_maps"
OUT_PLOT_DIR <- "time_series"
OUT_PLOT_COUNTRY_DIR <- "country"
dir_create(c(paste0("R_ignore/R_scripts/outputs/plots/", OUT_PLOT_DIR),
             paste0("R_ignore/R_scripts/outputs/plots/", OUT_PLOT_DIR, "/", OUT_PLOT_COUNTRY_DIR)))

# --------------------------- Load & filter data ----------------------------
# read in prevalence data
dat <- read.csv("R_ignore/R_scripts/data/all_who_get_prevalence.csv") |>
  mutate(collection_day = as.Date(collection_day))

# Read in shape files
africa_shp_admin1 <- readRDS(file = "R_ignore/R_scripts/data/sf_admin1_africa.rds")
target_crs <- 4326
africa_shp_admin1 <- st_transform(africa_shp_admin1, crs = target_crs)

# Define East Africa box and crop Admin1 shapefile
bbox_east_africa <- st_bbox(c(
  xmin = 28.48,
  xmax = 48.43,
  ymin = -4.6,
  ymax = 15.29
), crs = target_crs)

# Using the cropped Admin 1 data for East Africa
admin1_regions <- st_make_valid(africa_shp_admin1) |> st_crop(bbox_east_africa)

# which mutation we are focusing on
# all_who_mutations <- c("k13:622:I", "k13:469:Y", "k13:675:V", "k13:446:I", "k13:458:Y", "k13:476:I",   "k13:493:H",   "k13:539:T",
#                        "k13:543:T",  "k13:553:L",   "k13:561:H",   "k13:574:L",  "k13:580:Y",
#                        "k13:441:L", "k13:449:A",   "k13:469:F",   "k13:481:V",
#                        "k13:515:K", "k13:527:H",  "k13:537:I", "k13:537:D", "k13:538:V",  "k13:568:G")

all_who_mutations <- c("k13:622:I", "k13:675:V", "k13:561:H")

fig5_plots <- list()

for (mut in all_who_mutations){
  print(paste0("Processing ", mut, "........."))
  clean_mut <- gsub("^k13:(\\d+):([A-Za-z])$", "k13 \\1\\2", mut)
  
  # subset data for this mutation and AOI
  dat_sub <- dat |>
    filter(mutation == mut,
           year > min(plot_times), year < max(plot_times)+1)
  # quick skip if no positive prevalence
  if (!any(dat_sub$prevalence > 0, na.rm = TRUE)) next
  
  # load cached model output (skip cleanly if missing)
  cached <- try(make_or_load(mut, ell_km, tau2), silent = TRUE)
  xs <- cached$xs; ys <- cached$ys
  df_posterior_timeseries <- cached$time_series
    
  # Compute summary statistics (median, CI)
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
    filter(t > min(plot_times), t < max(plot_times)+1)
  
  #Plot1: All Prev plots by Country
  annual_prev_plot_country <- ggplot(admin1_time_series_final, aes(x = t, y = mean_prevalence, group = id_1, color = name_0, fill = name_0)) +
    geom_ribbon(aes(ymin = lower_ci_bound, ymax = upper_ci_bound), 
                alpha = 0.10, 
                linetype = "blank") +
    geom_line(linewidth = 0.7, alpha =1) +
    facet_wrap(~ name_0) +
    labs(
      x = "Year", 
      y = "Mean Prevalence", 
      title = paste("Annual Prevalence Trends by Country for", clean_mut),
      subtitle = "Each line represents an Admin 1 region"
    ) +
     geom_text_repel(data = (admin1_time_series_final %>%
                              group_by(name_0, name_1) %>%
                              slice_max(order_by = t, n = 1) %>%
                              ungroup()),
      aes(label = name_1),
      color = "black",
      hjust = 1,
      nudge_x = 0.5,
      size = 2.5,
      show.legend = FALSE,
      segment.color = "grey60",
      max.overlaps = 10
    ) +
    #theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  save_figs(file.path(OUT_PLOT_DIR, paste0(mut, "country_prev_draws")), annual_prev_plot_country, width = 10)

  #Plot2: Prev plots by Country for prevalences > XX%
  #filter out at this number 
  prev_cutoff <- 0.05
  #filter the data
  select_districts <- admin1_time_series_final %>%
    filter(mean_prevalence > prev_cutoff) %>%
    distinct(name_0)
  adm1_timeseries_bycountry_data <- admin1_time_series_final %>%
    semi_join(select_districts, by = c("name_0"))  %>%
    filter(!(name_1 == "Lake Kivu" & (mut == "k13:561:H"|mut =="k13:675:V"))) %>%
    filter(!(name_1 == "Lake Victoria" & (mut == "k13:561:H"|mut =="k13:675:V")))
  # If we do timeseries for all mutations I would make a new color pal based on num countries for each mut
  if(mut == "k13:561:H" | mut == "k13:675:V"){
    my_colors <- c(
      "#332288","#117733","#44AA99"
    )
  }
  if(mut == "k13:622:I"){
    my_colors <- c(
      "#88CCEE","#DDCC77",
      "#CC6677","#882255","#AA4499","#DDDDDD","#000000",
      "#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
      "#D55E00","#CC79A7","#9A6324","#800000"
    )
  }
  #Plot with ADM1 labels
  annual_prev_larger_1perc_plot_country <- adm1_timeseries_bycountry_data |>
    ggplot(aes(x = t, y = mean_prevalence, 
               group = id_1,
               color = name_0,
               fill = name_0)) +
    geom_ribbon(aes(ymin = lower_ci_bound, ymax = upper_ci_bound), 
                alpha = 0.10, 
                linetype = "blank") +
    geom_line(linewidth = 0.7, alpha =1) +
    geom_text_repel(
      data = (adm1_timeseries_bycountry_data %>%
                group_by(name_0, name_1) %>%
                slice_max(order_by = t, n = 1) %>%
                ungroup()),
      aes(label = name_1),
      color = "black",
      hjust = 0,
      nudge_x = 0.5,
      size = 2.5,
      show.legend = FALSE,
      segment.color = "grey60",
      max.overlaps = 10
    ) +
    facet_wrap(~ name_0) +
    #scale_y_continuous(limits = c(0,0.65))+
    scale_color_manual(values = my_colors) +
    scale_fill_manual(values = my_colors)+
    scale_x_continuous(breaks = seq(min(adm1_timeseries_bycountry_data$t), max(adm1_timeseries_bycountry_data$t), by = 2)) +
#   scale_x_continuous(breaks = unique(adm1_timeseries_bycountry_data$t))+
    labs(
      x = "Year", 
      y = paste("Prevalence"), 
      title = clean_mut 
      #paste("Annual Prevalence Trends by Country for", mut, "in ADM 1 where prev is >", prev_cutoff),
    ) +
    theme_bw()+
    theme(
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      title = element_text(size = 12),
      axis.text.x = element_text(size = 8, angle = 30, hjust = 1),  # <-- key tweak
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.05),
      #axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  save_figs(file.path(OUT_PLOT_DIR, paste0(mut, "country_prev_larger_", prev_cutoff, "_draws")), annual_prev_larger_1perc_plot_country, width = 10)
  
  fig5_plots[[mut]] <- annual_prev_larger_1perc_plot_country
}

##### Combine all TimeSeries plots into Figure5 with Cow Plot 
panelA <- fig5_plots[["k13:675:V"]]
panelB <- fig5_plots[["k13:561:H"]]
panelC <- fig5_plots[["k13:622:I"]]
Fig5 <- plot_grid(panelA, panelB, panelC, ncol = 1, align = "v",
                  labels = c("A","B","C"), label_size = 14)

save_figs(file.path(OUT_PLOT_DIR, paste0("Fig5_combined_plot_", prev_cutoff, "_draws")), Fig5, width = 8, height = 10)

save_figs(file.path(OUT_PLOT_DIR, paste0("Fig5_combined_plot_LARGE_", prev_cutoff, "_draws")), Fig5, width = 15, height = 15)


#################PLOTTING INDIVIDUAL Countries with and without datapoints #####
#CURRENTLY NOT SAVING ANYTHING FROM HERE
ADM1greater01 <- unique(admin1_time_series_final[admin1_time_series_final$mean_prevalence >prev_cutoff ,]$name_0)
for (admin in ADM1greater01) {
  time_series_plot_admin <- admin1_time_series_final |> filter(admin == name_0) |>
    ggplot(aes(x = t, y = mean_prevalence, 
               group = id_1,
               color = name_1,
               fill = name_1)) +
    geom_ribbon(aes(ymin = lower_ci_bound, ymax = upper_ci_bound), 
                alpha = 0.10, 
                linetype = "blank") +
    geom_line(linewidth = 0.7, alpha =1) +
    facet_wrap(~ name_1) +
    labs(
      x = "Year", 
      y = "Mean Prevalence", 
      title = paste("Annual Prevalence Trends for", admin, mut),
      subtitle = "Each graph is Admin 1 Region"
    ) +
    #theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  #save_figs(file.path(OUT_PLOT_DIR, OUT_PLOT_COUNTRY_DIR, mut, paste0(mut, "_", admin ,"draws", n_post_draws, "prev")), time_series_plot_admin, width = 10)
  
  sub_with_geom <- dat_sub %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(admin1_regions), remove = FALSE) 
  
  # Admin 1 - First-level administrative region
  coords_with_admin1 <- st_join(sub_with_geom, africa_shp_admin1, join = st_within) %>%
    mutate(
      lon = st_coordinates(.)[, 1],  # Extract longitude
      lat = st_coordinates(.)[, 2]   # Extract latitude
    ) 
  
  coords_with_admin1 <- rename(coords_with_admin1, id_1 = id_1) 
  coords_with_admin1 <- coords_with_admin1 %>% filter(name_0 == admin)
  
  time_series_plot_admin_with_data <- admin1_time_series_final |> filter(admin == name_0) |>
    ggplot(aes(x = t, y = mean_prevalence, 
               group = id_1,
               color = name_1,
               fill = name_1)) +
    geom_ribbon(aes(ymin = lower_ci_bound, ymax = upper_ci_bound), 
                alpha = 0.10, 
                linetype = "blank") +
    geom_pointrange(aes(x = year, y = (prevalence/100), ymin = (prevalence_lower/100), ymax = (prevalence_upper/100)), data = coords_with_admin1) +
    geom_line(linewidth = 0.7, alpha =1) +
    facet_wrap(~ name_1) +
    labs(
      x = "Year", 
      y = "Mean Prevalence", 
      title = paste("Annual Prevalence Trends for", admin, mut, "with real data points"),
      subtitle = "Each graph is Admin 1 Region"
    ) +
   # theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  #save_figs(file.path(OUT_PLOT_DIR, OUT_PLOT_COUNTRY_DIR, mut, paste0(mut, "_", admin ,"draws", n_post_draws, "prev_with_data")), time_series_plot_admin_with_data, width = 10)
}
