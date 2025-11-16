
library(sf)
library(tidyverse)
library(plotly)
library(rstan)
library(epitools)
sf_use_s2(FALSE)

#------------------------------------------------------------------
# AGGREGATE FOR ADMIN 1 TIME SERIES ANALYSIS
#------------------------------------------------------------------

# 1. LOAD ADMIN BOUNDARIES
#rds_file_admin0 <- "R_ignore/R_scripts/data/sf_admin0_africa.rds"
rds_file_admin1 <- "R_ignore/R_scripts/data/sf_admin1_africa.rds"

#africa_shp_admin0 <- readRDS(file = rds_file_admin0)
africa_shp_admin1 <- readRDS(file = rds_file_admin1)

target_crs <- 4326
#africa_shp_admin0 <- st_transform(africa_shp_admin0, crs = target_crs)
africa_shp_admin1 <- st_transform(africa_shp_admin1, crs = target_crs)

# Define East Africa box and crop Admin1 shapefile
bbox_east_africa <- st_bbox(c(
  xmin = 28.48,
  xmax = 48.43,
  ymin = -4.6,
  ymax = 15.29
), crs = target_crs)

# Using the cropped Admin 1 data for East Africa
#admin0_regions <- st_make_valid(africa_shp_admin0) |> st_crop(bbox_east_africa)
admin1_regions <- st_make_valid(africa_shp_admin1) |> st_crop(bbox_east_africa)

dat <- read.csv("R_ignore/R_scripts/data/all_who_get_prevalence_africa.csv") |>
  mutate(collection_day = as.Date(collection_day),
         collection_year = year(collection_day))

# which mutation we are focusing on
all_who_mutations <- c("k13:622:I", "k13:469:Y", "k13:675:V", "k13:446:I", "k13:458:Y", "k13:476:I",   "k13:493:H",   "k13:539:T",
                       "k13:543:T",  "k13:553:L",   "k13:561:H",   "k13:574:L",  "k13:580:Y",
                       "k13:441:L", "k13:449:A",   "k13:469:F",   "k13:481:V",
                       "k13:515:K", "k13:527:H",  "k13:537:I", "k13:537:D", "k13:538:V",  "k13:568:G")

focal_mut <- "k13:675:V"
# ------------------------------------------------------------------
# RUN MODEL - Following 19_ implementation
# ------------------------------------------------------------------
# subset data
# define map range (longitude and latitude)
x_range <- c(-20, 55)  # Covers westernmost to easternmost Africa
y_range <- c(-35, 38)  # Covers southernmost to northernmost Africa

#ea map range - used to speed up model on prelim runs 
#x_range <- c(28.48, 48.43)
#y_range <- c(-4.6, 15.29)

# subset data
dat_sub <- dat |>
  filter(mutation == focal_mut) |>
  filter((longitude > x_range[1]) & (longitude < x_range[2]))

for (length_space_val in c(75, 110)){
  length_space <- length_space_val / 110 # dividing by 110 because the grid is in units of lat/lon, but numerator is in units of km. This is approximately correct near to the equator
  print(length_space)
  length_time <- 3
  pi_nugget <- 0.1 # the size of the nugget variance as a proportion of the GP variance
  
  # draw features
  n_features <- 100
  omega_space <- rbind(rnorm(n_features),
                       rnorm(n_features))
  omega_time <- rbind(rnorm(n_features))
  
  # make feature maps
  X_space <- dat_sub |>
    dplyr::select(longitude, latitude) |>
    as.matrix()
  X_time <- dat_sub$collection_year
  
  feat <- cbind(cos(X_space %*% omega_space / length_space + X_time %*% omega_time / length_time),
                sin(X_space %*% omega_space / length_space + X_time %*% omega_time / length_time)) / sqrt(n_features)
  
  # prepare data list for Stan
  data_list <- list(
    N = nrow(dat_sub),
    n_features = n_features,
    feat = feat,
    k = dat_sub$numerator,
    n_trials = dat_sub$denominator,
    pi_nugget = pi_nugget
  )
  
  # run stan model (takes a long time to compile the FIRST time, but then should be faster)
  fit <- stan(
    file = "R_ignore/R_scripts/stan/RFF_v1.stan",
    data = data_list,
    iter = 1e3,
    chains = 1
  )
  
  # extract MCMC samples
  draws <- extract(fit)
  n_draws <- nrow(draws$beta)
  
  beta_samples <- draws$beta
  mu_samples <- draws$mu
  sigma_samples <- draws$sigma
  
  #hist(mu_samples, breaks = 100)
  #hist(sigma_samples, breaks = 100)
  
  # ------------------------------------------------
  
  # make feature maps for prediction
  nx = 300
  ny = 300
  X_pred <- expand.grid(x = seq(x_range[1], x_range[2], l = nx),
                        y = seq(y_range[1], y_range[2], l = ny)) |>
    as.matrix()
  
  # choose plotting time slices based on years present in data)
  t_plot <- sort(unique(dat_sub$collection_year))
}


# 1. PREP SPATIAL DATA
admin1_regions <- admin1_regions %>%
  mutate(region_id = id_1)

# Maps the grid points (X_pred) to their admin regions.
points_to_regions <- st_as_sf(as.data.frame(X_pred), coords = c("x", "y"), crs = st_crs(admin1_regions)) %>%
  st_join(admin1_regions %>% dplyr::select(region_id, name_0, name_1)) %>%
  st_drop_geometry() %>%
  mutate(grid_index = 1:n()) %>% # Keep track of the original row index
  drop_na(region_id)

# 2. PREDICT regional average FOR ALL YEARS at EACH posterior draw 
t_pred_vec <- seq(min(dat_sub$collection_year), max(dat_sub$collection_year), by = 1)

l <- list()
for (i in seq_along(t_pred_vec)) {
  message(sprintf("Processing year %s (%s of %s)", t_pred_vec[i], i, length(t_pred_vec)))
  
  # Generate prediction features for this time point
  t_pred <- rep(t_pred_vec[i], nrow(X_pred))
  feat_pred <- cbind(cos(X_pred %*% omega_space / length_space + t_pred %*% omega_time / length_time),
                     sin(X_pred %*% omega_space / length_space + t_pred %*% omega_time / length_time)) / sqrt(n_features)
  
  # Get prevalence predictions over all posterior draws for all grid points
  post_all_points <- feat_pred %*% t(beta_samples) %>%
    sweep(MARGIN = 2, STATS = mu_samples, FUN = "+") %>%
    plogis()
  
  # For EACH posterior draw, calculate the mean prevalence per admin region
  regional_means_for_draws <- post_all_points[points_to_regions$grid_index, ] %>%
    as.data.frame() %>%
    bind_cols(points_to_regions %>% dplyr::select(region_id)) %>%
    pivot_longer(-region_id, names_to = "post_draw_raw", values_to = "p") %>%
    mutate(post_draw = as.integer(gsub("V", "", post_draw_raw))) %>%
    group_by(region_id, post_draw) %>%
    summarise(mean_p = mean(p), .groups = 'drop')
  
  #Store
  l[[i]] <- regional_means_for_draws %>%
    mutate(t = t_pred_vec[i])
}

# Combine all results into a single data frame
df_posterior_timeseries <- bind_rows(l)


# 3. Calculate Quantiles over posterior draws
admin1_time_series_final <- df_posterior_timeseries %>%
  group_by(t, region_id) %>%
  summarise(
    lower_ci_bound = quantile(mean_p, probs = 0.025),
    mean_prevalence = quantile(mean_p, probs = 0.5),
    upper_ci_bound = quantile(mean_p, probs = 0.975),
    .groups = 'drop'
  ) %>%
  left_join(admin1_regions %>% st_drop_geometry() %>% distinct(region_id, name_0, name_1), by = "region_id")



# 4. Plotting
time_series_plot <- ggplot(admin1_time_series_final, aes(x = t, y = mean_prevalence, 
                                                         group = region_id,
                                                         color = name_0,
                                                         fill = name_0)) +
  geom_ribbon(aes(ymin = lower_ci_bound, ymax = upper_ci_bound), 
              alpha = 0.10, 
              linetype = "blank") +
  geom_line(linewidth = 0.7, alpha =1) +
  facet_wrap(~ name_0, scales = "free_y") +
  labs(
    x = "Year", 
    y = "Mean Prevalence", 
    title = paste("Annual Prevalence Trends by Country for", focal_mut),
    subtitle = "Each line represents an Admin 1 region"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

time_series_plot_name = paste0("R_ignore/R_scripts/outputs/time_series/", focal_mut, "_All_Africa.png")
ggsave(filename = time_series_plot_name, plot = time_series_plot, width = 7, height = 10)

time_series_plot_large <- admin1_time_series_final |> filter(admin1_time_series_final$name_1 %in% unique(admin1_time_series_final[admin1_time_series_final$mean_prevalence >0.01 ,]$name_1)) |>
  ggplot(aes(x = t, y = mean_prevalence, 
             group = region_id,
             color = name_0,
             fill = name_0)) +
  geom_ribbon(aes(ymin = lower_ci_bound, ymax = upper_ci_bound), 
              alpha = 0.10, 
              linetype = "blank") +
  geom_line(linewidth = 0.7, alpha =1) +
  facet_wrap(~ name_0, scales = "free_y") +
  labs(
    x = "Year", 
    y = "Mean Prevalence", 
    title = paste("Annual Prevalence Trends by Country for", focal_mut, "in ADM 1 where prev is > 0.01"),
    subtitle = "Each line represents an Admin 1 region"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

#ADM1greater01 <- unique(admin1_time_series_final[admin1_time_series_final$mean_prevalence >0.01 ,]$name_0)
admin = "Uganda" 
for (admin in ADM1greater01) {
  time_series_plot_admin <- admin1_time_series_final |> filter(admin == name_0) |>
    ggplot(aes(x = t, y = mean_prevalence, 
               group = region_id,
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
      title = paste("Annual Prevalence Trends for", admin, focal_mut),
      subtitle = "Each graph is Admin 1 Region"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
}



sub_with_geom <- dat_sub %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(admin1_regions), remove = FALSE) 

# Admin 1 - First-level administrative region
coords_with_admin1 <- st_join(sub_with_geom, africa_shp_admin1, join = st_within) %>%
  mutate(
    lon = st_coordinates(.)[, 1],  # Extract longitude
    lat = st_coordinates(.)[, 2]   # Extract latitude
  ) 

coords_with_admin1 <- rename(coords_with_admin1, region_id = id_1) 
coords_with_admin1 <- coords_with_admin1 %>% filter(name_0 == "Uganda")

time_series_plot_admin_with_data <- admin1_time_series_final |> filter(admin == name_0) |>
  ggplot(aes(x = t, y = mean_prevalence, 
             group = region_id,
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
    title = paste("Annual Prevalence Trends for", admin, focal_mut, "with real data points"),
    subtitle = "Each graph is Admin 1 Region"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
