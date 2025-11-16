# GRFF_working_example.R
#
# Author: Bob Verity
# Date: 2025-08-01
#
# Inputs:
#   • data/shapefiles/sf_admin0_africa.rds   – SF object of African country boundaries
#   • data/shapefiles/africa_water_bodies.shp – Shapefile of major water bodies in Africa
#   • data/validated_and_candidate_get_prevalence.csv – CSV of marker prevalence
#
# Outputs:
#   • (Optional) PNG files in outputs/ illustrating prevalence surfaces for the focal mutation
#
# Purpose:
#   This script implements a Random Fourier Feature (RFF) approximation to a Gaussian Process
#   for modelling the prevalence of a specified genetic marker (e.g., k13 mutations) over
#   space and time across Africa. It performs the following steps:
#     1. Loads geographic shapefiles for plotting country borders and water bodies.
#     2. Reads and pre‐processes sample prevalence data (longitude, latitude, date, counts).
#     3. Constructs RFF feature maps in space and time to approximate the GP kernel.
#     4. Prepares and runs a Stan model to infer latent GP parameters and nugget variance.
#     5. Extracts posterior samples and computes prevalence predictions on a regular grid
#        for each observed year.
#     6. Summarizes posterior distributions (median, 95% credible intervals) and calculates
#        margin of error.
#     7. Generates and (optionally) saves spatial‐temporal raster plots of predicted prevalence,
#        faceted by year, for the focal mutation.
#
# ------------------------------------------------------------------

library(tidyverse)
library(sf)
library(rstan)
library(ggplot2)
library(fs)
library(arrow)
library(withr)

CACHE_DIR <- "R_ignore/R_scripts/outputs/GRFF_cache_group_years"
dir_create(CACHE_DIR)

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

# build consistent filenames per (mutation, length_space, length_time, n_features)
cache_path <- function(mut, lenS, lenT, what, ext = c("rds","parquet","rds.gz")) {
  ext <- match.arg(ext)
  fn <- sprintf("%s_lenS%s_lenT%s_%s.%s",
                gsub("[:/ ]", "_", mut), lenS, lenT, what, ext)
  file.path(CACHE_DIR, fn)
}

# read in shapefiles for plotting
shape_Africa <- readRDS("R_ignore/R_scripts/data/shapefiles/sf_admin0_africa.rds")
shape_water <- st_read("R_ignore/R_scripts/data/shapefiles/africa_water_bodies.shp", quiet = TRUE)

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

# which mutation we are focusing on
all_who_mutations <- c("k13:622:I", "k13:469:Y", "k13:675:V", "k13:446:I", "k13:458:Y", "k13:476:I",   "k13:493:H",   "k13:539:T",
                       "k13:543:T",  "k13:553:L",   "k13:561:H",   "k13:574:L",  "k13:580:Y",
                       "k13:441:L", "k13:449:A",   "k13:469:F",   "k13:481:V",
                       "k13:515:K", "k13:527:H",  "k13:537:I", "k13:537:D", "k13:538:V",  "k13:568:G")

all_who_mutations <- c("k13:622:I", "k13:469:Y", "k13:675:V", "k13:446:I", 
                       "k13:441:L")

all_who_mutations <- c("k13:622:I", "k13:469:Y", "k13:675:V", "k13:446:I", "k13:458:Y", "k13:476:I",   "k13:493:H",   "k13:539:T",
                       "k13:543:T",  "k13:553:L",   "k13:561:H",   "k13:574:L",  "k13:580:Y",
                       "k13:441:L", "k13:449:A",   "k13:469:F")

EA_long_lim <- c(15, 52)
EA_lat_lim <- c(-12, 18)

for (focal_mut in all_who_mutations){
  print(focal_mut)
  df <- k13_grouped %>% arrange(prevalence) |>
    filter(mutation == focal_mut)
  if(nrow(df) != 0){
    # option to plot raw data. Useful for identifying spatial range of interest
    if (TRUE) {
      k13_grouped %>% arrange(prevalence) |>
        filter(mutation == focal_mut) |>
        ggplot() + theme_bw() +
        #geom_sf(fill = NA, data = shape_Africa) +
        geom_point(aes(x = longitude, y = latitude, col = prevalence), size = 2) +
        scale_color_viridis_c(option = "magma") +
        coord_sf(xlim = c(-20, 55), ylim = c(-35, 38), expand = FALSE) +
        facet_wrap(~year_group_numeric) +
        ggtitle(focal_mut)
    }
    
    # define map range (longitude and latitude)
    x_range <- c(-20, 55)  # Covers westernmost to easternmost Africa
    y_range <- c(-35, 38)  # Covers southernmost to northernmost Africa
    
    # subset data
    dat_sub <- k13_grouped |>
      filter(mutation == focal_mut) |>
      filter((longitude > x_range[1]) & (longitude < x_range[2]))
    
    # ------------------------------------------------------------------
    
    # define model parameters
    length_space <- 75 / 110
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
    X_time <- dat_sub$year_group_numeric
    
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
    
    # Save the fit
    fit_path <- cache_path(focal_mut, length_space, length_time, "fit", ext = "rds")
    saveRDS(fit, fit_path)
    
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
    t_plot <- sort(unique(dat_sub$year_group_numeric))
    
    l <- list()
    for (i in seq_along(t_plot)) {
      message(sprintf("%s of %s", i, length(t_plot)))
      
      # fix prediction at this time
      t_pred <- rep(t_plot[i], nrow(X_pred))
      
      # generate prediction feature maps
      feat_pred <- cbind(cos(X_pred %*% omega_space / length_space + t_pred %*% omega_time / length_time),
                         sin(X_pred %*% omega_space / length_space + t_pred %*% omega_time / length_time)) / sqrt(n_features)
      
      # get predictions over all posterior draws
      post_all <- feat_pred %*% t(beta_samples) |>
        sweep(MARGIN = 2, STATS = mu_samples, FUN = "+") |>
        plogis()
      
      # Compute exceedance probabilities
      exceed_probs <- data.frame(
        exceed_prob_1 = apply(post_all, 1, function(x) mean(x > 0.01)),
        exceed_prob_5 = apply(post_all, 1, function(x) mean(x > 0.05)),
        exceed_prob_10 = apply(post_all, 1, function(x) mean(x > 0.10))
      )
      
      # summarise posterior predictions
      pred_quant <- apply(post_all, 1, function(x) quantile(x, probs = c(0.025, 0.5, 0.975))) |>
        t() |>
        as.data.frame()
      
      # get combined data.frame from inferred distribution and truth
      l[[i]] <- data.frame(x = X_pred[,1],
                           y = X_pred[,2],
                           t = t_pred) |>
        bind_cols(pred_quant) |>
        mutate(MOE = `97.5%` - `2.5%`,
               p_pred = `50%`) |>
        bind_cols(exceed_probs)
    }
    df_pred <- bind_rows(l)
    
    # R-native (simple)
    pred_path <- cache_path(focal_mut, length_space, length_time, "df_pred", ext = "rds")
    saveRDS(df_pred, pred_path)
    
    # OR fast cross-language columnar (great for big grids)
    pred_parquet <- cache_path(focal_mut, length_space, length_time, "df_pred", ext = "parquet")
    arrow::write_parquet(df_pred, pred_parquet)
    
    # Save metadata
    meta <- list(
      mutation = focal_mut,
      length_space = length_space,
      length_time  = length_time,
      n_features   = n_features,
      x_range = x_range,
      y_range = y_range,
      EA_long_lim = EA_long_lim,
      EA_lat_lim  = EA_lat_lim,
      years = sort(unique(dat_sub$year_group_numeric))
    )
    meta_path <- cache_path(focal_mut, length_space, length_time, "meta", ext = "rds")
    saveRDS(meta, meta_path)
  }
}

