# ── Packages ────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(sf)
  library(ggplot2)
  library(units)
  library(rlang)  # for .data pronoun
  library(fs)
  library(devtools)
  library(purrr)
})

load_all()

# --- Paths & constants -----------------------------------------------------------
CACHE_DIR      <- "R_ignore/R_scripts/outputs/GRFF_kalman_cache_annual_2012_2025"
OUT_BASE       <- "GRFF_updated_maps"
OUT_PLOT_DIR <- "total_speed_spread"
dir_create(c(paste0("R_ignore/R_scripts/outputs/plots/", OUT_PLOT_DIR)))

# --- Cache helpers -----------------------------------------------------------
cache_path <- function(mut, lenS, lenT, what, ext = c("rds","parquet")) {
  ext <- match.arg(ext)
  fn  <- sprintf("%s_lenS%s_lenT%s_%s.%s",
                 gsub("[:/ ]", "_", mut), lenS, lenT, what, ext)
  file.path(CACHE_DIR, fn)
}

make_or_load <- function(mut, ell_km, tau2) {
  cached_rds <- cache_path(mut, ell_km, tau2, "full_model_output", "rds") 
  stopifnot(file_exists(cached_rds))
  readRDS(cached_rds)
}

make_raster_long <- function(arr3, xs, ys, times) {
  nx <- length(xs); ny <- length(ys)
  
  df <- do.call(rbind, lapply(seq_along(times), function(k) {
    data.frame(
      x = rep(xs, times = ny),
      y = rep(ys, each  = nx),
      p = as.vector(arr3[ , , k]),
      t = factor(times[k], levels = times)
    )
  }))
  
  df %>%
    dplyr::mutate(t = as.numeric(as.character(t))) %>%
    dplyr::filter(t %in% plot_times)
}

# --- Data I/O --------------------------------------------------------
shape_Africa <- readRDS("R_ignore/R_scripts/data/shapefiles/sf_admin0_africa.rds")
shape_water  <- sf::st_read("R_ignore/R_scripts/data/shapefiles/africa_water_bodies.shp", quiet = TRUE)

dat <- read.csv("R_ignore/R_scripts/data/all_who_get_prevalence.csv") |>
  mutate(
    collection_day  = as.Date(collection_day),
    collection_year = lubridate::year(collection_day)
  )

# --------------------------- Add K13 overall data ----------------------------
k13_any_site_year <- dat %>%
  # keep only K13 mutations; adjust this filter to your naming convention
  filter(grepl("^k13", mutation, ignore.case = TRUE)) %>%
  group_by(year, longitude, latitude) %>%
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
    prevalence = prevalence
  )

dat_with_k13 <- dat %>%
  bind_rows(k13_any_site_year) %>%
  arrange(year, longitude, latitude, mutation)

# --- Parameters ---------------------------------------------------------------
# choose a single run to (re)plot from cache
ell_km <- 80          # RFF length-scale in **kilometres**
tau2   <- 0.1         # RW1 variance in feature space
plot_times <- seq(2012, 2023, by = 1)

# prediction parameters
nx <- 200
ny <- 200

# Area parameters
q_thr <- 0.80      # exceedance-probability threshold
D     <- 1000       # resamples for ribbon

EA_long_lim <- c(15, 52)
EA_lat_lim  <- c(-12, 18)
x_range     <- c(-20, 55)  # lon extent used to subset points
y_range     <- c(-35, 38)  # not used below, but kept for reference
plot_crs = 4326

all_who_mutations <- c(
  "k13:comb","k13:622:I", "k13:469:Y", "k13:675:V", "k13:446:I", "k13:458:Y", 
  "k13:476:I",   "k13:493:H", "k13:539:T", "k13:543:T",  "k13:553:L", 
  "k13:561:H",   "k13:574:L", "k13:580:Y", "k13:441:L", "k13:449:A", 
  "k13:469:F",   "k13:481:V", "k13:515:K", "k13:527:H",  "k13:537:I", 
  "k13:537:D", "k13:538:V",  "k13:568:G")
area_ts_list <- list()
for (mut in all_who_mutations){
  print(paste0("Processing ", mut, "........."))
  # subset data for this mutation and AOI
  dat_sub <- dat_with_k13 |>
    filter(mutation == mut,
           longitude > x_range[1], longitude < x_range[2],
           year > 2011)
  
  # quick skip if no positive prevalence
  if (!any(dat_sub$prevalence > 0, na.rm = TRUE)) next
  
  # load cached model output (skip cleanly if missing)
  cached <- try(make_or_load(mut, ell_km, tau2), silent = TRUE)
  
  exceed_prob_draws <- cached$exceed_post_draw_list
  exceed_prob_draws_1 <- exceed_prob_draws$`1`
  exceed_prob_draws_5 <- exceed_prob_draws$`5`
  exceed_prob_draws_10 <- exceed_prob_draws$`10`
  
  lat_min <- cached$lat_min
  lat_max <- cached$lat_max
  lon_min <- cached$lon_min
  lon_max <- cached$lon_max
  xs <- cached$xs
  ys <- cached$ys
  plot_times <- cached$plot_times

  posterior_draws <-  cached$posterior_draws
  
  exceed_prob_draws_long_10 <- make_raster_long(exceed_prob_draws_10, xs, ys, plot_times)
  exceed_prob_draws_long_5 <- make_raster_long(exceed_prob_draws_5, xs, ys, plot_times)
  exceed_prob_draws_long_1 <- make_raster_long(exceed_prob_draws_1, xs, ys, plot_times)
  
  # Combine into a single data frame
  df_pixel_exceed_prob <- exceed_prob_draws_long_1 %>%
    rename(P_ex_p1 = p) %>%
    left_join(exceed_prob_draws_long_5 %>% rename(P_ex_p5 = p),
              by = c("x", "y", "t")) %>%
    left_join(exceed_prob_draws_long_10 %>% rename(P_ex_p10 = p),
              by = c("x", "y", "t"))
  
  # Calculate area km2 per cell
  lat_mid <- (lat_min + lat_max) / 2
  # 1 degree lat ≈ 111.32 km and 1 degree lon ≈ 111.32 * cos(latitude)
  km_per_deg_lat <- 111.32
  km_per_deg_lon <- 111.32 * cos(lat_mid * pi / 180)
  dx_deg <- (lon_max - lon_min) / (nx - 1)  # cell width in degrees lon
  dy_deg <- (lat_max - lat_min) / (ny - 1)  # cell height in degrees lat
  dx_km <- dx_deg * km_per_deg_lon
  dy_km <- dy_deg * km_per_deg_lat
  cell_area_km2 <- dx_km * dy_km
  
  #--- Area calculation based on exceedance prob --------------------------------
  df_area_over_time_exceed <- df_pixel_exceed_prob %>%
    # Identify which pixels are "high-risk" at each time point
    mutate(p1 = P_ex_p1 > q_thr,
           p5 = P_ex_p5 > q_thr,
           p10 = P_ex_p10 > q_thr) %>%
    group_by(t) %>%
    # Calculate the total high-risk area for that time point
    summarise(
      # Count the number of high-risk pixels
      n_p1 = sum(p1),
      n_p5 = sum(p5),
      n_p10 = sum(p10),
      # Calculate the total area
      total_area_p1 = n_p1 * cell_area_km2,
      total_area_p5 = n_p5 * cell_area_km2,
      total_area_p10 = n_p10 * cell_area_km2,
      .groups = "drop"
    ) %>%
    mutate(mutation = mut,
           year = as.numeric(as.character(t)))
  
  area_ts_list[[mut]] <- df_area_over_time_exceed
}

muts_of_interest <- names(area_ts_list)
my_colors <- c(
  "#332288","#117733","#44AA99","#88CCEE","#DDCC77",
  "#CC6677","#882255","#AA4499","#DDDDDD","#000000",
  "#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
  "#D55E00","#CC79A7","#9A6324","#800000"
)


area_ts_all <- bind_rows(area_ts_list) %>%
  filter(mutation %in% muts_of_interest) %>%
  mutate(
    mutation = factor(mutation, levels = muts_of_interest),
    mutation_clean = mutation %>%
      stringr::str_remove("^k13:") %>%
      stringr::str_replace_all(":", ""),
    
    mutation_clean = case_when(
      mutation == "k13:comb" ~ "k13 combined",
      TRUE ~ mutation_clean
    ))

p_area_p5 <- ggplot(area_ts_all,
                    aes(x = year, y = total_area_p5, color = mutation_clean)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = sort(unique(area_ts_all$year))) +
  labs(
    x = "Year",
    y = expression("Area with Pr(prev > 5%) ≥ 80% (km"^2*")"),
    color = "Mutation",
    title = "Total high-prevalence area over time (5% prevalence threshold)"
  ) +
  scale_color_manual(values = my_colors) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

save_figs(file.path(OUT_PLOT_DIR , "all_mutatations_total_speed"), p_area_p5)
saveRDS(area_ts_list, file.path("R_ignore/R_scripts/outputs/plots/", OUT_PLOT_DIR , "all_mutatations_total_speed.rds"))

