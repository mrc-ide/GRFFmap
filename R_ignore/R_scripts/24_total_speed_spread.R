# ── Packages ────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(ggplot2)
  library(units)
  library(rlang)  # for .data pronoun
  library(fs)
  library(devtools)
})

load_all()

# --- Paths & constants -----------------------------------------------------------
CACHE_DIR      <- "R_ignore/R_scripts/outputs/GRFF_kalman_cache_annual_2012_2025"
OUT_BASE       <- "GRFF_updated_maps"
OUT_SPEED   <- file.path("speed_estimate")
OUT_PLOT_DIR <- "total_speed_spread"
OUT_PLOT_DIR_DRAWS <- "total_speed_spread_draws"
dir_create(c(paste0("R_ignore/R_scripts/outputs/plots/", OUT_PLOT_DIR),
             paste0("R_ignore/R_scripts/outputs/plots/", OUT_PLOT_DIR_DRAWS)))

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

# --- Helpfer function --------------------------------------------------------
# Build total-union polygons per year (plus area) 
total_union_by_year <- function(polys_3857, flag_col, link_km = 0) {
  stopifnot(sf::st_crs(polys_3857)$epsg == 3857)
  link_m <- link_km * 1000
  
  polys <- polys_3857 |>
    dplyr::filter(.data[[flag_col]]) |>
    sf::st_set_precision(1)             # kill tiny slivers
  
  # optional: close small gaps to merge near-touching cells
  if (link_m > 0) polys <- polys |> sf::st_buffer(link_m)
  
  unions <- polys |>
    dplyr::group_by(t) |>
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") |>
    sf::st_make_valid()
  
  if (link_m > 0) unions <- unions |> sf::st_buffer(-link_m)
  
  unions |>
    dplyr::mutate(area_km2 = units::set_units(sf::st_area(geometry), "km^2") |> as.numeric())
}

# --- Plot in the same "blob" style ------------------------------------------
# Compute area by year, trend (km²/yr), and speed (km/yr)
process_exced <- function(mut, col_flag, label_tag){
  area_by_year <- combined_polys %>%
    filter(.data[[col_flag]]) %>%
    group_by(t) %>%
    summarise(
      geom = st_union(geometry) |> st_make_valid() |> st_buffer(0),
      .groups = "drop"
    ) %>%
    mutate(area_km2 = as.numeric(set_units(st_area(geom), "km^2")))
  
  # linear trend: km^2 per year
  fit <- lm(area_km2 ~ t, data = area_by_year)
  slope_km2_per_year <- coef(fit)[["t"]]
  
  # speed (km/year) via circular approx
  A_mean <- mean(area_by_year$area_km2, na.rm = TRUE)
  r_mean <- sqrt(A_mean / pi)
  speed_km_per_year <- slope_km2_per_year / (2 * sqrt(pi * A_mean))
  
  return(list(
    model = fit,
    slope_km2_per_year = slope_km2_per_year,
    mean_area_km2 = A_mean,
    mean_radius_km = r_mean,
    speed_km_per_year = speed_km_per_year,
    area_by_year = area_by_year
  ))
}

process_exced_draws <- function(polys_3857,           # sf with CRS=3857, has columns t, geometry, and a logical flag column
                          flag_col,             # e.g. "high_p10" or "theta_hat_ge_080"
                          link_km = 0) {        # optional: close small gaps between adjacent cells
  stopifnot(sf::st_crs(polys_3857)$epsg == 3857)
  
  link_m <- link_km * 1000
  
  # keep only flagged cells
  polys <- polys_3857 |>
    dplyr::filter(.data[[flag_col]]) |>
    sf::st_set_precision(1)
  
  # union per year, fix topology, (optionally) erode back if we dilated
  area_by_year <- polys |>
    dplyr::group_by(t) |>
    dplyr::summarise(geom = sf::st_union(geometry), .groups = "drop") |>
    sf::st_make_valid() |>
    (\(x) if (link_m > 0) sf::st_buffer(x, -link_m) else x)() |>
    dplyr::mutate(area_km2 = units::set_units(sf::st_area(geom), "km^2") |> as.numeric())
  
  # ensure numeric year
  area_by_year$t <- as.numeric(area_by_year$t)
  
  # linear trend in area (km^2 / year)
  fit <- lm(area_km2 ~ t, data = area_by_year)
  slope <- unname(coef(fit)["t"])
  
  # effective front speed (km / year) using A = π r^2 -> dr/dt = (1/(2√(π A))) dA/dt
  A_mean <- mean(area_by_year$area_km2, na.rm = TRUE)
  r_mean <- sqrt(A_mean / pi)
  speed  <- slope / (2 * sqrt(pi * A_mean))
  
  # optional: 95% CI for slope and speed (delta method on speed)
  se_slope <- sqrt(vcov(fit)["t","t"])
  ci_slope <- slope + c(-1,1) * 1.96 * se_slope
  # d speed / d slope = 1 / (2 * sqrt(pi*A_mean))
  ds_dβ <- 1 / (2 * sqrt(pi * A_mean))
  se_speed <- abs(ds_dβ) * se_slope
  ci_speed <- speed + c(-1,1) * 1.96 * se_speed
  
  list(
    model               = fit,
    slope_km2_per_year  = slope,
    slope_ci_km2_per_yr = ci_slope,
    mean_area_km2       = A_mean,
    mean_radius_km      = r_mean,
    speed_km_per_year   = speed,
    speed_ci_km_per_yr  = ci_speed,
    area_by_year        = area_by_year
  )
}

# --- area from per-cell exceedance probabilities θ̂_k -------------------------
# point area: sum over cells of 1{θ̂_k >= q_thr} * area_k
area_from_theta_hat_point <- function(theta_hat, q_thr, cell_area_km2) {
  sum((theta_hat >= q_thr) * cell_area_km2, na.rm = TRUE)
}

# ribbon via Beta posterior for each cell's θ_k (conjugate w/ D draws)
# D = number of posterior draws used to form θ̂_k = succ_k / D
area_from_theta_hat_beta <- function(theta_hat, D, q_thr, cell_area_km2, B = 500) {
  alpha <- pmax(theta_hat * D + 1, 1e-8)  # Beta(1+succ, 1+fail)
  beta  <- pmax((1 - theta_hat) * D + 1, 1e-8)
  area_b <- numeric(B)
  for (b in seq_len(B)) {
    theta_b <- rbeta(length(theta_hat), alpha, beta)
    area_b[b] <- sum((theta_b >= q_thr) * cell_area_km2, na.rm = TRUE)
  }
  c(med = median(area_b, na.rm = TRUE),
    lo  = as.numeric(quantile(area_b, 0.025, names = FALSE, na.rm = TRUE)),
    hi  = as.numeric(quantile(area_b, 0.975, names = FALSE, na.rm = TRUE)))
}

# produce unions and a faceted map for any polygon sf with a logical flag column
total_union_by_year <- function(polys_3857, flag_col, link_km = 0) {
  stopifnot(sf::st_crs(polys_3857)$epsg == 3857)
  link_m <- link_km * 1000
  polys <- polys_3857 |>
    dplyr::filter(.data[[flag_col]]) |>
    sf::st_set_precision(1)
  if (link_m > 0) polys <- polys |> sf::st_buffer(link_m)
  
  unions <- polys |>
    dplyr::group_by(t) |>
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") |>
    sf::st_make_valid()
  
  if (link_m > 0) unions <- unions |> sf::st_buffer(-link_m)
  unions |>
    dplyr::mutate(area_km2 = units::set_units(sf::st_area(geometry), "km^2") |> as.numeric())
}

plot_unions_facets <- function(unions_3857, cells_3857 = NULL, plot_crs = 4326,
                               title = "Total exceedance area over time") {
  unions_ll <- sf::st_transform(unions_3857, plot_crs)
  bb        <- sf::st_bbox(unions_ll)
  xlim      <- as.numeric(bb[c("xmin","xmax")])
  ylim      <- as.numeric(bb[c("ymin","ymax")])
  p <- ggplot()
  if (!is.null(cells_3857)) {
    p <- p + geom_sf(data = sf::st_transform(cells_3857, plot_crs),
                     fill = "grey90", color = "grey70", linewidth = 0.15)
  }
  p +
    geom_sf(data = unions_ll, fill = NA, color = "steelblue4", linewidth = 0.6) +
    geom_sf(data = unions_ll, fill = "lightsteelblue1", color = NA, alpha = 0.4) +
    geom_sf_text(data = unions_ll |> dplyr::mutate(lbl = paste0(round(area_km2), " km²")),
                 aes(label = lbl), size = 2.6) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = sf::st_crs(plot_crs)) +
    facet_wrap(~ t, nrow = 2) +
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme_bw()
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

q_thr <- 0.80      # exceedance-probability threshold
B     <- 500       # resamples for ribbon

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

cell_size_km <- 5
half_m <- (cell_size_km * 1000) / 2

results_df <- tibble(
  mutation = character(),
  exceedance = character(),
  slope_km2_per_year = numeric(),
  r_mean_km = numeric(),
  speed_km_per_year = numeric(),
  slope_km2_per_year_draws = numeric(),
  r_mean_km_draws = numeric(),
  speed_km_per_year_draws = numeric()
)

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
  
  xs <- cached$xs; ys <- cached$ys
  
  exceed_prob <- cached$exceed_post_gaussian_approx_list
  exceed_prob_draws <- cached$exceed_post_draw_list
  
  exceed_prob_1 <- exceed_prob$`1`
  exceed_prob_5 <- exceed_prob$`5`
  exceed_prob_10 <- exceed_prob$`10`
  
  exceed_prob_draws_1 <- exceed_prob_draws$`1`
  exceed_prob_draws_5 <- exceed_prob_draws$`5`
  exceed_prob_draws_10 <- exceed_prob_draws$`10`
  
  # long form exceedance (for union polygons)
  exceed_prob_long_10 <- make_raster_long(exceed_prob_10, xs, ys, plot_times)
  exceed_prob_long_5 <- make_raster_long(exceed_prob_5, xs, ys, plot_times)
  exceed_prob_long_1 <- make_raster_long(exceed_prob_1, xs, ys, plot_times)
  
  exceed_prob_draws_long_10 <- make_raster_long(exceed_prob_draws_10, xs, ys, plot_times)
  exceed_prob_draws_long_5 <- make_raster_long(exceed_prob_draws_5, xs, ys, plot_times)
  exceed_prob_draws_long_1 <- make_raster_long(exceed_prob_draws_1, xs, ys, plot_times)
  
  # Combine exceedance prob data
  combined_polys <- exceed_prob_long_10 %>%
    select(x, y, t, p) %>% rename(p10 = p) %>%
    full_join(exceed_prob_long_5  %>% select(x, y, t, p) %>% rename(p5  = p),
              by = c("x","y","t")) %>%
    full_join(exceed_prob_long_1  %>% select(x, y, t, p) %>% rename(p1  = p),
              by = c("x","y","t")) |>
    st_as_sf(coords = c("x","y"), crs = 4326, remove = FALSE) |>
    st_transform(3857) |>
  mutate(
    high_p10 = p10 > 0.80,
    high_p5  = p5  > 0.80,
    high_p1  = p1  > 0.80
  ) |>
    st_buffer(dist = half_m, endCapStyle = "SQUARE")
  
  combined_draws_polys <- exceed_prob_draws_long_10 %>%
    select(x, y, t, p) %>% rename(p10 = p) %>%
    full_join(exceed_prob_draws_long_5  %>% select(x, y, t, p) %>% rename(p5  = p),
              by = c("x","y","t")) %>%
    full_join(exceed_prob_draws_long_1  %>% select(x, y, t, p) %>% rename(p1  = p),
              by = c("x","y","t")) |>
    st_as_sf(coords = c("x","y"), crs = 4326, remove = FALSE) |>
    st_transform(3857)  |>
  mutate(
    high_p10 = p10 > 0.80,
    high_p5  = p5  > 0.80,
    high_p1  = p1  > 0.80
  ) |>
    st_buffer(dist = half_m, endCapStyle = "SQUARE")

  # precompute grid cell polygons + areas (used for ribbons)
  grid_ll  <- expand.grid(lon = xs, lat = ys)
  cells_sf <- st_as_sf(grid_ll, coords = c("lon","lat"), crs = 4326) |>
    st_transform(3857) |>
    st_buffer(dist = half_m, endCapStyle = "SQUARE")
  cell_area_km2 <- as.numeric(units::set_units(st_area(cells_sf), "km^2"))

  # Run for each exceedance flag
  area_ts_per_mut <- tibble(
    t = numeric(),
    area_pt = numeric(),
    area_med = numeric(),
    area_lo = numeric(),
    area_hi = numeric(),
    p_thr = numeric()
  )
  for (flag in c("high_p1","high_p5","high_p10")) {
    flag_plot <- switch(flag,
                        high_p1  = "Pr(Prev > 1%) ≥ 80%",
                        high_p5  = "Pr(Prev > 5%) ≥ 80%",
                        high_p10 = "Pr(Prev > 10%) ≥ 80%")
    
    p_thr <- switch(flag,
                    high_p1  = 0.01,
                    high_p5  = 0.05,
                    high_p10 = 0.1)
    
    print(paste0("...... Processing prevalence threshold ", p_thr, "........."))
    
    if (any(combined_polys[[flag]], na.rm = TRUE)) {
      # --- Area over time based on the gaussian approx exceedance prob -----------
      # compute area trend/speed plot
      total_area_time_list <- process_exced(mut, flag, flag_plot)
      area_by_year = total_area_time_list$area_by_year
      # plot (use p, and save p — not an undefined object)
      plot_area_time_gaussian <- ggplot(area_by_year, aes(x = t, y = area_km2)) +
        theme_bw() +
        geom_point() +
        geom_line() +
        geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
        scale_x_continuous(breaks = seq(min(area_by_year$t), max(area_by_year$t), by = 1)) +
        labs(
          x = "Year",
          y = "Total Area (km²)"
        ) + 
        theme(
          strip.background = element_rect(fill = "white", color = NA),
          strip.text = element_text(color = "black", face = "plain"),
          panel.border = element_rect(color = "grey80", fill = NA),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          title = element_text(size = 14),
          axis.text.x = element_text(size = 8),  # <-- key tweak
          axis.text.y = element_text(size = 8))
      save_figs(file.path(OUT_PLOT_DIR, paste0(mut, "_total_spread_exceedance_", flag)), plot_area_time_gaussian, width = 6)
      
      # --- Plot Area over time in lat long based on the gaussian approx exceedance prob -----------
      # total unions by year
      unions <- total_union_by_year(combined_polys, flag_col = flag, link_km = 0)
      cells_flagged <- combined_polys |> dplyr::filter(.data[[flag]]) |> dplyr::select(t)
      
      unions_ll <- sf::st_transform(unions, plot_crs)
      bb        <- sf::st_bbox(unions_ll)
      xlim      <- as.numeric(bb[c("xmin","xmax")])
      ylim      <- as.numeric(bb[c("ymin","ymax")])
      
      # Union (total area) layer – outline + subtle fill
      cells_ll <- sf::st_transform(cells_flagged, plot_crs)
      plot_area_time_lat_long_gaussian <- ggplot() +
        geom_sf(data = cells_ll,  fill = "grey90",        color = "grey70",      linewidth = 0.15) +
        geom_sf(data = unions_ll, fill = NA,               color = "steelblue4",  linewidth = 0.6)  +
        geom_sf(data = unions_ll, fill = "lightsteelblue1", color = NA,           alpha = 0.4)     +
        geom_sf_text(
          data = unions_ll |> dplyr::mutate(lbl = paste0(round(area_km2), " km²")),
          aes(label = lbl), size = 2.6, color = "black"
        ) +
        coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = sf::st_crs(plot_crs)) +
        facet_wrap(~ t, nrow = 2) +
        labs(title = "Total exceedance area over time", x = "Longitude", y = "Latitude") +
        theme_bw()
      
      save_figs(file.path(OUT_PLOT_DIR, paste0(mut, "_total_area_facets_", flag)),
                plot_area_time_lat_long_gaussian, width = 10)
    }
    
    if (any(combined_draws_polys[[flag]], na.rm = TRUE)) {
      # --- Draws-based area-vs-time from per-cell θ̂ (combined_draws_polys)------
      # Pull the per-cell exceedance probabilities for this threshold, per year:
      # p10/p5/p1 columns are θ̂(s,t) in [0,1]; choose the right one:
      plot_crs <- 4326
      q_thr    <- 0.80
      D_draws  <- cached$num_posterior_draws   # e.g., 1000
      B        <- 500   
      theta_col <- switch(flag, high_p1 = "p1", high_p5 = "p5", high_p10 = "p10")
      
      area_ts <- combined_draws_polys |>
        dplyr::as_tibble() |>
        dplyr::select(t, !!rlang::sym(theta_col)) |>
        dplyr::rename(theta_hat = !!rlang::sym(theta_col)) |>
        dplyr::group_by(t) |>
        dplyr::summarise(
          area_pt  = area_from_theta_hat_point(theta_hat, q_thr, cell_area_km2),
          {
            rib <- area_from_theta_hat_beta(theta_hat, D_draws, q_thr, cell_area_km2, B)
            area_med <- rib[["med"]]
            area_lo  <- rib[["lo"]]
            area_hi  <- rib[["hi"]]
            tibble(area_med, area_lo, area_hi)
          },
          .groups = "drop"
        ) |>
        dplyr::arrange(t)
      area_ts_per_mut <- rbind(area_ts_per_mut, area_ts %>% mutate(p_thr = p_thr))
      area_ts_list[[mut]] <- area_ts_per_mut
      
      # Ribbon plot (draws-based)
      area_over_time_ribbon <- ggplot(area_ts, aes(x = t)) +
        geom_ribbon(aes(ymin = area_lo, ymax = area_hi), alpha = 0.25) +
        geom_line(aes(y = area_med)) +
        geom_point(aes(y = area_med)) +
        labs(x = "Year",
             y = paste0("Area (km²) where ", flag_plot)) + 
        scale_x_continuous(breaks = seq(min(area_ts$t), max(area_ts$t), by = 1)) +
        theme_bw() + 
        theme(
          strip.background = element_rect(fill = "white", color = NA),
          strip.text = element_text(color = "black", face = "plain"),
          panel.border = element_rect(color = "grey80", fill = NA),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          title = element_text(size = 14),
          axis.text.x = element_text(size = 8),  # <-- key tweak
          axis.text.y = element_text(size = 8))
      save_figs(file.path(OUT_PLOT_DIR_DRAWS, paste0(mut, "_total_spread_exceedance_ribbon_", flag)), area_over_time_ribbon, width = 10)
      
      total_area_time_draws_list <- process_exced_draws(combined_draws_polys, flag)
      
      results_df <- rbind(results_df,
                           tibble(
        mutation = mut,
        exceedance = flag,
        slope_km2_per_year = total_area_time_list$speed_km_per_year,
        r_mean_km = total_area_time_list$r_mean_km,
        speed_km_per_year = total_area_time_list$speed_km_per_year,
        slope_km2_per_year_draws = total_area_time_draws_list$slope_km2_per_year,
        r_mean_km_draws = total_area_time_draws_list$mean_radius_km,
        speed_km_per_year_draws = total_area_time_draws_list$speed_km_per_year
      ))
    }
  }
}

write.csv(results_df, file = file.path("R_ignore/R_scripts/outputs/plots", OUT_PLOT_DIR_DRAWS, "total_speed_spread_info.csv"))
saveRDS(area_ts_list, file = file.path("R_ignore/R_scripts/outputs/plots", OUT_PLOT_DIR_DRAWS, "area_year_mutations.rds"))

mutations_to_plot <- c("k13:622:I", "k13:675:V", "k13:561:H", "k13:comb")

# choose which prevalence exceedance threshold (on p) to show
p_threshold <- 0.01     # 0.01, 0.05, or 0.10
thr_label   <- if (p_threshold==0.01) "1%" else if (p_threshold==0.05) "5%" else if (p_threshold==0.10) "10%" else paste0(100*p_threshold,"%")

# bind the three mutations into one data frame
df_plot <- do.call(dplyr::bind_rows, lapply(mutations_to_plot, function(m) {
  if (is.null(area_ts_list[[m]])) return(NULL)
  area_ts_list[[m]] %>%
    mutate(mutation = m)
})) %>%
  filter(!is.na(area_med), p_thr == p_threshold) %>%
  arrange(mutation, t)

stopifnot(nrow(df_plot) > 0)

# combined ribbon plot
p <- ggplot(df_plot, aes(x = t, group = mutation, color = mutation, fill = mutation)) +
  geom_ribbon(aes(ymin = area_lo, ymax = area_hi), alpha = 0.18, linewidth = 0, show.legend = FALSE) +
  geom_line(aes(y = area_med), linewidth = 1) +
  geom_point(aes(y = area_med), size = 1.4) +
  scale_x_continuous(breaks = seq(min(df_plot$t, na.rm = TRUE), max(df_plot$t, na.rm = TRUE), by = 1)) +
  labs(
    x = "Year",
    y = paste0("Area (km²) where Pr(Prevalence > ", thr_label, ") ≥ 80%"),
    color = "Mutation",
    fill  = "Mutation"
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "grey80", fill = NA),
    legend.position = "bottom"
  )

p


