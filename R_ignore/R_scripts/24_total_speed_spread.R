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
# Total area where Pr(p > p_thr) ≥ q_thr (point estimate)
area_highconf_point <- function(Pk, p_thr, q_thr, cell_area_km2) {
  succ <- rowSums(Pk > p_thr)        # successes per cell
  Dk   <- ncol(Pk)
  theta_hat <- succ / Dk             # per-cell exceedance probability estimate
  sum((theta_hat >= q_thr) * cell_area_km2)
}

# Uncertainty via bootstrap resampling of posterior draws
area_highconf_boot <- function(Pk, p_thr, q_thr, cell_area_km2, B = 500) {
  Dk <- ncol(Pk)
  area_b <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(Dk, Dk, replace = TRUE)
    succ_b  <- rowSums(Pk[, idx, drop = FALSE] > p_thr, na.rm = TRUE)
    theta_b <- succ_b / Dk
    area_b[b] <- sum((theta_b >= q_thr) * cell_area_km2, na.rm = TRUE)
  }
  c(
    med = median(area_b, na.rm = TRUE),
    lo  = as.numeric(quantile(area_b, 0.025, names = FALSE, na.rm = TRUE)),
    hi  = as.numeric(quantile(area_b, 0.975, names = FALSE, na.rm = TRUE))
  )
}

# Faceted map of total exceedance area by year
plot_total_area_facets <- function(union_sf_3857,
                                   all_cells_3857 = NULL,    # optional light-grey context
                                   plot_crs = 4326,
                                   title_prefix = "Total exceedance area over time") {
  
  unions_ll <- sf::st_transform(union_sf_3857, plot_crs)
  bb        <- sf::st_bbox(unions_ll)
  xlim      <- as.numeric(bb[c("xmin","xmax")])
  ylim      <- as.numeric(bb[c("ymin","ymax")])
  
  p <- ggplot()
  
  if (!is.null(all_cells_3857)) {
    cells_ll <- sf::st_transform(all_cells_3857, plot_crs)
    p <- p + geom_sf(data = cells_ll, fill = "grey90", color = "grey70", linewidth = 0.15)
  }
  
  # Union (total area) layer – outline + subtle fill
  p <- p +
    geom_sf(data = unions_ll, fill = NA, color = "steelblue4", linewidth = 0.6) +
    geom_sf(data = unions_ll, fill = "lightsteelblue1", color = NA, alpha = 0.4) +
    # Area label inside each panel
    geom_sf_text(
      data = unions_ll |> dplyr::mutate(lbl = paste0(round(area_km2), " km²")),
      aes(label = lbl),
      size = 2.6, color = "black"
    ) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = sf::st_crs(plot_crs)) +
    facet_wrap(~ t, nrow = 2) +
    labs(title = title_prefix, x = "Longitude", y = "Latitude") +
    theme_bw()
  
  p
}

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
  
  # plot (use p, and save p — not an undefined object)
  p <- ggplot(area_by_year, aes(x = t, y = area_km2)) +
    theme_bw() +
    geom_point() +
    geom_line() +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    labs(
      title = paste0("Total Area where ", label_tag),
      subtitle = paste0("Trend: ", round(slope_km2_per_year, 1), " km²/yr;  Speed: ",
                        round(speed_km_per_year, 3), " km/yr"),
      x = "Year",
      y = "Total Area (km²)"
    )
  
  # append to results_df
  results_df <- bind_rows(
    results_df,
    tibble(
      mutation = mut,
      exceedance = col_flag,
      slope_km2_per_year = slope_km2_per_year,
      r_mean_km = r_mean,
      speed_km_per_year = speed_km_per_year
    )
  )
  p
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

all_who_mutations <- c(
  "k13:comb","k13:622:I", "k13:469:Y", "k13:675:V", "k13:446:I", "k13:458:Y", 
  "k13:476:I",   "k13:493:H", "k13:539:T", "k13:543:T",  "k13:553:L", 
  "k13:561:H",   "k13:574:L", "k13:580:Y", "k13:441:L", "k13:449:A", 
  "k13:469:F",   "k13:481:V", "k13:515:K", "k13:527:H",  "k13:537:I", 
  "k13:537:D", "k13:538:V",  "k13:568:G")

results_df <- tibble(
  mutation = character(),
  exceedance = character(),
  slope_km2_per_year = numeric(),
  r_mean_km = numeric(),
  speed_km_per_year = numeric()
)

area_ts_df = tibble(
  mutation = character(),
  exceedance = character(),
  t = numeric(),
  area_pt = numeric(),
  area_med = numeric(),
  area_lo = numeric(),
  area_hi = numeric()
)

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
  
  posterior_draws_list <- cached$posterior_draws
  exceed_prob_10 <- cached$exceed_prob_10_perc
  exceed_prob_5  <- cached$exceed_prob_5_perc
  exceed_prob_1  <- cached$exceed_prob_1_perc
  xs <- cached$xs; ys <- cached$ys
  
  # long form exceedance (for union polygons)
  exceed_prob_long_10 <- make_raster_long(exceed_prob_10, xs, ys, plot_times)
  exceed_prob_long_5  <- make_raster_long(exceed_prob_5,  xs, ys, plot_times)
  exceed_prob_long_1  <- make_raster_long(exceed_prob_1,  xs, ys, plot_times)
  
  # Combine exceedance prob data
  combined <- exceed_prob_long_10 %>%
    select(x, y, t, p) %>% rename(p10 = p) %>%
    full_join(exceed_prob_long_5  %>% select(x, y, t, p) %>% rename(p5  = p),
              by = c("x","y","t")) %>%
    full_join(exceed_prob_long_1  %>% select(x, y, t, p) %>% rename(p1  = p),
              by = c("x","y","t"))
  
  # Convert the data to sf squares in planar meters
  cell_size_km <- 5
  half_m <- (cell_size_km * 1000) / 2
  grid_ll <- expand.grid(lon = xs, lat = ys)
  
  combined_sf <- combined |>
    st_as_sf(coords = c("x","y"), crs = 4326, remove = FALSE) |>
    st_transform(3857)  # Web Mercator (meters)
  
  combined_polys <- combined_sf |>
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
  for (flag in c("high_p1","high_p5","high_p10")) {
    flag_plot <- switch(flag,
                        high_p1  = "Pr(Prev > 1%) \u2265 80%",
                        high_p5  = "Pr(Prev > 5%) \u2265 80%",
                        high_p10 = "Pr(Prev > 10%) \u2265 80%")
    
    p_thr <- switch(flag,
                    high_p1  = 0.01,
                    high_p5  = 0.05,
                    high_p10 = 0.1)
    
    print(paste0("...... Processing prevalence threshold ", p_thr, "........."))
    
    # compute + print trend/speed plot (process_exced updates results_df)
    area_over_time_plot <- process_exced(mut, flag, flag_plot)
    save_figs(file.path(OUT_PLOT_DIR, paste0(mut, "_total_spread_exceedance_", flag)), area_over_time_plot, width = 10)
    
    # total unions by year
    unions <- total_union_by_year(combined_polys, flag_col = flag, link_km = 0)
    cells_flagged <- combined_polys |> dplyr::filter(.data[[flag]]) |> dplyr::select(t)
    
    p_total_facets <- plot_total_area_facets(
      union_sf_3857  = unions,
      all_cells_3857 = cells_flagged,
      title_prefix   = paste0("Total exceedance area over time — ", flag_plot)
    )
    save_figs(file.path(OUT_PLOT_DIR, paste0(mut, "_total_area_facets_", flag)),
              p_total_facets, width = 10)
    
    # area vs time ribbon from posterior draws
    half_m <- (cell_size_km * 1000) / 2
    grid_ll <- expand.grid(lon = xs, lat = ys)
    cells_sf <- st_as_sf(grid_ll, coords = c("lon","lat"), crs = 4326) |>
      st_transform(3857) |>
      st_buffer(dist = half_m, endCapStyle = "SQUARE")
    cell_area_km2 <- as.numeric(set_units(st_area(cells_sf), "km^2"))
    
    # Choose one method for the ribbon: "beta" (recommended) or "bootstrap"
    ribbon_method <- "beta"
    
    area_ts <- lapply(seq_along(plot_times), function(k) {
      Pk  <- posterior_draws_list[[k]]                 # [M x D]
      pt  <- area_highconf_point(Pk, p_thr, q_thr, cell_area_km2)
      rib <- switch(ribbon_method,
                    beta      = area_highconf_beta(Pk, p_thr, q_thr, cell_area_km2, B),
                    bootstrap = area_highconf_boot(Pk, p_thr, q_thr, cell_area_km2, B),
                    stop("Unknown ribbon_method"))
      
      tibble(
        t        = plot_times[k],
        area_pt  = pt,
        area_med = rib[["med"]],
        area_lo  = rib[["lo"]],
        area_hi  = rib[["hi"]]
      )
    }) |> dplyr::bind_rows() |>
      dplyr::mutate(exceedance_prob = p_thr)
    
    area_ts_df <- bind_rows(
      area_ts_df,
      tibble(
        mutation = mut,
        exceedance = area_ts_df$exceedance_prob,
        t = area_ts_df$t,
        area_pt = area_ts_df$area_pt,
        area_med = area_ts_df$area_med,
        area_lo = area_ts_df$area_lo,
        area_hi = area_ts_df$area_hi
    ))
    
    area_over_time_ribbon <- ggplot(area_ts, aes(x = t)) +
      geom_ribbon(aes(ymin = area_lo, ymax = area_hi), alpha = 0.25) +
      geom_line(aes(y = area_med)) +
      geom_point(aes(y = area_med)) +
      labs(
        x = "Year",
        y = paste0("Area (km²) where ", flag_plot),
        title = "High-confidence exceedance area (median and 95% ribbon)"
      ) +
      theme_bw()
    save_figs(file.path(OUT_PLOT_DIR, paste0(mut, "_total_spread_exceedance_ribbon", flag)), area_over_time_ribbon, width = 10)
  }
}

write.csv(results_df, file = file.path("R_ignore/R_scripts/outputs/plots", OUT_PLOT_DIR, "total_speed_spread_info.csv"))
write.csv(area_ts_df, file = file.path("R_ignore/R_scripts/outputs/plots", OUT_PLOT_DIR, "total_speed_spread_info_draws.csv"))


