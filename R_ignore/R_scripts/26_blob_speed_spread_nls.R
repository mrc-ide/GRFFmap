# ── Packages ────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(sf)
  library(ggplot2)
  library(units)
  library(rlang)  # for .data pronoun
  library(fs)
})

# load_all()

# --- Paths & constants -----------------------------------------------------------
CACHE_DIR      <- "R_ignore/R_scripts/outputs/GRFF_kalman_cache_annual_2012_2025"
OUT_PLOT_DIR       <- "largest_blob_speed_nlss"
OUT_PLOT_BLOB <- "blob_plots"
dir_create(c(paste0("R_ignore/R_scripts/outputs/plots/", OUT_PLOT_DIR),
             paste0("R_ignore/R_scripts/outputs/plots/", OUT_PLOT_DIR, "/", OUT_PLOT_BLOB)))

EA_long_lim <- c(15, 52)
EA_lat_lim  <- c(-12, 18)
x_range     <- c(-20, 55)  # lon extent used to subset points
y_range     <- c(-35, 38)  # not used below, but kept for reference

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
# From flagged squares -> dissolve per year -> split into components with areas
components_by_year <- function(polys_3857, flag_col, link_km = 0) {
  stopifnot(sf::st_crs(polys_3857)$epsg == 3857)  # meters
  link_m <- link_km * 1000
  
  polys <- polys_3857 |>
    dplyr::filter(.data[[flag_col]]) |>
    sf::st_set_precision(1)  # 1 m precision to kill micro slivers
  
  # optional “closing” to bridge gaps < link_km
  if (link_m > 0) polys <- polys |> sf::st_buffer(link_m)
  
  comps <- polys |>
    dplyr::group_by(t) |>
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") |>
    sf::st_make_valid()
  
  if (link_m > 0) comps <- comps |> sf::st_buffer(-link_m)
  
  comps |>
    sf::st_cast("MULTIPOLYGON", warn = FALSE) |>
    sf::st_cast("POLYGON", warn = FALSE) |>
    dplyr::mutate(component_id = dplyr::row_number(),
                  area_km2 = units::set_units(sf::st_area(geometry), "km^2") |> as.numeric())
}

# 2) Choose the global anchor blob = largest component across all years
pick_anchor_blob <- function(components_sf) {
  # returns: list(anchor_centroid (POINT), anchor_t (year), anchor_area_km2)
  anchor_row <- components_sf |>
    arrange(desc(area_km2)) |>
    slice(1)
  
  list(
    anchor_centroid = st_centroid(anchor_row$geometry),
    anchor_t        = anchor_row$t[[1]],
    anchor_area_km2 = anchor_row$area_km2[[1]]
  )
}

# 3) For each year, pick the component that contains the anchor centroid
extract_series_for_anchor <- function(components_sf, anchor_pt) {
  years <- sort(unique(components_sf$t))
  rows <- lapply(years, function(tt) {
    comps_t <- filter(components_sf, t == tt)
    if (nrow(comps_t) == 0) {
      return(tibble(t = tt, area_km2 = NA_real_))
    }
    ix <- lengths(st_contains(comps_t, anchor_pt)) > 0
    if (!any(ix)) {
      # no component contains the anchor point this year
      tibble(t = tt, area_km2 = NA_real_)
    } else {
      # if multiple contain (rare), take the largest
      comps_pick <- comps_t[ix, ] |> arrange(desc(area_km2)) |> slice(1)
      tibble(t = tt, area_km2 = comps_pick$area_km2)
    }
  })
  bind_rows(rows)
}

# 4) Fit slope (km^2/yr) and convert to radial speed (km/yr)
fit_speed <- function(area_by_year) {
  df <- filter(area_by_year, is.finite(area_km2))
  if (nrow(df) < 2) return(NULL)
  fit <- lm(area_km2 ~ t, data = df)
  slope_km2_per_year <- unname(coef(fit)[["t"]])
  A_mean <- mean(df$area_km2, na.rm = TRUE)
  r_mean <- sqrt(A_mean / pi)
  speed_km_per_year <- slope_km2_per_year / (2 * sqrt(pi * A_mean))
  list(fit = fit,
       slope_km2_per_year = slope_km2_per_year,
       r_mean_km = r_mean,
       speed_km_per_year = speed_km_per_year)
}

# --- Data I/O --------------------------------------------------------
shape_Africa <- readRDS("R_ignore/R_scripts/data/shapefiles/sf_admin0_africa.rds")
shape_water  <- sf::st_read("R_ignore/R_scripts/data/shapefiles/africa_water_bodies.shp", quiet = TRUE)

dat <- read.csv("R_ignore/R_scripts/data/all_who_get_prevalence.csv") |>
  mutate(
    collection_day  = as.Date(collection_day),
    collection_year = lubridate::year(collection_day)
  )

# choose a single run to (re)plot from cache
ell_km <- 80          # RFF length-scale in **kilometres**
tau2   <- 0.1         # RW1 variance in feature space
plot_times <- seq(2012, 2023, by = 1) # should be within t_vec

results_df <- tibble(
  mutation = character(),
  slope_km2_per_year = numeric(),
  r_mean_km = numeric(),
  speed_km_per_year = numeric()
)

all_who_mutations <- c("k13:622:I", "k13:469:Y", "k13:675:V", "k13:446:I", "k13:458:Y", "k13:476:I",   "k13:493:H",   "k13:539:T",
                       "k13:543:T",  "k13:553:L",   "k13:561:H",   "k13:574:L",  "k13:580:Y",
                       "k13:441:L", "k13:449:A",   "k13:469:F",   "k13:481:V",
                       "k13:515:K", "k13:527:H",  "k13:537:I", "k13:537:D", "k13:538:V",  "k13:568:G")

for (mut in all_who_mutations){
  print(mut)
  dat_sub <- dat |>
    filter(mutation == mut,
           longitude > x_range[1], longitude < x_range[2],
           collection_year > 2011)
  
  # ── Load cached model outputs (df_pred + meta) ──────────────────────────────────
  # (If you still need to run the model, do it in a separate script and save RDS.)
  cached <- make_or_load(mut, ell_km, tau2)
  
  exceed_prob_draws <- cached$exceed_post_draw_list
  exceed_prob_draws_1 <- exceed_prob_draws$`1`
  exceed_prob_draws_5 <- exceed_prob_draws$`5`
  exceed_prob_draws_10 <- exceed_prob_draws$`10`
  
  xs <- cached$xs
  ys <- cached$ys
  
  exceed_prob_long_10 <- make_raster_long(exceed_prob_draws_10, xs, ys, plot_times)
  exceed_prob_long_5 <- make_raster_long(exceed_prob_draws_5, xs, ys, plot_times)
  exceed_prob_long_1 <- make_raster_long(exceed_prob_draws_1, xs, ys, plot_times)
  
  # Combine exceedance prob data
  combined <- exceed_prob_long_10 %>%
    select(x, y, t, p) %>% rename(p10 = p) %>%
    full_join(exceed_prob_long_5  %>% select(x, y, t, p) %>% rename(p5  = p),
              by = c("x","y","t")) %>%
    full_join(exceed_prob_long_1  %>% select(x, y, t, p) %>% rename(p1  = p),
              by = c("x","y","t"))
  
  # Convert the data in polygons
  cell_size_km <- 5
  half_m <- (cell_size_km * 1000) / 2
  
  combined_polys <- combined |>
    st_as_sf(coords = c("x","y"), crs = 4326, remove = FALSE) |>
    st_transform(3857) |>
    mutate(
      high_p10 = p10 > 0.80,
      high_p5  = p5  > 0.80,
      high_p1  = p1  > 0.80
    ) |>
    st_buffer(dist = half_m, endCapStyle = "SQUARE")  # square cells
  
  # --- Fit NLS model per threshold on area of largest blob per year ------------
  total_area_year   <- list()
  model_fits_draws  <- list()
  min_t_per_thr     <- list() 
  for (p_thr in c(0.01, 0.05, 0.10)) {
    print(paste0("    ", p_thr))
    
    flag <- switch(as.character(p_thr),
                        "0.01" = "high_p1",
                        "0.05" = "high_p5",
                        "0.1"  = "high_p10")
    p_thr_col <- switch(as.character(p_thr),
                        "0.01" = "area_p1",
                        "0.05" = "area_p5",
                        "0.1"  = "area_p10")
    flag_plot <- switch(as.character(p_thr),
                        "0.01" = "Pr(Prev > 1%) >= 80%",
                        "0.05" = "Pr(Prev > 5%) >= 80%",
                        "0.1"  = "Pr(Prev > 10%) >= 80%")
    
    # --- Identify largest blob over time --------------------------------------
    # Components per year with areas
    comps <- components_by_year(combined_polys, flag, link_km = 8)
    link_km = 20
    link_m <- link_km * 1000
    # Combine individual polygons per year
    polys <- combined_polys |>
      dplyr::filter(.data[[flag]]) |>
      sf::st_set_precision(1)  # 1 m precision to kill micro slivers
    
    # “closing” to bridge gaps < link_km
    if (link_m > 0) polys <- polys |> sf::st_buffer(link_m)
    
    # Calculate the area for each large comp
    comps <- polys |>
      group_by(t) |>
      summarise(geometry = sf::st_union(geometry), .groups = "drop") |>
      st_make_valid() |>
      st_cast("MULTIPOLYGON") |>
      st_cast("POLYGON") |>
      mutate(
        component_id = row_number(),
        area_km2 = set_units(st_area(geometry), "km^2") |> as.numeric()
      )
    
    # Find largest anchor blob
    anchor_row <- comps |>
      arrange(desc(area_km2)) |>
      slice(1)
    
    anchor_centroid = st_centroid(anchor_row$geometry)
    anchor_t        = anchor_row$t[[1]]
    anchor_area_km2 = anchor_row$area_km2[[1]] 
    
    # Label which polygon in each year contains the anchor
    tracked <- comps |>
      dplyr::group_by(t) |>
      dplyr::mutate(is_tracked = lengths(sf::st_contains(geometry, anchor_centroid)) > 0) |>
      dplyr::ungroup() |>
      dplyr::filter(is_tracked)
    
    # Track that same blob over time via anchor centroid
    years <- sort(unique(comps$t))
    rows <- lapply(years, function(tt) {
      comps_t <- filter(comps, t == tt)
      if (nrow(comps_t) == 0) {
        return(tibble(t = tt, area_km2 = NA_real_))
      }
      ix <- lengths(st_contains(comps_t, anchor_centroid)) > 0
      if (!any(ix)) {
        # no component contains the anchor point this year
        tibble(t = tt, area_km2 = NA_real_)
      } else {
        # if multiple contain (rare), take the largest
        comps_pick <- comps_t[ix, ] |> arrange(desc(area_km2)) |> slice(1)
        tibble(t = tt, area_km2 = comps_pick$area_km2)
      }
    })
    key <- as.character(p_thr)
    total_area_year[[key]] <- bind_rows(rows)
    
    # --- Plot the largest blob over time with the area --------------------------
    # Use all-year extent of the components as the plotting window
    comps_ll   <- sf::st_transform(comps, 4326)
    tracked_ll <- sf::st_transform(tracked, 4326)
    anchor_ll  <- sf::st_transform(anchor_centroid, 4326)
    
    bb  <- sf::st_bbox(comps_ll)
    xlim <- as.numeric(bb[c("xmin","xmax")])
    ylim <- as.numeric(bb[c("ymin","ymax")])
    
    p_blob <- ggplot() +
      # all components (context)
      geom_sf(data = comps_ll, fill = "grey90", color = "grey70", linewidth = 0.15) +
      # tracked blob outline + subtle fill
      geom_sf(data = tracked_ll, fill = NA, color = "firebrick", linewidth = 0.6) +
      geom_sf(data = tracked_ll, fill = "mistyrose", color = NA, alpha = 0.4) +
      # anchor centroid
      geom_sf(data = anchor_ll, shape = 4, size = 2, linewidth = 0.5, color = "firebrick") +
      # (optional) show area label inside the tracked blob
      geom_sf_text(
        data = tracked_ll |> dplyr::mutate(lbl = paste0(round(area_km2), " km²")),
        aes(label = lbl),
        size = 2.6, color = "black"
      ) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = sf::st_crs(4326)) +
      facet_wrap(~ t, nrow = 2) +
      labs(title = paste0("Tracked largest blob over time for ", flag_plot),
           x = "Longitude", y = "Latitude") +
      theme_bw()
    
    print(p_blob)
    
    #save_figs(file.path(OUT_PLOT_DIR, OUT_PLOT_BLOB, paste0(mut, flag, "_speed_spread_time_largest_blob")), p_blob, width = 10)
    
    # --- Fit model --------------------------
    y <- total_area_year[[key]]$area_km2
    
    # per-threshold start year (first year with finite area)
    t_start <- min(total_area_year[[key]]$t[is.finite(y)], na.rm = TRUE)
    t_centered <- total_area_year[[key]]$t - t_start
    min_t_per_thr[[key]] <- t_start
    
    A_max <- max(y, na.rm = TRUE)
    t_centered <- total_area_year[[key]]$t - min(total_area_year[[key]]$t, na.rm = TRUE)
    df_tmp <- data.frame(t = total_area_year[[key]]$t, t_centered = t_centered, y = y)
    start_list <- list(
      a  = sqrt(A_max / pi) * 1.05,  # radius s.t. pi*a^2 ≈ max area
      b  = 0.1,
      t0 = 0
    )
    
    # Build formula: total_area_pX ~ (pi * a^2) * (1 - exp(-b * (t - t0)))^2
    fit_formula <- as.formula(
      "y ~ pi * a^2 * (1 - exp(-b * (t_centered - t0)))^2"
    )
    
    fit_model <- nls(
      fit_formula,
      data    = df_tmp,
      start   = start_list,
      algorithm = "port",
      lower   = c(a = 0,   b = 0,   t0 = -5),
      upper   = c(a = Inf, b = 10,  t0 = max(t_centered) + 5),
      control = nls.control(maxiter = 500, minFactor = 1e-10, warnOnly = TRUE)
    )
    
    model_fits_draws[[key]] <- fit_model
    }
    
    # Combine total are over time for each threshold into one df
    combined_area_year <- total_area_year |>
      imap(~ mutate(.x, p_thr = .y)) |>
      list_rbind() |>
      mutate(p_thr = paste0("area_", p_thr)) |>
      pivot_wider(
        names_from  = p_thr,
        values_from = area_km2
      ) |>
      arrange(t)
    
    # Long format for observed areas
    df_area_exceed_long <- combined_area_year %>%
      select(t, area_0.01, area_0.05, area_0.1) %>%
      pivot_longer(
        cols      = starts_with("area_"),
        names_to  = "threshold",
        values_to = "total_area"
      ) %>%
      mutate(
        threshold = recode(
          threshold,
          "area_0.01" = "p > 0.01",
          "area_0.05" = "p > 0.05",
          "area_0.1"  = "p > 0.10"
        )
      )
    
    # Prediction from model fit
    predictions_exceed_long <- imap_dfr(
      model_fits_draws,
      ~ {
        key <- .y              # "0.01", "0.05", "0.1"
        fit <- .x
        
        # observed years for this threshold
        obs <- total_area_year[[key]]
        t_start <- min_t_per_thr[[key]]
        t_end   <- max(obs$t[is.finite(obs$area_km2)], na.rm = TRUE)
        
        t_seq <- seq(t_start, t_end, length.out = 100)
        
        newdat <- data.frame(
          t          = t_seq,
          t_centered = t_seq - t_start
        )
        
        tibble(
          t        = t_seq,
          area_fit = as.numeric(predict(fit, newdata = newdat)),
          threshold = dplyr::recode(
            key,
            "0.01" = "p > 0.01",
            "0.05" = "p > 0.05",
            "0.1"  = "p > 0.10"
          )
        )
      }
    )
    
    nls_model_fit <- ggplot() +
      # observed points + lines
      geom_point(
        data = df_area_exceed_long,
        aes(x = t, y = total_area, color = threshold),
        size = 1.5
      ) +
      geom_line(
        data = df_area_exceed_long,
        aes(x = t, y = total_area, color = threshold),
        linewidth = 0.8
      ) +
      # fitted curves (dashed)
      geom_line(
        data = predictions_exceed_long,
        aes(x = t, y = area_fit, color = threshold),
        linewidth = 0.9,
        linetype = "dashed"
      ) +
      labs(
        x     = "Year",
        y     = "Total high-risk area (km²)",
        color = "Prevalence threshold",
        title = "NLS Model fit for largest blob across time"
      ) +
      theme_bw() +
      theme(
        panel.border    = element_rect(color = "grey80", fill = NA),
        legend.position = "bottom"
      )
    save_figs(file.path(OUT_PLOT_DIR, paste0(mut, flag, "_speed_spread_time_largest_blob_nls")), nls_model_fit, width = 10)
}
