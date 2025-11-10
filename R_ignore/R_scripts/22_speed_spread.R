# ── Packages ────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(fs)
  library(arrow)       # optional but nice for big cached tables
  library(devtools)
  library(dplyr)
  library(fs)
})

load_all()

# ── Paths & constants ───────────────────────────────────────────────────────────
CACHE_DIR      <- "R_ignore/R_scripts/outputs/GRFF_kalman_cache_annual"
OUT_BASE       <- "GRFF_updated_maps"
OUT_SPEED   <- file.path("speed_estimate")

EA_long_lim <- c(15, 52)
EA_lat_lim  <- c(-12, 18)
x_range     <- c(-20, 55)  # lon extent used to subset points
y_range     <- c(-35, 38)  # not used below, but kept for reference

# ── Cache helpers ───────────────────────────────────────────────────────────────
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
    dplyr::filter(t > 2013 & t < 2026)
}

# ---------- 0) Utils ----------
.ensure_metric <- function(r, metric_crs) {
  if (is.character(metric_crs)) terra::project(r, metric_crs) else terra::project(r, crs(metric_crs))
}

# Align a raster to a template grid (nearest neighbor)
.align_to <- function(r, template) terra::project(r, template, method = "near")

# Compute area of TRUE/1 cells (km^2) using raster mask in metric CRS
.area_km2 <- function(mask_r) {
  # expanse sums areas of non-NA cells; ensure background is NA and blob cells 1
  terra::expanse(mask_r, unit = "km")
}

# Compute “circular radius” from area
.radius_km <- function(area_km2) sqrt(area_km2 / pi)

# Compute IoU (Jaccard) between two boolean masks on the SAME grid
.iou <- function(a, b) {
  # a, b: SpatRaster with {NA,0/1}, aligned
  inter <- terra::ifel( (a == 1) & (b == 1), 1, NA )
  union <- terra::ifel( (a == 1) | (b == 1), 1, NA )
  ia <- sum(terra::expanse(inter, unit="km"), na.rm=TRUE)
  ua <- sum(terra::expanse(union, unit="km"), na.rm=TRUE)
  if (ua == 0) 0 else ia / ua
}

# Compute/attach area_km2 no matter what
.ensure_area_km2 <- function(polys, metric_crs = "EPSG:6933", info = NULL) {
  # Ensure sf
  polys_sf <- if (inherits(polys, "SpatVector")) sf::st_as_sf(polys) else polys
  
  # Join blob attributes if provided
  if (!is.null(info) && all(c("blob_id","area_km2") %in% names(info))) {
    polys_sf <- dplyr::left_join(polys_sf, info, by = "blob_id")
  }
  
  # If area_km2 still missing, compute from geometry
  if (!"area_km2" %in% names(polys_sf)) {
    p_metric <- sf::st_transform(polys_sf, metric_crs)
    polys_sf$area_km2 <- as.numeric(sf::st_area(p_metric)) / 1e6
  }
  
  polys_sf
}

# ── Data I/O ────────────────────────────────────────────────────────────────────
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
plot_times <- seq(2014, 2025, by = 1) # should be within t_vec

all_who_mutations <- c("k13:622:I", "k13:469:Y", "k13:675:V", "k13:446:I", "k13:458:Y", "k13:476:I",   "k13:493:H",   "k13:539:T",
                       "k13:543:T",  "k13:553:L",   "k13:561:H",   "k13:574:L",  "k13:580:Y",
                       "k13:441:L", "k13:449:A",   "k13:469:F",   "k13:481:V",
                       "k13:515:K", "k13:527:H",  "k13:537:I", "k13:537:D", "k13:538:V",  "k13:568:G")

focal_mut <- "k13:675:V"

for (focal_mut in all_who_mutations){
  # subset raw points shown on the plots
  dat_sub <- dat |>
    filter(mutation == focal_mut,
           longitude > x_range[1], longitude < x_range[2],
           collection_year > 2013)
  
  # ── Load cached model outputs (df_pred + meta) ──────────────────────────────────
  # (If you still need to run the model, do it in a separate script and save RDS.)
  cached <- make_or_load(focal_mut, ell_km, tau2)
  exceed_prob_10 <- cached$exceed_prob_10_perc
  exceed_prob_5 <- cached$exceed_prob_5_perc
  exceed_prob_1 <- cached$exceed_prob_1_perc
  
  xs <- cached$xs
  ys <- cached$ys
  
  exceed_prob_long_10 <- make_raster_long(exceed_prob_10, xs, ys, plot_times)
  exceed_prob_long_5 <- make_raster_long(exceed_prob_5, xs, ys, plot_times)
  exceed_prob_long_1 <- make_raster_long(exceed_prob_1, xs, ys, plot_times)
  
  combined <- exceed_prob_long_10 %>%
    select(x, y, t, p) %>% rename(p10 = p) %>%
    full_join(exceed_prob_long_5  %>% select(x, y, t, p) %>% rename(p5  = p),
              by = c("x","y","t")) %>%
    full_join(exceed_prob_long_1  %>% select(x, y, t, p) %>% rename(p1  = p),
              by = c("x","y","t"))
  
  # Build binary exceedance mask for a given year, with optional land/water masking
  mask_exceed_for_year <- function(
    df, year, col_exceed = "p5", prob_cut = 0.80,
    mask_land = NULL, mask_water = NULL,
    metric_crs = "EPSG:6933"
  ){
    stopifnot(all(c("x","y","t", col_exceed) %in% names(df)))
    dy <- df[df$t == year & is.finite(df[[col_exceed]]), c("x","y", col_exceed), drop = FALSE]
    if (!nrow(dy)) return(NULL)
    
    # build raster in lon/lat
    r_ll <- terra::rast(dy, type = "xyz", crs = "EPSG:4326"); names(r_ll) <- "val"
    
    # mask in lon/lat (CRS-safe)
    if (!is.null(mask_land)) {
      land_ll <- sf::st_transform(sf::st_as_sf(mask_land), "EPSG:4326")
      r_ll <- terra::mask(r_ll, terra::vect(land_ll))
    }
    if (!is.null(mask_water)) {
      water_ll <- sf::st_transform(sf::st_as_sf(mask_water), "EPSG:4326")
      r_ll <- terra::mask(r_ll, terra::vect(water_ll), inverse = TRUE)
    }
    
    # threshold in lon/lat → boolean mask (1/NA)
    hot_ll <- terra::ifel(r_ll >= prob_cut, 1, NA); names(hot_ll) <- "hot"
    
    # project boolean mask with nearest-neighbor to avoid interpolation
    hot_m <- terra::project(hot_ll, metric_crs, method = "near")
    hot_m
  }
  
  # Label contiguous blobs (4/8-neighborhood) and return:
  # - labeled raster
  # - table with per-blob area, centroid, etc.
  label_blobs <- function(mask, connectivity = 8, min_cells = 1) {
    stopifnot(inherits(mask, "SpatRaster"))
    
    # normalize connectivity to 4 or 8
    con <- suppressWarnings(as.integer(connectivity))
    if (is.na(con) || !(con %in% c(4L, 8L))) {
      warning(sprintf("connectivity=%s invalid; coercing to 8", as.character(connectivity)))
      con <- 8L
    }
    
    # ensure single layer
    if (terra::nlyr(mask) > 1) mask <- mask[[1]]
    
    # clean binary: 1 foreground, NA background
    r_bin <- terra::ifel(mask > 0.5, 1L, NA)
    
    # connected components
    cl <- terra::patches(r_bin, directions = con, zeroAsNA = TRUE)
    
    # ---- FIXED: freq() call (no 'useNA='); drop NA classes manually
    fq <- as.data.frame(terra::freq(cl))
    if (!nrow(fq)) return(list(lab = cl, info = NULL, polys = NULL, centroids = NULL))
    # terra::freq usually returns columns like: layer, value, count
    fq <- fq[!is.na(fq$value), c("value","count")]
    names(fq) <- c("blob_id","n_cells")
    
    # size filter
    fq <- fq[fq$n_cells >= min_cells, , drop = FALSE]
    if (!nrow(fq)) return(list(lab = cl, info = NULL, polys = NULL, centroids = NULL))
    
    # area per blob (km^2)
    cs_km2   <- terra::cellSize(cl, unit = "km")
    area_tbl <- as.data.frame(terra::zonal(cs_km2, cl, fun = "sum", na.rm = TRUE))
    names(area_tbl) <- c("blob_id","area_km2")
    
    info <- dplyr::inner_join(fq, area_tbl, by = "blob_id") |>
      dplyr::arrange(dplyr::desc(area_km2))
    
    # polygons + centroids
    polys <- terra::as.polygons(cl, dissolve = TRUE, na.out = FALSE)
    polys <- polys[!is.na(polys$patches), ]
    names(polys)[names(polys) == "patches"] <- "blob_id"
    
    centroids <- sf::st_as_sf(polys) |>
      dplyr::mutate(centroid = sf::st_centroid(geometry)) |>
      sf::st_set_geometry("centroid")
    
    list(lab = cl, info = info, polys = polys, centroids = centroids)
  }
  
  # Extract a single-blob mask given a labeled raster and the chosen blob_id
  mask_for_blob <- function(lab, blob_id) terra::ifel(lab == blob_id, 1, NA)
  
  # ---------- 1) Build aligned yearly masks & labeled blobs ----------
  prepare_yearly_blobs <- function(
    df, years, col_exceed = "p5", prob_cut = 0.80,
    mask_land = NULL, mask_water = NULL,
    metric_crs = "EPSG:6933", connectivity = 8, min_cells = 10
  ){
    # Build masks per year
    masks <- purrr::map(years, ~ mask_exceed_for_year(df, .x, col_exceed, prob_cut, mask_land, mask_water, metric_crs))
    names(masks) <- years
    
    # Per-year mask stats (how many cells & how many TRUEs?)
    mask_stats <- purrr::imap_dfr(masks, function(m, yr) {
      if (is.null(m)) return(tibble::tibble(year = as.integer(yr), n_cells = 0L, n_true = 0L))
      # Works for terra SpatRaster or matrix/logical vector — adjust if needed
      vals <- if (inherits(m, "SpatRaster")) terra::values(m, mat = FALSE) else as.logical(m[])
      tibble::tibble(
        year   = as.integer(yr),
        n_cells = length(vals),
        n_true  = sum(vals, na.rm = TRUE)
      )
    })
    
    # Template from first non-NULL mask that has any TRUE
    first_idx <- which(!purrr::map_lgl(masks, is.null))[1]
    if (is.na(first_idx)) stop("No data after thresholding in any year.")
    template <- masks[[first_idx]]
    
    # Align all masks to the template grid
    masks_aligned <- purrr::map(masks, ~ if (is.null(.x)) NULL else .align_to(.x, template))
    
    # Label blobs & collect blob counts
    blob_objs <- purrr::map(masks_aligned, ~ if (is.null(.x)) {
      list(lab=NULL, info=NULL, polys=NULL, centroids=NULL)
    } else label_blobs(.x, connectivity, min_cells))
    names(blob_objs) <- years
    
    blob_stats <- purrr::imap_dfr(blob_objs, function(b, yr) {
      tibble::tibble(
        year    = as.integer(yr),
        n_blobs = if (is.null(b$info)) 0L else nrow(b$info)
      )
    })
    
    # Quick message table
    msg <- dplyr::left_join(mask_stats, blob_stats, by = "year")
    print(msg)
    
    list(
      template = template,
      masks = masks_aligned,
      blobs = blob_objs,
      years = years,
      metric_crs = metric_crs,
      stats = msg
    )
  }
  
  # ---------- 2) Pick the single largest blob at any time ----------
  pick_master_blob <- function(prep) {
    yr <- prep$years
    areas_by_year <- purrr::map2(yr, prep$blobs, ~ if (is.null(.y$info)) NULL else dplyr::mutate(.y$info, year=.x))
    areas_df <- dplyr::bind_rows(areas_by_year)
    if (is.null(areas_df) || nrow(areas_df)==0) {
      message("No blobs found in any year (after min_cells filter).")
      return(NULL)
    }
    i <- which.max(areas_df$area_km2)
    list(master_year = areas_df$year[i], master_blob_id = areas_df$blob_id[i], master_area_km2 = areas_df$area_km2[i])
  }
  
  # ---------- 3) Link the master blob across all years ----------
  # Greedy, transparent linking by maximum IoU with the previous year's chosen mask (forward and backward).
  # Fallback: if IoU = 0 for all blobs in a year, leave NA (manual fill if desired).
  link_blob_across_years <- function(prep, master,
                                     min_iou = 0.01,     # IoU cutoff after dilation
                                     dilation_km = 25,   # how much to grow previous year's mask
                                     max_jump_km = 250)  # centroid fallback if no IoU
  {
    yrs <- prep$years
    chosen <- setNames(rep(NA_integer_, length(yrs)), yrs)
    
    y0 <- master$master_year
    chosen[as.character(y0)] <- master$master_blob_id
    last_mask <- mask_for_blob(prep$blobs[[as.character(y0)]]$lab, master$master_blob_id)
    
    pick_by_overlap_or_centroid <- function(y, ref_mask) {
      b <- prep$blobs[[as.character(y)]]
      if (is.null(b$lab) || is.null(b$info) || nrow(b$info) == 0) return(list(id = NA_integer_, mask = ref_mask))
      
      cand <- b$info$blob_id
      # IoU with dilated reference
      ref_dil <- .dilate_mask_km(ref_mask, dilation_km)
      ious <- vapply(cand, function(id) .iou(ref_dil, mask_for_blob(b$lab, id)), numeric(1))
      
      if (length(ious) && max(ious, na.rm = TRUE) >= min_iou) {
        best <- cand[which.max(ious)]
        return(list(id = best, mask = mask_for_blob(b$lab, best)))
      }
      
      # fallback: nearest centroid, within max_jump_km
      # previous centroid:
      prev_poly <- suppressWarnings(terra::as.polygons(ref_mask, dissolve = TRUE, values = FALSE)) |>
        sf::st_as_sf()
      if (!is.null(prev_poly) && nrow(prev_poly)) {
        cprev <- sf::st_centroid(prev_poly)
        ccur  <- sf::st_centroid(b$polys)
        dkm   <- as.numeric(sf::st_distance(cprev, ccur)) / 1000
        j     <- which.min(dkm)
        if (length(j) && is.finite(dkm[j]) && dkm[j] <= max_jump_km) {
          best <- b$polys$blob_id[j]
          return(list(id = best, mask = mask_for_blob(b$lab, best)))
        }
      }
      list(id = NA_integer_, mask = ref_mask)
    }
    
    # forward pass
    for (y in yrs[yrs > y0]) {
      pick <- pick_by_overlap_or_centroid(y, last_mask)
      chosen[as.character(y)] <- pick$id
      if (!is.na(pick$id)) last_mask <- pick$mask
    }
    
    # backward pass
    last_mask <- mask_for_blob(prep$blobs[[as.character(y0)]]$lab, master$master_blob_id)
    for (y in rev(yrs[yrs < y0])) {
      pick <- pick_by_overlap_or_centroid(y, last_mask)
      chosen[as.character(y)] <- pick$id
      if (!is.na(pick$id)) last_mask <- pick$mask
    }
    
    chosen
  }
  
  # ---------- 4) Build the single-blob mask series, areas & radii ----------
  single_blob_series <- function(prep, chosen_ids) {
    yr <- prep$years
    out <- map(yr, function(y) {
      b <- prep$blobs[[as.character(y)]]
      if (is.null(b$lab) || is.na(chosen_ids[[as.character(y)]])) {
        list(year=y, mask=NULL, area_km2=NA_real_, radius_km=NA_real_)
      } else {
        m <- mask_for_blob(b$lab, chosen_ids[[as.character(y)]])
        area_km2 <- sum(terra::expanse(m, unit="km"), na.rm=TRUE)
        list(year=y, mask=m, area_km2=area_km2, radius_km=.radius_km(area_km2))
      }
    })
    df <- tibble(
      year      = map_int(out, "year"),
      area_km2  = map_dbl(out, "area_km2"),
      radius_km = map_dbl(out, "radius_km")
    ) |>
      arrange(year) |>
      mutate(speed_km_per_year = c(NA, diff(radius_km) / diff(year)))
    list(series = out, table = df)
  }
  
  # Inputs
  years_vec <- 2014:2023
  
  prep <- prepare_yearly_blobs(
    df          = combined,
    years       = years_vec,
    col_exceed  = "p5",   # your exceedance column
    prob_cut    = 0.80,              # >= 80% probability
    mask_land   = shape_Africa,
    mask_water  = shape_water,
    metric_crs  = "EPSG:6933",
    connectivity= 4,
    min_cells   = 1                 # adjust speckle filter as needed
  )
  
  master <- pick_master_blob(prep)
  # master$master_year; master$master_blob_id; master$master_area_km2
  
  chosen_ids <- link_blob_across_years(
    prep, master,
    min_iou     = 1e-7,   # or 0.05
    dilation_km = 25,     # grow last-year mask by ~25 km before computing IoU
    max_jump_km = 250     # centroid fallback if still no overlap
  )
  # chosen_ids is a named integer vector: year -> blob_id (NA if no match that year)
  
  one_blob <- single_blob_series(prep, chosen_ids)
  blob_table <- one_blob$table
  save_csv(file.path(OUT_SPEED, paste0(focal_mut, "_speed_spread")), blob_table)
}




##### TESTING ######

library(ggplot2)
library(sf)
plot_year_blobs <- function(prep, year, chosen_ids = NULL, basemap = NULL, metric_crs = "EPSG:6933") {
  by <- prep$blobs[[as.character(year)]]
  if (is.null(by) || is.null(by$polys) || nrow(by$polys) == 0) {
    return(ggplot() + ggtitle(paste("No blobs —", year)))
  }
  
  # Ensure area_km2 exists
  polys <- .ensure_area_km2(by$polys, metric_crs = metric_crs, info = by$info)
  
  p <- ggplot()
  if (!is.null(basemap)) {
    p <- p + geom_sf(data = st_transform(basemap, st_crs(polys)), fill = NA, color = "grey60", linewidth = 0.3)
  }
  
  p <- p +
    geom_sf(data = polys, aes(fill = as.factor(blob_id)), alpha = 0.45, color = "white", linewidth = 0.25) +
    geom_sf_text(data = polys,
                 aes(label = paste0("ID ", blob_id, "\n", round(area_km2, 1), " km²")),
                 size = 3)
  
  # Highlight chosen blob (if any)
  if (!is.null(chosen_ids)) {
    id <- chosen_ids[[as.character(year)]]
    if (!is.na(id)) {
      sel <- polys[polys$blob_id == id, , drop = FALSE]
      if (nrow(sel)) p <- p + geom_sf(data = sel, fill = NA, color = "black", linewidth = 1.0)
    }
  }
  
  p + guides(fill = "none") + ggtitle(paste("Blobs —", year))
}

# example
plot_year_blobs(prep, 2019, chosen_ids, basemap = shape_Africa)
plot_year_blobs(prep, 2020, chosen_ids, basemap = shape_Africa)


df_ex <- df_pred %>% filter(exceed_prob_5 > 0.8)
ggplot() +
  geom_sf(data = shape_Africa, fill = NA, color = "black") +
  geom_tile(data = df_ex %>% filter(t >= 2014), aes(x = x, y = y, fill = exceed_prob_5)) +
  scale_fill_viridis_c(name = "P(>5%)", limits = c(0.8, 1)) +
  coord_sf() +
  facet_wrap(~t) +
  theme_minimal()

