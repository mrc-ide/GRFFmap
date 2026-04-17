#' Build a buffered Africa bbox, cropped shapefiles, and mask
#'
#' This helper constructs a bounding box around the African continent
#' (optionally buffered by a given distance in km), crops the supplied
#' Africa and water `sf` layers to that box, and returns a union mask
#' plus the x/y limits for plotting.
#'
#' @param africa_sf An `sf` object containing Africa land polygons
#'   (e.g. admin0 or admin1 for the African continent).
#' @param water_sf An `sf` object containing water bodies for Africa.
#' @param buf_km Numeric. Buffer (in kilometers) by which to expand the
#'   Africa bounding box in all directions. Default is 200.
#' @param crs_ll Integer or `crs` object. The longitude/latitude CRS
#'   used for plotting (default: EPSG 4326).
#' @param crs_metric Integer or `crs` object. A planar CRS used for
#'   buffering and intersections (default: EPSG 3857).
#'
#' @return A list with components:
#' \describe{
#'   \item{xlim}{Length-2 numeric vector of longitude limits (xmin, xmax).}
#'   \item{ylim}{Length-2 numeric vector of latitude limits (ymin, ymax).}
#'   \item{bbox_sfc_ll}{An `sfc` polygon giving the buffered bbox in lon/lat.}
#'   \item{africa_crop}{Cropped `sf` object of Africa polygons.}
#'   \item{water_crop}{Cropped `sf` object of water bodies.}
#'   \item{africa_mask}{Unioned `sf` polygon mask for Africa (for raster masking).}
#' }
#'
#' @examples
#' \dontrun{
#'   bbox_info <- build_africa_bbox_and_crop(shape_Africa, shape_water)
#'   xlim <- bbox_info$xlim
#'   ylim <- bbox_info$ylim
#'   africa_mask <- bbox_info$africa_mask
#' }
#'
#' @export
#' @importFrom sf st_crs st_transform st_union st_bbox st_as_sfc st_buffer
#' @importFrom sf sf_use_s2 st_make_valid st_intersection
build_africa_bbox_and_crop <- function(africa_sf,
                                       water_sf,
                                       buf_km = 200,
                                       crs_ll = 4326,
                                       crs_metric = 3857) {
  
  CRS_LL     <- sf::st_crs(crs_ll)
  CRS_METRIC <- sf::st_crs(crs_metric)
  
  # 1) transform to lon/lat only for final outputs / plotting
  africa_ll <- sf::st_transform(africa_sf, CRS_LL)
  water_ll  <- sf::st_transform(water_sf, CRS_LL)
  
  # 2) get bbox in metric CRS (no planar warning)
  africa_bbox_metric <- africa_sf |>
    sf::st_transform(CRS_METRIC) |>
    sf::st_bbox() |>
    sf::st_as_sfc() |>
    sf::st_buffer(buf_km * 1000)   # buffer in meters
  
  # back to lon/lat
  bbox_ll <- africa_bbox_metric |>
    sf::st_transform(CRS_LL)
  
  bb <- sf::st_bbox(bbox_ll)
  stopifnot(all(!is.na(bb)))
  
  xlim <- sort(as.numeric(bb[c("xmin", "xmax")]))
  ylim <- sort(as.numeric(bb[c("ymin", "ymax")]))
  
  bbox_sfc_ll <- sf::st_as_sfc(
    sf::st_bbox(
      c(
        xmin = xlim[1], xmax = xlim[2],
        ymin = ylim[1], ymax = ylim[2]
      ),
      crs = CRS_LL
    )
  )
  
  # 3) crop Africa & water within that bbox (in metric CRS, then back)
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  bbox_planar <- sf::st_transform(bbox_sfc_ll, CRS_METRIC)
  
  africa_planar <- africa_ll |>
    sf::st_transform(CRS_METRIC) |>
    sf::st_make_valid() |>
    sf::st_buffer(0) |>
    sf::st_intersection(bbox_planar)
  
  water_planar <- water_ll |>
    sf::st_transform(CRS_METRIC) |>
    sf::st_make_valid() |>
    sf::st_buffer(0) |>
    sf::st_intersection(bbox_planar)
  
  africa_crop <- sf::st_transform(africa_planar, CRS_LL)
  water_crop  <- sf::st_transform(water_planar,  CRS_LL)
  
  africa_mask <- africa_crop |>
    sf::st_union() |>
    sf::st_make_valid()
  
  list(
    xlim        = xlim,
    ylim        = ylim,
    bbox_sfc_ll = bbox_sfc_ll,
    africa_crop = africa_crop,
    water_crop  = water_crop,
    africa_mask = africa_mask
  )
}

#' Project points to a local LAEA system centered on the data
#'
#' @param df A data frame with longitude/latitude columns.
#' @param lon_col Name of longitude column (string), default "longitude".
#' @param lat_col Name of latitude column (string), default "latitude".
#' @param crs_ll CRS for lon/lat, default EPSG 4326.
#'
#' @return A list with:
#' \describe{
#'   \item{xy_km}{Matrix of projected coordinates in km (columns Xkm, Ykm).}
#'   \item{crs_laea}{The LAEA proj4 string used.}
#'   \item{lon0}{Center longitude.}
#'   \item{lat0}{Center latitude.}
#' }
#' @export
#' @importFrom sf st_as_sf st_crs st_transform st_coordinates
project_to_local_laea <- function(df,
                                  lon_col = "longitude",
                                  lat_col = "latitude",
                                  crs_ll = 4326) {
  stopifnot(lon_col %in% names(df), lat_col %in% names(df))
  
  lon <- df[[lon_col]]
  lat <- df[[lat_col]]
  
  lon0 <- stats::median(lon, na.rm = TRUE)
  lat0 <- stats::median(lat, na.rm = TRUE)
  
  crs_laea <- paste0(
    "+proj=laea +lat_0=", lat0,
    " +lon_0=", lon0,
    " +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  )
  
  pts_ll <- sf::st_as_sf(
    df,
    coords = c(lon_col, lat_col),
    crs = sf::st_crs(crs_ll)
  )
  pts_xy <- sf::st_transform(pts_ll, crs_laea)
  xy_m   <- sf::st_coordinates(pts_xy)
  xy_km  <- xy_m / 1000
  
  list(
    xy_km   = xy_km,
    crs_laea = crs_laea,
    lon0     = lon0,
    lat0     = lat0
  )
}
