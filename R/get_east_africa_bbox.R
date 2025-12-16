#' Get bounding box for East Africa
#'
#' Returns a bounding box defining a standard East Africa spatial extent,
#' suitable for cropping or plotting maps focused on the East African region.
#'
#' The bounding box spans longitudes 28.48–44.50 and latitudes −4.60–16.00,
#' and is returned as an \code{sf::st_bbox} object with the specified
#' coordinate reference system.
#'
#' @param target_crs A coordinate reference system passed to
#'   \code{sf::st_bbox()}, either as an EPSG code, a \code{crs} object
#'   (e.g., from \code{sf::st_crs()}), or any valid input accepted by \code{sf}.
#'
#' @return
#' An \code{sf::st_bbox} object representing the East Africa bounding box
#' in the specified coordinate reference system.
#'
#' @examples
#' \dontrun{
#' # Bounding box in WGS84 (longitude/latitude)
#' bbox_ea <- get_east_africa_bbox(sf::st_crs(4326))
#'
#' # Use for cropping an sf object
#' africa_ea <- sf::st_crop(africa_admin0, bbox_ea)
#' }
#'
#' @importFrom sf st_bbox st_crs
#'
#' @export
get_east_africa_bbox <- function(target_crs) {
  sf::st_bbox(
    c(xmin = 28.48, xmax = 44.5, ymin = -4.60, ymax = 16.00),
    crs = target_crs
  )
}
