#' Mask prevalence values to the African continent
#'
#' Sets prevalence values to `NA` for observations whose spatial coordinates
#' fall outside a supplied Africa polygon mask.
#'
#' The function converts input coordinates to an `sf` point object and uses
#' spatial intersection to determine whether each point lies within the
#' African continent (or any polygon provided in `africa_mask`). Values of
#' `p` corresponding to points outside the mask are replaced with `NA`.
#'
#' @param df A data frame containing at least the columns:
#'   \describe{
#'     \item{`x`}{Longitude in decimal degrees (EPSG:4326).}
#'     \item{`y`}{Latitude in decimal degrees (EPSG:4326).}
#'     \item{`p`}{Numeric prevalence value to be masked.}
#'   }
#'
#' @param africa_mask An `sf` object containing polygon geometry defining the
#'   African continent (e.g. admin0 boundaries). The geometry must be in the
#'   same coordinate reference system as `df` (typically EPSG:4326).
#'
#' @details
#' Spatial containment is evaluated using `sf::st_intersects()` with
#' `sparse = TRUE` (the default). A point is considered inside Africa if it
#' intersects at least one polygon in `africa_mask`.
#'
#' Only the `p` column is modified; all other columns are returned unchanged.
#' The number and order of rows in `df` are preserved.
#'
#' @return
#' A data frame identical to `df`, except that `p` is set to `NA_real_` for
#' rows whose coordinates fall outside `africa_mask`.
#'
#' @examples
#' \dontrun{
#' masked_df <- mask_to_africa(
#'   df = prevalence_grid,
#'   africa_mask = africa_admin0
#' )
#'
#' summary(masked_df$p)
#' }
#'
#' @importFrom sf st_as_sf st_intersects
#'
#' @export
mask_to_africa <- function(df, africa_mask) {
  old_s2 <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  
  pts <- sf::st_as_sf(df, coords = c("x","y"), crs = 4326, remove = FALSE)
  inside <- lengths(sf::st_intersects(pts, africa_mask)) > 0  # sparse = TRUE by default
  df$p[!inside] <- NA_real_
  df
}
