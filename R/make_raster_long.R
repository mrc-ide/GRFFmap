#' Convert a 3D raster array into long-format data
#'
#' Converts a 3D array representing spatial raster values across multiple
#' time points into a long-format data frame suitable for use with
#' \code{ggplot2::geom_raster()} or similar geoms.
#'
#' The array is assumed to be indexed as \code{[nx, ny, t]}, where
#' \code{nx = length(xs)}, \code{ny = length(ys)}, and \code{t = length(times)}.
#'
#' @param arr3 A numeric 3D array of dimensions
#'   \code{[nx, ny, length(times)]}, typically representing predicted
#'   prevalence or probability surfaces over space and time.
#'
#' @param xs Numeric vector giving the x-coordinates (e.g. longitudes)
#'   corresponding to the first dimension of \code{arr3}.
#'
#' @param ys Numeric vector giving the y-coordinates (e.g. latitudes)
#'   corresponding to the second dimension of \code{arr3}.
#'
#' @param times Vector of time labels corresponding to the third dimension
#'   of \code{arr3}. These are converted into a factor to control facet
#'   ordering in plots.
#'
#' @details
#' The function expands the raster grid using \code{xs} and \code{ys},
#' flattens each spatial slice \code{arr3[, , k]}, and binds all time
#' slices into a single long data frame.
#'
#' The returned data frame contains the following columns:
#' \itemize{
#'   \item \code{x}: x-coordinate (repeated for each row)
#'   \item \code{y}: y-coordinate
#'   \item \code{p}: raster value at \code{(x, y, t)}
#'   \item \code{t}: factor indicating the time slice
#' }
#'
#' @return
#' A data frame in long format with \code{nx * ny * length(times)} rows.
#'
#' @examples
#' \dontrun{
#' p_long <- make_raster_long(
#'   arr3  = p_post_mean,
#'   xs    = xs,
#'   ys    = ys,
#'   times = plot_times
#' )
#'
#' head(p_long)
#' }
#'
#' @note
#'
#' @export
make_raster_long <- function(arr3, xs, ys, times) {
  nx <- length(xs)
  ny <- length(ys)
  
  do.call(rbind, lapply(seq_along(times), function(k) {
    data.frame(
      x = rep(xs, times = ny),
      y = rep(ys, each  = nx),
      p = as.vector(arr3[, , k]),
      t = factor(times[k], levels = times)
    )
  }))
}
