#' Crop a long-format spatial data frame to x/y limits
#'
#' Filters a long-format data frame containing spatial coordinates to retain
#' only rows whose `x` and `y` values fall within specified bounds.
#'
#' This helper is typically used to crop raster-like long data frames
#' (e.g. expanded grids for plotting) prior to visualization or aggregation.
#'
#' @param df A data frame containing numeric columns `x` and `y`, representing
#'   spatial coordinates (e.g. longitude and latitude in decimal degrees, or
#'   projected coordinates).
#'
#' @param xlim Numeric vector of length two specifying the minimum and maximum
#'   values of `x` to retain, in the form `c(xmin, xmax)`.
#'
#' @param ylim Numeric vector of length two specifying the minimum and maximum
#'   values of `y` to retain, in the form `c(ymin, ymax)`.
#'
#' @details
#' Rows are kept if and only if
#' \code{x >= xlim[1]}, \code{x <= xlim[2]},
#' \code{y >= ylim[1]}, and \code{y <= ylim[2]}.
#'
#' The function assumes that `xlim` and `ylim` are ordered as
#' `(min, max)`; no reordering or validation is performed internally.
#'
#' @return
#' A data frame containing only rows of `df` that fall within the specified
#' spatial bounds. All columns are preserved.
#'
#' @examples
#' \dontrun{
#' cropped_df <- crop_long_df(
#'   df   = p_long,
#'   xlim = c(28, 45),
#'   ylim = c(-5, 15)
#' )
#'
#' nrow(cropped_df)
#' }
#'
#' @importFrom dplyr filter
#'
#' @export
crop_long_df <- function(df, xlim, ylim) {
  dplyr::filter(df, x >= xlim[1], x <= xlim[2], y >= ylim[1], y <= ylim[2])
}
