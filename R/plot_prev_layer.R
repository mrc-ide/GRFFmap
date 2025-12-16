#' Plot a faceted prevalence raster with sf overlays and optional point observations
#'
#' Creates a raster map of predicted prevalence (in percent) faceted by `t`,
#' overlaid with an `sf` land outline and a water layer. The plot is cropped to
#' a user-supplied bounding box and uses axis tick breaks computed from that
#' bounding box. Optionally overlays point observations (e.g., measured site
#' prevalence) and optionally shows a legend.
#'
#' @param p_long_df A data frame in long format containing raster/grid values to
#'   plot. Must contain at least:
#'   \describe{
#'     \item{`x`}{Longitude (degrees; typically EPSG:4326).}
#'     \item{`y`}{Latitude (degrees; typically EPSG:4326).}
#'     \item{`p`}{Prevalence on the grid as a proportion in \eqn{[0,1]}; multiplied by 100 for plotting.}
#'     \item{`t`}{Facet variable (e.g., year or year-block); used in `facet_wrap(~t)`.}
#'   }
#' @param title_text Character string; plot title. If `""`, no title is added.
#' @param shp An `sf` object providing land boundaries (outline only). Defaults
#'   to `shape_Africa`.
#' @param shp_water An `sf` object providing water polygons/background. Defaults
#'   to `shape_water`.
#' @param lims A list defining the spatial extent. Expected structure:
#'   \code{list(xlim = c(xmin, xmax), ylim = c(ymin, ymax))}. The function uses
#'   \code{lims[[1]]} and \code{lims[[2]]} to extract x/y limits, respectively.
#' @param facet_n_row Integer; number of rows used in `facet_wrap()` for the
#'   `t` panels. Default is `1`.
#' @param x_axis_break Numeric; spacing (in degrees) between longitude axis tick
#'   marks. Default is `5`.
#' @param y_axis_break Numeric; spacing (in degrees) between latitude axis tick
#'   marks. Default is `5`.
#' @param points_df Optional data frame of point observations to overlay. If not
#'   `NULL` and has at least one row, points are added as filled circles.
#'   Must include columns `x`, `y`, and `p` where `p` is in percent \eqn{[0,100]}
#'   (to match the fill scale). Default is `NULL`.
#' @param add_legend Logical; if `TRUE`, show the fill legend at the bottom. If
#'   `FALSE` (default), hide the legend.
#'
#' @details
#' The raster layer maps `fill` to `p * 100` and uses a continuous gradient
#' defined by the global objects `pp_cols` and `pp_vals`, with limits fixed to
#' `0`--`100` and `na.value = "white"`.
#'
#' Axis tick locations are computed from the bounding box (`lims`) by rounding
#' down/up to the nearest multiples of `x_axis_break` and `y_axis_break`.
#'
#' Spatial cropping is performed via `ggplot2::coord_sf()` using the provided
#' `lims`. The coordinate reference system is treated as EPSG:4326 for both the
#' raster coordinates and the `sf` overlays.
#'
#' If `points_df` is provided, points are sorted by `p` to improve overplotting.
#' Points are drawn with `shape = 21` so that `fill` is mapped to `p`.
#'
#' Final styling is applied by `plot_theme()`, expected to be available in the
#' package namespace or session.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' mut_lims <- list(xlim = c(28, 37), ylim = c(-4, 5))
#' p <- plot_prev_layer(
#'   p_long_df   = p_long_mean_avg_2y,
#'   title_text  = "k13 675V",
#'   shp         = shape_Africa_crop,
#'   shp_water   = shape_water_crop,
#'   lims        = mut_lims,
#'   x_axis_break = 2,
#'   y_axis_break = 2,
#'   points_df   = points_df_2y,
#'   add_legend  = TRUE
#' )
#' p
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_raster geom_sf geom_point facet_wrap
#' @importFrom ggplot2 coord_sf scale_fill_gradientn scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 labs theme
#' @importFrom grid unit
#'
#' @export
plot_prev_layer <- function(p_long_df, 
                            title_text, 
                            shp = shape_Africa, 
                            shp_water = shape_water, 
                            lims,
                            facet_n_row = 1,
                            x_axis_break = 5,
                            y_axis_break = 5,
                            points_df = NULL, 
                            add_legend = FALSE) {
  
  xmin <- lims[[1]][1]; xmax <- lims[[1]][2]
  ymin <- lims[[2]][1]; ymax <- lims[[2]][2]
  
  x_breaks <- seq(floor(xmin / x_axis_break) * x_axis_break,
                  ceiling(xmax / x_axis_break) * x_axis_break,
                  by = x_axis_break)
  
  y_breaks <- seq(floor(ymin / y_axis_break) * y_axis_break,
                  ceiling(ymax / y_axis_break) * y_axis_break,
                  by = y_axis_break)
  
  p <- ggplot() +
    geom_raster(aes(x = x, y = y, fill = p*100), data = p_long_df) +
    geom_sf(data = shp, linewidth = 0.2, fill = NA, color = "white") +
    geom_sf(data = shp_water, fill = "white", colour = NA) +
    facet_wrap(~t, nrow = facet_n_row) +
    coord_sf(
      xlim = c(xmin, xmax),
      ylim = c(ymin, ymax),
      expand = FALSE,
      crs = sf::st_crs(4326),
      default_crs = sf::st_crs(4326),
      clip = "on"
    ) +
    scale_fill_gradientn(colours = pp_cols,
                         values  = pp_vals,
                         limits = c(0, 100),
                         name = "Prevalence (%)",
                         breaks = seq(0, 100, by = 10),
                         na.value = "white") +
    labs(x = "Longitude", y = "Latitude") +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks)
  
  # Optional add data points
  if (!is.null(points_df) && nrow(points_df) > 0) {
    p <- p + geom_point(
      data = points_df %>% arrange(p),
      aes(x = x, y = y, fill = p),
      shape = 21, colour = "darkgrey", size = 1.5, stroke = 0.2
    )
  }
  
  # Optional add legend
  if (!add_legend) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + theme(legend.position = "bottom",
                   legend.key.width = unit(2, "cm"),
                   legend.key.height = unit(0.2, "cm"))
  }
  
  # Optional add plot title
  if (title_text != "") p <- p + labs(title = title_text)
  
  plot_theme(p)
}
