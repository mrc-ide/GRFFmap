#' Plot an exceedance probability raster faceted by year block
#'
#' Creates a faceted raster map of exceedance probability (in percent) over a
#' spatial grid, overlaid with Africa boundaries and a water layer. The raster
#' is cropped to the current global `xlim`/`ylim` extent and faceted by
#' `t`. A viridis color scale (mako) is used for the exceedance
#' probability surface.
#'
#' @param exceed_prob A data frame in long format containing raster/grid values
#'   to plot. Must include at least:
#'   \describe{
#'     \item{`x`}{Longitude (or x-coordinate) in degrees (EPSG:4326).}
#'     \item{`y`}{Latitude (or y-coordinate) in degrees (EPSG:4326).}
#'     \item{`p`}{Exceedance probability on the grid (proportion in
#'       \eqn{[0,1]}), which is multiplied by 100 for plotting.}
#'     \item{`t`}{Facet variable (e.g., `"2012-2013"`, …).}
#'   }
#'
#' @param title_text Character string; plot title. If `""`, no title is added.
#'
#' @param legend_title Character string; label for the fill legend.
#'
#' @param shp An `sf` object providing Africa land boundaries (outline only).
#'   Defaults to `shape_Africa`.
#'
#' @param shp_water An `sf` object providing water polygons/background.
#'   Defaults to `shape_water`.
#'
#' @param add_points Logical; reserved for future use. Currently not used.
#'
#' @param add_legend Logical; reserved for future use. Currently not used (the
#'   legend is always shown at the bottom).
#'
#' @details
#' The plot uses `ggplot2::geom_raster()` with `fill = p * 100` and
#' `ggplot2::scale_fill_viridis_c(option = "mako")`, with fixed limits of
#' `0`–`100` percent and `na.value = "white"`.
#'
#' Spatial limits are set via `ggplot2::coord_sf()` using global objects `xlim`
#' and `ylim`. The CRS is treated as EPSG:4326 for both raster coordinates and
#' `sf` overlays.
#'
#' If the global object `mut` equals `"k13:561:H"`, x-axis tick density is
#' increased to every 2 degrees.
#'
#' A final styling helper `plot_theme()` is applied and is expected to exist in
#' the package namespace or session.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' p <- plot_exceedance(
#'   exceed_prob  = exceed_long_df,
#'   title_text   = "Exceedance probability",
#'   legend_title = "Pr(p > threshold) (%)",
#'   shp          = shape_Africa,
#'   shp_water    = shape_water
#' )
#' p
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_raster geom_sf coord_sf facet_wrap
#' @importFrom ggplot2 annotate theme_bw labs scale_fill_viridis_c theme
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom grid unit
#' @importFrom sf st_crs
#'
#' @export
plot_exceedance <- function(exceed_prob,
                            title_text,
                            legend_title,
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
    theme_bw() +
    annotate(
      "rect",
      xmin = xmin, xmax = xmax,
      ymin = ymin, ymax = ymax,
      fill = "white", colour = NA
    ) +
    geom_raster(aes(x = x, y = y, fill = p * 100), data = exceed_prob) +
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
    scale_fill_viridis_c(
      option = "mako",
      direction = 1,
      limits = c(0, 100),
      breaks = seq(0, 100, by = 10),
      name = legend_title,
      na.value = "white"
    ) +
    labs(x = "Longitude", y = "Latitude") +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks)
  
  if (title_text != "") p <- p + labs(title = title_text)
  
  if (!add_legend) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + theme(
      legend.position = "bottom",
      legend.key.width = unit(2, "cm"),
      legend.key.height = unit(0.2, "cm")
    )
  }
  
  plot_theme(p)
}
