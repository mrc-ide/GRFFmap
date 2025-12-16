#' Apply a standardized ggplot theme for spatial prevalence figures
#'
#' Adds a consistent set of theme customizations to a \code{ggplot} object,
#' tailored for faceted spatial prevalence and exceedance maps. This helper
#' centralizes styling choices to ensure visual consistency across figures.
#'
#' @param p A \code{ggplot} object to which the theme customizations
#'   will be applied.
#'
#' @details
#' The theme adjustments include:
#' \itemize{
#'   \item White facet strip backgrounds with no border
#'   \item Reduced strip, axis, legend, and title font sizes
#'   \item Light grey panel borders
#'   \item Angled x-axis text for improved readability
#'   \item Centered legend justification
#' }
#'
#' This function is intended to be used as the final step in a plotting
#' pipeline, e.g. \code{plot_theme(p)}.
#'
#' @return
#' A \code{ggplot} object with the standardized theme applied.
#'
#' @examples
#' \dontrun{
#' p <- ggplot(df, aes(x, y)) +
#'   geom_point()
#'
#' plot_theme(p)
#' }
#'
#' @importFrom ggplot2 theme element_rect element_text
#'
#' @export
plot_theme <- function(p) {
  p + theme(
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 9, color = "black", face = "plain"),
    panel.border = element_rect(color = "grey80", fill = NA),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    title = element_text(size = 10),
    axis.text.x = element_text(size = 8, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 9),
    plot.title = element_text(hjust = 0),
    legend.justification = "center"
  )
}
