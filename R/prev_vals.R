#' Rescaled prevalence break values for plotting
#'
#' Returns rescaled values in \eqn{[0,1]} corresponding to prevalence
#' percentage breakpoints, suitable for use with
#' \code{scale_fill_gradientn(values = ...)}.
#'
#' @param breaks Numeric vector of prevalence percentages.
#'   Defaults to \code{c(0, 1, 5, 10, 20, 30, 50, 70, 90, 95, 100)}.
#'
#' @return A numeric vector of rescaled values in \eqn{[0,1]}.
#'
#' @examples
#' pp_vals()
#'
#' @importFrom scales rescale
#'
#' @export
prev_vals <- function(
    breaks = c(0, 1, 5, 10, 20, 30, 50, 70, 90, 95, 100)
) {
  scales::rescale(breaks)
}
