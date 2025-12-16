#' Continuous prevalence color palette
#'
#' Returns the color palette used for continuous prevalence maps.
#'
#' @return A character vector of hex color codes and named colors defining
#'   the prevalence color scale, ordered from low to high prevalence.
#'
#' @examples
#' pp_cols()
#'
#' @export
prev_colors <- function() {
  c(
    "#5E3A9B",   # dark purple (0)
    "#8cc4e0",   # medium-dark blue
    "#a8ecbf",   # mint-teal
    "palegreen3",
    "khaki2",
    "#f3c86b",
    "orange",
    "hotpink",
    "red",
    "darkred"
  )
}
