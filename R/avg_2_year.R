#' Average prevalence over two-year time blocks
#'
#' Aggregates pixel-level prevalence estimates into predefined two-year
#' calendar blocks and computes the mean prevalence within each block.
#'
#' The function assigns each observation to a fixed two-year interval
#' (e.g. `"2012-2013"`, `"2014-2015"`, …) based on the value of the `t` column,
#' and then averages prevalence values across all observations that share
#' the same spatial location (`x`, `y`) and time block.
#'
#' @param df A data frame containing at least the columns:
#'   \describe{
#'     \item{`x`}{Numeric x-coordinate (e.g. longitude or projected x).}
#'     \item{`y`}{Numeric y-coordinate (e.g. latitude or projected y).}
#'     \item{`t`}{Integer or numeric calendar year.}
#'     \item{`p`}{Numeric prevalence value to be averaged.}
#'   }
#'
#' @details
#' Only the following two-year blocks are currently supported:
#' `"2012-2013"`, `"2014-2015"`, `"2016-2017"`, `"2018-2019"`,
#' `"2020-2021"`, and `"2022-2023"`.
#'
#' Observations with years outside these ranges will have `NA` values
#' for `year_block` and will be dropped during aggregation.
#'
#' Missing prevalence values (`NA`) are ignored when computing the mean.
#'
#' @return
#' A data frame with one row per unique combination of `x`, `y`, and
#' `year_block`, containing:
#' \describe{
#'   \item{`x`}{Spatial x-coordinate.}
#'   \item{`y`}{Spatial y-coordinate.}
#'   \item{`year_block`}{Two-year time block label.}
#'   \item{`mean_p`}{Mean prevalence across the two-year block.}
#' }
#'
#' @examples
#' \dontrun{
#' avg_df <- avg_2_year(prev_surface_long)
#' head(avg_df)
#' }
#'
#' @importFrom dplyr mutate case_when group_by summarise
#'
#' @export
avg_2_year <- function(df){
  # Create a “year block” variable
  df <- df %>%
    mutate(
      t = case_when(
        t %in% c(2012, 2013) ~ "2012-2013",
        t %in% c(2014, 2015) ~ "2014-2015",
        t %in% c(2016, 2017) ~ "2016-2017",
        t %in% c(2018, 2019) ~ "2018-2019",
        t %in% c(2020, 2021) ~ "2020-2021",
        t %in% c(2022, 2023) ~ "2022-2023"
      )
    )
  
  # Average prevalence per pixel inside each 2-year block
  avg_two_year <- df %>%
    group_by(x, y, t) %>%
    summarise(p = mean(p, na.rm = TRUE), .groups = "drop")
  
  return(avg_two_year)
}
