#' Compute 2-year averaged prevalence surfaces using partner-drug year cutoffs
#'
#' Aggregates yearly prevalence estimates into non-overlapping 2-year blocks
#' (e.g., 2001–2002, 2003–2004, …, 2023–2024) defined to align with
#' partner-drug resistance surveillance periods, and computes the mean
#' prevalence within each block at each spatial location.
#'
#' The function expects a long-format data frame with spatial coordinates,
#' a year variable, and a prevalence value. Years outside the supported
#' range (2001–2024) are dropped.
#'
#' @param df A data frame containing at least the following columns:
#'   \describe{
#'     \item{x}{Numeric x-coordinate (e.g., longitude).}
#'     \item{y}{Numeric y-coordinate (e.g., latitude).}
#'     \item{t}{Integer year (e.g., 2012).}
#'     \item{p}{Numeric prevalence or probability value.}
#'   }
#'
#' @return A data frame with one row per spatial location and 2-year block,
#'   containing:
#'   \describe{
#'     \item{x}{x-coordinate.}
#'     \item{y}{y-coordinate.}
#'     \item{year_block}{Character label of the 2-year interval
#'       (e.g., "2011-2012").}
#'     \item{mean_p}{Mean prevalence across the 2-year block.}
#'   }
#'
#' @details
#' The 2-year aggregation windows correspond to predefined partner-drug
#' surveillance cutoffs used in resistance monitoring analyses. Prevalence
#' values are averaged within each window using
#' \code{mean(p, na.rm = TRUE)}. Spatial locations are defined by unique
#' combinations of \code{x} and \code{y}. Only complete 2-year blocks
#' between 2001 and 2024 are retained.
#'
#' @examples
#' \dontrun{
#' avg_df <- avg_2_year(prev_long_df)
#' }
#'
#' @export
avg_2_year_pd <- function(df) {
  df <- df %>%
    mutate(
      t = case_when(
        t %in% c(2001, 2002) ~ "2001-2002",
        t %in% c(2003, 2004) ~ "2003-2004",
        t %in% c(2005, 2006) ~ "2005-2006",
        t %in% c(2007, 2008) ~ "2007-2008",
        t %in% c(2009, 2010) ~ "2009-2010",
        t %in% c(2011, 2012) ~ "2011-2012",
        t %in% c(2013, 2014) ~ "2013-2014",
        t %in% c(2015, 2016) ~ "2015-2016",
        t %in% c(2017, 2018) ~ "2017-2018",
        t %in% c(2019, 2020) ~ "2019-2020",
        t %in% c(2021, 2022) ~ "2021-2022",
        t %in% c(2023, 2024) ~ "2023-2024",
        TRUE ~ NA_character_
      )
    ) 
  
  # Average prevalence per pixel inside each 2-year block
  avg_two_year <- df %>%
    group_by(x, y, t) %>%
    summarise(p = mean(p, na.rm = TRUE), .groups = "drop")
  
  return(avg_two_year)
}
