#' Add partner-drug–aligned 2-year surveillance blocks
#'
#' Adds a categorical 2-year time block variable (\code{year_block}) to a data
#' frame based on predefined partner-drug resistance surveillance cutoffs.
#' The resulting year blocks are non-overlapping 2-year intervals
#' (e.g., 2001–2002, 2003–2004, …, 2023–2024).
#'
#' Years outside the supported range (2001–2024) are dropped.
#'
#' @param df A data frame containing a column \code{t} representing the year
#'   of observation as an integer (e.g., 2015).
#'
#' @return A data frame identical to the input but with an additional character
#'   column:
#'   \describe{
#'     \item{year_block}{2-year partner-drug surveillance interval label
#'       (e.g., "2011-2012").}
#'   }
#'
#' @details
#' The year blocks correspond to predefined partner-drug surveillance periods
#' used in antimalarial resistance monitoring. These cutoffs are applied to
#' ensure temporal consistency between prevalence estimates and partner-drug
#' reporting intervals. Rows with years falling outside the defined intervals
#' are removed.
#'
#' @examples
#' \dontrun{
#' df_with_blocks <- add_year_block_pd(prev_df)
#' }
#'
#' @export
add_year_group_pd <- function(df, year_col = year) {
  dplyr::mutate(
    df,
    year_group = dplyr::case_when(
      {{ year_col }} %in% 2001:2002 ~ "2001-2002",
      {{ year_col }} %in% 2003:2004 ~ "2003-2004",
      {{ year_col }} %in% 2005:2006 ~ "2005-2006",
      {{ year_col }} %in% 2007:2008 ~ "2007-2008",
      {{ year_col }} %in% 2009:2010 ~ "2009-2010",
      {{ year_col }} %in% 2011:2012 ~ "2011-2012",
      {{ year_col }} %in% 2013:2014 ~ "2013-2014",
      {{ year_col }} %in% 2015:2016 ~ "2015-2016",
      {{ year_col }} %in% 2017:2018 ~ "2017-2018",
      {{ year_col }} %in% 2019:2020 ~ "2019-2020",
      {{ year_col }} %in% 2021:2022 ~ "2021-2022",
      {{ year_col }} %in% 2023:2024 ~ "2023-2024",
      TRUE ~ NA_character_
    )
  )
}
