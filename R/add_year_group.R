#' Add year groups to a data frame
#'
#' Adds a \code{year_group} column based on fixed ranges.
#'
#' @param df A data frame.
#' @param year_col Unquoted column that contains the year (numeric/integer).
#' @return \code{df} with an added character column \code{year_group}.
#' @examples
#' df <- data.frame(year = 2012:2023)
#' add_year_group(df, year)
#' @export
add_year_group <- function(df, year_col = year) {
  dplyr::mutate(
    df,
    year_group = dplyr::case_when(
      {{ year_col }} %in% 2012:2013 ~ "2012-2013",
      {{ year_col }} %in% 2014:2015 ~ "2014-2015",
      {{ year_col }} %in% 2016:2017 ~ "2016-2017",
      {{ year_col }} %in% 2018:2019 ~ "2018-2019",
      {{ year_col }} %in% 2020:2021 ~ "2020-2021",
      {{ year_col }} %in% 2022:2023 ~ "2022-2023",
      TRUE ~ NA_character_
    )
  )
}
