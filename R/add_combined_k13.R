#' Add combined K13 prevalence per site-year to a dataset
#'
#' Aggregates all K13 mutations within each site and year to create a
#' combined K13 prevalence ("k13:comb"), then appends these records to
#' the original data set.
#'
#' @param dat A data frame containing site-level mutation prevalence data.
#'   Must include the columns `year`, `longitude`, `latitude`,
#'   `mutation`, `numerator`, `denominator`, and `country_name`.
#'
#' @param mutation_regex Character string; regular expression used to
#'   identify K13 mutations. Default is `"^k13"`.
#'
#' @param combined_label Character string; mutation label assigned to the
#'   aggregated records. Default is `"k13:comb"`.
#'
#' @details
#' For each unique combination of year and site (longitude/latitude),
#' the function:
#' \enumerate{
#'   \item Filters mutations matching `mutation_regex`.
#'   \item Sums numerators across mutations.
#'   \item Uses the maximum denominator across mutations.
#'   \item Computes prevalence as `numerator / denominator * 100`.
#' }
#'
#' The aggregated rows are appended to the original data and sorted by
#' year, longitude, latitude, and mutation.
#'
#' @return
#' A data frame containing the original records plus additional rows
#' representing combined K13 prevalence per site-year.
#'
#' @examples
#' \dontrun{
#' dat_with_k13 <- add_combined_k13(dat)
#' }
#'
#' @importFrom dplyr filter group_by summarise transmute bind_rows arrange
#'
#' @export
add_combined_k13 <- function(
    dat,
    mutation_regex = "^k13",
    combined_label = "k13:comb"
) {
  
  k13_any_site_year <- dat %>%
    dplyr::filter(grepl(mutation_regex, mutation, ignore.case = TRUE)) %>%
    dplyr::group_by(study_id, survey_id, year, longitude, latitude, country_name) %>%
    dplyr::summarise(
      numerator   = sum(numerator, na.rm = TRUE),
      denominator = max(denominator, na.rm = TRUE),
      prevalence  = numerator / denominator * 100,
      .groups     = "drop"
    ) %>%
    dplyr::transmute(
      year,
      study_id,
      survey_id,
      longitude,
      latitude,
      mutation    = combined_label,
      numerator,
      denominator,
      prevalence,
      country_name
    )
  
  dat %>%
    dplyr::bind_rows(k13_any_site_year) %>%
    dplyr::arrange(year, longitude, latitude, mutation)
}
