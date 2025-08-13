
#------------------------------------------------
#' @title Empirical logit function
#'   
#' @description The empirical logit is a version of the logit transform that is
#'   shrunk towards a value, in our case 0.5. This avoids issues with infinite
#'   numbers when prevalence is 0% or 100%.
#'   
#' @param x TODO
#'
#' @export

emplogit <- function(x, N) {
  log(x + 0.5) - log(N - x + 0.5)
}

#------------------------------------------------
#' @title Convert Latitude/Longitude to Local Cartesian Coordinates
#'
#' @description Converts geographic coordinates (latitude, longitude) to
#' Cartesian x and y coordinates in kilometres using a local tangent plane 
#' (equirectangular projection) centred on a specified point. 
#' This assumes a flat Earth approximation.
#'
#' @param lon Numeric vector of longitudes in decimal degrees.
#' @param lat Numeric vector of latitudes in decimal degrees.
#' @param lon_centre Longitude (in degrees) for the centre of the projection.
#' @param lat_centre Latitude (in degrees) for the centre of the projection.
#'
#' @return A list with elements `x` and `y` representing Cartesian coordinates
#'   in kilometres. These elements will mirror the structure of the inputs, for
#'   example if `lon` and `lat` are matrices then these will be matrices
#'
#' @export

lonlat_to_cartesian <- function(lon, lat,
                                lon_centre = mean(range(lon)),
                                lat_centre = mean(range(lat))) {
  
  # Earth's radius in kilometres (mean radius)
  R <- 6371
  
  # Convert degrees to radians
  deg_to_rad <- function(deg) deg * pi / 180
  
  lat_rad <- deg_to_rad(lat)
  lon_rad <- deg_to_rad(lon)
  lat_c_rad <- deg_to_rad(lat_centre)
  lon_c_rad <- deg_to_rad(lon_centre)
  
  # Equirectangular projection: simple x/y projection assuming local flatness
  x <- R * (lon_rad - lon_c_rad) * cos(lat_c_rad)
  y <- R * (lat_rad - lat_c_rad)
  
  list(x = x, y = y)
}

#------------------------------------------------
# Convert decimal value to date, e.g. 2010.5 would be the year 2010 plus half
# way through the year in days
#' @noRd

decimal_to_date <- function(decimal_year) {
  year <- floor(decimal_year)
  remainder <- decimal_year - year
  start_date <- as.Date(paste0(year, "-01-01"))
  end_date <- as.Date(paste0(year + 1, "-01-01"))
  start_date + as.integer((end_date - start_date) * remainder)
}
